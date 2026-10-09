"""Paired first-hit E-values under the existing independent-chain background model.

The nearest radius is max(alpha cost, beta cost). Historical scans count paired
reference observations; ``include_identity=True`` additionally exposes paired junction
identity so callers can count distinct pairs and use the same unit for target size N.
"""
from __future__ import annotations

import polars as pl
from seqtree import Index, SearchParams
from seqtree.evalue import evalue_result

from .first_hit import _cost_lists, scope


def build_paired_ref(df: pl.DataFrame) -> pl.DataFrame:
    """Validate TRA/TRB linkage and return one row per complete paired observation.

    Integer and rich string IDs work. Missing partners remain unpaired; multiple
    distinct chains at one locus or inconsistent peptide/MHC annotations raise.
    Exact duplicated chain rows do not create extra paired observations.
    """
    required = {"complex_id", "gene", "cdr3", "epitope"}
    if missing := required - set(df.columns):
        raise ValueError(f"paired reference requires {sorted(missing)}")
    nz = df.filter(~pl.col("complex_id").cast(pl.String).fill_null("0").is_in(["", "0"]))
    if nz.filter(~pl.col("gene").is_in(["TRA", "TRB"]) | pl.col("gene").is_null()).height:
        raise ValueError("paired reference supports TRA/TRB loci only")
    if nz.filter(pl.col("epitope").is_null() | (pl.col("epitope") == "")).height:
        raise ValueError("paired reference requires peptide annotations")
    annotations = [c for c in ("epitope", "mhc_a", "mhc_b", "mhc_class", "species") if c in nz.columns]
    if nz.group_by("complex_id").agg(*(pl.col(c).n_unique().alias(c) for c in annotations))\
            .filter(pl.any_horizontal(*(pl.col(c) > 1 for c in annotations))).height:
        raise ValueError("paired complex has conflicting peptide/MHC/species annotations")
    key = ["complex_id", "gene", "cdr3", *[c for c in ("v", "j") if c in nz.columns]]
    chains = nz.unique(subset=key, maintain_order=True)
    if chains.group_by("complex_id", "gene").len().filter(pl.col("len") > 1).height:
        raise ValueError("paired complex has multiple distinct chains at one locus")
    # Keep unsearchable observations in db.load; a paired search needs both junctions.
    chains = chains.filter(pl.col("cdr3").is_not_null() & pl.col("cdr3").str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$"))
    a = chains.filter(pl.col("gene") == "TRA").select("complex_id", pl.col("cdr3").alias("alpha"), *annotations)
    b = chains.filter(pl.col("gene") == "TRB").select("complex_id", pl.col("cdr3").alias("beta"))
    return a.join(b, on="complex_id", how="inner", validate="1:1").sort("complex_id")


def _hits(a_cost_q, b_cost_q, epi, exclude_exact, identities=None):
    """Paired hits sorted deterministically by radius, label and optional pair identity."""
    a_map, b_map = {}, {}
    for c, r in a_cost_q:
        a_map[r] = min(c, a_map.get(r, c))
    for c, r in b_cost_q:
        b_map[r] = min(c, b_map.get(r, c))
    hits = []
    for r in a_map.keys() & b_map.keys():
        ca, cb = a_map[r], b_map[r]
        if exclude_exact and ca == 0 and cb == 0:
            continue
        hit = (max(ca, cb), epi[r])
        hits.append((*hit, identities[r]) if identities is not None else hit)
    return sorted(hits)


def paired_scan(ref: pl.DataFrame, control_a: Index, control_b: Index, pairs, *,
                params: SearchParams | None = None, threads: int = 0,
                exclude_exact: bool = False, include_identity: bool = False):
    """One batched search per chain/reference/control.

    Default hits are historical ``(radius, epitope)`` observation associations.
    ``include_identity`` emits ``(radius, epitope, (alpha, beta))`` for distinct-pair
    counting by :func:`pvalue`; pass a distinct-pair target size N in that mode.
    """
    if not len(control_a) or not len(control_b):
        raise ValueError("paired calibration requires nonempty alpha and beta controls")
    params = params or scope()
    pairs = list(pairs)
    a_idx = Index.build(ref["alpha"].to_list(), "aa")
    b_idx = Index.build(ref["beta"].to_list(), "aa")
    epi = ref["epitope"].to_list()
    identities = list(zip(ref["alpha"], ref["beta"])) if include_identity else None
    qa, qb = [a for a, _ in pairs], [b for _, b in pairs]
    a_cost = _cost_lists(a_idx, qa, params, threads, False, 10000, "paired: alpha", False)
    b_cost = _cost_lists(b_idx, qb, params, threads, False, 10000, "paired: beta", False)
    ca = _cost_lists(control_a, qa, params, threads, False, 10000, "paired: ctrl-a", False)
    cb = _cost_lists(control_b, qb, params, threads, False, 10000, "paired: ctrl-b", False)
    hits = [_hits(a_cost[i], b_cost[i], epi, exclude_exact, identities) for i in range(len(pairs))]
    return hits, [[c for c, _ in x] for x in ca], [[c for c, _ in x] for x in cb]


def _joint_result(n_pair: int, n_ca: int, n_cb: int, N: int, Ma: int, Mb: int) -> dict:
    """Apply each chain's existing rule of three, then the public Poisson helper."""
    if Ma <= 0 or Mb <= 0:
        raise ValueError("paired calibration requires nonempty alpha and beta controls")
    if N < 0 or not 0 <= n_pair <= N or not 0 <= n_ca <= Ma or not 0 <= n_cb <= Mb:
        raise ValueError("paired calibration counts must lie within their reference/control sizes")
    # The product is the existing independent-chain model, with finite-control corrections.
    result = evalue_result(n_pair, (n_ca or 3) * (n_cb or 3), N, Ma * Mb)
    return {"E": result["E"], "p_enrichment": result["p_enrichment"],
            "rule_of_three_alpha": n_ca == 0, "rule_of_three_beta": n_cb == 0}


def pvalue(paired_hits, ctrl_a_costs, ctrl_b_costs, N: int, Ma: int, Mb: int,
           epitope: str | None = None) -> dict:
    """Conditional first-hit calibration at the nearest paired radius.

    N uses observation associations for historical two-tuples, or distinct paired
    junctions for identity-bearing three-tuples; when restricting epitope, use its N.
    An empty control errors. Zero neighbours use 3/M separately for each chain.
    """
    if Ma <= 0 or Mb <= 0:
        raise ValueError("paired calibration requires nonempty alpha and beta controls")
    if N < 0:
        raise ValueError("paired calibration target size must be nonnegative")
    ph = [h for h in paired_hits if epitope is None or h[1] == epitope]
    if not ph:
        return {"radius": None, "n_pair": 0, "E": 0.0, "p_enrichment": 1.0}
    R = min(h[0] for h in ph)
    nearest = [h for h in ph if h[0] <= R]
    if any(len(h) != len(ph[0]) for h in ph) or len(ph[0]) not in (2, 3):
        raise ValueError("paired hits must consistently include or omit pair identity")
    n_p = len({h[2] for h in nearest}) if len(ph[0]) == 3 else len(nearest)
    n_ca = sum(c <= R for c in ctrl_a_costs)
    n_cb = sum(c <= R for c in ctrl_b_costs)
    result = _joint_result(n_p, n_ca, n_cb, N, Ma, Mb)
    return {"radius": R, "n_pair": n_p, "n_control_alpha": n_ca,
            "n_control_beta": n_cb, **result}


def _demo():
    """Self-check: an exact paired self is highly enriched; a random pair is not."""
    from .control import background
    ref = pl.DataFrame({"complex_id": [1, 2, 3], "epitope": ["E1", "E1", "E2"],
                        "alpha": ["CAASYGGSQGNLIF", "CAVRDSNYQLIW", "CAGHTGNQFYF"],
                        "beta": ["CASSLAPGATNEKLFF", "CASSPGQGAYEQYF", "CASSIRSSYEQYF"]})
    ca, cb = background("TRA"), background("TRB")
    pairs = [("CAASYGGSQGNLIF", "CASSLAPGATNEKLFF"),       # exact pair of complex 1
             ("CGGGGGGGGGGGGF", "CHHHHHHHHHHHHF")]         # nonsense -> no hit
    hits, cca, ccb = paired_scan(ref, ca, cb, pairs, exclude_exact=False)
    N = ref.height
    p_self = pvalue(hits[0], cca[0], ccb[0], N, len(ca), len(cb))["p_enrichment"]
    p_rand = pvalue(hits[1], cca[1], ccb[1], N, len(ca), len(cb))["p_enrichment"]
    assert hits[0] and hits[0][0][0] == 0, f"exact pair should be a radius-0 hit: {hits[0]}"
    assert p_self < 1e-3, f"exact pair should be enriched, got p={p_self}"
    assert p_rand == 1.0, f"nonsense pair should not hit, got p={p_rand}"
    print(f"OK  exact-pair p_enrichment={p_self:.2e}  random-pair p={p_rand}")


if __name__ == "__main__":
    _demo()
