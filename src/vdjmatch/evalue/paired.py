"""Paired-chain (α/β) first-hit E-value.

A paired hit is a reference clonotype matching the query on **both** chains within budget; its radius is
``R = max(cost_α, cost_β)`` (both chains within ``R``). The first paired hit is the reference minimising
that max-cost. Under chain independence the joint null factorises — ``π0^{αβ} ≈ π0^α · π0^β`` — so the
expected paired matches within ``R`` is ``E = N · (n_ctrl_α(R)/M_α) · (n_ctrl_β(R)/M_β)`` and
``p_enrichment = Poisson P(X ≥ n_pair(R) | E)``. Mirrors :mod:`vdjmatch.evalue.first_hit` per chain.
"""
from __future__ import annotations

import polars as pl
from seqtree import Index, SearchParams

from seqtree.evalue import _poisson_sf

from .first_hit import _cost_lists, scope


def build_paired_ref(df: pl.DataFrame) -> pl.DataFrame:
    """VDJdb long frame (needs ``complex_id`` linking chains) -> one row per paired complex.

    Returns ``complex_id, alpha, beta, epitope`` for every complex carrying both a TRA and a TRB.
    """
    nz = df.filter(pl.col("complex_id") != 0)
    a = (nz.filter(pl.col("gene") == "TRA").group_by("complex_id")
           .agg(pl.col("cdr3").first().alias("alpha"), pl.col("epitope").first().alias("epitope")))
    b = (nz.filter(pl.col("gene") == "TRB").group_by("complex_id")
           .agg(pl.col("cdr3").first().alias("beta")))
    return a.join(b, on="complex_id", how="inner").select("complex_id", "alpha", "beta", "epitope")


def _hits(a_cost_q, b_cost_q, epi, exclude_exact):
    """Per-query paired hits: refs within budget on BOTH chains -> cost-sorted [(R=max-cost, epitope)]."""
    a_map: dict[int, int] = {}
    for c, r in a_cost_q:
        if r not in a_map or c < a_map[r]:
            a_map[r] = c
    b_map: dict[int, int] = {}
    for c, r in b_cost_q:
        if r not in b_map or c < b_map[r]:
            b_map[r] = c
    hits = []
    for r in a_map.keys() & b_map.keys():
        ca, cb = a_map[r], b_map[r]
        if exclude_exact and ca == 0 and cb == 0:               # exact paired self -> leakage
            continue
        hits.append((max(ca, cb), epi[r]))
    return sorted(hits, key=lambda x: x[0])


def paired_scan(ref: pl.DataFrame, control_a: Index, control_b: Index, pairs, *,
                params: SearchParams | None = None, threads: int = 0, exclude_exact: bool = False):
    """One wide search per chain. ``pairs`` = list of ``(alpha, beta)``. Returns
    ``(paired_hits, ctrl_a_costs, ctrl_b_costs)``: ``paired_hits[q]`` cost-sorted ``[(R, epitope)]``,
    ``ctrl_*_costs[q]`` the per-chain control edit-cost lists."""
    params = params or scope()
    a_idx = Index.build(ref["alpha"].to_list(), "aa")
    b_idx = Index.build(ref["beta"].to_list(), "aa")
    epi = ref["epitope"].to_list()
    qa = [a for a, _ in pairs]
    qb = [b for _, b in pairs]
    a_cost = _cost_lists(a_idx, qa, params, threads, False, 10000, "paired: alpha", False)
    b_cost = _cost_lists(b_idx, qb, params, threads, False, 10000, "paired: beta", False)
    ca = _cost_lists(control_a, qa, params, threads, False, 10000, "paired: ctrl-a", False)
    cb = _cost_lists(control_b, qb, params, threads, False, 10000, "paired: ctrl-b", False)
    hits = [_hits(a_cost[i], b_cost[i], epi, exclude_exact) for i in range(len(pairs))]
    return (hits, [[c for c, _ in x] for x in ca], [[c for c, _ in x] for x in cb])


def pvalue(paired_hits, ctrl_a_costs, ctrl_b_costs, N: int, Ma: int, Mb: int,
           epitope: str | None = None) -> dict:
    """Joint first-hit E-value at the nearest paired (optionally ``epitope``-restricted) hit.
    ``N`` = number of paired references (epitope's count when ``epitope`` set); ``Ma``/``Mb`` control
    sizes. ``E = N·(n_ctrl_α(R)/Ma)·(n_ctrl_β(R)/Mb)``, ``p = P(X ≥ n_pair(R) | E)``."""
    ph = paired_hits if epitope is None else [(r, e) for r, e in paired_hits if e == epitope]
    if not ph:
        return {"radius": None, "n_pair": 0, "E": 0.0, "p_enrichment": 1.0}
    R = ph[0][0]
    n_p = sum(1 for r, _ in ph if r <= R)
    n_ca = sum(1 for c in ctrl_a_costs if c <= R)
    n_cb = sum(1 for c in ctrl_b_costs if c <= R)
    E = N * (n_ca / Ma) * (n_cb / Mb)
    p = _poisson_sf(n_p, E) if E > 0 else (0.0 if n_p > 0 else 1.0)
    return {"radius": R, "n_pair": n_p, "E": E, "p_enrichment": p}


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
