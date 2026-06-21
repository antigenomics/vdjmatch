"""Paired α/β annotation with a joint control-calibrated E-value (ROADMAP §2.4).

A query pair (cdr3a, cdr3b) *matches* a VDJdb **complex** when its α is within the search ball of
the complex's α CDR3 *and* its β is within the ball of the complex's β CDR3. Under chain
independence the joint null factorizes, ``π0^{αβ} ≈ π0^α · π0^β``, so among the ``N`` paired VDJdb
complexes the expected number of chance joint matches is ``E = N · π0^α · π0^β`` and the joint
enrichment is the Poisson tail ``P(Poisson(E) ≥ n_joint)``. We also report the Fisher-combined
per-chain enrichment. The joint null is far tinier than either chain's, so true pairs get dramatically
smaller E-values than single-chain matching.
"""
from __future__ import annotations

import math
from collections import defaultdict

import polars as pl
from seqtree import Index, SearchParams

try:
    from seqtree.evalue import _poisson_sf
except ImportError:  # pragma: no cover - fallback if the private helper moves
    def _poisson_sf(k: int, lam: float) -> float:
        if k <= 0:
            return 1.0
        term, cum = math.exp(-lam), math.exp(-lam)
        for i in range(1, k):
            term *= lam / i
            cum += term
        return max(0.0, 1.0 - cum)


def _fisher(p_a: float, p_b: float) -> float:
    """Fisher's method combining two independent p-values (chi-square, 4 dof)."""
    stat = -2.0 * (math.log(max(p_a, 1e-300)) + math.log(max(p_b, 1e-300)))
    # survival of chi-square with 4 dof = (1 + x/2) e^{-x/2}
    x = stat / 2.0
    return min(1.0, (1.0 + x) * math.exp(-x))


class PairedVdjdbIndex:
    """Per-chain seqtree indices over the paired VDJdb complexes, with α/β CDR3 → complex maps."""

    def __init__(self, a_idx, b_idx, a_to_cplx, b_to_cplx, cplx_epitope, n_pairs):
        self._a_idx, self._b_idx = a_idx, b_idx
        self._a_to_cplx, self._b_to_cplx = a_to_cplx, b_to_cplx
        self._cplx_epitope = cplx_epitope
        self.n_pairs = n_pairs

    @classmethod
    def build(cls, vdjdb: pl.DataFrame, species: str | None = None) -> "PairedVdjdbIndex":
        """Build from a full VDJdb frame (needs ``complex_id`` pairing; TRA/TRB rows)."""
        if species is not None:
            vdjdb = vdjdb.filter(pl.col("species") == species)
        paired = vdjdb.filter(pl.col("complex_id") != 0)
        a = paired.filter(pl.col("gene") == "TRA").select("complex_id", "cdr3", "epitope")
        b = paired.filter(pl.col("gene") == "TRB").select("complex_id", "cdr3", "epitope")
        complexes = a.join(b, on="complex_id", suffix="_b")  # one row per complex with both chains
        a_uc = complexes.select("cdr3").unique(maintain_order=True)["cdr3"].to_list()
        b_uc = complexes.select(pl.col("cdr3_b")).unique(maintain_order=True)["cdr3_b"].to_list()
        a_to_cplx, b_to_cplx, cplx_epitope = defaultdict(set), defaultdict(set), {}
        for cid, ca, cb, epi in zip(complexes["complex_id"], complexes["cdr3"],
                                    complexes["cdr3_b"], complexes["epitope"]):
            a_to_cplx[ca].add(cid)
            b_to_cplx[cb].add(cid)
            cplx_epitope[cid] = epi
        return cls(Index.build(a_uc, "aa"), Index.build(b_uc, "aa"),
                   dict(a_to_cplx), dict(b_to_cplx), cplx_epitope, complexes.height)

    def _matched_complexes(self, idx: Index, refs: list[str], cdr3_to_cplx: dict,
                           queries: list[str], params: SearchParams, threads: int) -> list[set]:
        """For each query CDR3, the set of VDJdb complex ids whose chain it matches."""
        res = idx.search_batch(queries, params, threads)
        out = []
        for hl in res:
            cplx = set()
            for h in hl:
                cplx |= cdr3_to_cplx.get(refs[h.ref_id], set())
            out.append(cplx)
        return out

    def annotate_pairs(self, pairs: pl.DataFrame, control_a: Index, control_b: Index,
                       params: SearchParams, threads: int = 0) -> pl.DataFrame:
        """Joint E-value per query pair. ``pairs`` needs ``cdr3a, cdr3b`` (+ optional ``epitope``
        ground truth). Returns per-pair n_joint, E, p_joint (Poisson), p_fisher, and the predicted
        epitope (modal among joint matches)."""
        a_refs = [self._a_idx.ref_seq(i) for i in range(len(self._a_idx))]
        b_refs = [self._b_idx.ref_seq(i) for i in range(len(self._b_idx))]
        qa, qb = pairs["cdr3a"].to_list(), pairs["cdr3b"].to_list()
        ca = self._matched_complexes(self._a_idx, a_refs, self._a_to_cplx, qa, params, threads)
        cb = self._matched_complexes(self._b_idx, b_refs, self._b_to_cplx, qb, params, threads)

        Ma, Mb, N = max(1, len(control_a)), max(1, len(control_b)), self.n_pairs
        nca = [sum(1 for h in hl if h.score >= 0) for hl in control_a.search_batch(qa, params, threads)]
        ncb = [sum(1 for h in hl if h.score >= 0) for hl in control_b.search_batch(qb, params, threads)]

        rows = []
        for i in range(pairs.height):
            joint = ca[i] & cb[i]
            n_joint = len(joint)
            pi_a, pi_b = nca[i] / Ma, ncb[i] / Mb
            E = N * pi_a * pi_b
            p_joint = _poisson_sf(n_joint, E if E > 0 else 3.0 / N)
            p_a = _poisson_sf(len(ca[i]), max(N * pi_a, 3.0 / N))
            p_b = _poisson_sf(len(cb[i]), max(N * pi_b, 3.0 / N))
            epis = [self._cplx_epitope[c] for c in joint]
            top = max(set(epis), key=epis.count) if epis else None
            rows.append((qa[i], qb[i], n_joint, len(ca[i]), len(cb[i]), E, p_joint,
                         _fisher(p_a, p_b), top))
        return pl.DataFrame(rows, orient="row", schema=[
            ("cdr3a", pl.Utf8), ("cdr3b", pl.Utf8), ("n_joint", pl.Int64),
            ("n_alpha", pl.Int64), ("n_beta", pl.Int64), ("E", pl.Float64),
            ("p_joint", pl.Float64), ("p_fisher", pl.Float64), ("epitope", pl.Utf8)])
