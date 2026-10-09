"""Paired annotation under the existing independent-chain control model.

Calibration counts distinct (alpha junction, beta junction) keys. Independent
reference observations remain available separately in ``reference`` and do not
inflate the null target size or matched-pair counts.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict

import polars as pl
from seqtree import Index, SearchParams
from seqtree.evalue import evalue_result

from ..evalue.paired import _joint_result, build_paired_ref


def _fisher(p_a: float, p_b: float) -> float:
    """Existing Fisher combination of two independent chain enrichment p-values."""
    x = -(math.log(max(p_a, 1e-300)) + math.log(max(p_b, 1e-300)))
    return min(1.0, (1.0 + x) * math.exp(-x))


class PairedVdjdbIndex:
    """Per-chain indices linked to unique paired junctions and their source observations."""

    def __init__(self, a_idx, b_idx, a_to_cplx, b_to_cplx, cplx_epitope, n_pairs):
        self._a_idx, self._b_idx = a_idx, b_idx
        self._a_to_cplx, self._b_to_cplx = a_to_cplx, b_to_cplx
        self._cplx_epitope = cplx_epitope
        self.n_pairs = n_pairs

    @classmethod
    def build(cls, vdjdb: pl.DataFrame, species: str | None = None) -> "PairedVdjdbIndex":
        """Validate paired cardinality; retain observations and index unique paired keys."""
        if species is not None:
            if "species" not in vdjdb.columns or vdjdb.height and vdjdb["species"].null_count() == vdjdb.height:
                raise ValueError("reference has no species metadata; use species=None explicitly")
            vdjdb = vdjdb.filter(pl.col("species") == species)
        observations = build_paired_ref(vdjdb)
        keys = observations.select("alpha", "beta").unique().sort("alpha", "beta").with_row_index("pair_id")
        complexes = observations.join(keys, on=["alpha", "beta"], validate="m:1")
        a_uc, b_uc = sorted(keys["alpha"].unique()), sorted(keys["beta"].unique())
        a_to_cplx, b_to_cplx, labels = defaultdict(set), defaultdict(set), defaultdict(set)
        for pid, ca, cb, epi in complexes.select("pair_id", "alpha", "beta", "epitope").iter_rows():
            a_to_cplx[ca].add(pid)
            b_to_cplx[cb].add(pid)
            labels[pid].add(epi)
        out = cls(Index.build(a_uc, "aa"), Index.build(b_uc, "aa"),
                  dict(a_to_cplx), dict(b_to_cplx), dict(labels), keys.height)
        out.reference = vdjdb
        out.paired_observations = observations
        out.n_observations = observations.height
        return out

    def _matched_complexes(self, idx: Index, refs: list[str], cdr3_to_cplx: dict,
                           queries: list[str], params: SearchParams, threads: int) -> list[set]:
        """For each query, the unique paired keys with an accepted chain match."""
        res = idx.search_batch(queries, params, threads)
        return [set().union(*(cdr3_to_cplx[refs[h.ref_id]] for h in hl)) for hl in res]

    def annotate_pairs(self, pairs: pl.DataFrame, control_a: Index, control_b: Index,
                       params: SearchParams, threads: int = 0) -> pl.DataFrame:
        """Fixed-ball pair/chain counts and calibrated E/p values; ties sort by epitope."""
        Ma, Mb, N = len(control_a), len(control_b), self.n_pairs
        if not Ma or not Mb:
            raise ValueError("paired calibration requires nonempty alpha and beta controls")
        if not {"cdr3a", "cdr3b"}.issubset(pairs.columns):
            raise ValueError("query pairs require cdr3a and cdr3b")
        a_refs = [self._a_idx.ref_seq(i) for i in range(len(self._a_idx))]
        b_refs = [self._b_idx.ref_seq(i) for i in range(len(self._b_idx))]
        qa, qb = pairs["cdr3a"].to_list(), pairs["cdr3b"].to_list()
        ca = self._matched_complexes(self._a_idx, a_refs, self._a_to_cplx, qa, params, threads)
        cb = self._matched_complexes(self._b_idx, b_refs, self._b_to_cplx, qb, params, threads)
        nca = [len(hl) for hl in control_a.search_batch(qa, params, threads)]
        ncb = [len(hl) for hl in control_b.search_batch(qb, params, threads)]
        rows = []
        for i in range(pairs.height):
            joint = ca[i] & cb[i]
            result = _joint_result(len(joint), nca[i], ncb[i], N, Ma, Mb)
            p_a = evalue_result(len(ca[i]), nca[i], N, Ma)["p_enrichment"]
            p_b = evalue_result(len(cb[i]), ncb[i], N, Mb)["p_enrichment"]
            votes = Counter(epi for pid in joint for epi in self._cplx_epitope[pid])
            top = min(votes, key=lambda e: (-votes[e], e)) if votes else None
            rows.append((qa[i], qb[i], len(joint), len(ca[i]), len(cb[i]), result["E"],
                         result["p_enrichment"], _fisher(p_a, p_b), top,
                         result["rule_of_three_alpha"], result["rule_of_three_beta"]))
        return pl.DataFrame(rows, orient="row", schema=[
            ("cdr3a", pl.String), ("cdr3b", pl.String), ("n_joint", pl.Int64),
            ("n_alpha", pl.Int64), ("n_beta", pl.Int64), ("E", pl.Float64),
            ("p_joint", pl.Float64), ("p_fisher", pl.Float64), ("epitope", pl.String),
            ("rule_of_three_alpha", pl.Boolean), ("rule_of_three_beta", pl.Boolean)])
