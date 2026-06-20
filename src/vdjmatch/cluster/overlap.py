"""Pairwise sample overlap via fuzzy CDR3 matching (ports the legacy ``cluster`` routine)."""
from __future__ import annotations

import polars as pl
from seqtree import pairwise_batch

from ..match.scope import search_params


def overlap(a: list[str], b: list[str] | None = None, scope: str = "1,0,0,1",
            matrix=None, threads: int = 0) -> pl.DataFrame:
    """Fuzzy-matching pairs between two CDR3 sets (or within one if ``b`` is None).

    Returns a long frame ``a_idx, a_cdr3, b_idx, b_cdr3, score, n_subs`` of all within-scope
    pairs. Within-set mode (``b is None``) drops the trivial self-pairs (i == j).
    """
    within = b is None
    bb = a if within else b
    params = search_params(scope, engine="seqtm", matrix=matrix or "")
    res = pairwise_batch(a, bb, params, "aa", threads)
    rows = [(i, a[i], h.ref_id, bb[h.ref_id], h.score, h.n_subs)
            for i, hl in enumerate(res) for h in hl if not (within and i == h.ref_id)]
    return pl.DataFrame(rows, orient="row",
                        schema=["a_idx", "a_cdr3", "b_idx", "b_cdr3", "score", "n_subs"])


def overlap_metrics(a: list[str], b: list[str], scope: str = "1,0,0,1",
                    threads: int = 0) -> dict[str, float]:
    """Summary overlap metrics between two repertoires: number of matched pairs, and the
    fraction of each set with at least one fuzzy match in the other."""
    pairs = overlap(a, b, scope=scope, threads=threads)
    a_hit = pairs["a_idx"].n_unique() if pairs.height else 0
    b_hit = pairs["b_idx"].n_unique() if pairs.height else 0
    return {"pairs": float(pairs.height),
            "frac_a_matched": a_hit / len(a) if a else 0.0,
            "frac_b_matched": b_hit / len(b) if b else 0.0}
