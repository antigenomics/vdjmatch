"""Annotate one or more query samples against a VDJdb index (built once, reused)."""
from __future__ import annotations

import sys

import polars as pl

from .. import aggregate, evalue, io, match
from ..match.scoring import DEFAULT_SCALE


def annotate_sample(index: "match.VdjdbIndex", sample_path: str, *, scope: str = "1,0,0,1",
                    matrix=None, species: str = "human", with_evalue: bool = True,
                    match_v: bool = False, match_j: bool = False, align: bool = True,
                    threads: int = 0) -> dict[str, pl.DataFrame]:
    """Annotate a single-chain rearrangement sample. Returns {hits, summary, calls} frames.

    Per gene (locus) present in both the sample and the index: fuzzy-search VDJdb, optionally
    compute per-query control-calibrated E-values, and aggregate. Genes without a control just
    skip the E-value (a warning is printed).
    """
    queries = io.read_rearrangement(sample_path)
    gap = DEFAULT_SCALE if matrix is not None else 1
    params = match.search_params(scope, engine="seqtm", matrix=matrix or "",
                                 gap_open=gap, gap_extend=gap)
    hit_frames, eval_frames = [], []
    for gene in index.genes:
        gq = queries.filter(pl.col("locus") == gene)
        if gq.height == 0:
            continue
        h = index.annotate(gq, params, gene=gene, threads=threads,
                           match_v=match_v, match_j=match_j, align=align)
        if h.height:
            hit_frames.append(h)
        if with_evalue:
            try:
                ctrl = evalue.background(gene, species)
                ev = evalue.query_evalues(index.index_for(gene), ctrl,
                                          gq["cdr3"].to_list(), params, threads=threads)
                eval_frames.append(ev)
            except Exception as e:  # control unavailable (e.g. offline HF) -> skip gracefully
                print(f"  [warn] no E-value for {gene}: {e}", file=sys.stderr)

    hits = pl.concat(hit_frames) if hit_frames else index.annotate(queries.head(0),
                                                                    params, gene=index.genes[0])
    evals = pl.concat(eval_frames) if eval_frames else None
    return {"hits": hits, "summary": aggregate.epitope_summary(hits),
            "calls": aggregate.best_call(hits, evals)}
