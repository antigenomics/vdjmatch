"""Shared CLI/API candidate reduction, one sample and one native batch per locus."""

from __future__ import annotations

import polars as pl

from .. import aggregate, io, match
from ..api import Annotator, _append_calls, _prepare
from ..aggregate.candidates import PMHC
from ..match.scoring import DEFAULT_SCALE


def annotate_sample(
    index,
    sample_path,
    *,
    scope="1,0,0,1",
    matrix=None,
    species="human",
    with_evalue=True,
    match_v=False,
    match_j=False,
    align=True,
    threads=1,
    progress=False,
    source="auto",
    sequence_convention=None,
    control=None,
    paired=False,
    link=None,
):
    """Return hit/candidate/call/count tables plus explicit ingestion diagnostics.

    Requested calibration errors propagate. Summary counts are descriptive match counts;
    sample-level statistical inference is outside this operation.
    """
    ann = Annotator(index)
    gap = DEFAULT_SCALE if matrix is not None else 1
    params = match.search_params(
        scope, engine="seqtm", matrix=matrix or "", gap_open=gap, gap_extend=gap
    )
    if paired:
        queries, report = io.read_cell(
            sample_path,
            link=link,
            source=source,
            sequence_convention=sequence_convention,
            return_report=True,
        )
        # Joint calibration is handled explicitly by the paired API.
        c = ann.paired_candidates(
            queries,
            cdr3a="cdr3a",
            cdr3b="cdr3b",
            scope=params,
            threads=threads,
            control=control,
            calibrate=with_evalue,
            species=species,
        )
        _, qa = _prepare(queries, "cdr3a", locus="TRA")
        _, qb = _prepare(queries, "cdr3b", locus="TRB")
        q = qa.join(qb.select("query_id"), on="query_id")
        calls = _append_calls(queries, q, c, index.genes, "vdjmatch_", paired=True)
        return {
            "hits": c,
            "candidates": c,
            "calls": calls,
            "summary": c.group_by(PMHC).agg(
                pl.col("query_id").n_unique().alias("n_query_pairs")
            ),
            "ingestion": pl.DataFrame([report]),
        }
    queries, report = io.read_rearrangement(
        sample_path,
        source=source,
        sequence_convention=sequence_convention,
        return_report=True,
    )
    hits, c = ann._evidence(
        queries,
        params,
        threads=threads,
        match_v=match_v,
        match_j=match_j,
        align=align,
        calibrate=with_evalue,
        species=species,
        control=control,
        score_scale=400.0 if matrix is not None else 1.0,
    )
    calls = _append_calls(queries, queries, c, index.genes, "vdjmatch_")
    return {
        "hits": hits,
        "candidates": c,
        "calls": calls,
        "summary": aggregate.epitope_summary(hits),
        "ingestion": pl.DataFrame([report]),
    }
