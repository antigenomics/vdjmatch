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
    exclude_exact=False,
    search_mode="fixed",
    unified_distance="edit",
):
    """Return hit/candidate/call/count tables plus explicit ingestion diagnostics.

    Requested calibration errors propagate. Summary counts are descriptive match counts;
    sample-level statistical inference is outside this operation.
    """
    if search_mode not in {"fixed", "ball", "unified"}:
        raise ValueError("search_mode must be fixed, ball or unified")
    if search_mode == "ball" and paired:
        raise ValueError("paired graded balls are not supported; use fixed mode with an explicit scope")
    if search_mode == "unified" and (paired or match_v or match_j or not with_evalue or matrix is not None):
        raise ValueError("unified mode requires single-chain controls and its own geometry; paired, V/J restrictions, matrix overrides and --no-evalue are unsupported")
    ann = Annotator(index)
    gap = DEFAULT_SCALE if matrix is not None else 1
    params = match.search_params(
        scope, engine="seqtm", matrix=matrix or "", gap_open=gap, gap_extend=gap
    )
    if paired:
        if match_v or match_j:
            raise ValueError(
                "paired matching does not support match_v/match_j; "
                "use single-chain matching or disable these flags"
            )
        queries, report = io.read_cell(
            sample_path,
            link=link,
            source=source,
            sequence_convention=sequence_convention,
            return_report=True,
        )
        # Joint calibration is handled explicitly by the paired API.
        hits, c = ann.paired_candidates(
            queries,
            cdr3a="cdr3a",
            cdr3b="cdr3b",
            scope=params,
            threads=threads,
            control=control,
            calibrate=with_evalue,
            species=species,
            align=align,
            score_scale=400.0 if matrix is not None else 1.0,
            return_hits=True,
            progress=progress,
            exclude_exact=exclude_exact,
        )
        _, qa = _prepare(queries, "cdr3a", locus="TRA")
        _, qb = _prepare(queries, "cdr3b", locus="TRB")
        q = qa.join(qb.select("query_id"), on="query_id")
        calls = _append_calls(queries, q, c, index.genes, "vdjmatch_", paired=True)
        return {
            "hits": hits,
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
    if search_mode == "unified":
        hits, c = ann.unified_candidates(queries, threads=threads, control=control,
            species=species, exclude_exact=exclude_exact, return_hits=True, distance=unified_distance)
        return {"hits": hits, "candidates": c,
                "calls": _append_calls(queries, queries, c, index.genes, "vdjmatch_"),
                "summary": aggregate.epitope_summary(hits),
                "ingestion": pl.DataFrame([report])}
    if search_mode == "ball":
        hits, c = ann.graded_candidates(
            queries, threads=threads, matrix=matrix, control=control,
            calibrate=with_evalue, species=species, exclude_exact=exclude_exact,
            return_hits=True, match_v=match_v, match_j=match_j,
        )
        calls = pl.concat([
            _append_calls(queries, queries, c.filter(pl.col("radius") == radius),
                          index.genes, "vdjmatch_").with_columns(pl.lit(radius).alias("radius"))
            for radius in range(1, 6)
        ], how="diagonal_relaxed")
        return {"hits": hits, "candidates": c, "calls": calls,
                "summary": c.group_by("radius", *PMHC).agg(pl.col("query_id").n_unique().alias("n_queries")),
                "ingestion": pl.DataFrame([report])}
    hits, c = ann._evidence(
        queries,
        params,
        threads=threads,
        match_v=match_v,
        match_j=match_j,
        align=align,
        progress=progress,
        calibrate=with_evalue,
        species=species,
        control=control,
        score_scale=400.0 if matrix is not None else 1.0,
        exclude_exact=exclude_exact,
    )
    calls = _append_calls(queries, queries, c, index.genes, "vdjmatch_")
    return {
        "hits": hits,
        "candidates": c,
        "calls": calls,
        "summary": aggregate.epitope_summary(hits),
        "ingestion": pl.DataFrame([report]),
    }
