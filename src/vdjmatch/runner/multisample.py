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
    top_k=10,
):
    """Return hit/candidate/call/count tables plus explicit ingestion diagnostics.

    Requested calibration errors propagate. Summary counts are descriptive match counts;
    sample-level statistical inference is outside this operation.
    """
    if search_mode not in {"fixed", "ball", "ranked"}:
        raise ValueError("search_mode must be fixed, ball or ranked")
    if search_mode != "fixed" and (match_v or match_j):
        raise ValueError("graded/ranked modes use sequence-only predicates; filter the reference explicitly")
    if search_mode == "ball" and paired:
        raise ValueError("paired graded balls are not supported; use ranked paired mode")
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
        if search_mode == "ranked":
            if with_evalue and control is None:
                from ..evalue.control import background
                organism = ann._control_species(["TRA", "TRB"], species)
                control = {gene: background(gene, organism) for gene in ("TRA", "TRB")}
            h, c, stats = ann.ranked_candidates(
                queries, k=top_k, threads=threads, matrix=matrix,
                control=control if with_evalue else None, exclude_exact=exclude_exact, paired=True,
            )
            _, qa = _prepare(queries, "cdr3a", locus="TRA")
            _, qb = _prepare(queries, "cdr3b", locus="TRB")
            q = qa.join(qb.select("query_id"), on="query_id")
            return _ranked_tables(index, queries, q, h, c, stats, report, paired=True)
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
    if search_mode == "ranked":
        if with_evalue and control is None:
            from ..evalue.control import background
            active = sorted(set(queries["locus"]) & set(index.genes))
            organism = ann._control_species(active, species) if active else species
            control = {gene: background(gene, organism) for gene in active}
        h, c, stats = ann.ranked_candidates(
            queries, k=top_k, threads=threads, matrix=matrix,
            control=control if with_evalue else None, exclude_exact=exclude_exact,
        )
        return _ranked_tables(index, queries, queries, h, c, stats, report)
    if search_mode == "ball":
        hits, c = ann.graded_candidates(
            queries, threads=threads, matrix=matrix, control=control,
            calibrate=with_evalue, species=species, exclude_exact=exclude_exact,
            return_hits=True,
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


def _ranked_tables(index, queries, q, hits, candidates, statistics, report, *, paired=False):
    calls = _append_calls(queries, q, candidates, index.genes, "vdjmatch_", paired=paired)
    calls = calls.join(statistics, on="query_id", how="left", validate="1:1")
    return {"hits": hits, "candidates": candidates, "calls": calls,
            "global_statistics": statistics,
            "summary": candidates.group_by(PMHC).agg(pl.col("query_id").n_unique().alias("n_queries")),
            "ingestion": pl.DataFrame([report])}
