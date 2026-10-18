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
    retain_hits=True,
    paired_marginals=False,
):
    """Return hit/candidate/call/count tables plus explicit ingestion diagnostics.

    Requested calibration errors propagate. Summary counts are descriptive match counts;
    sample-level statistical inference is outside this operation.
    """
    if search_mode not in {"fixed", "ball", "unified"}:
        raise ValueError("search_mode must be fixed, ball or unified")
    if paired_marginals and (not paired or search_mode == 'unified'):
        raise ValueError('paired marginal scopes require paired fixed or ball mode')
    if search_mode == "ball" and paired and not paired_marginals:
        raise ValueError("paired graded balls are not supported; use fixed mode with an explicit scope")
    if search_mode == "unified" and (paired or match_v or match_j or not with_evalue or matrix is not None):
        raise ValueError("unified mode requires single-chain controls and its own geometry; paired, V/J restrictions, matrix overrides and --no-evalue are unsupported")
    ann = Annotator(index)
    gap = DEFAULT_SCALE if matrix is not None else 1
    params = match.search_params(
        scope, engine="seqtm", matrix=matrix or "", gap_open=gap, gap_extend=gap
    )
    if paired and paired_marginals:
        queries,report = io.read_cell(sample_path,link=link,source=source,
            sequence_convention=sequence_convention,return_report=True)
        axes,prepared = ann.paired_scope_candidates(queries,scope=params,
            ball=search_mode=='ball',matrix=matrix,match_v=match_v,match_j=match_j,
            control=control,calibrate=with_evalue,species=species,threads=threads,
            progress=progress,exclude_exact=exclude_exact)
        calls = _paired_scope_calls(queries,axes,species=species,radii=range(1,6) if search_mode=='ball' else [None])
        for name,q in zip(('alpha','beta'),prepared):
            flags = q.select('query_id',
                (pl.col('cdr3').is_not_null() & pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$')).fill_null(False).alias(name+'_query_valid'),
                *[pl.col(field).alias(name+'_'+field) for kind in ('v','j')
                    for field in ('original_'+kind,kind+'_imputed') if field in q.columns])
            calls = calls.join(flags,on='query_id',how='left',validate='m:1')
            for kind in ('v','j'):
                if kind+'_imputed' in q.columns:report[kind+'_imputed_'+name+'_rows']=int(q[kind+'_imputed'].sum())
        summary = pl.concat([frame.group_by('radius',*PMHC).agg(
            pl.col('query_id').n_unique().alias('n_queries')).with_columns(pl.lit(axis).alias('axis'))
            for axis,frame in axes.items()],how='diagonal_relaxed')
        return {'candidates':axes['linked'],'TRA.candidates':axes['TRA'],
            'TRB.candidates':axes['TRB'],'calls':calls,'summary':summary,
            'ingestion':pl.DataFrame([report])}
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
        result = ann.paired_candidates(
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
            return_hits=retain_hits,
            progress=progress,
            exclude_exact=exclude_exact,
        )
        hits,c = result if retain_hits else (index.empty_hits(queries.schema['query_id']),result)
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
    if match_v or match_j:
        queries = ann._matching_queries(queries, match_v=match_v, match_j=match_j, species=species)
        for kind,flag in [('v',match_v),('j',match_j)]:
            if flag:report[kind+'_imputed_rows']=int(queries[kind+'_imputed'].sum())
    if search_mode == "unified":
        hits, c = ann.unified_candidates(queries, threads=threads, control=control,
            species=species, exclude_exact=exclude_exact, return_hits=True, distance=unified_distance)
        return {"hits": hits, "candidates": c,
                "calls": _append_calls(queries, queries, c, index.genes, "vdjmatch_"),
                "summary": aggregate.epitope_summary(hits),
                "ingestion": pl.DataFrame([report])}
    if search_mode == "ball":
        result = ann.graded_candidates(
            queries, threads=threads, matrix=matrix, control=control,
            calibrate=with_evalue, species=species, exclude_exact=exclude_exact,
            return_hits=retain_hits, match_v=match_v, match_j=match_j, progress=progress,
        )
        hits, c = result if retain_hits else (index.empty_hits(queries.schema['query_id']), result)
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


def _paired_scope_calls(queries, axes, *, species, radii=None):
    """Adapt independent scope axes to the existing linked-first label reducer."""
    from ..aggregate.paired_labels import paired_labels, AXES
    organism = {'human':'HomoSapiens','mouse':'MusMusculus'}.get(species,species)
    if radii is None:
        radii = sorted(set(value for frame in axes.values() for value in frame['radius'].drop_nulls())) or [None]
    rows = []
    for radius in radii:
        selection = queries.select('query_id')
        evidence = []
        for source,(_,column,maximum,count,label) in zip(('linked','TRA','TRB'),AXES):
            frame = axes[source]
            frame = frame.filter(pl.col('radius')==radius) if radius is not None else frame
            positive = frame.filter(pl.col('ned_score').is_finite() & (pl.col('ned_score')>0))
            maxima = positive.group_by('query_id').agg(pl.col('ned_score').max().alias(maximum))
            won = positive.join(maxima,left_on=['query_id','ned_score'],right_on=['query_id',maximum])
            stats = won.group_by('query_id').agg(pl.col('epitope').n_unique().cast(pl.UInt64).alias(count),
                pl.col('epitope').min().alias(label),pl.len().alias('_pmhc_ties'))
            status = 'linked_selection_status' if source=='linked' else ('alpha' if source=='TRA' else 'beta')+'_selection_status'
            selection = selection.join(maxima,on='query_id',how='left').join(stats,on='query_id',how='left').with_columns(
                pl.col(count).fill_null(0),pl.when(pl.col('_pmhc_ties')>1).then(pl.lit('ambiguous'))
                .when(pl.col('_pmhc_ties')==1).then(pl.lit('ranked_evidence')).otherwise(pl.lit('no_support')).alias(status)).drop('_pmhc_ties')
            evidence.append(won.select('query_id',*PMHC,pl.col('ned_score').alias(column)))
        evidence = pl.concat(evidence,how='diagonal_relaxed')
        # The kernel field names are only the existing reducer's internal schema.
        # Export scope names and no inherited kernel/final-selection probability.
        calls,_ = paired_labels(evidence,selection,species=organism)
        calls = calls.rename({'kernel_score':'linked_score','kernel_alpha_score':'alpha_score',
            'kernel_beta_score':'beta_score'}).with_columns(
                pl.lit(radius,dtype=pl.UInt16).alias('radius'),
                pl.lit('uncalibrated_paired_scope_selection').alias('calibration'),
                pl.lit(None,dtype=pl.Float64).alias('E'),pl.lit(None,dtype=pl.Float64).alias('p_enrichment'),
                pl.col('paired_epitope').alias('vdjmatch_epitope'),
                pl.col('paired_label_status').alias('vdjmatch_status'))
        rows.append(queries.join(calls,on='query_id',how='left',validate='1:m',maintain_order='left'))
    return pl.concat(rows,how='diagonal_relaxed')
