"""Compact requested detail with exact complete-label randomization maxima."""
from dataclasses import replace

import numpy as np
import polars as pl

from ..aggregate.candidates import PMHC
from ..runner.component_bank import ComponentBank, _mark

_KEYS=['query_id',*PMHC]
_QUERY_NULL_LIMIT_BYTES=2*1024**3


def _query_null_bytes(queries,draws,targets,finite=0):
    # Include coverage, maxima/bin indexes, finite records and sort/copy scratch.
    return ((queries*targets+7)//8+(draws*targets+7)//8+
        32*(queries+draws)+48*finite+16*1024**2)


def _complete(seen,cases):
    if len(seen) and (not np.all(seen[:-1]==255) or
            seen[-1] != ((1<<(cases%8))-1 if cases%8 else 255)):
        raise RuntimeError('randomized scoring did not preserve the complete observed case grid')


def _counts(scores,distribution,draws):
    finite=np.sort(distribution[np.isfinite(distribution)])
    counts=len(finite)-np.searchsorted(finite,scores,side='left')
    counts[np.isnan(scores)]=draws
    return counts,len(finite)


def _tails(frame,reference,query,reference_hits,query_hits,draws):
    return frame.with_columns(pl.Series('reference_exceedances',reference),pl.Series('query_exceedances',query),
        pl.lit(draws).alias('reference_draws'),pl.lit(draws).alias('query_draws'),
        pl.Series('reference_hit_draws',reference_hits),pl.Series('query_hit_draws',query_hits)).with_columns(
        (pl.col('reference_hit_draws')/draws).alias('null_reference_match_rate'),
        (pl.col('query_hit_draws')/draws).alias('null_query_match_rate'),
        ((pl.col('reference_exceedances')+1)/(draws+1)).alias('p_reference'),
        ((pl.col('query_exceedances')+1)/(draws+1)).alias('p_query')).with_columns(
        pl.max_horizontal('p_reference','p_query').alias('p_both'))


def calibrate_compact(tagged,normalizers,test_queries,reference,*,components,combine,models,
                      detail_targets,replicates,report_p,pgen_columns,combine_options):
    """Retain declarations and every reportable label; calibrate all labels/axes.

    Reported-label tails compare that label's score with complete-label maxima
    in both directions. Query-null screening only discards labels that cannot
    pass the final max-of-directions criterion. Fixed-target diagnostic tails
    remain separate. Only finite query-null values are stored, keyed by Pgen
    stratum and target; no observations or null no-hits leave any denominator.
    """
    from ..runner.combined import _check_case_budget,_CASE_LIMIT
    if not np.isfinite(report_p) or not 0<report_p<1:
        raise ValueError('report_p must be finite and between zero and one')
    if not set(_KEYS)<=set(detail_targets.columns):
        raise ValueError('detail targets require query identities and full pMHC')
    requested=detail_targets.select(pl.col('query_id').cast(pl.String),pl.col(PMHC).cast(pl.String)).unique()
    if requested.select('query_id').unique().join(tagged.select('query_id'),on='query_id',how='anti').height:
        raise ValueError('detail targets contain unknown query identities')
    _check_case_budget(requested.height)
    targets=reference.select(PMHC).unique(maintain_order=True).with_row_index('_target_row')
    nq,nt=tagged.height,targets.height
    if not nt:
        raise ValueError('compact calibration requires reference targets')
    test_count=test_queries.height
    if _query_null_bytes(nq,test_count,nt)>_QUERY_NULL_LIMIT_BYTES:
        raise ValueError('compact query-null buffer budget exceeded')
    query_rows=tagged.select('query_id').with_row_index('_source_row')
    bins=tagged.select('pgen_bin').unique(maintain_order=True).with_row_index('_bin')
    qbins=tagged.select('pgen_bin').join(bins,on='pgen_bin',validate='m:1',maintain_order='left')['_bin'].to_numpy()
    test_bins=test_queries.select('pgen_bin').join(bins,on='pgen_bin',validate='m:1',maintain_order='left')['_bin'].to_numpy()
    bin_rows=[np.flatnonzero(test_bins==i) for i in range(bins.height)]
    if any(len(rows)!=replicates for rows in bin_rows):
        raise RuntimeError('query-null scoring did not preserve every draw in every testing group')
    declared=requested.join(query_rows,on='query_id',validate='m:1').join(targets,on=PMHC,
        how='left',nulls_equal=True,validate='m:1')
    missing=declared.filter(pl.col('_target_row').is_null())
    declared=declared.filter(pl.col('_target_row').is_not_null())
    declared_codes=np.sort(declared['_source_row'].to_numpy().astype(np.uint64)*nt+declared['_target_row'].to_numpy())
    test_max=np.full(test_count,np.nan)
    null_parts=[];finite_nulls=0
    observed_max=np.full(nq,np.nan)
    parts=[];retained=0;template=None;query_max=None
    def bank(q,model,plan):
        raw=components(q,model)
        if not isinstance(raw,ComponentBank) or not raw.targets.equals(targets.drop('_target_row')):
            raise RuntimeError('compact scoring requires the complete fixed target population')
        norm=pl.concat(raw.frames(normalizers))
        mapping=plan.join(raw.mapping,on='query_id',how='left',validate='1:1',maintain_order='left')
        if mapping['_query_row'].null_count():
            raise RuntimeError('compact scoring lost source query identities')
        return replace(raw,mapping=mapping,targets=targets),norm
    real_plan=pl.concat([test_queries.select('query_id').with_row_index('_source_row').with_columns(
        pl.lit(False).alias('_observed')),query_rows.with_columns(pl.lit(True).alias('_observed'))])
    all_queries=pl.concat([normalizers.drop('replicate'),test_queries.drop('replicate'),tagged],how='diagonal_relaxed')
    raw,norm=bank(all_queries,reference,real_plan)
    seen_test=np.zeros((test_count*nt+7)//8,dtype=np.uint8)
    seen_observed=np.zeros((nq*nt+7)//8,dtype=np.uint8)
    def consume_real(frame):
        nonlocal query_max,retained,template,finite_nulls
        template=frame.head(0)
        null=frame.filter(~pl.col('_observed'))
        rows=null['_source_row'].to_numpy();groups=null['_target_row'].to_numpy()
        _mark(seen_test,rows,groups,test_count,nt)
        null_scores=null['call_score'].to_numpy()
        np.fmax.at(test_max,rows,null_scores)
        finite=np.isfinite(null_scores)
        finite_nulls+=int(finite.sum())
        if _query_null_bytes(nq,test_count,nt,finite_nulls)>_QUERY_NULL_LIMIT_BYTES:
            raise ValueError('compact query-null buffer budget exceeded')
        if finite.any():
            null_parts.append(pl.DataFrame({'_null_group':test_bins[rows[finite]].astype(np.uint64)*nt+groups[finite],
                '_null_score':null_scores[finite]}))
        observed=frame.filter(pl.col('_observed'))
        if not observed.height:
            return
        if query_max is None:
            _complete(seen_test,test_count*nt)
            query_max=[test_max[rows] for rows in bin_rows]
        rows=observed['_source_row'].to_numpy();groups=observed['_target_row'].to_numpy()
        _mark(seen_observed,rows,groups,nq,nt)
        scores=observed['call_score'].to_numpy()
        np.fmax.at(observed_max,rows,scores)
        counts=np.zeros(len(rows),dtype=np.int64);hits=np.zeros(len(rows),dtype=np.int64)
        for i in np.unique(qbins[rows]):
            at=np.flatnonzero(qbins[rows]==i)
            counts[at],hits[at]=_counts(scores[at],query_max[i],replicates)
        codes=rows.astype(np.uint64)*nt+groups
        pos=np.searchsorted(declared_codes,codes)
        is_declared=np.zeros(len(rows),dtype=bool)
        valid=pos<len(declared_codes)
        is_declared[valid]=declared_codes[pos[valid]]==codes[valid]
        keep=is_declared | (np.isfinite(scores)&((counts+1)/(replicates+1)<=report_p))
        retained+=int(keep.sum())
        if retained>_CASE_LIMIT:
            raise ValueError('compact retained detail budget exceeded')
        if keep.any():
            parts.append(observed.filter(keep).with_columns(pl.Series('declared_target',is_declared[keep]),
                pl.Series('_query_any',counts[keep]),pl.Series('_query_any_hits',hits[keep])))
    combine(raw,norm,reference,all_queries,_consume=consume_real,**combine_options)
    _complete(seen_test,test_count*nt);_complete(seen_observed,nq*nt)
    if query_max is None:
        query_max=[test_max[rows] for rows in bin_rows]
    if missing.height:
        absent=missing.select(_KEYS).with_columns(*[pl.lit(None,dtype=dtype).alias(c)
            for c,dtype in template.schema.items() if c not in _KEYS]).join(query_rows,
                on='query_id',suffix='_known',validate='m:1').with_columns(
            pl.col('_source_row_known').alias('_source_row'),pl.lit(-1).alias('_target_row'),
            pl.lit(True).alias('_observed'),pl.lit(True).alias('declared_target'),
            pl.lit(replicates,dtype=pl.Int64).alias('_query_any'),
            pl.lit(0,dtype=pl.Int64).alias('_query_any_hits'),pl.lit('no_reference_target').alias('combined_status')).drop('_source_row_known')
        parts.append(absent)
    if parts:
        detail=pl.concat(parts,how='diagonal_relaxed').with_row_index('_case')
    else:
        detail=template.with_columns(pl.lit(False).alias('declared_target'),
            pl.lit(0,dtype=pl.Int64).alias('_query_any'),pl.lit(0,dtype=pl.Int64).alias('_query_any_hits')).with_row_index('_case')
    _check_case_budget(detail.height)
    del raw,norm,all_queries,parts,consume_real,seen_test,seen_observed
    sparse=(pl.concat(null_parts,rechunk=True) if null_parts else
        pl.DataFrame(schema={'_null_group':pl.UInt64,'_null_score':pl.Float64}))
    null_parts.clear()
    sparse=sparse.sort('_null_group')
    null_groups=sparse['_null_group'].to_numpy();null_scores=sparse['_null_score'].to_numpy()
    obs_scores=detail['call_score'].to_numpy();obs_rows=detail['_source_row'].to_numpy()
    groups=detail['_target_row'].to_numpy();valid=groups>=0
    codes=obs_rows[valid].astype(np.uint64)*nt+groups[valid].astype(np.uint64)
    order=np.argsort(codes);sorted_codes=codes[order];case_indices=np.flatnonzero(valid)[order]
    query_fixed=np.full(detail.height,replicates,dtype=np.int64);query_hits=np.zeros(detail.height,dtype=np.int64)
    work=detail.filter(pl.col('_target_row')>=0).with_columns(pl.Series('_bin',qbins[obs_rows[valid]]))
    for (b,g),part in work.partition_by(['_bin','_target_row'],as_dict=True).items():
        indices=part['_case'].to_numpy()
        key=b*nt+g
        lo,hi=np.searchsorted(null_groups,key,side='left'),np.searchsorted(null_groups,key,side='right')
        query_fixed[indices],query_hits[indices]=_counts(obs_scores[indices],null_scores[lo:hi],replicates)
    database_query=np.zeros(nq,dtype=np.int64);database_query_hits=np.zeros(nq,dtype=np.int64)
    for i,distribution in enumerate(query_max):
        indices=np.flatnonzero(qbins==i)
        database_query[indices],database_query_hits[indices]=_counts(observed_max[indices],distribution,replicates)
    del sparse,null_groups,null_scores,null_parts,test_max,work,query_max
    fixed=np.where(valid,0,replicates).astype(np.int64);fixed_hits=np.zeros(detail.height,dtype=np.int64)
    any_counts=np.zeros(detail.height,dtype=np.int64)
    database_reference=np.zeros(nq,dtype=np.int64);database_reference_hits=np.zeros(nq,dtype=np.int64)
    observed_queries=pl.concat([normalizers.drop('replicate'),tagged],how='diagonal_relaxed')
    plan=query_rows.with_columns(pl.lit(True).alias('_observed'))
    draws=0
    for model in models:
        raw,norm=bank(observed_queries,model,plan)
        seen=np.zeros((nq*nt+7)//8,dtype=np.uint8);maxima=np.full(nq,np.nan)
        def consume_mock(frame):
            rows=frame['_source_row'].to_numpy();target=frame['_target_row'].to_numpy()
            _mark(seen,rows,target,nq,nt)
            scores=frame['call_score'].to_numpy();np.fmax.at(maxima,rows,scores)
            code=rows.astype(np.uint64)*nt+target
            positions=np.searchsorted(sorted_codes,code);inside=positions<len(sorted_codes)
            matched=np.flatnonzero(inside)
            matched=matched[sorted_codes[positions[matched]]==code[matched]]
            indices=case_indices[positions[matched]];values=scores[matched]
            fixed_hits[indices]+=np.isfinite(values)
            fixed[indices]+=np.isnan(obs_scores[indices]) | (values>=obs_scores[indices])
        combine(raw,norm,model,observed_queries,_consume=consume_mock,**combine_options)
        _complete(seen,nq*nt)
        any_counts+=np.isnan(obs_scores) | (maxima[obs_rows]>=obs_scores)
        database_reference+=np.isnan(observed_max) | (maxima>=observed_max)
        database_reference_hits+=np.isfinite(maxima)
        draws+=1
        del raw,norm,consume_mock,seen,maxima
    if draws!=replicates:
        raise RuntimeError('reference null did not preserve every draw')
    detail=detail.with_columns(pl.Series('p_reference_any',(any_counts+1)/(replicates+1)),
        ((pl.col('_query_any')+1)/(replicates+1)).alias('p_query_any')).with_columns(
        pl.max_horizontal('p_reference_any','p_query_any').alias('p_any')).with_columns(
        (pl.col('call_score').is_not_null()&(pl.col('p_any')<=report_p)).alias('reported_specificity'))
    tails=_tails(detail.select(*_KEYS,pl.col('call_score').alias('score')),fixed,query_fixed,fixed_hits,query_hits,
        replicates).with_columns(pl.lit('epitope').alias('test_scope'))
    keep=detail.select(pl.col('declared_target')|pl.col('reported_specificity')).to_series()
    detail=detail.filter(keep).drop('_case','_source_row','_target_row','_observed','_query_any','_query_any_hits').join(
        tagged.select('query_id','pgen_bin',*pgen_columns),on='query_id',validate='m:1')
    tails=tails.filter(keep)
    database=tagged.select('query_id').with_columns(pl.Series('score',observed_max,nan_to_null=True),
        *[pl.lit(None,dtype=reference.schema[c]).alias(c) for c in PMHC])
    database=_tails(database,database_reference,database_query,database_reference_hits,database_query_hits,
        replicates).with_columns(pl.lit('database').alias('test_scope'))
    return detail,pl.concat([tails,database.select(tails.columns)],how='vertical_relaxed')
