"""Two directional empirical tails for a frozen score and testing scope."""
from __future__ import annotations

import polars as pl

_KEYS = ['query_id', 'group_id']


def _validate_scores(frame: pl.DataFrame) -> pl.DataFrame:
    if not {'score', *_KEYS}.issubset(frame.columns):
        raise ValueError('tail inputs require query_id, group_id and score')
    if any(frame[key].null_count() for key in _KEYS):
        raise ValueError('tail case keys cannot be missing')
    if not frame.schema['score'].is_numeric() and frame.schema['score'] != pl.Null:
        raise ValueError('tail scores must be numeric or null for no-hit')
    frame = frame.with_columns(pl.col('score').cast(pl.Float64))
    if frame.filter(pl.col('score').is_not_null() & ~pl.col('score').is_finite()).height:
        raise ValueError('tail scores must be finite or null for no-hit')
    return frame


def _tail(observed: pl.DataFrame, null: pl.DataFrame, direction: str) -> pl.DataFrame:
    null = _validate_scores(null)
    if 'replicate' not in null.columns or not null.height or null['replicate'].null_count():
        raise ValueError(f'{direction} null requires nonempty, non-null replicate identities')
    if null.select(*_KEYS, 'replicate').is_duplicated().any():
        raise ValueError(f'{direction} null contains duplicate case/replicate rows')
    cases = observed.select(_KEYS)
    supplied = null.select(_KEYS).unique()
    if supplied.join(cases, on=_KEYS, how='anti').height or cases.join(supplied, on=_KEYS, how='anti').height:
        raise ValueError(f'{direction} null case keys must match observed cases exactly')
    if null.group_by('replicate').len().filter(pl.col('len') != observed.height).height:
        raise ValueError(f'{direction} null must include every case for every replicate')
    joined = null.select(*_KEYS, 'score').join(
        observed.select(*_KEYS, pl.col('score').alias('_observed')), on=_KEYS, validate='m:1')
    counts = joined.group_by(_KEYS).agg(
        pl.len().cast(pl.Int64).alias(f'{direction}_draws'),
        (pl.col('_observed').is_null() | (pl.col('score') >= pl.col('_observed')).fill_null(False))
        .sum().cast(pl.Int64).alias(f'{direction}_exceedances'))
    return counts.with_columns(
        ((pl.col(f'{direction}_exceedances') + 1) / (pl.col(f'{direction}_draws') + 1)).alias(f'p_{direction}'))


def directional_tails(observed: pl.DataFrame, reference_null: pl.DataFrame,
                      query_null: pl.DataFrame) -> pl.DataFrame:
    """Return reference/query randomization tails and their maximum per case.

    Frames have ``query_id,group_id,score``; null frames additionally have a
    ``replicate`` identity and a complete row for every observed case on every
    draw. Higher scores are better, ties count as exceedances, and null no-hits
    (``score=None``) stay in the denominator. An observed no-hit has tail one.
    Each tail is ``(1 + exceedances)/(1 + draws)``. The returned ``p_both`` is
    their maximum, a test against the union of the two null hypotheses whenever
    each applicable marginal test is valid; it requires no independence.

    The caller must freeze the same statistic and scope across observed and
    randomized scores and own the preselection event. A component tail cannot
    be reused for a selected winning group or another score. These tails are
    conditional randomization tests under exchangeability, not binding evidence
    probabilities or posterior probabilities.
    """
    reserved = {'p_reference', 'p_query', 'p_both', 'reference_exceedances',
                'query_exceedances', 'reference_draws', 'query_draws'}
    if reserved.intersection(observed.columns):
        raise ValueError('observed contains reserved tail output columns')
    normalized = _validate_scores(observed)
    if not normalized.height or normalized.select(_KEYS).is_duplicated().any():
        raise ValueError('observed requires nonempty unique case keys')
    reference = _tail(normalized, reference_null, 'reference')
    query = _tail(normalized, query_null, 'query')
    return observed.join(reference, on=_KEYS, validate='1:1', maintain_order='left').join(
        query, on=_KEYS, validate='1:1', maintain_order='left').with_columns(
            pl.max_horizontal('p_reference', 'p_query').alias('p_both'))


def _null_group_axes(frame, nulls, columns):
    """Align nullable metadata and lossless scalar widths before key joins."""
    axes=[part.select(columns) for part in (frame,nulls)]
    for column in columns:
        left,right=(part.schema[column] for part in axes)
        if left==right:
            continue
        if left==pl.Null or right==pl.Null:
            dtype=right if left==pl.Null else left
        elif left.is_integer() and right.is_integer():
            dtype=pl.Int128
        elif left.is_float() and right.is_float():
            dtype=pl.Float64
        elif all(t.base_type() in (pl.String,pl.Categorical,pl.Enum) for t in (left,right)):
            dtype=pl.String
        else:
            raise ValueError('null group column types must match')
        axes=[part.with_columns(pl.col(column).cast(dtype)) for part in axes]
    return axes


def null_tail_probabilities(frame, nulls, *, higher_is_better=True, score_column="score",
                            group_columns=(), _details=True):
    """Tail of the existing selected score against each frozen null population.

    P=(1+number of null scores at least as extreme)/(1+all null queries),
    retaining ties and no-hits. It is a conservative rank P-value if query and
    null scores are exchangeable under the same complete pipeline. E=k/N is
    the observed expected false calls per query at that score threshold, not a
    confidence bound. Maxima retain the most conservative represented source.
    Unsupported queries have P=1/E=0 because they make no call.
    Optional groups keep their own null populations, including nullable pMHC
    keys. Internal normalization can omit unused per-source detail tables.
    """
    import numpy as np
    import polars as pl
    if not isinstance(score_column,str) or not score_column.strip():
        raise ValueError('score_column must name a numeric input column')
    if {'p_null','E_null'} & set(frame.columns):
        raise ValueError('null-tail output columns already exist')
    if (isinstance(group_columns, str) or len(set(group_columns)) != len(group_columns) or
            any(not isinstance(c, str) or not c.strip() or c in
                {'query_id', score_column, 'null_source', '_tail_group'} for c in group_columns)):
        raise ValueError('distinct group column names required')
    group_columns = list(group_columns)
    for part in (frame, nulls):
        if missing := {'query_id',score_column,*group_columns} - set(part.columns):
            raise ValueError('missing columns: '+', '.join(sorted(missing)))
        if part['query_id'].null_count() or part.select('query_id',*group_columns).n_unique()!=part.height:
            raise ValueError('query_id must be unique and non-null within each group; reduce candidate labels first')
        if not part[score_column].cast(pl.Float64).drop_nulls().is_finite().all():
            raise ValueError('non-null scores must be finite')
    if not nulls.height or not {'outcome','null_source'} <= set(nulls.columns):
        raise ValueError('nonempty null cases require outcome and null_source')
    if nulls['outcome'].null_count() or not nulls['outcome'].eq('null').all():
        raise ValueError('only null outcomes define null tails')
    if nulls['null_source'].null_count() or nulls['null_source'].cast(pl.String).str.strip_chars().eq('').any():
        raise ValueError('every null query requires a nonempty null_source')
    sign=1 if higher_is_better else -1
    values=sign*frame[score_column].cast(pl.Float64).to_numpy()
    missing=np.isnan(values)
    if group_columns:
        axes=_null_group_axes(frame,nulls,group_columns)
        groups=pl.concat(axes).unique(maintain_order=True).with_row_index('_tail_group')
        query_group,null_group=[part.join(groups,on=group_columns,how='left',nulls_equal=True,
            validate='m:1',maintain_order='left')['_tail_group'].to_numpy() for part in axes]
        group_count=groups.height
    else:
        query_group=np.zeros(frame.height,dtype=np.uint32)
        null_group=np.zeros(nulls.height,dtype=np.uint32)
        group_count=1
    covered=np.bincount(null_group,minlength=group_count)[query_group]>0
    if np.any(~covered & ~missing):
        raise ValueError('missing null calibration group')
    queries=np.rec.fromarrays([query_group,values],names=['group','score']) if group_columns else None
    p=np.zeros(frame.height); e=np.zeros(frame.height); rows=[]
    for source in sorted(nulls['null_source'].unique().to_list()):
        mask=(nulls['null_source']==source).to_numpy()
        source_groups=null_group[mask]
        scores=sign*nulls.filter(pl.Series(mask))[score_column].cast(pl.Float64).to_numpy()
        valid=~np.isnan(scores)
        if group_columns:
            ordered=np.sort(np.rec.fromarrays([source_groups[valid],scores[valid]],names=['group','score']),
                order=['group','score'])
            # Lexicographic lower bounds count strictly smaller scores within
            # each group. Its end minus that bound retains ties as exceedances.
            counts=np.searchsorted(ordered['group'],query_group,side='right')-np.searchsorted(ordered,queries,side='left')
        else:
            ordered=np.sort(scores[valid])
            counts=len(ordered)-np.searchsorted(ordered,values,side='left')
        counts[missing]=0
        draws=np.bincount(source_groups,minlength=group_count)[query_group]
        present=draws>0
        tail=np.where(missing,1.,(counts+1)/(draws+1))
        rate=np.divide(counts,draws,out=np.zeros(frame.height),where=present)
        p[present]=np.maximum(p[present],tail[present]); e[present]=np.maximum(e[present],rate[present])
        if _details:
            rows.append(frame.select('query_id',*group_columns).filter(pl.Series(present)).with_columns(
                pl.lit(source).alias('null_source'),
                pl.Series('null_queries',draws[present],dtype=pl.Int32 if int(mask.sum())<2**31 else pl.Int64),
                pl.Series('null_tail_count',counts[present]),pl.Series('p_null',tail[present]),pl.Series('E_null',rate[present])))
    p[~covered]=1.
    return frame.with_columns(pl.Series('p_null',p),pl.Series('E_null',e)),pl.concat(rows) if _details else None
