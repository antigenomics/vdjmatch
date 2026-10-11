"""Experimental J usage conditional on V gene and full junction length.

Population frequency evidence only: no V/length marginal multiplier, transferred
P-value or biological confidence. Linked pairs use joint J and stratum categories.
"""
import polars as pl


def conditional_j_usage(queries, reference, background, *, paired=False, exclude_exact=False):
    """Return keyed residual J ratio using one pooled pseudo-observation per stratum.

    Frames contain query_id (queries only) and seq0/v0/j0 [seq1/v1/j1].
    Distinct full junction/resolved V keys are the observation unit. Conflicting
    or absent J annotations stay explicit categories; query missing J is neutral.
    """
    axes=range(2 if paired else 1)
    seq=[f'seq{i}' for i in axes];v=[f'v{i}' for i in axes];j=[f'j{i}' for i in axes]
    strata=[*v,*[f'length{i}' for i in axes]]
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('unique non-null query IDs required')
    def features(frame):
        return frame.with_columns(*[pl.col(c).alias('_identity_'+c) for c in v]).with_columns(
            *[pl.col(c).str.replace(r'\*.*$', '') for c in v],
            *[pl.when(pl.col(c).str.contains(r'^TR[AB]J[0-9]+(?:-[0-9]+)?(?:\*[0-9]+)?$').fill_null(False))
              .then(pl.col(c).str.replace(r'\*.*$', '')).otherwise(pl.lit('__missing__')).alias(c) for c in j],
            *[pl.col(c).str.len_chars().alias(f'length{i}') for i,c in enumerate(seq)])
    q=features(queries)
    def population(frame):
        # Never select a first J annotation for a conflicting search key.
        return features(frame).group_by([*seq,*['_identity_'+c for c in v]],maintain_order=True).agg(
            *[pl.col(c).first() for c in v],
            *[pl.when(pl.col(c).n_unique()==1).then(pl.col(c).first())
              .otherwise(pl.lit('__ambiguous__')).alias(c) for c in j],
            *[pl.col(c).first() for c in strata if c not in v])
    ref,bg=population(reference),population(background)
    # One unseen-J category keeps predictions independent of query batch contents.
    categories=pl.concat([ref.select(j),bg.select(j)]).unique().height+1
    def predictive(pop,prefix):
        out=q.select('query_id',*seq,*strata,*j)
        for keys,col in [(strata,'stratum'),([*strata,*j],'joint'),(j,'global')]:
            counts=pop.group_by(keys).len().rename({'len':col})
            out=out.join(counts,on=keys,how='left',validate='m:1',maintain_order='left',nulls_equal=True)
        out=out.with_columns(pl.col('stratum','joint','global').fill_null(0).cast(pl.Int64),pl.lit(pop.height).alias('population'))
        if exclude_exact and pop.height:
            exact=q.select('query_id',*seq,*strata,*j).join(pop,on=seq,how='inner',suffix='_ref',nulls_equal=True)
            sx=pl.all_horizontal([pl.col(c).eq_missing(pl.col(c+'_ref')) for c in strata if c not in [f'length{i}' for i in axes]])
            jx=pl.all_horizontal([pl.col(c)==pl.col(c+'_ref') for c in j])
            removed=exact.group_by('query_id').agg(pl.len().alias('_total'),sx.sum().alias('_stratum'),(sx&jx).sum().alias('_joint'),jx.sum().alias('_global'))
            out=out.join(removed,on='query_id',how='left',validate='1:1',maintain_order='left').with_columns(
                *[(pl.col(c)-pl.col('_'+('total' if c=='population' else c)).fill_null(0)).alias(c) for c in ['population','stratum','joint','global']])
        out=out.with_columns(((pl.col('joint')+(pl.col('global')+.5)/(pl.col('population')+.5*categories))/(pl.col('stratum')+1)).alias(prefix+'_probability'))
        return out.select('query_id',pl.col(prefix+'_probability'),*[pl.col(c).alias(prefix+'_'+c) for c in ['population','stratum','joint','global']])
    result=predictive(ref,'j_reference').join(predictive(bg,'j_control'),on='query_id',validate='1:1',maintain_order='left')
    missing=pl.any_horizontal([pl.col(c)=='__missing__' for c in j])
    status=q.select('query_id',missing.alias('_missing'))
    result=result.join(status,on='query_id',validate='1:1',maintain_order='left').with_columns(
        pl.when(pl.col('_missing')).then(pl.lit(1.0)).otherwise(pl.col('j_reference_probability')/pl.col('j_control_probability')).alias('conditional_j_usage_ratio'),
        pl.when(pl.col('_missing')).then(pl.lit('missing_query_j_neutral')).otherwise(pl.lit('empirical_conditional_j_ranking_only')).alias('conditional_j_status'))
    return result.drop('_missing')
