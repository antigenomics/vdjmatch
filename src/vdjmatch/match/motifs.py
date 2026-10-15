"""masked-motif V-conditioned MF×IKF evidence; not a calibrated likelihood or P-value."""
import polars as pl

from ..aggregate.candidates import PMHC
from ..evalue.epitope_rules import _tokens


def masked_motif_scores(queries, reference, *, exclude_exact=False, target_reference=None):
    """Return supported V+masked-4-mer scores per query/full pMHC.

    Inputs use full junction ``cdr3`` and ``v``; queries have unique query_id,
    reference has PMHC columns. Call separately for each species and locus.
    Follow masked_motif: IMGT classes, all four masks, distinct junction presence,
    x_length,V,target>=2, x_V,target>=3, x_target>=5, and sum(max(MF_V×IKF))
    per motif divided by query length. All model support gates and global motif
    frequencies remain fixed. Benchmark exclusion masks exact-junction target/V
    contributions only in the hit numerator, without refitting the model.
    Missing/unsupported motifs emit no row, never a prior-only specificity.
    Optional target_reference defines compatible-restriction unions labelled with
    each query's declared PMHC; global motif frequencies still use reference.
    """
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('unique non-null query IDs required')
    for frame in (queries,reference):
        if any(frame[c].null_count() for c in ('cdr3','v')):
            raise ValueError('complete junction and V required')
    ref=reference.select(*PMHC,'cdr3','v').with_columns(pl.col('v').str.replace(r'\*.*$',''))
    ref=ref.unique(maintain_order=True).rename({'cdr3':'seq0','v':'v0'}).with_row_index('_id')
    ref=ref.with_columns(pl.col('seq0').str.len_chars().alias('_length'))
    q=queries.select('query_id','cdr3','v',*(PMHC if target_reference is not None else [])).rename({'cdr3':'seq0','v':'v0'}).with_row_index('_id')
    q=q.with_columns(pl.col('v0').str.replace(r'\*.*$',''))
    def tokens(frame):
        return _tokens(frame,[0],'imgt_wild4',full_junction=True,wildcard_positions=(0,1,2,3))
    rt=tokens(ref).join(ref.drop('v0'),on='_id',validate='m:1')
    target_rt=rt
    if target_reference is not None:
        target=target_reference.select(*PMHC,'cdr3','v').with_columns(pl.col('v').str.replace(r'\*.*$',''))
        target=target.unique(maintain_order=True).rename({'cdr3':'seq0','v':'v0'}).with_row_index('_id')
        target=target.with_columns(pl.col('seq0').str.len_chars().alias('_length'))
        target_rt=tokens(target).join(target.drop('v0'),on='_id',validate='m:1')
    strata=[*PMHC,'v0','_length','t0']
    model=target_rt.group_by(strata).len().rename({'len':'x_lva'}).with_columns(
        pl.col('x_lva').sum().over([*PMHC,'v0','t0']).alias('x_va'),
        pl.col('x_lva').sum().over([*PMHC,'t0']).alias('x_a'),
        pl.lit(ref.height).alias('X')).join(rt.group_by('t0').len().rename({'len':'x'}),on='t0',validate='m:1')
    # Model support gates are fixed before matching any query.
    model=model.filter((pl.col('x_lva')>=2)&(pl.col('x_va')>=3)&(pl.col('x_a')>=5))
    hits=tokens(q).join(q.drop('v0'),on='_id',validate='m:1').join(model,
        on=['v0','t0',*(PMHC if target_reference is not None else [])],how='inner')
    hits=hits.with_columns(pl.col('x_va').alias('hit_x_va'))
    if exclude_exact:
        removed=target_rt.group_by('seq0',*PMHC,'v0','t0').len().rename({'len':'_removed'})
        hits=hits.join(removed,on=['seq0',*PMHC,'v0','t0'],how='left',validate='m:1',nulls_equal=True).with_columns(
            (pl.col('hit_x_va').cast(pl.Int64)-pl.col('_removed').fill_null(0)).alias('hit_x_va')).drop('_removed')
    hits=hits.filter(pl.col('hit_x_va')>0).with_columns(
        (pl.col('hit_x_va')/pl.col('x')*(pl.col('X')/pl.col('x')).log()).alias('_mfikf'))
    return (hits.group_by('query_id',*PMHC,'t0').agg(pl.col('_mfikf').max(),pl.col('seq0').first())
        .group_by('query_id',*PMHC).agg((pl.col('_mfikf').sum()/pl.col('seq0').first().str.len_chars()).alias('motif_score'),
            pl.len().alias('supported_motifs')))
