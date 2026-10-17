"""Linked-first paired label decisions, separate from any-call detection."""
import polars as pl

AXES = (('linked', 'kernel_linked_score', 'kernel_score', 'tied_peptides', 'selected_linked_epitope'),
        ('alpha', 'kernel_alpha_score', 'kernel_alpha_score', 'alpha_tied_peptides', 'selected_alpha_epitope'),
        ('beta', 'kernel_beta_score', 'kernel_beta_score', 'beta_tied_peptides', 'selected_beta_epitope'))


def paired_labels(evidence, selection, *, species):
    """Return flat calls and selected-label indicators from native paired evidence.

    Positive finite linked evidence retains every linked maximum peptide. Only
    absent linked evidence permits the union of alpha/beta maximum peptides.
    Full-pMHC duplicate peptides collapse within each axis. Native maxima, tie
    counts and lexical winners must match exactly; all selection metadata stays
    unchanged. Output scores are decision indicators, never continuous evidence.
    """
    if species not in ('HomoSapiens', 'MusMusculus'):
        raise ValueError('canonical human or mouse species required')
    if selection['query_id'].null_count() or selection['query_id'].n_unique()!=selection.height:
        raise ValueError('native selection requires unique non-null query IDs')
    if evidence['query_id'].null_count() or evidence.select('query_id').unique().join(
            selection.select('query_id'),on='query_id',how='anti').height:
        raise ValueError('paired evidence contains an unknown query ID')
    if any(c in selection.columns for c in ('paired_epitope','paired_label_source','paired_label_status','paired_label_count','_label_order')):
        raise ValueError('native selection contains reserved label output columns')
    winners=[]
    for axis,column,maximum,count,label in AXES:
        native=selection.select('query_id',pl.col(maximum).cast(pl.Float64).alias('_expected'),
            pl.col(count).cast(pl.UInt64).alias('_ties'),pl.col(label).alias('_label'))
        positive=pl.col(column).cast(pl.Float64).is_finite() & (pl.col(column).cast(pl.Float64)>0)
        values=evidence.filter(positive).select('query_id','epitope',pl.col(column).cast(pl.Float64).alias('_value'))
        if values.filter(pl.col('epitope').is_null() | (pl.col('epitope').str.strip_chars()=='')).height:
            raise ValueError('positive paired evidence requires an epitope label')
        maxima=values.group_by('query_id').agg(pl.col('_value').max())
        won=values.join(maxima,on=['query_id','_value'],how='inner').select('query_id','epitope').unique()
        observed=won.group_by('query_id').agg(pl.len().alias('_observed_ties'),pl.col('epitope').min().alias('_observed_label'))
        check=native.join(maxima,on='query_id',how='left').join(observed,on='query_id',how='left')
        supported=pl.col('_expected').is_finite() & (pl.col('_expected')>0)
        drift=(supported & (~pl.col('_value').eq_missing(pl.col('_expected')) |
            ~pl.col('_ties').eq_missing(pl.col('_observed_ties')) |
            ~pl.col('_label').eq_missing(pl.col('_observed_label')))) | (~supported.fill_null(False)&pl.col('_value').is_not_null())
        if check.filter(drift.fill_null(False)).height:
            raise ValueError(f'{axis} evidence disagrees with native maxima, ties or lexical winner')
        winners.append(won.with_columns(pl.lit(axis).alias('_axis')))
    linked=winners[0].drop('_axis').with_columns(pl.lit('linked').alias('paired_label_source'))
    marginal=pl.concat(winners[1:]).select('query_id','epitope').unique().join(
        linked.select('query_id').unique(),on='query_id',how='anti').with_columns(pl.lit('marginal_fallback').alias('paired_label_source'))
    chosen=pl.concat([linked,marginal]).rename({'epitope':'paired_epitope'})
    calls=selection.with_row_index('_label_order').join(chosen,on='query_id',how='left').with_columns(
        pl.col('paired_label_source').fill_null('none'),
        pl.col('paired_epitope').count().over('query_id').cast(pl.UInt32).alias('paired_label_count'))
    native_ambiguity=pl.lit(False)
    for axis,_,maximum,_,_ in AXES:
        status='linked_selection_status' if axis=='linked' else axis+'_selection_status'
        if status not in calls.columns: continue
        source='linked' if axis=='linked' else 'marginal_fallback'
        value=pl.col(maximum).cast(pl.Float64)
        native_ambiguity=native_ambiguity | ((pl.col('paired_label_source')==source) &
            value.is_finite() & (value>0) & (pl.col(status)=='ambiguous')).fill_null(False)
    calls=calls.with_columns(pl.when(pl.col('paired_label_count')==0).then(pl.lit('no_support'))
        .when((pl.col('paired_label_count')>1) | native_ambiguity).then(pl.lit('ambiguous'))
        .otherwise(pl.lit('ranked_evidence')).alias('paired_label_status'))
    calls=calls.sort('_label_order','paired_epitope').drop('_label_order')
    scores=calls.select(pl.lit(species).alias('species'),pl.lit('PAIRED').alias('locus'),'query_id',
        pl.col('paired_epitope').alias('epitope'),pl.lit('stem4_linked_first').alias('scorer'),
        pl.when(pl.col('paired_label_count')>0).then(pl.lit(1.,dtype=pl.Float64)).otherwise(None).alias('score'),
        (pl.col('paired_label_count')>0).alias('supported'))
    return calls,scores
