"""One-vs-other annotation retrieval from complete scorer candidate outputs."""
import polars as pl

CONTEXT = ['species', 'locus']
QUERY = [*CONTEXT, 'query_id']
TARGET = [*CONTEXT, 'epitope']


def reference_confusion(scores, truth, profiles):
    """Return flat confusion counts and tie-preserving winning-label calls.

    Truth retains every annotated epitope per query. Missing labels are
    other-class labels for retrieval, not experimentally demonstrated nonbinding.
    Every scorer must explicitly represent every truth query, including no-hits.
    Supported finite maxima win; all tied epitopes win. Models and exact-hit
    filtering belong to the score producer and are never altered here.
    """
    if (not {*QUERY, 'epitope', 'scorer', 'score', 'supported'} <= set(scores.columns)
            or not {*QUERY, 'epitope'} <= set(truth.columns)
            or not {*TARGET, 'valid_source_records', 'level'} <= set(profiles.columns)):
        raise ValueError('confusion requires complete scorer/query/target/support, truth and profile columns')
    if scores['supported'].dtype != pl.Boolean or scores['supported'].null_count():
        raise ValueError('supported must be an explicit non-null Boolean')
    if not truth.height or not scores.height:
        raise ValueError('confusion requires nonempty truth and explicit scorer query dispositions')
    scores = scores.with_columns(pl.col('score').cast(pl.Float64))
    if truth.select(pl.any_horizontal(pl.col([*QUERY, 'epitope']).is_null()).any()).item():
        raise ValueError('truth identities and annotated epitopes must be non-null')
    if scores.select(pl.any_horizontal(pl.col([*QUERY, 'scorer']).is_null()).any()).item():
        raise ValueError('score query and scorer identities must be non-null')
    truth = truth.select(*QUERY, 'epitope').unique()
    queries = truth.select(QUERY).unique()
    methods = queries.select(CONTEXT).unique().join(scores.select('scorer').unique(), how='cross')
    expected = queries.join(methods, on=CONTEXT, validate='m:m')
    actual = scores.select(*QUERY, 'scorer').unique()
    if (expected.height != actual.height or
            expected.join(actual, on=[*QUERY, 'scorer'], how='anti').height or
            queries.join(methods.select(CONTEXT).unique(), on=CONTEXT, how='anti').height):
        raise ValueError('every scorer must explicitly contain every truth query, including no-hits')
    eligible = profiles.filter((pl.col('level') == 'epitope') & (pl.col('valid_source_records') >= 30)).select(
        *TARGET, 'valid_source_records', *[c for c in ('profile_class', 'full20_profile_class') if c in profiles.columns])
    if eligible.select(TARGET).is_duplicated().any():
        raise ValueError('eligible pooled profiles must be unique by species/locus/epitope')
    candidates = scores.filter(pl.col('supported') & pl.col('score').is_finite().fill_null(False) & pl.col('epitope').is_not_null())
    winners = candidates.filter(pl.col('score') == pl.col('score').max().over(*QUERY, 'scorer')).select(
        *QUERY, 'scorer', 'epitope').unique()
    calls = expected.join(winners, on=[*QUERY, 'scorer'], how='left', validate='1:m').with_columns(
        pl.col('epitope').count().over(*QUERY, 'scorer').alias('winning_labels')).with_columns(
            pl.when(pl.col('winning_labels') == 0).then(pl.lit('no_supported_hit'))
            .when(pl.col('winning_labels') > 1).then(pl.lit('tied_winners')).otherwise(pl.lit('winner')).alias('status'))
    population = queries.group_by(CONTEXT).len(name='queries')
    positives = truth.group_by(TARGET).len(name='positives')
    predicted = winners.group_by(*TARGET, 'scorer').len(name='predicted_positives')
    tp = winners.join(truth, on=[*QUERY, 'epitope'], how='semi').group_by(*TARGET, 'scorer').len(name='tp')
    disposition = calls.unique([*QUERY, 'scorer']).group_by(*CONTEXT, 'scorer').agg(
        (pl.col('winning_labels') == 0).sum().alias('no_hit_queries'),
        (pl.col('winning_labels') > 1).sum().alias('tied_queries'))
    metrics = eligible.join(methods, on=CONTEXT, validate='m:m').join(population, on=CONTEXT, validate='m:1').join(
        positives, on=TARGET, how='left', validate='m:1').join(predicted, on=[*TARGET, 'scorer'], how='left', validate='1:1').join(
        tp, on=[*TARGET, 'scorer'], how='left', validate='1:1').join(disposition, on=[*CONTEXT, 'scorer'], validate='m:1').with_columns(
            pl.col('positives', 'predicted_positives', 'tp').fill_null(0).cast(pl.Int64)).with_columns(
                (pl.col('queries')-pl.col('positives')).alias('negatives'),
                (pl.col('positives')-pl.col('tp')).alias('fn'),
                (pl.col('predicted_positives')-pl.col('tp')).alias('fp')).with_columns(
                    (pl.col('negatives')-pl.col('fp')).alias('tn'))
    tp, tn, fp, fn = [pl.col(c).cast(pl.Float64) for c in ('tp', 'tn', 'fp', 'fn')]
    denominator = ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)).sqrt()
    if {'profile_class', 'classification_eligible', 'motif_supported'} <= set(profiles.columns):
        from .fast_profiles import joint_reference_profiles
        joint = joint_reference_profiles(profiles)
        metrics = metrics.join(joint, on=['species', 'epitope'], how='left', validate='m:1')
    return metrics.with_columns(
        ((tp+tn)/pl.col('queries')).alias('accuracy'),
        pl.when((tp+fn > 0) & (tn+fp > 0)).then((tp/(tp+fn)+tn/(tn+fp))/2).alias('balanced_accuracy'),
        pl.when(denominator > 0).then((tp*tn-fp*fn)/denominator).alias('mcc')
    ).sort(*TARGET, 'scorer'), calls.sort(*QUERY, 'scorer', 'epitope', nulls_last=True)
