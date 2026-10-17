"""One-vs-other annotation retrieval from complete scorer candidate outputs."""
import polars as pl

CONTEXT = ['species', 'locus']
QUERY = [*CONTEXT, 'query_id']
TARGET = [*CONTEXT, 'epitope']


def class_group_confusion(metrics, *, nboot=0, seed=0):
    """Equal-epitope macro summaries for alpha, beta and joint class heatmaps.

    Align marginal class axes to scored locus. Joint classes also describe
    single-chain metrics when both marginal chains are observed; the scored
    locus, species and scorer stay separate. No paired score is inferred.
    Report defined metric denominators and separate core3/full20 definitions.
    Individual epitope metrics remain the primary input.
    """
    import numpy as np
    if type(nboot) is not int or nboot < 0 or nboot == 1 or type(seed) is not int or seed < 0:
        raise ValueError('nboot must be zero or at least two; seed must be nonnegative')
    rng = np.random.default_rng(seed)
    joint_rng = np.random.default_rng(np.random.SeedSequence(seed).spawn(1)[0])
    schema = {c: pl.String for c in [*CONTEXT, 'scorer', 'class_axis', 'class_definition', 'epitope_class']}
    schema.update({c: pl.Int64 for c in ['epitopes', 'query_target_appearances', 'positives', 'negatives', 'tp', 'tn', 'fp', 'fn']})
    for metric in ('accuracy', 'balanced_accuracy', 'mcc'):
        schema[metric+'_defined_epitopes'] = pl.Int64
        schema.update({stat+'_'+metric: pl.Float64 for stat in ('mean', 'median', 'min', 'max')})
        if nboot:
            schema['bootstrap_valid_'+metric] = pl.Int64
            schema['bootstrap_status_'+metric] = pl.String
            schema.update({'bootstrap_'+stat+'_'+metric: pl.Float64 for stat in ('se', 'bias', 'lo', 'hi')})
    if metrics.select(*TARGET, 'scorer').is_duplicated().any():
        raise ValueError('class-group input must contain one row per pooled epitope/scorer/locus')
    parts = []
    for definition, prefix in [('core3', ''), ('full20', 'full20_')]:
        for axis, loci, field in [('alpha', ['TRA'], 'alpha_'+prefix+'profile_class'),
                ('beta', ['TRB'], 'beta_'+prefix+'profile_class'),
                ('joint', ['TRA', 'TRB', 'PAIRED', 'TRA+TRB'], prefix+'joint_class')]:
            if field not in metrics.columns:
                continue
            part = metrics.filter(pl.col('locus').str.to_uppercase().is_in(loci))
            if axis == 'joint':
                observed = (pl.col('both_chains_observed').fill_null(False)
                    if 'both_chains_observed' in part.columns else pl.lit(False))
                part = part.filter(~pl.col('locus').str.to_uppercase().is_in(['TRA', 'TRB']) | observed)
            if not part.height:
                continue
            part = part.with_columns(pl.lit(axis).alias('class_axis'),
                pl.lit(definition).alias('class_definition'), pl.col(field).fill_null('missing_class').alias('epitope_class'),
                *[pl.when(pl.col(c).is_finite()).then(pl.col(c)).alias(c) for c in ('accuracy', 'balanced_accuracy', 'mcc')])
            aggregates = [pl.col('epitope').n_unique().alias('epitopes'),
                pl.col('queries').sum().alias('query_target_appearances'),
                *[pl.col(c).sum() for c in ('positives', 'negatives', 'tp', 'tn', 'fp', 'fn')]]
            for metric in ('accuracy', 'balanced_accuracy', 'mcc'):
                aggregates += [pl.col(metric).count().alias(metric+'_defined_epitopes'),
                    *[getattr(pl.col(metric), stat)().alias(stat+'_'+metric) for stat in ('mean', 'median', 'min', 'max')]]
            groups = [*CONTEXT, 'class_axis', 'class_definition', 'epitope_class']
            summary = part.group_by(*groups, 'scorer').agg(aggregates)
            if nboot:
                if axis == 'joint':
                    single = pl.col('locus').str.to_uppercase().is_in(['TRA', 'TRB'])
                    # New descriptive groups must not consume legacy bootstrap draws.
                    boot = pl.concat([_bootstrap_class_means(subset, groups, nboot=nboot, rng=stream)
                        for subset, stream in [(part.filter(~single), rng), (part.filter(single), joint_rng)]
                        if subset.height], how='vertical_relaxed')
                else:
                    boot = _bootstrap_class_means(part, groups, nboot=nboot, rng=rng)
                summary = summary.join(boot, on=[*groups, 'scorer'], validate='1:1')
            parts.append(summary.select(*[pl.col(c).cast(dtype) for c, dtype in schema.items()]))
    return (pl.concat(parts)
        .sort(*CONTEXT, 'class_definition', 'class_axis', 'epitope_class', 'scorer') if parts else pl.DataFrame(schema=schema))


def _bootstrap_class_means(part, groups, *, nboot, rng):
    """Efron empirical epitope bootstrap, paired across scorers and metrics.

    Multinomial multiplicities equal sampling N epitope units with replacement.
    Fixed references and class memberships define the conditional estimand.
    Do not resample dependent graph edges, tokens or duplicated source records.
    """
    import numpy as np
    records, metrics = [], ('accuracy', 'balanced_accuracy', 'mcc')
    for key, frame in part.sort([*groups, 'epitope', 'scorer']).group_by(groups, maintain_order=True):
        epitopes = frame.select('epitope').unique().sort('epitope').with_row_index('_e')
        scorers = frame.select('scorer').unique().sort('scorer').with_row_index('_s')
        indexed = frame.join(epitopes, on='epitope', validate='m:1').join(scorers, on='scorer', validate='m:1')
        n, methods = epitopes.height, scorers.height
        values = np.full((n, methods, 3), np.nan)
        values[indexed['_e'].to_numpy(), indexed['_s'].to_numpy(), :] = indexed.select(metrics).to_numpy()
        values = values.reshape(n, -1)
        finite = np.isfinite(values)
        draws = np.full((nboot, methods*3), np.nan)
        for start in range(0, nboot, 32):
            weights = rng.multinomial(n, np.full(n, 1/n), size=min(32, nboot-start))
            denominator = weights @ finite.astype(np.float64)
            draws[start:start+len(weights)] = np.divide(weights @ np.nan_to_num(values), denominator,
                out=np.full(denominator.shape, np.nan), where=denominator > 0)
        for method, scorer in enumerate(scorers['scorer']):
            record = dict(zip(groups, key)) | {'scorer': scorer}
            for k, metric in enumerate(metrics):
                column = method*3+k
                sample = draws[:, column]
                valid = int(np.isfinite(sample).sum())
                sufficient = int(finite[:, column].sum()) >= 2
                usable = sufficient and valid == nboot
                record.update({'bootstrap_valid_'+metric: valid,
                    'bootstrap_status_'+metric: 'empirical_epitope_percentile' if usable else
                        'insufficient_epitopes' if not sufficient else 'undefined_resamples',
                    'bootstrap_se_'+metric: float(sample.std(ddof=1)) if usable else None,
                    'bootstrap_bias_'+metric: float(sample.mean()-np.nanmean(values[:, column])) if usable else None,
                    'bootstrap_lo_'+metric: float(np.quantile(sample, .025)) if usable else None,
                    'bootstrap_hi_'+metric: float(np.quantile(sample, .975)) if usable else None})
            records.append(record)
    return pl.DataFrame(records, infer_schema_length=None)


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
    pooled = profiles.filter((pl.col('level') == 'epitope') & pl.col('locus').is_in(['TRA', 'TRB']))
    if pooled.select(TARGET).is_duplicated().any():
        raise ValueError('eligible pooled profiles must be unique by species/locus/epitope')
    eligible = pooled.filter(pl.col('valid_source_records') >= 30).select(
        *TARGET, 'valid_source_records', *[c for c in ('profile_class', 'full20_profile_class') if c in profiles.columns])
    paired_loci = scores.filter(pl.col('locus').str.to_uppercase().is_in(['PAIRED', 'TRA+TRB'])).select('locus').unique()
    if paired_loci.height:
        paired = pooled.group_by('species', 'epitope').agg(
            pl.col('valid_source_records').sum(), pl.col('valid_source_records').max().alias('eligibility_chain_max_records'))
        paired = paired.filter(pl.col('eligibility_chain_max_records') >= 30).join(paired_loci, how='cross')
        eligible = pl.concat([eligible, paired], how='diagonal_relaxed')
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
