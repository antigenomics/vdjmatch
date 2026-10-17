"""One-vs-other counts retain no-hits, contradictory truth and score ties."""
import polars as pl
import pytest

from vdjmatch.aggregate.reference_confusion import reference_confusion


def test_retrieval_confusion_counts_and_complete_query_validation(tmp_path):
    truth = pl.DataFrame({'query_id': ['q1', 'q2', 'q3', 'q3', 'q4'], 'epitope': ['A', 'B', 'A', 'B', 'B']})
    scores = pl.DataFrame({'query_id': ['q1', 'q1', 'q2', 'q3', 'q3', 'q4'],
        'epitope': ['A', 'B', 'A', 'A', 'B', None], 'scorer': ['test']*6,
        'score': [1., 0., 1., 1., 1., None], 'supported': [True]*5+[False]})
    profiles = pl.DataFrame({'epitope': ['A', 'B', 'tail'], 'valid_source_records': [30, 31, 29], 'level': ['epitope']*3})
    def context(frame):
        return frame.with_columns(pl.lit('HomoSapiens').alias('species'), pl.lit('TRB').alias('locus'))
    scores, truth, profiles = map(context, (scores, truth, profiles))
    metrics, calls = reference_confusion(scores, truth, profiles)
    assert metrics['epitope'].to_list() == ['A', 'B']
    a, b = metrics.to_dicts()
    assert (a['tp'], a['fn'], a['fp'], a['tn']) == (2, 0, 1, 1)
    assert (b['tp'], b['fn'], b['fp'], b['tn']) == (1, 2, 0, 1)
    assert a['accuracy'] == a['balanced_accuracy'] == .75
    assert a['mcc'] == pytest.approx(1/3**.5)
    assert b['balanced_accuracy'] == pytest.approx(2/3)
    assert a['queries'] == b['queries'] == 4
    assert a['no_hit_queries'] == a['tied_queries'] == 1
    assert calls.filter(pl.col('query_id') == 'q3')['winning_labels'].to_list() == [2, 2]
    assert calls.filter(pl.col('query_id') == 'q4')['status'].to_list() == ['no_supported_hit']
    with pytest.raises(ValueError, match='every scorer'):
        reference_confusion(scores.filter(pl.col('query_id') != 'q4'), truth, profiles)
    with pytest.raises(ValueError, match='every scorer'):
        reference_confusion(pl.concat([scores, scores.with_columns(pl.lit('TRA').alias('locus'),
            pl.lit('second').alias('scorer'))]),
            pl.concat([truth, truth.with_columns(pl.lit('TRA').alias('locus'))]), profiles)
    absent, empty_calls = reference_confusion(scores.with_columns(pl.lit(None).alias('score'),
        pl.lit(False).alias('supported')), truth, profiles)
    assert absent['mcc'].null_count() == 2
    assert absent['balanced_accuracy'].to_list() == [.5, .5]
    assert absent['no_hit_queries'].to_list() == [4, 4]
    assert empty_calls.height == 4 and set(empty_calls['status']) == {'no_supported_hit'}
    from vdjmatch.cli.__main__ import main
    paths = [tmp_path/name for name in ('scores.tsv', 'truth.tsv', 'profiles.tsv')]
    for path, frame in zip(paths, (scores, truth, profiles)):
        frame.write_csv(path, separator='\t')
    assert main(['score-confusion', str(paths[0]), '--truth', str(paths[1]), '--profiles', str(paths[2]),
        '--output-prefix', str(tmp_path/'summary')]) == 0
    assert pl.read_csv(tmp_path/'summary.confusion.tsv', separator='\t')['tp'].to_list() == [2, 1]
    paired, _ = reference_confusion(scores.with_columns(pl.lit('paired').alias('locus')),
        truth.with_columns(pl.lit('paired').alias('locus')), profiles)
    assert paired['epitope'].to_list() == ['A', 'B']
    assert paired['valid_source_records'].to_list() == [30, 31]
    small = pl.concat([profiles, profiles.with_columns(pl.lit('TRA').alias('locus'))]).with_columns(
        pl.when(pl.col('epitope') == 'A').then(15).otherwise(pl.col('valid_source_records')).alias('valid_source_records'))
    paired, _ = reference_confusion(scores.with_columns(pl.lit('paired').alias('locus')),
        truth.with_columns(pl.lit('paired').alias('locus')), small)
    assert paired['epitope'].to_list() == ['B']  # 15+15 source records cannot manufacture eligibility
    assert paired['queries'].to_list() == [4]


def test_class_groups_keep_equal_epitope_weight_and_undefined_counts():
    from vdjmatch.aggregate.reference_confusion import class_group_confusion
    metrics = pl.DataFrame({'species': ['HomoSapiens']*2, 'locus': ['TRA']*2,
        'epitope': ['A', 'B'], 'scorer': ['test']*2, 'alpha_profile_class': ['v_only']*2,
        'alpha_full20_profile_class': ['v_and_motif']*2, 'beta_profile_class': ['motif_only']*2,
        'queries': [10, 100], 'positives': [5, 50], 'negatives': [5, 50],
        'tp': [4, 20], 'tn': [4, 20], 'fp': [1, 30], 'fn': [1, 30],
        'accuracy': [.8, .4], 'balanced_accuracy': [.8, .4], 'mcc': [.6, None]})
    groups = class_group_confusion(metrics)
    assert groups.height == 2  # beta metadata must not create a beta-scoring cohort
    assert groups['class_definition'].to_list() == ['core3', 'full20']
    assert groups['class_axis'].to_list() == ['alpha']*2
    assert groups['mean_accuracy'].to_list() == pytest.approx([.6]*2)
    assert groups['mean_mcc'].to_list() == [.6]*2
    assert groups['mcc_defined_epitopes'].to_list() == [1, 1]
    assert groups['epitopes'].to_list() == [2, 2]
    assert groups['query_target_appearances'].to_list() == [110, 110]
    with pytest.raises(ValueError, match='one row'):
        class_group_confusion(pl.concat([metrics, metrics]))
    other = metrics.with_columns(pl.lit('other').alias('scorer'),
        (1-pl.col('accuracy')).alias('accuracy'))
    boot = class_group_confusion(pl.concat([metrics, other]), nboot=1000, seed=0)
    original = boot.filter(pl.col('scorer') == 'test').row(0, named=True)
    opposite = boot.filter(pl.col('scorer') == 'other').row(0, named=True)
    assert original['bootstrap_se_accuracy'] == pytest.approx(.2/2**.5, abs=.02)
    assert original['bootstrap_lo_accuracy'] == .4 and original['bootstrap_hi_accuracy'] == .8
    assert original['bootstrap_lo_accuracy'] == pytest.approx(1-opposite['bootstrap_hi_accuracy'])
    assert original['bootstrap_se_accuracy'] == pytest.approx(opposite['bootstrap_se_accuracy'])
    assert original['bootstrap_valid_accuracy'] == 1000
    assert original['bootstrap_status_mcc'] == 'insufficient_epitopes'
    assert original['bootstrap_lo_mcc'] is None
    assert original['bootstrap_valid_mcc'] < 1000
    mixed = pl.concat([metrics, metrics.head(1).with_columns(pl.lit('TRB').alias('locus'),
        pl.lit(float('inf')).alias('mcc'), pl.lit(float('nan')).alias('balanced_accuracy'))])
    combined = class_group_confusion(mixed, nboot=1000, seed=0)
    singleton = combined.filter(pl.col('locus') == 'TRB').row(0, named=True)
    assert singleton['bootstrap_lo_accuracy'] is None
    assert singleton['mcc_defined_epitopes'] == singleton['balanced_accuracy_defined_epitopes'] == 0
    assert singleton['bootstrap_valid_mcc'] == singleton['bootstrap_valid_balanced_accuracy'] == 0


def test_joint_classes_describe_single_chain_scores_without_fusing_populations():
    from vdjmatch.aggregate.reference_confusion import class_group_confusion
    from polars.testing import assert_frame_equal

    base = pl.DataFrame({'species': ['HomoSapiens']*3, 'locus': ['TRA']*3,
        'epitope': ['A', 'B', 'C'], 'scorer': ['first']*3,
        'alpha_profile_class': ['v_only']*3, 'beta_profile_class': ['motif_only']*3,
        'alpha_full20_profile_class': ['v_and_motif']*3,
        'beta_full20_profile_class': ['v_and_motif']*3,
        'joint_class': ['alpha: v_only; beta: motif_only']*3,
        'full20_joint_class': ['alpha: v_and_motif; beta: v_and_motif']*3,
        'both_chains_observed': [True, False, None],
        'alpha_valid_source_records': [30]*3, 'beta_valid_source_records': [1, None, None],
        'both_chains_eligible': [False]*3,
        'queries': [10, 10, 10], 'positives': [5]*3, 'negatives': [5]*3,
        'tp': [4, 2, 1], 'tn': [4, 2, 1], 'fp': [1, 3, 4], 'fn': [1, 3, 4],
        'accuracy': [.8, .4, .2], 'balanced_accuracy': [.8, .4, .2], 'mcc': [.6, -.2, -.6]})
    # Locus, scorer and species remain independent metric populations.
    populations = []
    for species in ('HomoSapiens', 'MusMusculus'):
        for locus in ('TRA', 'TRB', 'paired'):
            for scorer in ('first', 'second'):
                populations.append(base.with_columns(pl.lit(species).alias('species'),
                    pl.lit(locus).alias('locus'), pl.lit(scorer).alias('scorer')))
    metrics = pl.concat(populations)
    grouped = class_group_confusion(metrics)
    joint = grouped.filter(pl.col('class_axis') == 'joint')
    assert joint.height == 3*2*2*2
    singles = joint.filter(pl.col('locus').is_in(['TRA', 'TRB']))
    assert singles['epitopes'].to_list() == [1]*16
    assert singles['query_target_appearances'].to_list() == [10]*16
    assert singles['mean_accuracy'].to_list() == pytest.approx([.8]*16)
    assert set(singles['class_definition']) == {'core3', 'full20'}
    assert set(singles['epitope_class']) == {base['joint_class'][0], base['full20_joint_class'][0]}
    assert singles.select('species', 'locus', 'scorer', 'class_definition').unique().height == 16
    assert joint.filter(pl.col('locus') == 'paired')['epitopes'].to_list() == [3]*8
    # Mixed joint groups combine real paired intervals with singleton null intervals.
    singleton_bootstrap = class_group_confusion(metrics, nboot=1000, seed=3)
    singleton_joint = singleton_bootstrap.filter((pl.col('class_axis') == 'joint') &
        pl.col('locus').is_in(['TRA', 'TRB']))
    paired_joint = singleton_bootstrap.filter((pl.col('class_axis') == 'joint') &
        (pl.col('locus') == 'paired'))
    assert singleton_joint['bootstrap_lo_accuracy'].null_count() == singleton_joint.height
    assert set(singleton_joint['bootstrap_status_accuracy']) == {'insufficient_epitopes'}
    assert paired_joint['bootstrap_lo_accuracy'].null_count() == 0
    assert singleton_bootstrap.schema['bootstrap_lo_accuracy'] == pl.Float64
    legacy_singleton = class_group_confusion(metrics.drop('both_chains_observed'), nboot=1000, seed=3)
    assert_frame_equal(singleton_bootstrap.filter(~((pl.col('class_axis') == 'joint') &
        pl.col('locus').is_in(['TRA', 'TRB']))), legacy_singleton, check_exact=True)
    reverse_metrics = metrics.filter((pl.col('locus') != 'paired') | (pl.col('epitope') == 'A')).with_columns(
        pl.lit(True).alias('both_chains_observed'))
    reverse = class_group_confusion(reverse_metrics, nboot=1000, seed=3)
    reverse_legacy = class_group_confusion(reverse_metrics.drop('both_chains_observed'), nboot=1000, seed=3)
    assert_frame_equal(reverse.filter(~((pl.col('class_axis') == 'joint') &
        pl.col('locus').is_in(['TRA', 'TRB']))), reverse_legacy, check_exact=True)
    # Missing joint observation metadata preserves the earlier paired behavior.
    bootstrap_metrics = pl.concat([metrics, metrics.filter(pl.col('epitope') == 'B').with_columns(
        pl.lit('D').alias('epitope'), pl.lit(True).alias('both_chains_observed'))])
    legacy_bootstrap = class_group_confusion(bootstrap_metrics.drop('both_chains_observed'), nboot=1000, seed=3)
    extended_bootstrap = class_group_confusion(bootstrap_metrics, nboot=1000, seed=3)
    new_rows = (pl.col('class_axis') == 'joint') & pl.col('locus').is_in(['TRA', 'TRB'])
    assert_frame_equal(extended_bootstrap.filter(~new_rows), legacy_bootstrap, check_exact=True)
    previous = class_group_confusion(metrics.drop('both_chains_observed'))
    assert_frame_equal(grouped.filter(pl.col('class_axis') != 'joint'),
        previous.filter(pl.col('class_axis') != 'joint'))
    assert_frame_equal(joint.filter(pl.col('locus') == 'paired'),
        previous.filter(pl.col('class_axis') == 'joint'))
    assert previous.filter((pl.col('class_axis') == 'joint') & pl.col('locus').is_in(['TRA', 'TRB'])).is_empty()
