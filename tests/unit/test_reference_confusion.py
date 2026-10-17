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
