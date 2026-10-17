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
