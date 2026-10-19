"""Synthetic V imputation must not rewrite the accepted resolved-V model."""
import polars as pl
import pytest
from polars.testing import assert_frame_equal
from vdjmatch.evalue.reference_rule_bank import stream_reference_kernels


@pytest.mark.parametrize('paired', [False, True])
def test_resolved_first_preserves_baseline_and_extended_query_scores(paired):
    queries = pl.DataFrame(dict(query_id=[0, 1], seq0=['CASSLGQAYEQYF']*2,
        v0=['TRBV19*01']*2, baseline_v0=['TRBV19*01', None], key0=['known', 'unknown']))
    reference = pl.DataFrame(dict(seq0=['CASSLGRAYEQYF', 'CASSLGQAYEQWF', 'CASSLGQAYEQYF'],
        v0=['TRBV19*01']*3, baseline_v0=['TRBV19*01', 'TRBV19*01', None],
        key0=['known', 'known', 'unknown'], epitope=['A', 'B', 'C']))
    if paired:
        queries = queries.rename({'seq0': 'seq1', 'v0': 'v1', 'baseline_v0': 'baseline_v1'})
        reference = reference.rename({'seq0': 'seq1', 'v0': 'v1', 'baseline_v0': 'baseline_v1'})
        for name, value in [('seq0', 'CAVRDSNYQLIW'), ('v0', 'TRAV12-2*01'), ('baseline_v0', 'TRAV12-2*01')]:
            queries = queries.with_columns(pl.lit(value).alias(name))
            reference = reference.with_columns(pl.lit(value).alias(name))
    def run(q, r, **options):
        frames = []
        stream_reference_kernels(q, r, paired=paired, identity_keys=['key0'],
            consume=lambda scorer, frame: frames.append(frame.with_columns(pl.lit(scorer).alias('scorer'))),
            **options)
        return pl.concat(frames).sort('query_id', 'epitope', 'scorer')
    extended = run(queries, reference)
    baseline = run(queries.head(1), reference.head(2))
    repaired = run(queries, reference, v_model_policy='resolved-first')
    assert_frame_equal(repaired.filter((pl.col('query_id') == 0) & pl.col('epitope').is_in(['A', 'B'])), baseline)
    assert_frame_equal(repaired.filter(pl.col('query_id') == 1), extended.filter(pl.col('query_id') == 1))
    absent = repaired.filter((pl.col('query_id') == 0) & (pl.col('epitope') == 'C'))
    assert absent.height == 6 and absent['score'].null_count() == 6 and not absent['supported'].any()
    assert repaired.height == 2*3*6
    # The fixed model is independent of query order, batch membership and repeats.
    for q in (queries.reverse(), queries.head(1), queries.tail(1)):
        assert_frame_equal(run(q, reference, v_model_policy='resolved-first'),
            repaired.filter(pl.col('query_id').is_in(q['query_id'].implode())))


def test_resolved_first_empty_trusted_model_is_explicit():
    q = pl.DataFrame(dict(query_id=[0], seq0=['CASSLGQAYEQYF'], v0=['TRBV19*01'],
        baseline_v0=['TRBV19*01']))
    r = pl.DataFrame(dict(seq0=['CASSLGRAYEQYF'], v0=['TRBV19*01'],
        baseline_v0=pl.Series([None], dtype=pl.String), epitope=['A']))
    frames = []
    info = stream_reference_kernels(q, r, v_model_policy='resolved-first',
        consume=lambda scorer, frame: frames.append(frame))
    assert info['native_batches'] == 0 and info['query_target_rows'] == 1
    scores = pl.concat(frames)
    assert scores.height == scores['score'].null_count() == 6
    assert not scores['supported'].any()
    with pytest.raises(ValueError, match='original resolved allele columns'):
        stream_reference_kernels(q.drop('baseline_v0'), r, v_model_policy='resolved-first', consume=lambda *args: None)
