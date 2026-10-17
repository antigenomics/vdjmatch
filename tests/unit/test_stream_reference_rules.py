"""All-label sparse projections agree with the fixed scalar rule oracle."""
import polars as pl
import pytest


@pytest.mark.parametrize('paired', [False, True])
def test_streamed_rules_match_every_scalar_model(paired):
    from vdjmatch.evalue.reference_rule_bank import stream_reference_rules
    from vdjmatch.evalue.epitope_rules import reference_epitope_rules
    reference = pl.DataFrame({'seq0': ['CASSLGQAYEQYF', 'CASSLGRAYEQYF', 'CASSLGRAYEQYF',
        'CASSWGRAYEQYF', 'CASSLGQAYEQYF', 'CASSAAAAYEQYF'],
        'v0': ['TRBV19*01', 'TRBV19*01', 'TRBV20-1*01', 'TRBV19*01', 'TRBV19*01', 'TRBV20-1*01'],
        'epitope': ['A', 'A', 'B', 'B', 'B', 'C']})
    if paired:
        reference = reference.with_columns(pl.col('seq0').alias('seq1'), pl.col('v0').alias('v1'),
            pl.lit('CAVRLGQAGKSTF').alias('seq0'), pl.lit('TRAV1-1*01').alias('v0'))
    keys = ['seq0', 'v0', *(['seq1', 'v1'] if paired else [])]
    queries = reference.select(keys).unique(maintain_order=True).with_row_index('query_id')
    extra = queries.head(1).with_columns(pl.lit(queries.height, dtype=pl.UInt32).alias('query_id'),
        pl.lit('TRBV6-5*01' if not paired else 'TRAV2*01').alias('v0'))
    queries = pl.concat([queries, extra])
    collected = {}
    def consume(scorer, frame):
        collected.setdefault(scorer, []).append(frame)
    stream_reference_rules(queries, reference, paired=paired, species='human', threads=1, consume=consume)
    expected = []
    for epitope in sorted(reference['epitope'].unique()):
        expected.append(reference_epitope_rules(queries,
            reference.filter(pl.col('epitope') == epitope).select(keys),
            reference.filter(pl.col('epitope') != epitope).select(keys),
            paired=paired, exclude_exact=True).with_columns(pl.lit(epitope).alias('epitope')))
    expected = pl.concat(expected).sort('query_id', 'epitope')
    score_columns = [c for c in expected.columns if c.endswith('_score')]
    assert len(collected) == len(score_columns) == 15
    for column in score_columns:
        actual = pl.concat(collected[column.removesuffix('_score')]).sort('query_id', 'epitope')
        assert actual.select('query_id', 'epitope').equals(expected.select('query_id', 'epitope'))
        assert actual['score'].to_list() == pytest.approx(expected[column].to_list(), abs=1e-14)
        assert actual['supported'].all()
    assert stream_reference_rules(queries.head(0), reference, paired=paired, consume=consume)['native_batches'] == 0
    empty_background = reference.with_columns(pl.lit('A').alias('epitope'))
    collected.clear()
    stream_reference_rules(queries, empty_background, paired=paired, consume=consume)
    assert all(not pl.concat(frames)['supported'].any() and pl.concat(frames)['score'].null_count() == queries.height
        for frames in collected.values())


def test_streamed_rules_with_no_core_motif_tokens():
    from vdjmatch.evalue.reference_rule_bank import stream_reference_rules
    from vdjmatch.evalue.epitope_rules import reference_epitope_rules
    reference = pl.DataFrame({'seq0': ['CASAAQYF', 'CASACQYF'], 'v0': ['TRBV19*01']*2, 'epitope': ['A', 'B']})
    queries = reference.select('seq0', 'v0').with_row_index('query_id')
    frames = []
    stream_reference_rules(queries, reference, consume=lambda name, frame: frames.append((name, frame)))
    for name, frame in frames:
        expected = []
        for epitope in ('A', 'B'):
            expected.append(reference_epitope_rules(queries,
                reference.filter(pl.col('epitope') == epitope).drop('epitope'),
                reference.filter(pl.col('epitope') != epitope).drop('epitope'), exclude_exact=True)
                .with_columns(pl.lit(epitope).alias('epitope')))
        expected = pl.concat(expected).sort('query_id', 'epitope')
        assert frame.sort('query_id', 'epitope')['score'].to_list() == pytest.approx(expected[name+'_score'].to_list())
