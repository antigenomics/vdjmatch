"""All-label sparse projections agree with the fixed scalar rule oracle."""
import polars as pl
import pytest
import numpy as np


@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('one_peptide', [False, True])
def test_native_peptide_channels_exactly_match_legacy_artificial_models(monkeypatch, paired, one_peptide):
    from types import SimpleNamespace
    from vdjmatch.evalue.reference_rule_bank import stream_reference_rules
    from vdjmatch.match import geometry

    reference = pl.DataFrame({'seq0': ['CASSLGQAYEQYF', 'CASSLGQAYEQYF',
        'CASSLGRAYEQYF', 'CASSWGRAYEQYF', 'CAS'+'A'*120+'QYF'],
        'v0': ['TRBV19*01', 'TRBV19*01', 'TRBV19*02', 'TRBV20-1*01', 'TRBV19*01'],
        'epitope': ['A', 'B', 'A', 'B', 'C']})
    if one_peptide:
        reference = reference.with_columns(pl.lit('A').alias('epitope'))
    if paired:
        reference = reference.rename({'seq0': 'seq1', 'v0': 'v1'}).with_columns(
            pl.lit('CAVRLGQAGKSTF').alias('seq0'), pl.lit('TRAV1-1*01').alias('v0'))
    keys = ['seq0', 'v0', *(['seq1', 'v1'] if paired else [])]
    queries = reference.select(keys).unique(maintain_order=True).with_row_index('query_id')
    queries = pl.concat([queries, queries.head(1).with_columns(
        pl.lit(queries.height, dtype=pl.UInt32).alias('query_id'),
        pl.lit('TRAV2*01' if paired else 'TRBV6-5*01').alias('v0'))])
    native = geometry.grouped_total_distance_count_batch
    calls = []
    def legacy_channels(*args, **kwargs):
        T, R = len(kwargs['target_peptide_ids']), args[1].height
        assert kwargs['target_peptide_ids'] == list(range(T))
        assert 'return_hits' not in kwargs
        assert kwargs['chunk_size']*4*((1+3*T)*(len(args[2])+1)+T) <= 128*1024**2
        ref = np.repeat(np.arange(R), np.diff(kwargs['reference_target_offsets']))
        labels = np.asarray(kwargs['reference_target_ids'])
        exclusive = np.asarray(kwargs['reference_exclusive_peptide_ids'])
        chosen = np.flatnonzero(exclusive >= 0)
        refs = np.r_[ref, np.arange(R), chosen]
        ids = np.r_[labels, np.full(R, T), T+1+exclusive[chosen]]
        order = np.lexsort((ids, refs))
        models = 2*T+1
        legacy = kwargs | dict(reference_target_offsets=np.r_[0, np.cumsum(np.bincount(refs, minlength=R))],
            reference_target_ids=ids[order], target_peptide_ids=[0]*models,
            reference_exclusive_peptide_ids=[-1]*R, chunk_size=1)
        calls.append(kwargs['same_v_gene'])
        for chunk in native(*args, **legacy):
            counts = np.asarray(chunk.counts)
            assert counts.shape[1] == 4*T+4
            assert np.array_equal(counts[:, 0], counts[:, 1+T])
            channels = np.r_[0, 1+np.arange(T), 1+models+np.arange(T), 1+T+1+np.arange(T)]
            yield SimpleNamespace(start=chunk.start, stop=chunk.stop, counts=counts[:, channels])
    def collect():
        frames = {}
        stream_reference_rules(queries, reference, paired=paired,
            consume=lambda name, frame: frames.setdefault(name, []).append(frame))
        return {name: pl.concat(parts) for name, parts in frames.items()}
    actual = collect()
    monkeypatch.setattr(geometry, 'grouped_total_distance_count_batch', legacy_channels)
    expected = collect()
    assert calls == [True, False]
    assert len(actual) == 15
    for name in actual:
        assert actual[name].equals(expected[name]), name


@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('radius', [.005, 90.005])
def test_streamed_rules_match_every_scalar_model(paired, radius):
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
    stream_reference_rules(queries, reference, paired=paired, species='human', threads=1, radius=radius, consume=consume)
    expected = []
    for epitope in sorted(reference['epitope'].unique()):
        expected.append(reference_epitope_rules(queries,
            reference.filter(pl.col('epitope') == epitope).select(keys),
            reference.filter(pl.col('epitope') != epitope).select(keys),
            paired=paired, exclude_exact=True, radius=radius).with_columns(pl.lit(epitope).alias('epitope')))
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
