"""Exact graph-partition mixture and matched three-channel sum checks."""
import numpy as np
import polars as pl
import pytest

from vdjmatch.aggregate.reference_combination import (
    graph_mixture_weights, _component_colour_squares, combine_reference_evidence)


def test_graph_law_and_articulation_puncture_are_exact():
    reference = pl.DataFrame({'group_id': [0] * 5 + [1],
        'cdr3': ['CAAF', 'CCAF', 'CCCF', 'CAAFF', 'CAAF', 'CWWF'],
        'v': ['TRBV1'] * 6})
    full, punctures = graph_mixture_weights(reference, threads=1)
    assert full['source_records'].to_list() == [5, 1]
    assert full['clonotype_keys'].to_list() == [4, 1]
    np.testing.assert_allclose(full.select('original_weight', 'gapblock_weight', 'stem4_weight').row(0), [.625, .375, 0])
    bridge = punctures.filter((pl.col('group_id') == 0) & (pl.col('excluded_junction') == 'CCAF'))
    np.testing.assert_allclose(bridge.select('original_weight', 'gapblock_weight', 'stem4_weight').row(0), [3/9, 4/9, 2/9])
    assert bridge['components'].item() == 3
    assert punctures.filter((pl.col('group_id') == 0) & (pl.col('excluded_junction') == 'CAAF'))['source_records'].item() == 3
    singleton = punctures.filter(pl.col('group_id') == 1)
    assert singleton['original_weight'].item() is None
    assert singleton['weight_status'].item() == 'empty_reference'
    reversed_full, reversed_punctures = graph_mixture_weights(reference.reverse(), threads=1)
    assert full.equals(reversed_full) and punctures.equals(reversed_punctures)


def test_coloured_punctures_match_direct_fragment_oracle():
    edges = [(0, 1), (1, 2), (1, 3), (3, 4)]
    colours = np.array([[1, 1], [1, 0], [0, 1], [1, 1], [1, 0], [0, 1]], dtype=np.int64)
    full, punctures = _component_colour_squares(6, edges, colours)
    def oracle(removed):
        pending = set(range(6)) - {removed}
        total = 0
        while pending:
            component = {pending.pop()}
            stack = list(component)
            while stack:
                u = stack.pop()
                for a, b in edges:
                    v = b if a == u else a if b == u else None
                    if v in pending:
                        pending.remove(v)
                        component.add(v)
                        stack.append(v)
            mass = colours[list(component)].sum(axis=0)
            total += mass @ mass
        return total
    assert full == oracle(None)
    assert punctures.tolist() == [oracle(i) for i in range(6)]


def test_actual_normalization_and_missing_component_no_redistribution():
    context = dict(species='human', locus='TRB', epitope='E', mhc_a='A', mhc_b=None, mhc_class='MHCI')
    frame = pl.DataFrame([dict(query_id=str(i), group_id=0, original_score=o, gapblock_score=g, stem4_score=s, **context)
        for i, (o, g, s) in enumerate([(3.0, 30.0, 300.0), (None, 30.0, None), (None, None, None)])])
    null = pl.DataFrame([dict(query_id=str(i), original_score=float(i), gapblock_score=float(i * 10),
        stem4_score=float(i * 100), outcome='null', null_source='raw-control', **context) for i in range(3)])
    weights = pl.DataFrame({'group_id': [0], 'original_weight': [0.75], 'gapblock_weight': [0.25], 'stem4_weight': [0.0]})
    result = combine_reference_evidence(frame, null, weights)
    np.testing.assert_allclose(result['normalized_score'].head(2), [np.log(4), 0.25 * np.log(4)])
    assert result['normalized_score'][2] is None
    assert result['query_id'].to_list() == ['0', '1', '2']
    assert all(c in result.columns for c in ['contribution_original', 'contribution_gapblock', 'contribution_stem4'])
    with pytest.raises(ValueError, match='missing reference weights'):
        combine_reference_evidence(frame, null, weights.with_columns(pl.lit(1).alias('group_id')))


def test_structure_cli_exports_mixture_weights_for_unmodeled_v_once(tmp_path, monkeypatch):
    from vdjmatch.cli.__main__ import main
    from vdjmatch import cluster
    reference = tmp_path / 'reference.tsv'
    pl.DataFrame({'cdr3': ['CASSLGQETQYF', 'CASSIGQETQYF'], 'v': ['TRBV999', 'TRBV1'],
        'j': ['TRBJ1-1'] * 2, 'gene': ['TRB'] * 2, 'species': ['HomoSapiens'] * 2,
        'epitope': ['E'] * 2, 'mhc_a': ['HLA-A*02:01'] * 2, 'mhc_b': ['B2M'] * 2,
        'mhc_class': ['MHCI'] * 2}).write_csv(reference, separator='\t')
    calls = []
    original = cluster.overlap
    def spy(*args, **kwargs):
        calls.append(len(args[0]))
        return original(*args, **kwargs)
    monkeypatch.setattr(cluster, 'overlap', spy)
    prefix = tmp_path / 'weights'
    assert main(['reference-structure', '--vdjdb', str(reference), '--species', 'human', '--locus', 'TRB',
        '--reference-population', 'all-valid', '--threads', '1', '--output-prefix', str(prefix)]) == 0
    full = pl.read_csv(str(prefix) + '.structures.tsv', separator='\t')
    punctures = pl.read_csv(str(prefix) + '.punctures.tsv', separator='\t')
    assert calls == [2]
    assert full['junctions'].item() == 1
    assert full['mixture_junctions'].item() == 2
    assert full['hamming_junctions'].item() == 2
    degrees = pl.read_csv(str(prefix) + '.degrees.tsv', separator='\t')
    assert degrees.select('degree', 'junction_count').rows() == [(1, 2)]
    assert sum(full.select('original_weight', 'gapblock_weight', 'stem4_weight').row(0)) == 1
    unknown = punctures.filter(pl.col('excluded_junction') == 'CASSLGQETQYF')
    assert unknown['junctions'].item() is None
    assert unknown['mixture_junctions'].item() == 1
    assert unknown['mixture_only_puncture'].item()
    assert punctures.height == 2


def test_automatic_normalization_strata_keep_matching_banks_separate(tmp_path):
    from vdjmatch.aggregate.combination import null_normalized_evidence
    from vdjmatch.cli.__main__ import main
    context = dict(species='human', locus='TRB', epitope='E', mhc_a='A', mhc_b='B', mhc_class='MHCI')
    observed = pl.DataFrame([dict(query_id='q' + s, normalization_stratum=s, original_score=3.0,
        gapblock_score=None, stem4_score=None, original_weight=1.0, gapblock_weight=0.0, stem4_weight=0.0, **context)
        for s in ('A', 'B')])
    null = pl.DataFrame([dict(query_id=str(i), normalization_stratum=s, original_score=value,
        gapblock_score=None, stem4_score=None, outcome='null', null_source='bank', **context)
        for i, (s, value) in enumerate([('A', 1.0), ('A', 2.0), ('B', 3.0), ('B', 4.0)])])
    scored = null_normalized_evidence(observed, null, three_channel=True)
    np.testing.assert_allclose(scored['p_component_original'], [1 / 3, 1.0])
    with pytest.raises(ValueError, match='both evidence'):
        null_normalized_evidence(observed, null.drop('normalization_stratum'), three_channel=True)
    with pytest.raises(ValueError, match='non-null'):
        null_normalized_evidence(observed.with_columns(pl.lit(None).alias('normalization_stratum')), null, three_channel=True)
    source, bank, output = [tmp_path / s for s in ('source.tsv', 'bank.tsv', 'score.tsv')]
    observed.write_csv(source, separator='\t')
    null.write_csv(bank, separator='\t')
    assert main(['combine-evidence', str(source), '--null-evidence', str(bank), '--three-channel',
        '--output', str(output)]) == 0
    assert pl.read_csv(output, separator='\t')['combined_score'].to_list() == pytest.approx([np.log(3), 0.0])


def test_shipping_sum_retains_raw_score_and_applies_declared_support(tmp_path):
    from vdjmatch.cli.__main__ import main
    context = dict(species='human', locus='TRB', epitope='E', mhc_a='A', mhc_b='B', mhc_class='MHCI')
    observed = pl.DataFrame([dict(query_id='q', original_score=3.0, original_supported=False,
        gapblock_score=None, stem4_score=None, original_weight=1.0, gapblock_weight=0.0,
        stem4_weight=0.0, **context)])
    null = pl.DataFrame([dict(query_id=str(i), original_score=float(i), gapblock_score=None,
        stem4_score=None, outcome='null', null_source='bank', **context) for i in range(2)])
    source, bank, output = [tmp_path / s for s in ('source.tsv', 'bank.tsv', 'score.tsv')]
    observed.write_csv(source, separator='\t'); null.write_csv(bank, separator='\t')
    command = ['combine-evidence', str(source), '--null-evidence', str(bank), '--three-channel', '--output', str(output)]
    assert main(command) == 0
    result = pl.read_csv(output, separator='\t')
    assert result['normalized_score'].item() == pytest.approx(np.log(3))
    assert result['combined_score'].item() is None
    assert result['combined_status'].item() == 'unsupported'
    observed.drop('original_supported').write_csv(source, separator='\t')
    assert main(command) == 0
    assert pl.read_csv(output, separator='\t')['combined_score'].item() == pytest.approx(np.log(3))
