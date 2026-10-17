"""Tiny raw-input checks of fresh control calibration and shared native search."""
from argparse import Namespace

import polars as pl
import pytest
from seqtree import Index

from vdjmatch.cli import control_calibration as calibration


def _inputs(tmp_path, monkeypatch, statistic='nearest'):
    reference = tmp_path / 'custom.tsv'
    pl.DataFrame({'cdr3': ['CASSF', 'CATSF', 'CASSF'], 'gene': ['TRB'] * 3,
        'species': ['HomoSapiens'] * 3, 'v': ['TRBV1'] * 3, 'j': ['TRBJ1-1'] * 3,
        'epitope': ['A', 'B', 'A'], 'evidence_type': ['binding', 'functional', 'binding']}).write_csv(reference, separator='\t')
    sample = tmp_path / 'query.tsv'
    pl.DataFrame({'query_id': ['exact', 'near', 'far', 'invalid'],
        'junction_aa': ['CASSF', 'CARSF', 'YYYYY', 'CA*F'],
        'v_call': ['TRBV1'] * 4, 'j_call': ['TRBJ1-1'] * 4, 'locus': ['TRB'] * 4,
        'pgen': [0.9] * 4}).write_csv(sample, separator='\t')
    control = tmp_path / 'control.tsv'
    pl.DataFrame({'junction_aa': ['CASSF', 'CATSF', 'CAAAF', 'CGGGF', 'WWWWW'],
        'occurrences': [2, 1, 3, 1, 2], 'reads': [100, 1, 200, 1, 1]}).write_csv(control, separator='\t')
    calls = []
    def pgen(reference, control, **kw):
        seqs = pl.concat([reference.select('junction_aa'), control.select('junction_aa')]).unique()
        calls.append(seqs.height)
        return seqs.with_columns(pl.lit(0.01).alias('pgen'))
    monkeypatch.setattr(calibration, 'compute_pgen', pgen)
    args = Namespace(sample=str(sample), vdjdb=str(reference), control=str(control),
        species='human', locus='TRB', threads=1, max_edits=1, scope=None, statistic=statistic,
        mock_replicates=3, seed=7, epitope=None, evidence_type=None, min_score=0,
        output_prefix=str(tmp_path / 'result'), input_format='auto', sequence_convention=None)
    return args, calls


def test_single_batch_raw_control_tails_and_filtered_recalibration(tmp_path, monkeypatch):
    args, pgen_calls = _inputs(tmp_path, monkeypatch)
    calls = []
    native = Index.search_batch
    def spy(self, queries, params, threads=0):
        calls.append((queries, threads))
        return native(self, queries, params, threads)
    monkeypatch.setattr(Index, 'search_batch', spy)
    from vdjmatch.cli.__main__ import main
    assert main(['database-evidence', args.sample, '--vdjdb', args.vdjdb,
        '--control', args.control, '--locus', 'TRB', '--max-edits', '1',
        '--mock-replicates', '3', '--seed', '7', '--threads', '1',
        '--output-prefix', args.output_prefix]) == 0
    result = pl.read_csv(args.output_prefix + '.database.tsv', separator='\t')
    assert len(calls) == len(pgen_calls) == 1
    assert result['query_id'].to_list() == ['exact', 'near', 'far', 'invalid']
    assert result['database_status'].to_list() == ['matched', 'matched', 'no_hit', 'invalid_sequence']
    assert result['database_nearest_edits'].to_list() == [0, 1, None, None]
    assert result['pgen'].to_list() == [0.01, 0.01, 0.01, None]
    assert result['p_both'][2] == result['p_query'][2] == result['p_reference'][2] == 1
    assert result['p_both'][3] is None
    assert result['database_exact_gene_clonotype'].to_list() == [True, False, False, None]
    assert all(result['p_both'][i] == max(result['p_query'][i], result['p_reference'][i]) for i in range(3))
    assert len(set(result['null_query_match_rate'].head(3))) == 1
    assert result['p_query'][0] <= result['p_query'][1]
    assert result['null_query_match_rate'][2] == result['query_hit_draws'][2] / args.mock_replicates
    groups = pl.read_csv(args.output_prefix + '.calibration.tsv', separator='\t')
    assert groups['reference_event_sampler'].unique().to_list() == ['uniform-event-slots-v1']
    assert groups['reference_source_records'].to_list() == [3, 2, 1]
    assert groups['native_batches'].to_list() == [1, 1, 1]
    assert groups['pgen_bin'].to_list() == [-7, -7, -7]
    assert groups['pgen_bin_reference_source_records'].to_list() == [3, 2, 1]
    assert groups['pgen_bin_control_unique_junctions'].to_list() == [5, 5, 5]
    assert groups['pgen_bin_control_events'].to_list() == [9, 9, 9]
    assert groups['pgen_model_source'].to_list() == ['olga'] * 3
    assert groups['seqtree_native_sha256'].str.len_chars().to_list() == [64] * 3
    candidates = pl.read_csv(args.output_prefix + '.epitopes.tsv', separator='\t')
    assert set(candidates['database_input_row']) == {0, 1}
    assert set(candidates['epitope']) == {'A', 'B'}
    assert not list(tmp_path.glob('result*.json'))
    args.evidence_type = ['binding']
    assert calibration.main(args) == 0
    groups = pl.read_csv(args.output_prefix + '.calibration.tsv', separator='\t')
    assert groups['reference_source_records'].to_list() == [2, 2]
    assert len(calls) == len(pgen_calls) == 2


def test_count_statistic_no_valid_queries_and_insufficient_pool(tmp_path, monkeypatch):
    args, _ = _inputs(tmp_path, monkeypatch, statistic='count')
    assert calibration.main(args) == 0
    result = pl.read_csv(args.output_prefix + '.database.tsv', separator='\t')
    assert result['database_hit_count'].to_list() == [2, 2, 0, None]
    assert result['database_score'].to_list() == [2.0, 2.0, None, None]
    pl.DataFrame({'junction_aa': ['CA*F']}).write_csv(args.sample, separator='\t')
    assert calibration.main(args) == 0
    assert pl.read_csv(args.output_prefix + '.database.tsv', separator='\t')['p_both'].to_list() == [None]
    pl.DataFrame({'junction_aa': ['CASSF']}).write_csv(args.control, separator='\t')
    with pytest.raises(ValueError, match='insufficient'):
        calibration.main(args)


def test_explicit_scope_freezes_both_null_directions(tmp_path, monkeypatch):
    args, _ = _inputs(tmp_path, monkeypatch)
    args.max_edits = 5
    args.scope = '0,0,0,0'
    assert calibration.main(args) == 0
    exact_scope = pl.read_csv(args.output_prefix + '.database.tsv', separator='\t')
    assert exact_scope['database_status'].to_list() == ['matched', 'no_hit', 'no_hit', 'invalid_sequence']
    assert exact_scope['p_both'][1] == 1
    args.scope = '1,0,0,1'
    assert calibration.main(args) == 0
    wider_scope = pl.read_csv(args.output_prefix + '.database.tsv', separator='\t')
    assert wider_scope['database_status'][1] == 'matched'
    assert wider_scope['null_query_match_rate'][0] >= exact_scope['null_query_match_rate'][0]
    assert wider_scope['null_reference_match_rate'][0] >= exact_scope['null_reference_match_rate'][0]


def test_query_control_hardlink_is_rejected(tmp_path, monkeypatch):
    import os
    args, _ = _inputs(tmp_path, monkeypatch)
    linked = tmp_path / 'same-input.tsv'
    os.link(args.sample, linked)
    args.control = str(linked)
    with pytest.raises(ValueError, match='cannot also define'):
        calibration.main(args)


def test_operation_specific_scopes_and_reserved_output_columns(tmp_path, monkeypatch):
    args, _ = _inputs(tmp_path, monkeypatch)
    args.scope = '0,0,0,1'
    assert calibration.main(args) == 0
    result = pl.read_csv(args.output_prefix + '.database.tsv', separator='\t')
    assert result['database_status'][1] == 'no_hit'
    args.scope = '1,0,0,1'
    assert calibration.main(args) == 0
    assert pl.read_csv(args.output_prefix + '.database.tsv', separator='\t')['database_status'][1] == 'matched'
    pl.DataFrame({'junction_aa': ['CASF', 'CASSAF', 'CARSF']}).write_csv(args.sample, separator='\t')
    args.scope = '0,1,0,1'
    assert calibration.main(args) == 0
    result = pl.read_csv(args.output_prefix + '.database.tsv', separator='\t')
    assert result['database_status'].to_list() == ['matched', 'no_hit', 'no_hit']
    assert result['database_nearest_edits'].to_list() == [1, None, None]
    for fields in ({'pgen': [0.1], 'source_pgen': [0.2]}, {'query_hit_draws': [3]}, {'null_query_match_rate': [0.1]}):
        pl.DataFrame({'junction_aa': ['CASSF'], **fields}).write_csv(args.sample, separator='\t')
        with pytest.raises(ValueError, match='reserved'):
            calibration.main(args)
