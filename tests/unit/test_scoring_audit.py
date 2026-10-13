"""Synthetic legacy reproduction checks; never run the scientific dataset."""
from collections import Counter
import importlib.util
from pathlib import Path
import sys

import numpy as np
import polars as pl
import pytest

from vdjmatch.match import scoring_audit as audit


def _legacy(name):
    # Oracle source is intentionally available only in the repository test.
    path = Path(__file__).resolve().parents[2] / 'bench' / (name+'.py')
    spec = importlib.util.spec_from_file_location('_legacy_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _frame():
    return pl.DataFrame({'cdr3': ['CASSAAAAAAAF', 'CASSACAAAAAF', 'CASSADAAAAAF',
                                 'CASSAAAAAAAF', 'CASSACAAAAAF', 'CASSAEAAAAAF'],
        'epitope': ['X']*3+['Y']*3, 'v': ['TRBV1']*6, 'j': ['TRBJ1']*6,
        'mhc_class': ['MHCI']*6, 'reference_id': ['S']*6})


def test_native_events_preserve_shared_epitope_multiplicity_and_threads():
    frame = _frame()
    old = _legacy('gen_vdjam')
    expected = []
    for epi in ['X', 'Y']:
        expected.extend((epi, *row) for row in old.hamming1_events(
            frame.filter(pl.col('epitope') == epi)['cdr3'].to_list(), 2))
    actual = audit.substitution_events(frame, threads=1).sort('epitope', 'a', 'b', 'pos', 'L')
    assert Counter(actual.iter_rows()) == Counter(expected)
    assert actual.equals(audit.substitution_events(frame, threads=2).sort('epitope', 'a', 'b', 'pos', 'L'))
    assert actual.height == 6  # Shared A/C pair is retained separately for X and Y.


def test_matrix_and_metric_match_legacy():
    old = _legacy('gen_vdjam')
    metrics = _legacy('metrics')
    seqs = _frame()['cdr3'].to_list()
    events = old.hamming1_events(seqs, 2)
    background = audit.position_background(seqs)
    assert audit.estimate_matrix(events, background) == old.estimate_matrix(events, background)
    for pairs in [[(1, 0), (0, 0)], [(0, 1), (1, 0), (0, 0), (1, -2)]]:
        assert audit.pr_auc_balanced(pairs) == metrics.pr_auc_balanced(pairs)
        assert audit.pr_auc_balanced(pairs[::-1]) == metrics.pr_auc_balanced(pairs)
    assert audit.pr_auc_balanced([(1, 0), (0, 0)]) == .75  # Accepted trapezoid, not AP.


def test_loo_scores_and_no_hits_preserve_legacy_contract():
    frame = _frame()
    ret, sig = {}, ([.9, .5, .2], [.9, .5, .2])
    tasks, summary, pairs, queries, refs = audit.loo_audit(frame, locus='TRB', retention=ret,
        significance=sig, structural=None, min_epi=2, top=2, max_queries=3, threads=1)
    assert set(tasks['method']) == {'unit', 'blosum', 'pam', 'vdjamr', 'bps', 'vdjam', 'region'}
    assert tasks['candidate_pairs'].to_list() == [9]*14
    assert set(refs.filter(pl.col('cdr3') == 'CASSAAAAAAAF')['epitope']) == {'X'}
    assert queries.height == 6
    assert summary.height == 7
    assert pairs.height == 126
    assert tasks.equals(audit.loo_audit(frame, locus='TRB', retention=ret,
        significance=sig, structural=None, min_epi=2, top=2, max_queries=3, threads=2)[0])
    nohit = audit.loo_audit(frame, locus='TRB', retention=ret, significance=sig, structural=None,
        min_epi=2, top=2, max_queries=3, subs=0)[0]
    assert set(nohit['candidate_pairs']) == {0}
    assert set(nohit['no_hit_queries']) == {3}
    assert all(np.isnan(nohit['balanced_pr_auc'].to_numpy()))


def test_vectorized_scores_match_legacy_pair_scores(monkeypatch):
    # Load the old scorer without importing its CLI-only module paths globally.
    monkeypatch.setitem(sys.modules, '_bench', _legacy('_bench'))
    monkeypatch.setitem(sys.modules, 'metrics', _legacy('metrics'))
    monkeypatch.setitem(sys.modules, 'gen_vdjam', _legacy('gen_vdjam'))
    old = _legacy('loo_vdjam')
    frame = _frame().unique('cdr3', maintain_order=True)
    cand = [[0, 1, 2, 3]] * frame.height
    from types import SimpleNamespace
    native = [[SimpleNamespace(ref_id=i) for i in hit] for hit in cand]
    dis = audit.named_dissim('BLOSUM62')
    ret = {('TRB', 'V', 'TRBV1'): [.9]*5}
    ids, scores = audit._score_candidates(frame, frame, native, locus='TRB', retention=ret,
        significance=([.9, .5, .2], [.9, .5, .2]), matrices={'flat': (dis, None), 'region': (dis, 'region')})
    for name, mode in [('flat', None), ('region', 'region')]:
        expected = old.score_pairs(frame['cdr3'].to_list(), frame['v'].to_list(), frame['j'].to_list(),
            cand, frame['cdr3'].to_list(), frame['epitope'].to_list(), 'X', 'TRB', ret, dis, mode)
        np.testing.assert_allclose(scores[name], [s for _, s in expected], rtol=0, atol=0)
    assert len(ids) == 12
    assert audit.codon_dissim() == old.codon_dissim()
    assert audit.named_dissim('BLOSUM62') == old.named_dissim('BLOSUM62')
    assert audit.named_dissim('PAM250') == old.named_dissim('PAM250')


def test_structural_malformed_fails(tmp_path):
    path = tmp_path/'structural.inc'
    path.write_text('kStructural[24 * 24] = {1, 2};')
    with pytest.raises(ValueError, match='576'):
        audit.structural_dissim(path)


@pytest.mark.parametrize('task', ['matrix', 'loo'])
def test_cli_rebuilds_flat_tables_and_provenance(tmp_path, task):
    import argparse
    import json
    from vdjmatch.cli.scoring_audit import register
    raw = tmp_path/'raw.tsv'
    _frame().with_columns(pl.lit('TRB').alias('gene'), pl.lit('HomoSapiens').alias('species')).write_csv(raw, separator='\t')
    prefix = tmp_path/task
    parser = argparse.ArgumentParser()
    register(parser.add_subparsers())
    args = parser.parse_args(['scoring-audit', '--task', task, '--pmhc', str(raw), '--locus', 'TRB',
        '--min-epi', '2', '--top', '2', '--max-queries', '3', '--output-prefix', str(prefix)])
    assert args.func(args) == 0
    manifest = json.loads(Path(str(prefix)+'.provenance.json').read_text())
    assert manifest['source_chain_records'] == 6
    assert manifest['reference']['input_sha256']
    assert manifest['structural_status'].startswith('omitted')
    assert manifest['options']['seed'] == 0
    for name, count in manifest['output_rows'].items():
        table = pl.read_csv(str(prefix)+'.'+name+'.tsv', separator='\t')
        assert table.height == count
        assert not any(isinstance(dtype, (pl.List, pl.Struct)) for dtype in table.dtypes)
    if task == 'matrix':
        corr = pl.read_csv(str(prefix)+'.correlations.tsv', separator='\t')
        assert set(corr['region']) == {'all', 'V', 'NDN', 'J'}
        assert corr.filter(pl.col('region') == 'V')['events'].item() == 0
