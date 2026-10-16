"""Small independent graph/occurrence and finite-sample Chao checks."""
import math
import numpy as np
import polars as pl
import pytest


def test_coverage_keeps_occurrences_and_distinguishes_unseen_mass_from_types():
    from vdjmatch.aggregate.reference_coverage import abundance_coverage
    summary, curve = abundance_coverage([1, 1, 2, 4], nboot=5, seed=7, threads=1)
    assert summary['observations'] == 8 and summary['observed_taxa'] == 4
    assert summary['singletons'] == 2 and summary['doubletons'] == 1
    assert summary['estimated_taxa'] == pytest.approx(4 + 7/8 * 4/2)
    assert summary['coverage'] == pytest.approx(1 - 2/8 * 7/8)
    assert summary['unseen_mass'] == pytest.approx(1 - summary['coverage'])
    assert summary['unseen_type_fraction'] == pytest.approx(1 - 4/summary['estimated_taxa'])
    assert summary['entropy_nats'] == pytest.approx(-sum(p*math.log(p) for p in [1/8, 1/8, 2/8, 4/8]))
    observed = curve.filter(pl.col('method') == 'observed').row(0, named=True)
    assert observed['sample_coverage'] == pytest.approx(summary['coverage'])
    assert observed['qD'] == 4 and observed['qD_lo'] <= 4 <= observed['qD_hi']
    assert curve['coverage_lo'].min() >= 0 and curve['coverage_hi'].max() <= 1
    assert curve.equals(abundance_coverage([1, 1, 2, 4], nboot=5, seed=7, threads=1)[1])
    zero, _ = abundance_coverage([1, 1, 3], nboot=5, seed=7, threads=1)
    assert zero['estimated_taxa'] == pytest.approx(3 + 4/5)
    single, _ = abundance_coverage([1], nboot=3, seed=7, threads=1)
    assert single['estimated_taxa'] == 1 and single['coverage_status'] == 'zero_doubleton_formula'
    for counts in ([], [0, 2], [1.5, 2], [True]):
        with pytest.raises(ValueError):
            abundance_coverage(counts, nboot=5, seed=7)
    with pytest.raises(ValueError, match='int32'):
        abundance_coverage([2**31], nboot=5, seed=7)


@pytest.mark.parametrize('geometry', ['levenshtein', 'sw'])
def test_graphs_batch_once_preserve_duplicate_occurrences_and_real_indels(geometry):
    pytest.importorskip('leidenalg')
    from vdjmatch.aggregate.reference_coverage import reference_graph, leiden_taxa
    records = pl.DataFrame({'cdr3': ['CASSLAAAF', 'CASSLAAAF', 'CASSLAAAAF', 'CAGGGGGGF'],
        'v': ['TRBV1'] * 4, 'j': ['TRBJ1-1'] * 4, 'group_id': [0] * 4})
    nodes, edges = reference_graph(records, geometry=geometry, radius=20,
        species='human', locus='TRB', threads=1)
    assert nodes.height == 3 and nodes['occurrences'].sum() == 4
    assert edges.height == 1
    membership = leiden_taxa(nodes, edges, resolution=1., seed=7)
    assert membership.group_by('community').agg(pl.col('occurrences').sum()).height == 2
    assert sorted(membership.group_by('community').agg(pl.col('occurrences').sum())['occurrences']) == [1, 3]
    assert membership.equals(leiden_taxa(nodes, edges, resolution=1., seed=7))
    with pytest.raises(ValueError, match='pair budget'):
        reference_graph(records, geometry=geometry, radius=20, species='human', locus='TRB', pair_budget=1)
    if geometry == 'levenshtein':
        with pytest.raises(ValueError, match='symmetric'):
            reference_graph(records, geometry=geometry, scope='1,0,1,1')


def test_full_sw_cdr_loop_distance_matches_independent_local_alignment():
    from vdjmatch.aggregate.reference_coverage import reference_graph
    from vdjmatch.match.geometry import load_v_loops
    from Bio.Align import PairwiseAligner, substitution_matrices
    records = pl.DataFrame({'cdr3': ['CASSLAAAF', 'CASSLAAAAF'],
        'v': ['TRBV1', 'TRBV19'], 'j': ['TRBJ1-1'] * 2, 'group_id': [0, 0]})
    aligner = PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'),
        open_gap_score=-11, extend_gap_score=-1)
    model = load_v_loops('human')
    expected = 0
    for a, b in [(records['cdr3'][0], records['cdr3'][1]),
                 *[(model['TRBV1*01'][c].replace('.', ''), model['TRBV19*01'][c].replace('.', ''))
                   for c in ('cdr1', 'cdr2', 'cdr25')]]:
        expected += aligner.score(a, a) + aligner.score(b, b) - 2*aligner.score(a, b)
    nodes, edges = reference_graph(records, geometry='sw', radius=int(expected), species='human', locus='TRB')
    assert edges.height == 1 and edges['distance_ceiling'].item() == expected
    assert reference_graph(records, geometry='sw', radius=int(expected)-1, species='human', locus='TRB')[1].is_empty()


@pytest.mark.parametrize('geometry', ['levenshtein', 'sw'])
def test_reference_coverage_cli_uses_one_batch_for_real_and_control_graphs(tmp_path, monkeypatch, geometry):
    pytest.importorskip('leidenalg')
    import json
    import vdjmatch.evalue.mock_reference as mocks
    import vdjmatch.cluster as cluster
    from seqtree import gapblock
    from vdjmatch.cli.__main__ import main
    reference, control, prefix = (tmp_path / p for p in ('reference.tsv', 'control.tsv', 'coverage'))
    junctions = ['CASSLAAAF', 'CASSLAAAF', 'CASSLAAAAF', 'CAGGGGGGF']
    pl.DataFrame({'cdr3': junctions, 'v': ['TRBV1'] * 4, 'j': ['TRBJ1-1'] * 4,
        'species': ['HomoSapiens'] * 4, 'gene': ['TRB'] * 4, 'epitope': ['DEMO'] * 4}).write_csv(reference, separator='\t')
    pl.DataFrame({'junction_aa': sorted(set(junctions)), 'v_call': ['TRBV1'] * 3,
        'j_call': ['TRBJ1-1'] * 3, 'occurrences': [10] * 3}).write_csv(control, separator='\t')
    calls = []
    def pgen(ref, pool, **options):
        calls.append('pgen')
        return pool.select('junction_aa').with_columns(pl.lit(.25).alias('pgen'))
    monkeypatch.setattr(mocks, 'compute_pgen', pgen)
    if geometry == 'levenshtein':
        native = cluster.overlap
        def batched(*args, **kwargs):
            calls.append('graph')
            return native(*args, **kwargs)
        monkeypatch.setattr(cluster, 'overlap', batched)
    else:
        from functools import wraps
        native = gapblock.grouped_count_batch
        @wraps(native)
        def batched(*args, **kwargs):
            calls.append('graph')
            return native(*args, **kwargs)
        monkeypatch.setattr(gapblock, 'grouped_count_batch', batched)
    args = ['reference-coverage', '--vdjdb', str(reference), '--control', str(control),
        '--locus', 'TRB', '--geometry', geometry, '--nboot', '3', '--mock-replicates', '2',
        '--output-prefix', str(prefix)]
    if geometry == 'sw':
        args += ['--sw-radius', '20']
    assert main(args) == 0
    assert calls == ['pgen', 'graph']
    summary = pl.read_csv(str(prefix)+'.coverage.tsv', separator='\t')
    assert summary.height == 6 and set(summary['population']) == {'reference', 'matched_control'}
    assert summary['observations'].to_list() == [4] * 6
    observed = summary.filter(pl.col('population') == 'reference')
    assert observed['observed_taxa'].to_list() == [3, 2]
    members = pl.read_csv(str(prefix)+'.communities.tsv', separator='\t')
    assert members.group_by('group_id').agg(pl.col('occurrences').sum())['occurrences'].to_list() == [4]*3
    manifest = json.loads((tmp_path/'coverage.manifest.json').read_text())
    assert manifest['search_batches'] == 1 and manifest['nboot'] == 3
    assert len(manifest['tables']) == 5 and manifest['control']['replicates'] == 2
