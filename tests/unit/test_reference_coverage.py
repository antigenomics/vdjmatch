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


@pytest.mark.parametrize('geometry', ['levenshtein', 'both'])
def test_combined_diagnostics_loads_control_and_pgen_once(tmp_path, monkeypatch, geometry):
    pytest.importorskip('leidenalg')
    import json
    import vdjmatch.evalue.mock_reference as mocks
    from vdjmatch.cli.__main__ import main
    reference, control, prefix = (tmp_path / p for p in ('ref.tsv', 'control.tsv', 'out'))
    seqs = ['CASSLAAAF', 'CASSLAAAAF', 'CAGGGGGGF']
    pl.DataFrame({'cdr3': seqs, 'v': ['TRBV1']*3, 'j': ['TRBJ1-1']*3,
        'species': ['HomoSapiens']*3, 'gene': ['TRB']*3, 'epitope': ['DEMO']*3}).write_csv(reference, separator='\t')
    pl.DataFrame({'junction_aa': seqs, 'v_call': ['TRBV1']*3, 'j_call': ['TRBJ1-1']*3,
        'occurrences': [10]*3}).write_csv(control, separator='\t')
    calls = []
    def pgen(ref, pool, **options):
        calls.append('pgen')
        return pool.select('junction_aa').with_columns(pl.lit(.25).alias('pgen'))
    monkeypatch.setattr(mocks, 'compute_pgen', pgen)
    prepare = mocks._prepare_mock_control
    def prepare_once(*args, **kwargs):
        calls.append('prepare')
        return prepare(*args, **kwargs)
    monkeypatch.setattr(mocks, '_prepare_mock_control', prepare_once)
    args = ['reference-diagnostics', '--vdjdb', str(reference), '--control', str(control),
        '--locus', 'TRB', '--geometry', geometry, '--scope', '2,1,1,2', '--mock-replicates', '2', '--nboot', '3',
        '--output-prefix', str(prefix)]
    if geometry == 'both':
        args += ['--sw-radius', '20']
    assert main(args) == 0
    assert calls == ['pgen', 'prepare']
    motifs = pl.read_csv(str(prefix)+'.motifs.tsv', separator='\t')
    assert set(motifs['motif_mode']) == {'exact3', 'wild4', 'imgt3', 'imgt_wild4'}
    manifest = json.loads((tmp_path/'out.levenshtein.manifest.json').read_text())
    assert manifest['shared_invocation_inputs'] is True


def test_sw_knee_uses_reference_nearest_within_pmhc_and_keeps_missing_models(monkeypatch):
    from Bio.Align import PairwiseAligner, substitution_matrices
    from vdjmatch.aggregate.reference_coverage import sw_knee, reference_graph
    from vdjmatch.match import geometry
    records = pl.DataFrame({'cdr3': ['CASSLAAAF', 'CASSLAAAAF', 'CAGGGGGGF', 'CASSLAAAF', 'CASSLAAAF'],
        'v': ['TRBV1', 'TRBV1', 'TRBV1', 'TRBV1', 'TRBV999'],
        'j': ['TRBJ1-1']*5, 'group_id': [0, 0, 0, 1, 0]})
    captured = []
    def knee(values):
        captured.extend(sorted(values))
        return min(values)
    monkeypatch.setattr(geometry, 'kneedle_radius', knee)
    radius, curve = sw_knee(records, species='human', locus='TRB',
        cdr_loops=('cdr1', 'cdr2', 'cdr25'))
    aligner = PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'),
        open_gap_score=-11, extend_gap_score=-1)
    seqs = records['cdr3'][:3].to_list()
    expected = sorted(min(aligner.score(a, a)+aligner.score(b, b)-2*aligner.score(a, b)
        for j, b in enumerate(seqs) if i != j) for i, a in enumerate(seqs))
    assert captured == expected and radius == min(expected)
    assert curve['nearest_distance'].null_count() == 2  # singleton pMHC and unresolved V
    nodes, edges = reference_graph(records, geometry='sw', radius=radius, species='human', locus='TRB')
    missing = nodes.filter(~pl.col('geometry_available'))
    assert missing.height == 1 and missing['occurrences'].sum() == 1
    assert not edges.filter(pl.col('a').is_in(missing['node'].implode()) | pl.col('b').is_in(missing['node'].implode())).height


def test_edit_search_deduplicates_geometry_only_and_bounds_actual_edges(monkeypatch):
    import vdjmatch.cluster as cluster
    from vdjmatch.aggregate.reference_coverage import reference_graph
    records = pl.DataFrame({'cdr3': ['CASSLAAAF', 'CASSLAAAF', 'CASSLAAAAF']*4,
        'v': ['TRBV1']*12, 'j': ['TRBJ1-1', 'TRBJ1-2', 'TRBJ1-1']*4,
        'group_id': np.repeat(np.arange(4), 3)})
    native = cluster.overlap
    calls = []
    def batched(seqs, **kwargs):
        calls.append(seqs)
        return native(seqs, **kwargs)
    monkeypatch.setattr(cluster, 'overlap', batched)
    nodes, edges = reference_graph(records, geometry='levenshtein', pair_budget=1)
    assert calls == [['CASSLAAAAF', 'CASSLAAAF']]
    assert nodes.height == 12 and edges.height == 12
    groups = nodes.select(pl.col('node').cast(pl.Int64).alias('a'), 'group_id')
    assert edges.join(groups, on='a').group_by('group_id').len()['len'].to_list() == [3]*4


def test_combined_invalid_graph_fails_before_loading_raw_inputs(tmp_path, capsys):
    from vdjmatch.cli.__main__ import main
    with pytest.raises(SystemExit, match='2'):
        main(['reference-diagnostics', '--vdjdb', str(tmp_path/'missing'), '--control', str(tmp_path/'missing-control'),
            '--locus', 'TRB', '--output-prefix', str(tmp_path/'out')])
    assert '--sw-radius or --sw-radius-knee' in capsys.readouterr().err


def test_native_bootstrap_coverage_matches_owning_python_curve():
    from vdjtools.stats import bootstrap_abundances, sample_coverage
    from vdjmatch.aggregate.reference_coverage import abundance_coverage
    for counts in ([1, 1, 2, 4], [1]*15000):
        _, curve = abundance_coverage(counts, nboot=5, seed=7, threads=2)
        depths = curve['m'].to_list()
        expected = np.asarray([sample_coverage(draw, m=depths)['sample_coverage'].to_numpy()
            for draw in bootstrap_abundances(counts, nboot=5, seed=7)])
        np.testing.assert_allclose(curve['coverage_lo'], np.quantile(expected, .025, axis=0), rtol=0, atol=1e-9)
        np.testing.assert_allclose(curve['coverage_hi'], np.quantile(expected, .975, axis=0), rtol=0, atol=1e-9)
    with pytest.raises(ValueError, match='working memory'):
        abundance_coverage([10_000_000], nboot=5, seed=7)
