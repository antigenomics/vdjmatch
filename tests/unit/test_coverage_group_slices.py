"""Independent population-filter oracle for sorted graph reduction and budgets."""
from argparse import Namespace
import hashlib
import json

import polars as pl
import pytest


@pytest.fixture
def coverage_fixture(tmp_path, monkeypatch):
    pytest.importorskip('leidenalg')
    from vdjmatch import db
    from vdjmatch.cli import coverage, epitope, structure
    from vdjmatch.aggregate import reference_coverage as owner

    sources = []
    for group in range(5):
        seqs = (['CASSF', 'CASSF', 'CATTF'] if group == 0 else
                ['CASGF']*3 if group == 2 else ['CAGGF'] if group == 4 else ['X'])
        gene = 'TRAV1-1' if group == 0 else 'TRAV2' if group == 2 else None
        sources.append(({'epitope': f'fixture{group}'}, pl.DataFrame(dict(
            cdr3=seqs, v=[gene]*len(seqs), j=[None]*len(seqs),
            reference_valid=[group in (0, 2, 4)]*len(seqs)))))
    state = dict(graph_calls=0, leiden=[], abundance=[], unsorted=False)

    def graph(records, **kwargs):
        state['graph_calls'] += 1
        nodes = records.group_by('group_id', 'cdr3', 'v', 'j').len(name='occurrences').sort(
            'group_id', 'cdr3', 'v', 'j', nulls_last=True).with_row_index('node').with_columns(
            pl.col('v').is_not_null().alias('geometry_available'))
        state['nodes'] = nodes
        if state['unsorted']:
            nodes = nodes.reverse()
        empty = pl.DataFrame(schema={'a': pl.Int64, 'b': pl.Int64, 'distance_ceiling': pl.Int64})
        return nodes, empty

    def leiden(nodes, edges, **kwargs):
        state['leiden'].append((nodes['node'].to_list(), kwargs['seed']))
        return nodes.with_columns(pl.lit(0, dtype=pl.Int64).alias('community'))

    def abundance(counts, **kwargs):
        state['abundance'].append((counts.tolist(), kwargs['seed']))
        return (dict(observations=int(counts.sum()), coverage_status='fixture'),
                pl.DataFrame(dict(m=[1], sample_coverage=[.5])))

    def mocks(selected, events, *, replicates, **kwargs):
        for replicate in range(replicates):
            yield selected.with_columns(pl.lit(replicate).alias('replicate'))

    from vdjmatch.evalue import mock_reference
    monkeypatch.setattr(mock_reference, 'iter_mock_references', mocks)
    monkeypatch.setattr(structure, 'reference_groups', lambda raw, options: sources)
    monkeypatch.setattr(owner, 'reference_graph', graph)
    monkeypatch.setattr(owner, 'leiden_taxa', leiden)
    monkeypatch.setattr(owner, 'abundance_coverage', abundance)
    monkeypatch.setattr(epitope, 'write_tables', lambda *args: [])
    monkeypatch.setattr(db, 'provenance', lambda path: {'fixture': True})
    (tmp_path/'reference').write_text('software fixture')
    (tmp_path/'control').write_text('software fixture')
    (tmp_path/'scratch').mkdir()
    options = Namespace(vdjdb=str(tmp_path/'reference'), species='human', locus='TRA',
        control=None, control_model=None, pgen_features=None, targets_from_sample=None,
        reference_population='curated-benchmark', epitope=None, geometry='sw', sw_radius=5,
        sw_radius_knee=False, scope='1,1,1,1', cdr_loops='cdr12plus', gap_open=11, gap_extend=1,
        threads=1, nboot=200, seed=7, mock_replicates=0, pair_budget=100,
        leiden_resolution=1., output_prefix=str(tmp_path/'tables'))
    inputs = dict(raw=sources[0][1], epitope_info={}, prepared_control={},
        pgen_features={'model': 'software fixture'}, events=pl.DataFrame({'event_weight': [1]}),
        pgen=None, control_info={})
    return coverage, options, inputs, state


def seed_for(identity, resolution):
    return int.from_bytes(hashlib.sha256(json.dumps(
        [7, 'HomoSapiens', 'TRA', identity, resolution], sort_keys=True).encode()
        ).digest()[:4], 'little') % 2**31


@pytest.mark.parametrize('replicates', [0, 2])
def test_group_slices_match_independent_filter_rows_weights_and_seeds(
        tmp_path, coverage_fixture, replicates):
    coverage, options, inputs, state = coverage_fixture
    options.mock_replicates = replicates
    if replicates:
        options.control = str(tmp_path/'control')
    assert coverage._main(options, inputs=inputs, scratch=tmp_path/'scratch') == 0
    assert state['graph_calls'] == 1
    populations = pl.read_csv(tmp_path/'tables.populations.tsv', separator='\t')
    expected_leiden, expected_abundance = [], []
    for context in populations.to_dicts():
        selected = state['nodes'].filter(pl.col('group_id') == context['group_id'])
        if selected.is_empty():
            continue
        identity = {key: value for key, value in context.items() if key not in
            ('group_id', 'source_records', 'retained_records', 'excluded_records')}
        modeled = selected.filter(pl.col('geometry_available'))
        expected_leiden.append((modeled['node'].to_list(), seed_for(identity, 'leiden')))
        expected_abundance.append((selected['occurrences'].to_list(), seed_for(identity, 'exact_clonotype')))
        if modeled.height:
            expected_abundance.append(([int(modeled['occurrences'].sum())],
                                       seed_for(identity, 'leiden_community')))
    assert state['leiden'] == expected_leiden
    assert state['abundance'] == expected_abundance
    communities = pl.read_csv(tmp_path/'tables.communities.tsv', separator='\t')
    assert communities.sort('node').select('node', 'occurrences').rows() == state['nodes'].select('node', 'occurrences').rows()
    assert communities.filter(~pl.col('geometry_available'))['community'].is_null().all()
    assert communities['j'].null_count() == communities.height
    assert communities['occurrences'].sum() == 7*(1+replicates)
    summary = pl.read_csv(tmp_path/'tables.coverage.tsv', separator='\t')
    assert summary.filter(pl.col('coverage_status') == 'no_valid_occurrences')['group_id'].to_list() == [1, 3]
    curves = pl.read_csv(tmp_path/'tables.rarefaction.tsv', separator='\t')
    assert curves.height == len(expected_abundance)
    assert curves.schema['m'] == pl.Int64 and curves.schema['sample_coverage'] == pl.Float64
    manifest = json.loads((tmp_path/'tables.manifest.json').read_text())
    assert manifest['membership_budget_bytes'] == 2*1024**3
    assert manifest['actual_membership_rows'] == manifest['estimated_membership_rows'] == 7*(1+replicates)
    assert manifest['estimated_membership_bytes'] == 7*(1+replicates)*1024


def test_unsorted_graph_groups_fail_before_reduction(tmp_path, coverage_fixture):
    coverage, options, inputs, state = coverage_fixture
    state['unsorted'] = True
    with pytest.raises(ValueError, match='sorted by group_id'):
        coverage._main(options, inputs=inputs, scratch=tmp_path/'scratch')
    assert not state['leiden'] and not state['abundance']


def test_membership_budget_rejects_before_mock_or_graph_allocation(tmp_path, coverage_fixture, monkeypatch):
    coverage, options, inputs, state = coverage_fixture
    from vdjmatch.evalue import mock_reference
    monkeypatch.setattr(mock_reference, 'iter_mock_references',
        lambda *args, **kwargs: pytest.fail('mock allocation before budget guard'))
    options.mock_replicates = 200
    options.control = str(tmp_path/'control')
    options.membership_budget_bytes = 7*201*1024-1
    with pytest.raises(ValueError, match='membership estimate .* explicit budget'):
        coverage._main(options, inputs=inputs, scratch=tmp_path/'scratch')
    assert state['graph_calls'] == 0
    for budget in (-1,):
        options.membership_budget_bytes = budget
        with pytest.raises(ValueError, match='membership budget'):
            coverage.validate_options(options)
    for budget in (0, 8*1024**3+1):
        options.membership_budget_bytes = budget
        coverage.validate_options(options)


def test_zero_membership_budget_runs_without_the_estimate_abort(tmp_path, coverage_fixture):
    coverage, options, inputs, state = coverage_fixture
    options.membership_budget_bytes = 0
    options.mock_replicates = 2
    options.control = str(tmp_path/'control')
    assert coverage._main(options, inputs=inputs, scratch=tmp_path/'scratch') == 0
    manifest = json.loads((tmp_path/'tables.manifest.json').read_text())
    assert manifest['membership_budget_bytes'] == 0
    assert manifest['actual_membership_rows'] == 21 and state['graph_calls'] == 1
