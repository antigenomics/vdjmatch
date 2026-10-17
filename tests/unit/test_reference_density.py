"""Independent short-sequence hand fixture for one-batch kth-neighbour export."""
import argparse
import json

import polars as pl
import pytest
from seqtree import SubstitutionMatrix, gapblock

from vdjmatch.aggregate import reference_coverage as coverage
from vdjmatch.cli import reference_density as cli


def similarity(q, r, matrix, go, ge):
    """Independent affine local DP, including separate insertion/deletion states."""
    h = [[0]*(len(r)+1) for _ in range(len(q)+1)]
    e = [[-10**9]*(len(r)+1) for _ in range(len(q)+1)]
    f = [[-10**9]*(len(r)+1) for _ in range(len(q)+1)]
    for i, a in enumerate(q, 1):
        for j, b in enumerate(r, 1):
            e[i][j] = max(h[i-1][j]-go, e[i-1][j]-ge)
            f[i][j] = max(h[i][j-1]-go, f[i][j-1]-ge)
            h[i][j] = max(0, h[i-1][j-1]+matrix.similarity(a, b), e[i][j], f[i][j])
    return max(map(max, h))


def test_full_key_mhc_k_self_gaps_nulls_and_one_native_batch(tmp_path, monkeypatch):
    from vdjmatch.match.geometry import kneedle_radius
    matrix = SubstitutionMatrix.blosum62()
    def geometry(frame, **kwargs):
        groups = frame['v'].replace_strict({'TRBV1': 0, 'TRBV2': 1}, default=None).cast(pl.Int64)
        return frame.with_columns(groups.alias('_vgroup')), [[0, 3], [3, 0]], matrix
    monkeypatch.setattr(coverage, '_sw_geometry', geometry)
    native = gapblock.grouped_count_batch
    calls = []
    def count(*args, **kwargs):
        calls.append(kwargs)
        return native(*args, **kwargs)
    monkeypatch.setattr(gapblock, 'grouped_count_batch', count)
    rows = []
    for seq, v, j, mhc in [
        ('CAAAAF', 'TRBV1', 'TRBJ1', 'HLA-A*01:01'),
        ('CAAAAF', 'TRBV1', 'TRBJ1', 'HLA-A*01:01'),  # occurrence duplicate, not a second neighbor
        ('CAAAAF', 'TRBV1', 'TRBJ2', 'HLA-A*01:01'),  # zero-distance distinct full clonotype
        ('CAAAAF', 'TRBV2', 'TRBJ1', 'HLA-A*01:01'),
        ('CAARAF', 'TRBV1', 'TRBJ1', 'HLA-A*01:01'),
        ('CAAARAF', 'TRBV1', 'TRBJ1', 'HLA-A*01:01'),  # true length-changing SW case
        ('CAAAAF', 'TRBV1', 'TRBJ1', 'HLA-A*02:01'),  # separate MHC singleton
        ('CAGGGF', 'TRBV999', 'TRBJ1', 'HLA-A*03:01'),  # unmodeled loop geometry
    ]:
        rows.append(dict(cdr3=seq, v=v, j=j, mhc_a=mhc, mhc_b='B2M', mhc_class='MHCI',
            epitope='DEMO', gene='TRB', species='HomoSapiens'))
    source = tmp_path/'reference.tsv'
    pl.DataFrame(rows).write_csv(source, separator='\t')
    parser = argparse.ArgumentParser()
    cli.register(parser.add_subparsers())
    args = parser.parse_args(['reference-density', '--vdjdb', str(source), '--locus', 'TRB',
        '--epitope', 'DEMO', '--gap-open', '3', '--gap-extend', '1', '--cdr-loops', 'cdr12plus',
        '--output-prefix', str(tmp_path/'out')])
    assert args.func(args) == 0
    assert len(calls) == 1
    assert calls[0]['alignment_mode'] == 'local' and calls[0]['gap_open'] == 3
    assert calls[0]['gap_extend'] == 1 and not calls[0].get('exclude_exact', False)
    curves = pl.read_csv(tmp_path/'out.neighbor-distances.tsv', separator='\t')
    knees = pl.read_csv(tmp_path/'out.neighbor-knees.tsv', separator='\t')
    manifest = json.loads((tmp_path/'out.manifest.json').read_text())
    assert manifest['neighbor_orders'] == [1, 2, 5, 10] and manifest['search_batches'] == 1
    keys = list(dict.fromkeys((r['mhc_a'], r['cdr3'], r['v'], r['j']) for r in rows))
    assert curves.height == len(keys)*4
    group = {'TRBV1': 0, 'TRBV2': 1}
    def distance(a, b, go=3):
        return max(0, sum(matrix.similarity(x, x) for x in a[1]+b[1])-
            2*similarity(a[1], b[1], matrix, go, 1)) + 3*(group[a[2]] != group[b[2]])
    assert distance(keys[0], keys[4], 3) != distance(keys[0], keys[4], 11)
    expected_by_k = {k: [] for k in [1, 2, 5, 10]}
    for row in curves.to_dicts():
        key = (row['mhc_a'], row['cdr3'], row['v'], row['j'])
        candidates = sorted(distance(key, other) for other in keys
            if key != other and key[0] == other[0] and key[2] in group and other[2] in group)
        k = row['neighbor_order']
        expected = candidates[k-1] if len(candidates) >= k else None
        assert row['neighbor_distance'] == expected
        if expected is not None: expected_by_k[k].append(expected)
        if row['mhc_a'] == 'HLA-A*02:01':
            assert row['distance_status'] == 'insufficient_neighbors'
        if row['mhc_a'] == 'HLA-A*03:01':
            assert row['distance_status'] == 'geometry_unavailable'
    for row in knees.to_dicts():
        k = row['neighbor_order']
        assert row['distances_count'] == len(expected_by_k[k])
        assert row['epsilon_knee'] == kneedle_radius(expected_by_k[k])
        assert row['reference_clonotypes'] == 7 and row['modeled_clonotypes'] == 6
        assert row['unmodeled_clonotypes'] == 1
    assert knees.filter(pl.col('neighbor_order').is_in([5, 10]))['epsilon_knee'].null_count() == 2
    # The original graph wrapper keeps the same radius, node order, columns and failure behavior.
    nodes = curves.filter(pl.col('neighbor_order') == 1).select('group_id', 'cdr3', 'v', 'j')
    options = dict(species='human', locus='TRB', cdr_loops=('cdr1', 'cdr2', 'cdr25'), gap_open=3, gap_extend=1)
    expected_radius = kneedle_radius(expected_by_k[1])
    if expected_radius is None:
        with pytest.raises(ValueError, match='reference nearest-distance curve has no confirmed Kneedle boundary; supply --sw-radius'):
            coverage.sw_knee(nodes, **options)
    else:
        radius, curve = coverage.sw_knee(nodes, **options)
        assert radius == int(expected_radius)
        expected_curve = curves.filter(pl.col('neighbor_order') == 1).select(
            'node', 'group_id', 'cdr3', 'v', 'j', pl.col('neighbor_distance').alias('nearest_distance'))
        assert curve.columns == ['node', 'group_id', 'cdr3', 'v', 'j', 'nearest_distance']
        assert curve.to_dicts() == expected_curve.to_dicts()
    before = len(calls)
    unresolved = nodes.filter(pl.col('v') == 'TRBV999')
    absent, summary = coverage.sw_neighbor_knees(unresolved, neighbor_orders=(1, 2), **options)
    assert len(calls) == before and absent['neighbor_distance'].null_count() == 2
    assert summary['knee_status'].to_list() == ['no_modeled_vertices']*2
    for bad in [(), (0,), (True,), (1.5,)]:
        with pytest.raises(ValueError, match='positive integers'):
            coverage.sw_neighbor_knees(nodes, neighbor_orders=bad, **options)

    # A second native chunk must exclude only its own global search index.
    many = pl.DataFrame({'group_id': [0]*17+[1, 2], 'cdr3': ['CAAAAF']*19,
        'v': ['TRBV1']*18+['TRBV999'], 'j': [f'TRBJ{i}' for i in range(19)]})
    before = len(calls)
    distances, _ = coverage.sw_neighbor_knees(many, neighbor_orders=(1, 2, 5, 10), **options)
    assert len(calls) == before+1 and calls[-1]['return_hits']
    assert 'query_reference_ids' not in calls[-1]
    assert distances.filter(pl.col('group_id') == 0)['neighbor_distance'].to_list() == [0]*68
    assert distances.filter(pl.col('group_id') == 1)['distance_status'].to_list() == ['insufficient_neighbors']*4
    assert distances.filter(pl.col('group_id') == 2)['distance_status'].to_list() == ['geometry_unavailable']*4
