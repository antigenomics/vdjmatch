"""Graph acquisition agrees with enumeration; degree bins are not taxon counts."""
import itertools

import numpy as np
import pytest

from vdjmatch.aggregate.graph_acquisition import graph_acquisition, degree_diversity, clonotype_graph_degrees


def test_graph_noise_ranks_triangle_support_and_keeps_connected_zero_ties(tmp_path):
    import polars as pl
    from vdjmatch.cli.__main__ import main

    nodes, edges, output = (tmp_path / name for name in ('nodes.tsv', 'edges.tsv', 'out.tsv'))
    nodes.write_text('node_id\tmetadata\n005\torphan\n001\ta\n002\tb\n003\tc\n004\tleaf\n')
    edges.write_text('a\tb\n001\t002\n001\t003\n002\t003\n003\t004\n002\t001\n004\t004\n')
    args = ['graph-noise', str(nodes), '--edges', str(edges), '--output', str(output)]
    assert main(args) == 0
    result = pl.read_csv(output, separator='\t', schema_overrides={'node_id': pl.String})
    assert result['node_id'].to_list() == ['005', '001', '002', '003', '004']
    assert result['degree'].to_list() == [0, 2, 2, 3, 1]
    assert result['triangle_degree'].to_list() == [0, 2, 2, 2, 0]
    assert result['noise_rank'].to_list() == [5, 1, 1, 1, 4]
    for budget, expected in [(0, [0.] * 5), (2, [0., 2/3, 2/3, 2/3, 0.]),
                             (4, [0., 1., 1., 1., 1.]), (5, [1.] * 5)]:
        retained = tmp_path / f'retained-{budget}.tsv'
        assert main(args[:-1] + [str(retained), '--retain', str(budget)]) == 0
        weighted = pl.read_csv(retained, separator='\t', schema_overrides={'node_id': pl.String})
        assert weighted.drop('retention_weight').equals(result)
        assert weighted['retention_weight'].to_list() == pytest.approx(expected)
    for budget in (-1, 6):
        with pytest.raises(SystemExit):
            main(args[:-1] + [str(tmp_path / 'invalid-budget.tsv'), '--retain', str(budget)])
    assert not (tmp_path / 'invalid-budget.tsv').exists()
    edges.write_text('a\tb\n')
    assert main(args[:-1] + [str(tmp_path / 'no-edges.tsv')]) == 0
    edges.write_text('a\tb\n001\t006\n')
    with pytest.raises(SystemExit):
        main(args[:-1] + [str(tmp_path / 'invalid.tsv')])


def test_graph_acquisition_matches_every_induced_subgraph():
    # Triangle, pendant vertex and orphan; includes degree one and degree zero.
    edges = {(0, 1), (0, 2), (1, 2), (2, 3)}
    degrees = np.bincount(np.array(sorted(edges)).ravel(), minlength=5)
    curves = graph_acquisition(degrees, np.arange(1, 6))
    for row in curves.to_dicts():
        m = row['sample_junctions']
        counts = []
        for sample in itertools.combinations(range(5), m):
            observed = [e for e in edges if set(e) <= set(sample)]
            counts.append((len(observed), len(set(itertools.chain.from_iterable(observed)))))
        assert row['expected_graph_edges'] == pytest.approx(np.mean([x[0] for x in counts]))
        assert row['expected_connected_junctions'] == pytest.approx(np.mean([x[1] for x in counts]))
        assert row['expected_orphan_junctions'] + row['expected_connected_junctions'] == pytest.approx(m)
    assert graph_acquisition([0], [1])['expected_connected_fraction'][0] == 0
    for degrees, sizes in [([-1], [1]), ([1], [1]), ([0], [0]), ([0], [2]), ([0.], [1])]:
        with pytest.raises(ValueError):
            graph_acquisition(degrees, sizes)


def test_similarity_hills_match_orphan_clique_and_two_family_limits():
    for degrees, expected in [([0]*4, 4), ([3]*4, 1), ([1]*4, 2)]:
        row = degree_diversity(degrees)
        for q in range(3):
            assert row[f'graph_hill{q}'] == pytest.approx(expected)
            assert row[f'graph_hill{q}_enrichment'] == pytest.approx(4/expected)
        assert row['degree_chao_status'] == 'degree_mass_proxy'
    row = degree_diversity([1, 2, 1, 0])
    assert row['degree_mass'] == 8
    assert row['degree_singletons'] == 1 and row['degree_doubletons'] == 2
    assert row['degree_chao_coverage'] == pytest.approx(1-1/8*7/11)
    assert row['graph_hill0'] >= row['graph_hill1'] >= row['graph_hill2']


def test_vj_distinct_clonotypes_keep_exact_junction_graph_links(tmp_path):
    import polars as pl
    from vdjmatch.cli.__main__ import main
    degrees = clonotype_graph_degrees([2, 1, 1], [[0, 1]])
    assert degrees.tolist() == [2, 2, 2, 0]
    reference, output = tmp_path/'reference.tsv', tmp_path/'profile'
    records = pl.DataFrame({'cdr3': ['CASSF', 'CASSF', 'CASAF', 'CWWWF'],
        'v': ['TRBV19', 'TRBV20-1', 'TRBV19', 'TRBV19'], 'j': ['TRBJ1-1']*4,
        'epitope': ['A']*4, 'species': ['HomoSapiens']*4, 'gene': ['TRB']*4})
    records = pl.concat([records, records.head(1).with_columns(pl.lit('TRBV19*02').alias('v'),
        pl.lit('TRBJ1-1*02').alias('j'))])
    pl.concat([records, records.reverse().with_columns(pl.lit('B').alias('epitope'))]).write_csv(reference, separator='\t')
    assert main(['reference-structure', '--vdjdb', str(reference), '--locus', 'TRB',
        '--threads', '1', '--output-prefix', str(output)]) == 0
    summary = pl.read_csv(str(output)+'.structures.tsv', separator='\t').row(0, named=True)
    assert summary['hamming_junctions'] == 3 and summary['hamming_graph_edges'] == 1
    assert summary['degree_clonotypes'] == 4 and summary['degree_clonotype_edges'] == 3
    histogram = pl.read_csv(str(output)+'.clonotype-degrees.tsv', separator='\t').filter(pl.col('epitope') == 'A')
    assert histogram['inclusive_degree'].to_list() == [1, 3]
    assert histogram['clonotype_count'].to_list() == [1, 3]
    curve = pl.read_csv(str(output)+'.graph-acquisition.tsv', separator='\t').row(-1, named=True)
    assert curve['sample_clonotypes'] == curve['graph_clonotypes'] == 4
    assert curve['expected_graph_edges'] == curve['expected_connected_clonotypes'] == 3
