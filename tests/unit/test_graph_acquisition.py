"""Graph acquisition agrees with enumeration; degree bins are not taxon counts."""
import itertools

import numpy as np
import pytest

from vdjmatch.aggregate.graph_acquisition import graph_acquisition, degree_diversity


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
