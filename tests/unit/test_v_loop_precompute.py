"""Exact run-local V-loop reuse; only synthetic junctions and shipped templates."""
import numpy as np
import pytest

from vdjmatch.match import geometry as td


@pytest.mark.parametrize('species', ['human', 'mouse'])
@pytest.mark.parametrize('alphabet,matrix', [('full', 'clipped'), ('imgt', 'clipped'), ('full', 'blosum62')])
def test_precomputed_loops_match_original_allele_subsets(species, alphabet, matrix, monkeypatch):
    table = td.precompute_v_loop_distances(species=species, alphabet=alphabet, matrix=matrix)
    costs = td._substitution_costs(alphabet, matrix)
    model = table['model']
    for locus in ['TRA', 'TRB']:
        alleles = sorted(a for a in model if model[a]['locus'] == locus)
        genes = [alleles[-1], alleles[0], alleles[0], alleles[len(alleles)//2]]
        old, old_ids = td._v_loop_groups(genes, model, costs)
        new, new_ids = td._v_loop_groups(genes, model, costs, table)
        np.testing.assert_array_equal(new, old)
        np.testing.assert_array_equal(new_ids, old_ids)
    assert all(not distances.flags.writeable for _, distances in table['tables'].values())
    monkeypatch.setattr(td, 'load_v_loops', lambda *args: pytest.fail('precomputed model must be reused'))
    genes = sorted(a for a in model if model[a]['locus'] == 'TRB')[:2]
    params = dict(species=species, alphabet=alphabet, matrix=matrix, v_loop_distances=table)
    td.distance_matrix(['CASSAAAAAAAF'], ['CASSACAAAAAF'], genes[:1], genes[1:], **params)


def test_precomputation_preserves_single_paired_counts_and_hard_v_matching():
    table = td.precompute_v_loop_distances()
    qs, rs = ['CASSAAAAAAAF', 'CASSACAAAAAF'], ['CASSADAAAAAF', 'CASSAAAAAAAF']
    model = table['model']
    beta = sorted(a for a in model if model[a]['locus'] == 'TRB')[:2]
    alpha = sorted(a for a in model if model[a]['locus'] == 'TRA')[:2]
    for same_v_gene in [False, True]:
        kwargs = dict(thresholds=[[0, 10000], [0, 10000]], same_v_gene=same_v_gene, linear_mass=True)
        original = td.total_distance_count_batch(qs, rs, beta, beta[::-1], **kwargs)
        reused = td.total_distance_count_batch(qs, rs, beta, beta[::-1], v_loop_distances=table, **kwargs)
        assert reused == original
        kwargs.pop('linear_mass')
        original = td.paired_total_distance_mass_batch(qs, qs, rs, rs, alpha, beta, alpha[::-1], beta[::-1], **kwargs)
        reused = td.paired_total_distance_mass_batch(qs, qs, rs, rs, alpha, beta, alpha[::-1], beta[::-1], v_loop_distances=table, **kwargs)
        assert reused == original
    np.testing.assert_array_equal(td.distance_matrix(qs, rs, beta, beta[::-1]),
                                  td.distance_matrix(qs, rs, beta, beta[::-1], v_loop_distances=table))


def test_precomputed_context_mismatch_fails():
    table = td.precompute_v_loop_distances()
    beta = sorted(a for a in table['model'] if table['model'][a]['locus'] == 'TRB')[0]
    with pytest.raises(ValueError, match='matching species'):
        td.distance_matrix(['CASSAAAAAAAF'], ['CASSAAAAAAAF'], [beta], [beta], species='mouse', v_loop_distances=table)
    with pytest.raises(ValueError, match='matching substitution costs'):
        td.distance_matrix(['CASSAAAAAAAF'], ['CASSAAAAAAAF'], [beta], [beta], matrix='blosum62', v_loop_distances=table)
