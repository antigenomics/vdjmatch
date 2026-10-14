"""Finite exchangeability oracle for ties, censoring and sparse rank evidence."""
from itertools import combinations, combinations_with_replacement

import numpy as np
import pytest

from vdjmatch.evalue.local_rank import local_rank_evidence


def test_all_small_binnings_preserve_null_with_data_dependent_far_branch():
    from vdjmatch.evalue.local_rank import nearest_conditioned_rank

    # Exhaust all tied pooled geometries, label assignments and population sizes.
    # Bins 2 and 3 are censored: the branch is selected by the observed labels.
    strict = False
    for distances in combinations_with_replacement(range(4), 5):
        for n in range(1, 5):
            refs, controls, nearest = [], [], []
            for target in combinations(range(5), n):
                background = set(range(5)).difference(target)
                refs.append([sum(distances[i] <= t for i in target) for t in (0, 1)])
                controls.append([sum(distances[i] <= t for i in background) for t in (0, 1)])
                first = min(distances[i] for i in target)
                nearest.append(sum(distances[i] <= first for i in background))
            N, M = np.full(len(refs), n), np.full(len(refs), 5-n)
            original = local_rank_evidence(refs, controls, N, M)['rank_bayes_factor']
            refined = nearest_conditioned_rank(original, np.array(refs)[:, -1], N, M, nearest)
            assert np.mean(original) == pytest.approx(1, abs=2e-14)
            assert np.mean(refined) <= 1 + 2e-14
            strict |= np.mean(refined) < 1 - 2e-14
    assert strict  # Upper ties can destroy equality, but not the e-value bound.


@pytest.mark.parametrize('distances', [(0, 1, 2, 3, 4), (0, 0, 1, 1, 3), (3, 3, 3, 3, 3)])
@pytest.mark.parametrize('n', [1, 2, 3])
def test_finite_label_null_and_bound(distances, n):
    scores, pvalues, nearest_pvalues, conditional_scores = [], [], [], []
    for target in combinations(range(5), n):
        ref = np.array([sum(distances[i] <= t for i in target) for t in range(3)])
        control = np.array([sum(distances[i] <= t for i in range(5) if i not in target) for t in range(3)])
        result = local_rank_evidence(ref[None, :], control[None, :], [n], [5-n], prior=.2)
        score, pvalue, posterior = (result[k][0] for k in ('rank_bayes_factor', 'p_rank_bound', 'posterior_rank_signal'))
        assert posterior == pytest.approx(.2*score/(.8+.2*score))
        scores.append(score); pvalues.append(pvalue)
        nearest_pvalues.append(result['p_nearest_rank'][0])
        from vdjmatch.evalue.local_rank import nearest_conditioned_rank
        first_distance = min(distances[i] for i in target)
        k = sum(distances[i] <= first_distance for i in range(5) if i not in target)
        conditional_scores.append(nearest_conditioned_rank([score], [ref[-1]], [n], [5-n], [k])[0])
    assert np.mean(scores) == pytest.approx(1, abs=2e-14)
    assert np.mean(conditional_scores) <= 1 + 2e-14
    for alpha in set(pvalues):
        assert np.mean(np.array(pvalues) <= alpha) <= alpha + 1e-14
    for alpha in set(nearest_pvalues):
        assert np.mean(np.array(nearest_pvalues) <= alpha) <= alpha + 1e-14


def test_empty_and_no_near_evidence_are_uninformative():
    result = local_rank_evidence(np.zeros((3, 2), dtype=int), np.zeros((3, 2), dtype=int), [1, 0, 4], [20, 20, 0], prior=.3)
    np.testing.assert_allclose(result['rank_bayes_factor'], 1)
    np.testing.assert_allclose(result['p_rank_bound'], 1)
    np.testing.assert_allclose(result['posterior_rank_signal'], .3)


def test_invalid_counts_and_prior_rejected():
    with pytest.raises(ValueError):
        local_rank_evidence([[2, 1]], [[0, 0]], [3], [4])
    with pytest.raises(ValueError):
        local_rank_evidence([[0]], [[0]], [1], [1], prior=1)


def test_large_population_singleton_near_evidence_is_finite():
    result = local_rank_evidence([[1, 1]], [[0, 17_954_000]], [1], [17_954_000])
    from scipy.special import digamma
    expected = 17_954_001 / (digamma(17_954_002) + np.euler_gamma)
    assert result['rank_bayes_factor'][0] == pytest.approx(expected)
    assert result['p_nearest_rank'][0] == pytest.approx(1/17_954_001)


def test_nearest_conditioned_far_mean_matches_finite_oracle():
    from vdjmatch.evalue.local_rank import nearest_conditioned_rank
    from collections import defaultdict
    groups = defaultdict(list)
    for target in combinations(range(7), 3):
        ranks = [1 + sum(i not in target for i in range(j)) for j in target]
        groups[ranks[0]-1].append(np.mean([1/r for r in ranks]))
    harmonic = sum(1/r for r in range(1, 6))
    for first, values in groups.items():
        result = nearest_conditioned_rank([1.], [0], [3], [4], [first])
        assert result[0] == pytest.approx(5*np.mean(values)/harmonic)
    # Existing multi-reference evidence remains exact whenever radius support exists.
    np.testing.assert_array_equal(nearest_conditioned_rank([2., 3.], [1, 2], [3, 4], [5, 6], [0, 0]), [2., 3.])
    with pytest.raises(ValueError):
        nearest_conditioned_rank([1.], [0], [1], [3], [4])
