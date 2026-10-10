"""Finite exchangeability oracle for ties, censoring and sparse rank evidence."""
from itertools import combinations

import numpy as np
import pytest

from vdjmatch.evalue.local_rank import local_rank_evidence


@pytest.mark.parametrize('distances', [(0, 1, 2, 3, 4), (0, 0, 1, 1, 3), (3, 3, 3, 3, 3)])
@pytest.mark.parametrize('n', [1, 2, 3])
def test_finite_label_null_and_bound(distances, n):
    scores, pvalues, nearest_pvalues = [], [], []
    for target in combinations(range(5), n):
        ref = np.array([sum(distances[i] <= t for i in target) for t in range(3)])
        control = np.array([sum(distances[i] <= t for i in range(5) if i not in target) for t in range(3)])
        result = local_rank_evidence(ref[None, :], control[None, :], [n], [5-n], prior=.2)
        score, pvalue, posterior = (result[k][0] for k in ('rank_bayes_factor', 'p_rank_bound', 'posterior_rank_signal'))
        assert posterior == pytest.approx(.2*score/(.8+.2*score))
        scores.append(score); pvalues.append(pvalue)
        nearest_pvalues.append(result['p_nearest_rank'][0])
    assert np.mean(scores) == pytest.approx(1, abs=2e-14)
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
