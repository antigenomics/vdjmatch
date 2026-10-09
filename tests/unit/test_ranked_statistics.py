"""Deterministic finite-count oracles for the global fixed-K test."""

import pytest
from itertools import combinations
from math import comb
from scipy.stats import beta, binom

from vdjmatch.evalue.ranked import order_statistic, paired_order_statistic
from vdjmatch.evalue.paired import _joint_punctured_results


def test_fixed_k_exact_stopping_count_tail():
    k, n, m = 5, 20, 50
    for nc in (0, 1, 10, m):
        result = order_statistic(k, n, nc, m)
        prefix = k + nc
        expected = sum(comb(n, i) * comb(m, prefix-i)
                       for i in range(k, min(n, prefix)+1) if 0 <= prefix-i <= m) / comb(n+m, prefix)
        assert result["p_global"] == pytest.approx(expected)
        assert result["n_control"] == nc and result["control_size"] == m
        assert result["E_raw"] == n * nc / m
        assert result["p_upper"] is result["E"] is result["finite_control_delta"] is None


def test_ties_and_symmetric_infinity_exclusion_are_conservative():
    n, m = 3, 4
    for scores in ([0, 0, 0, 1, 1, 2, 2], [0, 0, 1, 1, 2, float("inf"), float("inf")]):
        for k in range(1, n+1):
            probabilities = []
            for positions in combinations(range(n+m), n):
                target = sorted(scores[i] for i in positions)
                threshold = target[k-1]
                c = sum(score <= threshold for i, score in enumerate(scores) if i not in positions)
                probabilities.append(1.0 if threshold == float("inf") else order_statistic(k, n, c, m)["p_global"])
            for alpha in set(probabilities):
                assert sum(p <= alpha for p in probabilities) / len(probabilities) <= alpha + 1e-12


def test_full_reference_gate_is_resolvable_without_more_controls():
    result = order_statistic(10, 77681, 0, 250000)
    assert result["p_global"] == pytest.approx(comb(77681, 10) / comb(327681, 10))
    assert result["p_global"] < .001


def test_paired_fixed_k_uses_existing_punctured_bound_and_swaps():
    args = (5, 100, 12, 17, 1, 1, 1000, 2000)
    result = paired_order_statistic(*args, exclude_exact=True)
    old = _joint_punctured_results([5], [12], [17], [100], 1000, 2000, [1], [1])
    assert result["p_global"] == old["p_enrichment"][0]
    swapped = paired_order_statistic(5, 100, 17, 12, 1, 1, 2000, 1000, True)
    assert swapped == pytest.approx(result)
    unpunctured = paired_order_statistic(*args)
    upper = beta.ppf(1 - 5e-7, 13, 988) * beta.ppf(1 - 5e-7, 18, 1983)
    assert unpunctured["p_upper"] == pytest.approx(upper)


@pytest.mark.parametrize("args", [(0, 10, 0, 10), (11, 10, 0, 10), (5, 10, 11, 10), (5, 10, 0, 0), (True, 10, 0, 10)])
def test_invalid_counts_fail(args):
    with pytest.raises(ValueError):
        order_statistic(*args)
