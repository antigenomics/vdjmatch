"""Deterministic finite-count oracles for the global fixed-K test."""

import pytest
from scipy.stats import beta, binom

from vdjmatch.evalue.ranked import order_statistic, paired_order_statistic
from vdjmatch.evalue.paired import _joint_punctured_results


def test_fixed_k_confidence_bound_and_discrete_null():
    k, n, m, delta = 5, 20, 50, 1e-6
    for nc in (0, 1, 10, m):
        upper = 1.0 if nc == m else beta.ppf(1 - delta, nc + 1, m - nc)
        result = order_statistic(k, n, nc, m)
        assert result["p_upper"] == pytest.approx(upper)
        assert result["p_global"] == pytest.approx(min(1, delta + binom.sf(k - 1, n, upper)))
    # A two-score IID null: T_K=0 iff at least K targets score zero.
    # Exact enumeration verifies the combined target/control rejection mass.
    p_zero, alpha = 0.05, 0.2
    reject = sum(binom.pmf(c, m, p_zero)
                 for c in range(m + 1) if order_statistic(k, n, c, m)["p_global"] <= alpha)
    assert binom.sf(k - 1, n, p_zero) * reject <= alpha


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
