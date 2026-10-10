"""Distance-local rank evidence under an explicit exchangeable-label model.

These are conditional model probabilities, not calibrated biological recognition
probabilities. Reference/control ascertainment and observation units must support
exchangeability before interpreting the null bound operationally.
"""
from __future__ import annotations

import numpy as np
from scipy.special import digamma
from scipy.stats import hypergeom


def _harmonic_interval(lower, upper):
    """Mean reciprocal integer rank over inclusive [lower, upper]."""
    lower, upper = np.broadcast_arrays(lower, upper)
    width = upper - lower + 1
    value = (digamma(upper + 1) - digamma(lower)) / width
    # Avoid cancellation for large, closely spaced integer ranks.
    small = width <= 8
    if np.any(small):
        exact = np.zeros(np.count_nonzero(small))
        lo, size = lower[small], width[small]
        for offset in range(8):
            exact += np.where(offset < size, 1 / (lo + offset), 0)
        value[small] = exact / size
    return value


def local_rank_evidence(reference_counts, control_counts, reference_population,
                        control_population, *, prior=.5):
    """Reduce batched cumulative distance balls to finite-control rank evidence.

    Rows are queries; columns are identical increasing, fixed distance cutoffs.
    Distances between cutoffs are tied in the same bin on both sides. References
    beyond the last cutoff contribute their null-conditional mean rank evidence.
    No target density, minimum target size, random tie breaking or fitted scale is
    needed. Empty populations return neutral evidence and the declared prior;
    callers must retain an explicit unavailable-population disposition.

    ``rank_bayes_factor`` is the mean normalized reciprocal background rank.
    ``p_rank_bound=min(1,1/BF)`` is a conservative bound, not an exact tail.
    ``posterior_rank_signal`` is the rank-tilted model posterior, not Prob(TP).
    The model conditions on fixed pooled distances and exchangeable labels;
    generated-draw/unique-reference differences do not establish that assumption.
    """
    if not np.isfinite(prior) or not 0 < prior < 1:
        raise ValueError('rank signal prior must be strictly between zero and one')
    n, m = np.asarray(reference_counts, dtype=float), np.asarray(control_counts, dtype=float)
    N, M = np.asarray(reference_population, dtype=float), np.asarray(control_population, dtype=float)
    if n.ndim != 2 or n.shape != m.shape or n.shape[1] == 0 or N.shape != (len(n),) or M.shape != N.shape:
        raise ValueError('aligned cumulative count matrices and populations required')
    for x in (n, m, N, M):
        if not np.isfinite(x).all() or np.any(x < 0) or np.any(x != np.floor(x)) or np.any(x > 2**53-1):
            raise ValueError('finite nonnegative exactly representable integer counts required')
    if np.any(np.diff(n, axis=1) < 0) or np.any(np.diff(m, axis=1) < 0) or np.any(n > N[:, None]) or np.any(m > M[:, None]):
        raise ValueError('cumulative counts must increase and not exceed populations')
    if np.any(N + M > 2**53-1):
        raise ValueError('combined population exceeds exactly representable integer range')
    lower = np.concatenate((np.zeros((len(m), 1)), m[:, :-1]), axis=1)
    target_shells = np.diff(n, axis=1, prepend=0)
    near = np.sum(target_shells * _harmonic_interval(lower + 1, m + 1), axis=1)
    far = (N-n[:, -1]) * _harmonic_interval(m[:, -1] + 1, M + 1)
    harmonic = digamma(M + 2) + np.euler_gamma
    bf = np.ones(len(n))
    usable = (N > 0) & (M > 0)
    bf[usable] = (M[usable]+1) * (near[usable]+far[usable]) / (N[usable]*harmonic[usable])
    bf[(n[:, -1] == 0) & (m[:, -1] == 0)] = 1
    if not np.isfinite(bf).all() or np.any(bf <= 0):
        raise ArithmeticError('nonfinite or nonpositive rank evidence')
    nearest_p = np.ones(len(n))
    hit = usable & (n[:, -1] > 0)
    if np.any(hit):
        first = np.argmax(n[hit] > 0, axis=1)
        before = m[hit][np.arange(np.count_nonzero(hit)), first]
        # Exact first-target rank tail, conservative with upper-bin tie ranks.
        # This tests nearest-hit evidence, not the complete averaged statistic.
        nearest_p[hit] = hypergeom.sf(0, N[hit]+M[hit], N[hit], before+1)
    return {'rank_bayes_factor': bf, 'p_rank_bound': np.minimum(1, 1/bf),
            'p_nearest_rank': nearest_p,
            'posterior_rank_signal': prior*bf/(1-prior+prior*bf)}
