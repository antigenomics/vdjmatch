"""Frozen evidence on an aligned query-by-target-key distance rectangle."""

import numpy as np

from .competitive import reduced_hits


def fixed_target_scores(native, correction, exact, background, source_log2, rules, chains):
    """Return native and vX scores while masking exact hits only."""
    if not native.shape[1]:
        return tuple(np.full(len(native), -np.inf) for _ in range(2))
    eligible = np.where(exact, np.inf, native)
    native_score = -eligible.min(axis=1)
    supported = np.isfinite(native_score)
    target = eligible[supported] + correction[supported] - rules["source_bonus"] * source_log2
    width = min(2, target.shape[1])
    hits = np.partition(target, width - 1, axis=1)[:, :width]
    first = hits[:, 0]
    second = first if width == 1 else np.where(np.isfinite(hits[:, 1]), hits[:, 1], first)
    candidate = np.full(len(native), -np.inf)
    candidate[supported] = reduced_hits(first, second, background[supported], rules, chains)
    return native_score, candidate
