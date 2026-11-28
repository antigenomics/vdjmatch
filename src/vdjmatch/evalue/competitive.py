"""Evidence reduction for the frozen vX research score."""

import numpy as np

MAX_MIDPOINT_GAP_PER_CHAIN = 24
ASSIGNMENT_SIZE_SLOPE = 0.002832


def reduced_hits(first, second, background, rules, chains):
    """Score ordered finite target hits with the configured competitive margin.

    Callers repeat the first hit for singletons and handle unsupported targets
    before reduction. An infinite background distance yields zero excess.
    """
    # Keep the midpoint before subtraction to preserve the frozen binary64 scores.
    midpoint = (first + second) / 2
    gap = np.clip(midpoint - first, 0, MAX_MIDPOINT_GAP_PER_CHAIN * chains)
    excess = np.maximum(first + chains * rules["margin_per_chain"] - background, 0)
    return -first - rules["margin_penalty"] * excess + rules["second_hit_gap_weight"] * gap


def assignment_score(scores, target_key_counts):
    """Shift evidence by fixed target size for the score-only assay comparison."""
    return scores + ASSIGNMENT_SIZE_SLOPE * np.log1p(target_key_counts)


def assay_local_score(ordered_distances, width, rank_power):
    """Score each assay case from ascending distances to its same-target peers.

    The caller masks self-hits and any other exclusions with positive infinity.
    A row without the specified number of finite neighbors has no score.
    """
    distances = np.asarray(ordered_distances, dtype=float)
    if distances.ndim != 2 or not isinstance(width, (int, np.integer)) or not 1 <= width <= distances.shape[1]:
        raise ValueError("width must select at least one column from a distance matrix")
    if not np.isfinite(rank_power) or rank_power < 0:
        raise ValueError("rank power must be finite and nonnegative")
    weights = np.arange(width, 0, -1, dtype=float) ** rank_power
    return -(distances[:, :width] @ weights) / weights.sum()


def paired_assay_local_score(distances, *, exact_hits=None):
    """Rank one annotated assay group from a square paired-chain distance matrix.

    The diagonal is always masked. An optional boolean mask excludes additional
    exact hits without changing the neighborhood width or the assay population.
    The width and rank power were selected on the TCRvdb panel.
    """
    group = np.array(distances, dtype=float, copy=True)
    if group.ndim != 2 or group.shape[0] != group.shape[1] or not group.shape[0]:
        raise ValueError("expected a nonempty square same-target distance matrix")
    if exact_hits is not None:
        mask = np.asarray(exact_hits)
        if mask.shape != group.shape or mask.dtype != np.bool_:
            raise ValueError("exact_hits must be a boolean mask aligned to distances")
        group[mask] = np.inf
    np.fill_diagonal(group, np.inf)
    width = max(1, (group.shape[0] - 1) * 124 // 1000)
    order = np.sort(group, axis=1)
    return assay_local_score(order, width, 7 / 3)
