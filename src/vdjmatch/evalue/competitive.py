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
