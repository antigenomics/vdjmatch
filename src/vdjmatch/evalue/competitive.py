"""Evidence reduction for the frozen vX research score."""

from collections import Counter

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
    values = np.asarray(distances)
    if values.ndim != 2 or values.shape[0] != values.shape[1] or not values.shape[0]:
        raise ValueError("expected a nonempty square same-target distance matrix")
    if values.dtype.kind not in "iuf" or np.isnan(values).any() or np.any(values < 0):
        raise ValueError("distances must be real and nonnegative, or positive infinity")
    group = values.astype(float, copy=True)
    if exact_hits is not None:
        mask = np.asarray(exact_hits)
        if mask.shape != group.shape or mask.dtype != np.bool_:
            raise ValueError("exact_hits must be a boolean mask aligned to distances")
        group[mask] = np.inf
    np.fill_diagonal(group, np.inf)
    width = max(1, (group.shape[0] - 1) * 124 // 1000)
    order = np.sort(group, axis=1)
    return assay_local_score(order, width, 7 / 3)


def _anchored_trimers(junction):
    core = junction[3:-2]
    width = min(4, len(core) - 2)
    return (
        tuple(("N", i, core[i:i + 3]) for i in range(width))
        + tuple(("C", i, core[len(core) - i - 3:len(core) - i]) for i in range(width))
    )


def anchored_motif_scores(reference_keys, junctions, targets, *,
                          reference_pool=None, sibling_references=None):
    """Return fixed-catalogue TRB motif enrichment in the requested target order.

    Each reference item is a (target, full junction) observation. Callers
    deduplicate their declared reference keys before passing the retained
    items. Distinct keys may carry the same junction. Targets with no
    reference observations are unsupported. The score uses no query labels.
    An explicit pool can use physical keys counted once. Optional sibling
    junctions come from compatible same-peptide targets, with shared
    reference keys removed by the caller.
    """
    targets = tuple(targets)
    if not targets or len(set(targets)) != len(targets):
        raise ValueError("expected distinct target labels")
    target_counts = {target: Counter() for target in targets}
    sizes = Counter()
    pooled = Counter()
    for target, junction in reference_keys:
        counts = target_counts[target]
        tokens = _anchored_trimers(junction)
        counts.update(tokens)
        pooled.update(tokens)
        sizes[target] += 1
    denominator = 2 + sizes.total()
    if reference_pool is not None:
        pooled.clear()
        denominator = 2
        for junction in reference_pool:
            pooled.update(_anchored_trimers(junction))
            denominator += 1
    queries = [_anchored_trimers(junction) for junction in junctions]
    scores = np.full((len(queries), len(targets)), -np.inf)
    for column, target in enumerate(targets):
        if not sizes[target]:
            continue
        counts = target_counts[target]
        siblings = tuple(sibling_references.get(target, ())) if sibling_references is not None else ()
        sibling_counts = Counter(token for junction in siblings for token in _anchored_trimers(junction))
        borrowing = len(siblings) / (len(siblings) + 2)
        for row, tokens in enumerate(queries):
            if not tokens:
                continue
            prior = np.asarray([(1 + pooled[token]) / denominator for token in tokens])
            if siblings:
                sibling_prior = np.asarray([(1 + sibling_counts[token]) / (len(siblings) + 2)
                                            for token in tokens])
                prior = (1 - borrowing) * prior + borrowing * sibling_prior
            observed = np.asarray([counts[token] for token in tokens])
            enrichment = np.log(
                (observed + 2 * prior) / ((sizes[target] + 2) * prior)
            )
            if reference_pool is None and sibling_references is None:
                scores[row, column] = enrichment.mean()
            else:
                scores[row, column] = np.add.accumulate(enrichment)[-1] / len(tokens)
    return scores


def paired_motif_ranking_scores(base, tra, trb, target_key_counts, *,
                               support_weight, conflict_weight, lower_tier_base=None):
    """Refine lower paired-chain ranks while keeping the frozen winning tier.

    Inputs share the query-by-target axes. Use a complete compatible target
    context here, then select any reporting labels from the resulting scores.
    The weights are explicit inputs so species settings can be verified.
    An optional nearest-hit score matrix supplies the lower-tier baseline.
    A missing nearest-hit score falls back to the original score.
    """
    base = np.asarray(base, dtype=float)
    tra, trb = np.asarray(tra, dtype=float), np.asarray(trb, dtype=float)
    counts = np.asarray(target_key_counts)
    if (base.ndim != 2 or not base.shape[1] or tra.shape != base.shape
            or trb.shape != base.shape or counts.shape != (base.shape[1],)):
        raise ValueError("expected aligned query-by-target scores and one count per target")
    if np.isnan(base).any() or np.isposinf(base).any():
        raise ValueError("base scores must be finite or negative infinity")
    lower = base
    if lower_tier_base is not None:
        alternative = np.asarray(lower_tier_base, dtype=float)
        if alternative.shape != base.shape or np.isnan(alternative).any() or np.isposinf(alternative).any():
            raise ValueError("lower-tier scores must be aligned and finite or negative infinity")
        lower = np.where(np.isfinite(alternative), alternative, base)
    if (not np.isfinite([support_weight, conflict_weight]).all()
            or min(support_weight, conflict_weight) < 0):
        raise ValueError("motif weights must be finite and nonnegative")
    motif = np.minimum(np.where(np.isfinite(tra), tra, 0), np.where(np.isfinite(trb), trb, 0))
    maximum = base.max(axis=1, keepdims=True)
    cap = np.nextafter(maximum, -np.inf)
    adjust = np.isfinite(base) & (base < maximum)
    motif = np.where(adjust & (counts > 2), motif, 0)
    with np.errstate(over="ignore"):
        correction = support_weight * np.maximum(motif, 0) + conflict_weight * np.minimum(motif, 0)
        result = np.where(adjust, np.minimum(lower + correction, cap), base)
    if np.any(np.isfinite(base) & ~np.isfinite(result)):
        raise ValueError("motif refinement exceeds the finite score range")
    return result
