"""Evidence reduction for the frozen vX research score."""

from collections import Counter

import numpy as np

# The midpoint cap is two maximum-cost junction edits per chain.
# One edit costs four clipped units at junction weight three.
MAX_MIDPOINT_GAP_PER_CHAIN = 2 * (3 * 4)
ASSIGNMENT_SIZE_SLOPE = 1 / 400


def fr3_sw_distance_correction(context, locus):
    """Return the research FR3 correction in native distance units.

    This replaces 1/512 of the frozen FR3 cost with Smith–Waterman
    dissimilarity. Add it only where the frozen geometry is available
    for every observed chain. The allele order matches the context.
    """
    from seqtree import SubstitutionMatrix, gapblock
    from ..match.geometry import _v_loop_groups

    alleles = context["tables"][locus][0]
    model = context["model"]
    region = "cdr25"  # Historical key for the selected FR3 segment.
    only_fr3 = {gene: dict(entry, cdr1="", cdr2="") for gene, entry in model.items()}
    old, indices = _v_loop_groups(alleles.tolist(), only_fr3, context["costs"])
    old = old[np.ix_(indices, indices)]
    segments = [model[gene][region].replace(".", "").replace("-", "") for gene in alleles]
    keep = np.flatnonzero([bool(segment) for segment in segments])
    aligned = [segments[i] for i in keep]
    n = len(aligned)
    correction = np.zeros(old.shape, dtype=float)
    if not n:
        return alleles, correction

    native = gapblock.grouped_count_batch(
        aligned, aligned, [10000], reference_target_offsets=list(range(n + 1)),
        reference_target_ids=[0] * n, target_peptide_ids=[0],
        reference_exclusive_peptide_ids=[0] * n, query_gene_ids=[0] * n,
        reference_gene_ids=[0] * n, query_target_ids=[0] * n,
        return_hits=True, chunk_size=32, threads=1, matrix=SubstitutionMatrix.blosum62(),
        gap_open=2, gap_extend=2, alignment_mode="local",
    )
    distance = np.full((n, n), np.nan)
    for chunk in native:
        q, r, numerator = np.asarray(chunk.hits).T
        distance[q, r] = numerator / 2
    if not np.isfinite(distance).all():
        raise ValueError("FR3 alignment did not cover the resolved allele pairs")
    pair = np.ix_(keep, keep)
    frozen_fr3_weight = 3 / 2 + 1 / 256
    correction[pair] = (distance - old[pair]) * frozen_fr3_weight / 512
    return alleles, correction


def reduced_hits(first, second, background, rules, chains):
    """Score the nearest hit, adjusted for competition and second-hit support.

    The margin is in distance units per observed chain. The margin penalty
    and second-hit gap weight are dimensionless. Callers repeat the first
    distance when only one finite target hit remains. Infinite background
    distance removes the competition penalty.
    """
    # Keep the midpoint before subtraction to preserve the frozen binary64 scores.
    midpoint = (first + second) / 2
    midpoint_gap = np.clip(midpoint - first, 0, MAX_MIDPOINT_GAP_PER_CHAIN * chains)
    margin_shortfall = np.maximum(first + chains * rules["margin_per_chain"] - background, 0)
    return -first - rules["margin_penalty"] * margin_shortfall + rules["second_hit_gap_weight"] * midpoint_gap


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
    width = max(1, 31 * (group.shape[0] - 1) // 250)
    order = np.sort(group, axis=1)
    return assay_local_score(order, width, 7 / 3)


def _anchored_tokens(junction, motif_length):
    core = junction[3:-2]
    width = min(4, len(core) - motif_length + 1)
    return (
        tuple(("N", i, core[i:i + motif_length]) for i in range(width))
        + tuple(("C", i, core[len(core) - i - motif_length:len(core) - i]) for i in range(width))
    )


def anchored_motif_scores(reference_keys, junctions, targets, *,
                          reference_pool=None, sibling_references=None, motif_length=3):
    """Return fixed-catalogue TRB motif enrichment in the requested target order.

    Each reference item is a (target, full junction) observation. Callers
    deduplicate their declared reference keys before passing the retained
    items. Distinct keys may carry the same junction. Targets with no
    reference observations are unsupported. The score uses no query labels.
    An explicit pool can use physical keys counted once. Optional sibling
    junctions come from compatible same-peptide targets, with shared
    reference keys removed by the caller.
    Global and sibling frequencies start with one occurrence and one
    nonoccurrence. Both motif lengths borrow two prior observations
    when estimating the target frequency.
    """
    if not isinstance(motif_length, (int, np.integer)) or motif_length not in (3, 4):
        raise ValueError("motif_length must be 3 or 4")
    targets = tuple(targets)
    if not targets or len(set(targets)) != len(targets):
        raise ValueError("expected distinct target labels")
    target_counts = {target: Counter() for target in targets}
    sizes = Counter()
    pooled = Counter()
    prior_mass = 2
    for target, junction in reference_keys:
        counts = target_counts[target]
        tokens = _anchored_tokens(junction, motif_length)
        counts.update(tokens)
        pooled.update(tokens)
        sizes[target] += 1
    denominator = prior_mass + sizes.total()
    if reference_pool is not None:
        pooled.clear()
        denominator = prior_mass
        for junction in reference_pool:
            pooled.update(_anchored_tokens(junction, motif_length))
            denominator += 1
    queries = [_anchored_tokens(junction, motif_length) for junction in junctions]
    scores = np.full((len(queries), len(targets)), -np.inf)
    for column, target in enumerate(targets):
        if not sizes[target]:
            continue
        counts = target_counts[target]
        siblings = tuple(sibling_references.get(target, ())) if sibling_references is not None else ()
        sibling_counts = Counter(token for junction in siblings for token in _anchored_tokens(junction, motif_length))
        borrowing = len(siblings) / (len(siblings) + prior_mass)
        for row, tokens in enumerate(queries):
            if not tokens:
                continue
            prior = np.asarray([(1 + pooled[token]) / denominator for token in tokens])
            if siblings:
                sibling_prior = np.asarray([(1 + sibling_counts[token]) / (len(siblings) + prior_mass)
                                            for token in tokens])
                prior = (1 - borrowing) * prior + borrowing * sibling_prior
            observed = np.asarray([counts[token] for token in tokens])
            enrichment = np.log(
                (observed + prior_mass * prior) / ((sizes[target] + prior_mass) * prior)
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
