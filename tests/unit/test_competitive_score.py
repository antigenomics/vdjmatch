"""Protect the frozen midpoint arithmetic and target/background boundary cases."""

import numpy as np

from vdjmatch.evalue.competitive import (
    anchored_motif_scores, assay_local_score, assignment_score, paired_assay_local_score,
    reduced_hits,
)


def test_competitive_score_binary64():
    first = np.array([-0.125, 1.3, 4.25, 21.0, 9.0])
    second = np.array([-0.125, 1.4, 4.25, 222.0, 15.0])
    background = np.array([np.inf, 2.5, 12.0, 20.0, 9.0])
    for chains, gamma, margin, weight in (
        (1, 0.1, 12, -0.1),
        (1, 0.25, 4.25, 0.00625),
        (2, 0.025, 6, -0.1),
    ):
        rules = {
            "margin_penalty": gamma,
            "margin_per_chain": margin,
            "second_hit_gap_weight": weight,
        }
        gap = np.stack((first, second), axis=1).sum(axis=1) / 2 - first
        penalty = np.maximum(first + chains * margin - background, 0)
        frozen = -first - gamma * penalty + weight * np.clip(gap, 0, 24 * chains)
        np.testing.assert_array_equal(
            reduced_hits(first, second, background, rules, chains).view("u8"),
            frozen.view("u8"),
        )


def test_size_adjustment_uses_rounded_slope():
    scores = np.array([[-12.0, -12.0, -np.inf]])
    counts = np.array([1, 7, 15])
    adjusted = assignment_score(scores, counts)
    np.testing.assert_array_equal(adjusted, scores + 0.002832 * np.log1p(counts))
    assert adjusted.argmax(axis=1).item() == 1


def test_assay_local_score_keeps_equal_distances_and_masked_hits():
    ordered = [[2, 2, 8], [1, 4, 4], [1, np.inf, np.inf]]
    np.testing.assert_array_equal(assay_local_score(ordered, 2, 1), [-2, -2, -np.inf])


def test_paired_assay_local_score_preserves_membership_and_input():
    distances = np.array([[0, 1, 3], [1, 0, 2], [3, 2, 0]])
    exact_hits = distances == 1
    np.testing.assert_array_equal(paired_assay_local_score(distances), [-1, -1, -2])
    np.testing.assert_array_equal(paired_assay_local_score(distances, exact_hits=exact_hits), [-3, -2, -2])
    np.testing.assert_array_equal(distances.diagonal(), [0, 0, 0])
    assert paired_assay_local_score([[0]]).item() == -np.inf


def test_anchored_motif_scores_counts_distinct_reference_observations():
    refs = [("a", "CASAAAFF"), ("a", "CASAAAFF"), ("b", "CASRRRFF")]
    result = anchored_motif_scores(refs, ["CASAAAFF", "CAS"], ["a", "b", "empty"])
    prior = 3 / 5
    np.testing.assert_allclose(result[0, :2], np.log([(2 + 2 * prior) / (4 * prior), 2 / 3]))
    assert np.isneginf(result[1]).all()
    assert np.isneginf(result[:, 2]).all()
