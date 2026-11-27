"""Protect the frozen midpoint arithmetic and target/background boundary cases."""

import numpy as np

from vdjmatch.evalue.competitive import assignment_score, reduced_hits


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
