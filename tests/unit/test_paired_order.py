"""Checks for cohort-relative, label-free paired ordering."""
import subprocess
import sys

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from vdjmatch.aggregate.paired_order import paired_evidence_order, OUTPUTS, COUNTS


def sample():
    return pl.DataFrame({
        "query_id": list("abcdefg"), "cohort": ["x"] * 7,
        "species": ["HomoSapiens"] * 7, "epitope": ["PEPTIDE"] * 7,
        "mhc_a": ["A*02:01"] * 7, "mhc_b": [""] * 7, "mhc_class": ["MHCI"] * 7,
        "paired_score": [1., 1., 9., 8., 7., 6., None],
        "fallback_score": [1., 2., 3., 4., 5., 6., None],
        "alpha_score": [100., -100., 1., 2., 3., 0., None],
        "beta_score": [100., -100., 3., 2., 1., 0., None],
        "n_reference_pair": [1, 1, 0, 0, 0, 0, 0],
        "n_reference_alpha": [1, 1, 1, 1, 1, 0, 1],
        "n_reference_beta": [1, 1, 0, 0, 0, 0, 0],
        "paired_available": [True] * 6 + [False],
    })


def test_blocks_ties_and_permutation():
    frame = sample()
    out = paired_evidence_order(frame)
    assert out["evidence_block"].to_list() == [1, 1, 2, 2, 2, 3, 3]
    assert out["evidence_block"].dtype == pl.UInt8
    assert out["evidence_score"].to_list() == [2.5, 2.5, 1.5, 1.5, 1.5, 2/3, 1/3]
    assert_frame_equal(out.select(frame.columns), frame.with_columns(pl.col(COUNTS).cast(pl.UInt32)))
    assert_frame_equal(out, paired_evidence_order(frame.reverse()).sort("query_id"))
    scaled = frame.with_columns([(pl.col(c) * 3 + 9).alias(c) for c in
                                 ["paired_score", "fallback_score", "alpha_score", "beta_score"]])
    assert_frame_equal(out.select(OUTPUTS), paired_evidence_order(scaled).select(OUTPUTS))
    # Moving beta alone breaks the balanced rank sums; block-1 outliers have no effect.
    changed = frame.with_columns(pl.Series("beta_score", [100., -100., 1., 2., 3., 0., None]))
    assert paired_evidence_order(changed)["evidence_score"][2:5].n_unique() == 3


def test_grouping_and_unavailable_scores():
    frame = sample()
    other = frame.with_columns(pl.lit("A*01:01").alias("mhc_a"), -pl.col("alpha_score"))
    assert paired_evidence_order(frame.with_columns(pl.lit(None).alias("mhc_b"))).height == 7
    combined = paired_evidence_order(pl.concat([frame, other]))
    assert_frame_equal(combined.head(7), paired_evidence_order(frame))
    observed_no_hit = frame.with_columns(pl.when(pl.col("query_id") == "c").then(None)
                                         .otherwise(pl.col("beta_score")).alias("beta_score"))
    assert paired_evidence_order(observed_no_hit)["evidence_block"][2] == 2
    empty_evidence = frame.with_columns([pl.lit(None).alias(c) for c in
        ["paired_score", "fallback_score", "alpha_score", "beta_score",
         "n_reference_pair", "n_reference_alpha", "n_reference_beta"]])
    assert paired_evidence_order(empty_evidence)["evidence_block"].to_list() == [3] * 7
    assert paired_evidence_order(frame.head(0)).height == 0
    assert paired_evidence_order(frame.head(1))["evidence_score"].item() == 2.5
    assert paired_evidence_order(frame.slice(2, 3))["evidence_score"].to_list() == [1.5] * 3


@pytest.mark.parametrize("name,value", [("alpha_score", float("inf")),
    ("beta_score", float("nan")), ("n_reference_pair", -1), ("n_reference_alpha", .5)])
def test_invalid_values(name, value):
    with pytest.raises(ValueError):
        paired_evidence_order(sample().with_columns(pl.lit(value).alias(name)))


def test_invalid_identity_and_types():
    frame = sample()
    for invalid in [pl.concat([frame, frame.head(1)]), frame.drop("cohort"),
                    frame.with_columns(pl.lit("1").alias("alpha_score")),
                    frame.with_columns(pl.lit(1).alias("paired_available")),
                    frame.with_columns(pl.lit("").alias("query_id")),
                    frame.with_columns(pl.lit(2**32).alias("n_reference_pair"))]:
        with pytest.raises(ValueError):
            paired_evidence_order(invalid)


def test_cli(tmp_path):
    source, target = tmp_path / "input.tsv", tmp_path / "output.tsv"
    sample().write_csv(source, separator="\t")
    result = subprocess.run([sys.executable, "-m", "vdjmatch.cli", "paired-evidence-order",
                             str(source), "--output", str(target)],
                            capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stderr
    out = pl.read_csv(target, separator="\t")
    assert out["evidence_block"].to_list() == [1, 1, 2, 2, 2, 3, 3]
    assert out["evidence_score"].to_list() == paired_evidence_order(sample())["evidence_score"].to_list()


def test_support_oracle_and_balanced_noise():
    # Two balanced tied blocks: perfect support gives AP=1, independent support AP=1/2.
    base = sample().head(4).with_columns(
        pl.Series("n_reference_pair", [1, 1, 0, 0]),
        pl.lit(0).alias("n_reference_alpha"), pl.lit(0).alias("n_reference_beta"),
        pl.lit(1.).alias("paired_score"), pl.lit(1.).alias("fallback_score"),
    )
    scores = paired_evidence_order(base)["evidence_score"].to_list()

    def ap(labels):
        total, positives, result = 0, 0, 0.
        for score in sorted(set(scores), reverse=True):
            hits = sum(label for label, value in zip(labels, scores) if value == score)
            total += scores.count(score)
            positives += hits
            result += hits * positives / total
        return result / sum(labels)

    assert ap([1, 1, 0, 0]) == 1.
    assert ap([1, 0, 1, 0]) == .5
