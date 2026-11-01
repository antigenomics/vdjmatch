"""Exact count identities and origin/record accounting, without scientific reruns."""

import numpy as np
import polars as pl
import pytest
from vdjmatch.match.scoring_audit import column_pair_counts, origin_block_matrices


def test_origin_pair_counter_matches_explicit_distinct_pairs():
    encoded = np.array([[0, 0], [1, 0], [0, 0]], dtype=np.uint8)
    weights = np.array([[0.2, 0], [0.5, 0], [0.7, 0]])
    counts, _, stats = column_pair_counts(encoded, np.arange(3), weights)
    assert np.allclose(counts[:2, :2], [[0.28, 0.45], [0.45, 0]])
    assert np.isclose(stats["column_pairs"], 0.59)
    assert column_pair_counts(encoded, np.array([0, 0, 1]), weights)[0][
        0, 1
    ] == pytest.approx(0.175)
    with pytest.raises(ValueError):
        column_pair_counts(encoded, np.arange(3), weights * np.nan)
    assert (
        column_pair_counts(
            np.zeros((0, 2), dtype=int), np.array([], dtype=int), np.zeros((0, 2))
        )[0].sum()
        == 0
    )


def _inputs():
    seqs = ["AAA", "AAC", "CAA", "CAC", "DDD"]
    clones = pl.DataFrame(
        dict(
            query_id=list(range(5)),
            cdr3=seqs,
            v=["V"] * 5,
            j=["J"] * 5,
            epitope=["E"] * 5,
        )
    )
    rows = [
        dict(
            query_id=i,
            position=p,
            n_probability=(None if i == 4 else 0.0 if p == 1 else 1.0),
        )
        for i in range(5)
        for p in range(3)
    ]
    return clones, pl.DataFrame(rows)


def test_origin_blocks_remove_fixed_identities_and_retain_missing_cells():
    clones, origins = _inputs()
    matrix, blocks, support, pairs = origin_block_matrices(clones, origins)
    assert support.filter(~pl.col("origin_supported"))["query_id"].to_list() == [4]
    assert blocks.height == 1 and blocks["clonotypes"][0] == 4
    assert pairs["n_weighted_ordered_site_pairs"].sum() == 24
    assert origin_block_matrices(clones.reverse(), origins.reverse())[0].equals(matrix)
    pooled = matrix.filter(pl.col("weighting") == "pooled")
    # Two N columns, each with 2 A and 2 C: ordered AA=CC=4, AC=CA=8.
    assert (
        pooled.filter((pl.col("a") == "A") & (pl.col("b") == "A"))["pair_weight"][0]
        == 4
    )
    ac = pooled.filter((pl.col("a") == "A") & (pl.col("b") == "C")).row(0, named=True)
    assert ac["pair_weight"] == 8 and ac["supported"] and ac["distance_bits"] == 0
    assert (
        pooled.filter((pl.col("a") == "D") & (pl.col("b") == "D"))["distance_bits"][0]
        is None
    )
    assert (
        np.isclose(pooled["q"].sum(), 1.0)
        and blocks["n_weighted_unordered_site_pairs"][0] == 12
    )
    # A smaller epitope has AA=2, AC=CA=1: equal-epitope and pooled differ.
    second = clones.filter(pl.col("query_id") < 2).with_columns(
        pl.lit("F").alias("epitope")
    )
    repeated, _, _, _ = origin_block_matrices(pl.concat([clones, second]), origins)
    aa = repeated.filter((pl.col("a") == "A") & (pl.col("b") == "A"))
    assert aa.filter(pl.col("weighting") == "equal_epitope")["q"][0] == pytest.approx(1/3)
    assert aa.filter(pl.col("weighting") == "pooled")["q"][0] == pytest.approx(3/14)
    with pytest.raises(ValueError, match="every codon"):
        origin_block_matrices(clones, origins.slice(1))
    with pytest.raises(ValueError, match="unique"):
        origin_block_matrices(pl.concat([clones, clones.slice(0, 1)]), origins)


def test_origin_cli_preserves_source_units_and_zero_support(tmp_path, monkeypatch):
    import argparse
    import json
    import vdjtools.model
    from vdjmatch import db
    from vdjmatch.cli.scoring_audit import register

    seqs = [
        "CAVRDSNYQLIW",
        "CAVRDSNYQLIW",
        "CAVRESNYQLIW",
        "CAVRDSNYQLIF",
        "CAVRDSNYQLIW",
    ]
    frame = pl.DataFrame(
        dict(
            cdr3=seqs,
            v=["TRAV1-1*01"] * 4 + ["TRAV1-1*01,TRAV2*01"],
            j=["TRAJ1*01"] * 5,
            species=["HomoSapiens"] * 5,
            gene=["TRA"] * 5,
            epitope=["E"] * 5,
            mhc_a=["A"] * 5,
            mhc_b=["B"] * 5,
            mhc_class=["MHCI"] * 5,
            reference_valid=[True] * 5,
        )
    )
    monkeypatch.setattr(db, "load", lambda *args, **kwargs: frame)
    monkeypatch.setattr(db, "provenance", lambda p: dict(test_raw=str(p)))

    def posterior(model, cdr3_aas, **kwargs):
        assert (
            len(cdr3_aas) == 3
        )  # Duplicate source records do not become extra origin queries.
        assert all("," not in v for v in kwargs["v"])
        return pl.DataFrame(
            [
                dict(
                    row=i,
                    position=p,
                    n_probability=(
                        None if s.endswith("F") else 0.5 if 0 < p < len(s) - 1 else 0.0
                    ),
                    supported=not s.endswith("F"),
                )
                for i, s in enumerate(cdr3_aas)
                for p in range(len(s))
            ]
        )

    monkeypatch.setattr(
        vdjtools.model, "n_region_posterior_batch", posterior, raising=False
    )
    raw = tmp_path / "raw.tsv"
    raw.write_text("fixture source")
    parser = argparse.ArgumentParser()
    register(parser.add_subparsers())
    prefix = tmp_path / "origin"
    args = parser.parse_args(
        [
            "scoring-audit",
            "--task",
            "origin-matrix",
            "--pmhc",
            str(raw),
            "--locus",
            "TRA",
            "--min-records",
            "1",
            "--output-prefix",
            str(prefix),
        ]
    )
    assert args.func(args) == 0
    support = pl.read_csv(str(prefix) + ".support.tsv", separator="\t")
    assert support["source_records"].sum() == 5
    observed = {r["status"]: r["source_records"] for r in support.to_dicts()}
    assert observed == dict(
        origin_supported=3, zero_sequence_likelihood=1, unresolved_call=1
    )
    calls = pl.read_csv(str(prefix) + ".calls.tsv", separator="\t")
    assert calls.filter(pl.col("status") == "ambiguous_call").height == 1
    manifest = json.loads((tmp_path / "origin.provenance.json").read_text())
    assert manifest["origin_queries"] == 3 and manifest["source_records"] == 5
    assert pl.read_csv(str(prefix) + ".matrix.tsv", separator="\t").height == 800
