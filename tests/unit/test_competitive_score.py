"""Protect the frozen midpoint arithmetic and target/background boundary cases."""

import numpy as np
import polars as pl
import pytest

from vdjmatch.evalue.competitive import (
    anchored_motif_scores, assay_local_score, assignment_score, paired_assay_local_score,
    fr3_sw_distance_correction, paired_motif_ranking_scores, reduced_hits,
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
    np.testing.assert_array_equal(adjusted, scores + np.log1p(counts) * (1 / 400))
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
    for invalid in (np.nan, -np.inf, -1, 1j):
        with pytest.raises(ValueError, match="real and nonnegative"):
            paired_assay_local_score([[0, invalid], [invalid, 0]])


def test_anchored_motif_scores_counts_distinct_reference_observations():
    refs = [("a", "CASAAAFF"), ("a", "CASAAAFF"), ("b", "CASRRRFF")]
    result = anchored_motif_scores(refs, ["CASAAAFF", "CAS"], ["a", "b", "empty"])
    prior = 3 / 5
    np.testing.assert_allclose(result[0, :2], np.log([(2 + 2 * prior) / (4 * prior), 2 / 3]))
    assert np.isneginf(result[1]).all()
    assert np.isneginf(result[:, 2]).all()


def test_anchored_motif_scores_borrow_from_compatible_siblings():
    refs = [("a", "CASAAAFF"), ("b", "CASRRRFF")]
    result = anchored_motif_scores(
        refs, ["CASAAAFF"], ["a", "b"],
        reference_pool=["CASAAAFF", "CASRRRFF"],
        sibling_references={"a": ["CASAAAFA"]},
    )
    prior = (2 / 3) * (1 / 2) + (1 / 3) * (2 / 3)
    np.testing.assert_allclose(result[0], np.log([(1 + 2 * prior) / (3 * prior), 2 / 3]))


def test_anchored_tetramer_uses_longer_windows():
    refs = [("a", "CASAAAAFF"), ("b", "CASRRRRFF")]
    result = anchored_motif_scores(refs, ["CASAAAAFF", "CASAAAFF"], ["a", "b"], motif_length=4)
    np.testing.assert_allclose(result[0], np.log([4 / 3, 2 / 3]))
    assert np.isneginf(result[1]).all()


def test_research_score_cli_keeps_output_identity(tmp_path):
    from vdjmatch.cli.__main__ import main

    queries = tmp_path / "queries.tsv"
    queries.write_text("query_id\tjunction_aa\n01\tCASAAAFF\n02\tCASRRRFF\n")
    reference = tmp_path / "reference.tsv"
    reference.write_text(
        "target\tjunction_aa\tv_call\tj_call\nA\tCASAAAFF\tV1\tJ1\n"
        "A\tCASAAAFF\tV1\tJ1\nB\tCASRRRFF\tV1\tJ1\n"
    )
    output = tmp_path / "scores.tsv"
    assert main(["anchored-motif-score", "--reference", str(reference), "--queries",
                 str(queries), "--output", str(output)]) == 0
    calls = pl.read_csv(output, separator="\t", schema_overrides={"query_id": pl.String})
    assert calls.filter(pl.col("winner")).select("query_id", "target").rows() == [("01", "A"), ("02", "B")]
    distances = tmp_path / "distances.npy"
    np.save(distances, [[0, 3], [3, 0]])
    assert main(["paired-assay-score", "--distances", str(distances), "--queries",
                 str(queries), "--output", str(output)]) == 0
    assert output.read_text() == "query_id\tscore\n01\t-3.0\n02\t-3.0\n"
    targets = tmp_path / "targets.tsv"
    targets.write_text("target\nA\nB\n")
    scores = tmp_path / "rank-input.npz"
    np.savez(scores, base=[[-1, -3], [-1, -1]], tra=[[0, -1], [0, 0]],
             trb=[[0, -2], [0, 0]], target_key_counts=[5, 6])
    assert main(["paired-motif-rank", "--scores", str(scores), "--queries", str(queries),
                 "--targets", str(targets), "--support-weight", "0", "--conflict-weight",
                 "4", "--output", str(output)]) == 0
    assert output.read_text() == "query_id\ttarget\tscore\n01\tA\t-1.0\n01\tB\t-11.0\n02\tA\t-1.0\n02\tB\t-1.0\n"


def test_paired_motif_cli_rebuilds_fixed_class_memberships(tmp_path):
    import polars as pl

    from vdjmatch.cli.__main__ import main

    path = {stem: tmp_path / (stem + ".tsv") for stem in ("queries", "targets", "reference", "output")}
    junctions = ("junction_aa_alpha", "junction_aa_beta")
    queries = pl.DataFrame({"query_id": ["01", "02"], junctions[0]: ["CASAAAFF", "CASRRRFF"],
                            junctions[1]: ["CASRRRFF", "CAS"]})
    targets = pl.DataFrame({"target": ["B", "C", "A"], "mhc_class": ["MHCI", "MHCII", "MHCI"]})
    members = pl.DataFrame([
        ("p1", "A", "CASAAAFF", "CASRRRFF"),
        ("p1", "B", "CASAAAFF", "CASRRRFF"),
        ("p2", "A", "CASAAAFF", "CASAAAFF"),
        ("p3", "A", "CASAAAFF", "CASAAAFF"),
        ("p1", "C", "CASAAAFF", "CASRRRFF"),
        ("p5", "C", "CASAAAFF", "CASAAAFF"),
        ("p6", "C", "CASAAAFF", "CASAAAFF"),
    ], schema=["reference_key", "target", *junctions], orient="row")
    for name, frame in (
        ("queries", queries), ("targets", targets),
        ("reference", pl.concat([members, members.head(1)])),
    ):
        frame.write_csv(path[name], separator="\t")
    base = np.array([[-1, -3, -2], [-1, -1, -2]])
    archive = tmp_path / "base.npz"
    np.savez(archive, base=base)
    assert main(["paired-motif-rank", "--scores", str(archive), "--queries", str(path["queries"]),
                 "--targets", str(path["targets"]), "--reference", str(path["reference"]),
                 "--support-weight", "0", "--conflict-weight", "4", "--output", str(path["output"])]) == 0
    expected = []
    for chain in junctions:
        pairs = list(members.select("target", chain).iter_rows())
        class_i = anchored_motif_scores([pair for pair in pairs if pair[0] != "C"], queries[chain], ["B", "A"])
        class_ii = anchored_motif_scores([pair for pair in pairs if pair[0] == "C"], queries[chain], ["C"])
        expected.append(np.column_stack((class_i[:, 0], class_ii[:, 0], class_i[:, 1])))
    scores = paired_motif_ranking_scores(base, *expected, [1, 3, 3], support_weight=0, conflict_weight=4)
    actual = pl.read_csv(path["output"], separator="\t", schema_overrides={"query_id": pl.String})
    assert actual.select("query_id", "target").rows() == [
        (query, target) for query in queries["query_id"] for target in targets["target"]
    ]
    assert np.any(scores != base)
    np.testing.assert_array_equal(actual["score"].to_numpy().reshape(base.shape), scores)


def test_frozen_paired_base_masks_hits_in_both_lists_and_preserves_unsupported_rows(monkeypatch):
    import polars as pl

    from vdjmatch.cli import paired_vx_base
    from vdjmatch.evalue import equal_tail
    from vdjmatch.match import geometry

    fixed = pl.DataFrame({
        "cdr3a": ["CASAAAFF"] * 3, "cdr3b": ["CASRRRFF", "CASAAAFF", "CASYYYFF"],
        "nativea": ["a"] * 3, "nativeb": ["b"] * 3,
    })
    query = fixed.head(1)
    queries = pl.concat([query, query.with_columns(pl.lit(None, dtype=pl.String).alias("nativea"))]
                        ).with_row_index("query_id")
    # Equal-tail also imports the native function. Load it before stubbing so
    # other tests cannot inherit this test's one-row distance matrix.
    assert equal_tail.paired_distance_matrix is geometry.paired_distance_matrix
    monkeypatch.setattr(geometry, "paired_distance_matrix", lambda *args, **kwargs: np.array([[0, 4, 9]]))
    monkeypatch.setattr(paired_vx_base, "_regional_corrections", lambda _, locus: (
        np.array(["a" if locus == "TRA" else "b"]), np.zeros((1, 1))))
    groups, backgrounds = ([0, 1], [2]), ([2], [0, 1])
    for exclude_exact, expected in ((False, [-.275, -9.525]), (True, [-4.175, -9.425])):
        actual = paired_vx_base._base_scores(queries, fixed, groups, backgrounds, {},
                                             threads=1, exclude_exact=exclude_exact)
        np.testing.assert_allclose(actual[0], expected)
        np.testing.assert_array_equal(actual[1], [-np.inf, -np.inf])
        maxima = paired_vx_base._vx_maxima(queries, fixed, groups, backgrounds, {}, ["A", "B"], None,
                                           threads=1, exclude_exact=exclude_exact)
        assert maxima.equals(paired_vx_base._maxima_table(queries["query_id"], actual, ["A", "B"], None))


def test_frozen_paired_maxima_exports_only_unique_winners_as_callable_scores():
    from vdjmatch.cli.paired_vx_base import _maxima_table

    scores = np.array([[0, -1], [2, 2], [-np.inf, -np.inf]])
    expected = [("01", "A", 0., 0., False), ("02", None, 2., None, True),
                ("03", None, None, None, False)]
    table = _maxima_table(["01", "02", "03"], scores, ["A", "B"], None)
    assert table.rows() == expected
    nulls = _maxima_table(["01", "02", "03"], scores, ["A", "B"], "generated")
    assert nulls.rows() == [(*row, "null", "generated") for row in expected]


def test_paired_winner_bound_retains_exact_ties_and_checks_the_runner_up():
    from vdjmatch.cli.paired_vx_base import _competitive_winners

    native = np.array([[0, 0, 1, 100], [0, 1, 5, 100], [np.inf, np.inf, np.inf, np.inf]])
    hits = native.copy()
    maximum, winner, tied = _competitive_winners(
        native, hits, np.arange(4), np.arange(4), np.arange(4), ([2], [2], [3], [0]))
    np.testing.assert_array_equal(maximum, [-11 * (1 / 40), -7 * (1 / 40), -np.inf])
    np.testing.assert_array_equal(winner[:2], [0, 0])
    np.testing.assert_array_equal(tied, [True, False, False])

    # The closest hit loses to a hit with a slightly larger corrected distance
    # when the latter has no close other-peptide competitor.
    hits[1, 1] = .1
    maximum, winner, tied = _competitive_winners(
        native, hits, np.arange(4), np.arange(4), np.arange(4), ([2], [3], [3], [0]))
    assert (maximum[1], winner[1], tied[1]) == (-.1, 1, False)

    # The first target has two distinct keys at the same nearest distance.
    # Both targets attain the incumbent exactly at their upper bound.
    maximum, _, tied = _competitive_winners(
        native[:1], hits[:1], [0, 1, 2, 0], [0, 3], [0, 0, 0, 1], ([3], [3]))
    assert (maximum[0], tied[0]) == (0., True)


def test_single_chain_bound_includes_trb_bonus_and_tra_uses_fixed_rank():
    from vdjmatch.cli.paired_vx_base import _competitive_winners, _competitor

    native = np.array([[0, 1, 50, 100]])
    hits = np.array([[0., .4, 50., 100.]])
    maximum, winner, tied = _competitive_winners(
        native, hits, [0, 1, 2], [0, 1], [0, 1, 1], ([3], [3]),
        locus="TRB", source_adjustment=np.array([0., 9 / 32, 0.])
    )
    assert winner.tolist() == [1] and tied.tolist() == [False]
    np.testing.assert_allclose(maximum, [3 / 20 - (.4 - 9 / 32)])

    backgrounds = np.array([[2, 5, 4, 1, 3, 9, 7, np.inf, 6, np.inf]])
    np.testing.assert_array_equal(_competitor(backgrounds, 8), [9])
    np.testing.assert_array_equal(_competitor(backgrounds[:, :8], 8), [np.inf])


def test_paired_withholding_filters_target_labels_and_competitor_memberships(monkeypatch, tmp_path):
    import polars as pl

    from vdjmatch.cli import paired_vx_base

    holdouts = tmp_path / "withheld.tsv"
    holdouts.write_text("query_id\tepitope\n01\tH\n01\tH\n02\tP\n02\tH\n02\tQ\n")
    queries = pl.DataFrame({"query_id": ["01", "02"]})
    targets = pl.DataFrame({"epitope": ["H", "P", "H", "Q"]})
    folds = paired_vx_base._withholding_folds(holdouts, queries["query_id"], {"H", "P", "Q"})
    assert folds == {("H",): [0], ("H", "P", "Q"): [1]}
    groups = [[0], [1], [2], [3]]
    backgrounds = [
        (np.array([3]), np.array(["Q"])),
        (np.array([2, 3, 3]), np.array(["H", "H", "Q"])),
        (np.array([3]), np.array(["Q"])),
        (np.array([1]), np.array(["P"])),
    ]

    def score(part, fixed, members, others, context, *, threads, exclude_exact):
        assert part["query_id"].to_list() == ["01"]
        assert members == [[1], [3]]
        assert [x.tolist() for x in others] == [[3], [1]]
        assert exclude_exact
        return np.array([[-3, -5]])

    monkeypatch.setattr(paired_vx_base, "_base_scores", score)
    result = paired_vx_base._withheld_scores(
        queries, None, targets, groups, backgrounds, None, folds=folds, threads=1, exclude_exact=True)
    np.testing.assert_array_equal(result, [[-np.inf, -3, -np.inf, -5], [-np.inf] * 4])
    with pytest.raises(ValueError, match="every query ID"):
        paired_vx_base._withholding_folds(holdouts, ["01"], {"H", "P", "Q"})


def test_paired_motif_refinement_leaves_top_ties_and_unsupported_cells():
    frozen = [[-1, -1, -2, -np.inf], [-1, -2, -3, -4]]
    tra = [[4, -4, 5, 3], [0, -2, 4, 7]]
    trb = [[7, 7, 3, 0], [0, -1, 3, 6]]
    actual = paired_motif_ranking_scores(frozen, tra, trb, [9, 8, 7, 2],
                                         support_weight=1, conflict_weight=1/8)
    cap = np.nextafter(-1., -np.inf)
    np.testing.assert_array_equal(actual, [[-1, -1, cap, -np.inf], [-1, -2.25, cap, -4]])
    alternative = [[0, 0, -10, 4], [4, -np.inf, -6, -2]]
    refined = paired_motif_ranking_scores(frozen, tra, trb, [9, 8, 7, 2],
                                          support_weight=1, conflict_weight=1/8,
                                          lower_tier_base=alternative)
    np.testing.assert_array_equal(refined, [[-1, -1, -7, -np.inf], [-1, -2.25, -3, -2]])
    with pytest.raises(ValueError, match="finite score range"):
        paired_motif_ranking_scores([[0, -1]], [[0, -2]], [[0, -2]], [3, 3],
                                    support_weight=0, conflict_weight=1e308)


def test_fr3_sw_research_correction_matches_verified_native_units():
    from vdjmatch.match.geometry import precompute_v_loop_distances

    context = precompute_v_loop_distances()
    for locus, pair, residual in (
        ("TRA", ("TRAV1-1*01", "TRAV9-2*04"), 0.5),
        ("TRB", ("TRBV1*01", "TRBV9*03"), -4.5),
    ):
        alleles, correction = fr3_sw_distance_correction(context, locus)
        np.testing.assert_array_equal(alleles, context["tables"][locus][0])
        np.testing.assert_array_equal(correction, correction.T)
        np.testing.assert_array_equal(correction.diagonal(), 0)
        i, j = np.searchsorted(alleles, pair)
        assert correction[i, j] == residual * (3 / 2 + 1 / 256) / 512
