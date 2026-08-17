"""Precursor-frequency estimators.

The maths tests are regressions against an independent oracle, not golden numbers: the closed-form
ball is checked against brute-force enumeration of the same ball, and the exact union against
enumerating the union. Both oracles are slow and correct; the fast paths are what ships.
"""
import math

import pytest

from vdjmatch.precursor import (
    RecombinationEvent,
    ball_mass,
    check_junctions,
    closed_ball_mass,
    coverage_corrected_mass,
    event_ratio,
    expected_cells,
    load_model,
    motif_mass,
    observed_mass,
    p_at_least,
    paired_frequency,
    pgen,
    precursor_frequency,
    shell_profile,
    union_mass,
    unseen_junctions,
)

vdjtools = pytest.importorskip("vdjtools", reason="needs the optional vdjmatch[precursor] extra")

A = "CASSLAPGATNEKLFF"
B = "CASSLAPGATNEKLYF"      # Hamming-1 from A -> balls share 20 sequences
C = "CASSLAPGATNEKQYF"      # Hamming-2 from A -> balls share 2 sequences
FAR = "CASSQDRDTQYF"        # different length -> balls cannot overlap at any radius


@pytest.fixture(scope="module")
def model():
    return load_model("TRB")


# --------------------------------------------------------------------------- no model needed
def test_check_junctions_splits_on_the_conserved_anchors():
    ok, bad = check_junctions([A, "ASSLAPGATNEKL", B, "CASSLAPGATNEKLW"])
    assert ok == [A, B, "CASSLAPGATNEKLW"]
    assert bad == ["ASSLAPGATNEKL"]          # anchor-stripped IMGT CDR3 would score 0.0 silently


@pytest.mark.parametrize("bad,why", [
    (dict(donor="", v="V", j="J", junction_nt="ACG" * 3, junction_aa="CAF"), "pooled donor"),
    (dict(donor="d", v="V", j="J", junction_nt="CAF", junction_aa="CAF"), "aa as the nt key"),
    (dict(donor="d", v="V", j="J", junction_nt="ACG" * 2, junction_aa="CAF"), "length mismatch"),
])
def test_recombination_event_rejects_the_silent_mistakes(bad, why):
    with pytest.raises(ValueError):
        RecombinationEvent(**bad)


def test_p_at_least_matches_poisson_and_is_monotone_in_k():
    f, n = 1e-6, 1e7                          # lambda = 10
    assert p_at_least(f, n, 1) == pytest.approx(1 - math.exp(-10.0))
    assert p_at_least(f, n, 1) > p_at_least(f, n, 10) > p_at_least(f, n, 30)
    assert p_at_least(0.0, n, 1) == 0.0
    with pytest.raises(ValueError):
        p_at_least(f, n, 0)


def test_expected_cells_and_pairing_are_plain_products():
    assert expected_cells(1e-6, n_cells=1e11, compartment=0.3) == pytest.approx(3e4)
    assert paired_frequency(1e-3, 1e-3, kappa=0.5) == pytest.approx(5e-7)


# --------------------------------------------------------------------------- the closed ball
def test_closed_ball_r1_reproduces_the_native_hamming_1_mass(model):
    fast = closed_ball_mass(model, [A, B, C], r=1)
    native = pgen(model, [A, B, C], mismatches=1)
    for f, n in zip(fast, native):
        assert f == pytest.approx(n, rel=1e-12)


@pytest.mark.parametrize("r", [0, 1, 2])
def test_closed_ball_matches_brute_force_enumeration(model, r):
    from seqtree.distance import neighbourhood
    brute = sum(pgen(model, neighbourhood(A, r=r)))
    assert closed_ball_mass(model, [A], r=r)[0] == pytest.approx(brute, rel=1e-9)


def test_closed_ball_rejects_a_radius_past_the_sequence(model):
    with pytest.raises(ValueError):
        closed_ball_mass(model, ["CASF"], r=9)


# --------------------------------------------------------------------------- the union
@pytest.mark.parametrize("seqs", [[A, B], [A, C], [A, FAR], [A, B, C], [A, B, C, FAR]])
@pytest.mark.parametrize("r", [1, 2])
def test_union_mass_equals_the_enumerated_union(model, seqs, r):
    fast = union_mass(model, seqs, r=r)
    oracle = ball_mass(model, seqs, r=r)
    assert fast["union"] == pytest.approx(oracle["union"], rel=1e-9)
    assert fast["naive_sum"] == pytest.approx(oracle["naive_sum"], rel=1e-12)


def test_union_is_bounded_by_the_naive_sum_and_exceeds_the_observed_mass(model):
    u = union_mass(model, [A, B], r=1)
    assert observed_mass(model, [A, B]) < u["union"] < u["naive_sum"]
    assert u["overlap"] > 0
    assert u["n_multiply_covered"] == 20          # two centres at Hamming 1 share exactly 20


def test_disjoint_balls_have_nothing_to_correct(model):
    u = union_mass(model, [A, FAR], r=1)          # different lengths, substitutions only
    assert u["overlap"] == pytest.approx(0.0, abs=1e-12)
    assert u["n_multiply_covered"] == 0
    assert u["n_clustered"] == 0


def test_a_single_centre_union_is_its_closed_ball(model):
    u = union_mass(model, [A], r=1)
    assert u["union"] == pytest.approx(closed_ball_mass(model, [A], r=1)[0], rel=1e-12)
    assert u["overlap"] == 0.0


def test_union_of_an_empty_set_is_zero(model):
    assert union_mass(model, [], r=1)["union"] == 0.0
    assert union_mass(model, ["", None], r=1)["n_seqs"] == 0


def test_a_component_above_the_ceiling_raises_rather_than_thrashing(model):
    with pytest.raises(MemoryError, match="connected component"):
        union_mass(model, [A, B], r=2, max_members=10)


# --------------------------------------------------------------------------- shells
def test_shells_partition_the_union_and_retention_interpolates(model):
    prof = shell_profile(model, [A, B], r=2)
    assert sum(s["mass"] for s in prof["shells"]) == pytest.approx(prof["union"], rel=1e-9)
    assert all(s["mass"] >= 0 for s in prof["shells"])
    assert observed_mass(model, [A, B]) < prof["retained"] < prof["union"]
    assert prof["shells"][0]["mass"] == pytest.approx(observed_mass(model, [A, B]), rel=1e-9)


def test_alpha_1_is_the_raw_union_and_alpha_0_is_the_observed_bound(model):
    assert shell_profile(model, [A, B], r=1, alpha=1.0)["retained"] == pytest.approx(
        shell_profile(model, [A, B], r=1)["union"], rel=1e-12)
    assert shell_profile(model, [A, B], r=1, alpha=0.0)["retained"] == pytest.approx(
        observed_mass(model, [A, B]), rel=1e-9)


# --------------------------------------------------------------------------- motifs
def test_a_fully_pinned_motif_is_one_sequence_and_widening_only_adds_mass(model):
    assert motif_mass(model, list(A)) == pytest.approx(observed_mass(model, [A]), rel=1e-12)
    wide = list(A)
    wide[5] = ""
    assert motif_mass(model, wide) > motif_mass(model, list(A))


# --------------------------------------------------------------------------- events
def test_event_ratio_counts_donors_as_distinct_objects():
    nt = "TGT" * len(A)
    events = [RecombinationEvent("d1", "TRBV1", "TRBJ1", nt, A),
              RecombinationEvent("d2", "TRBV1", "TRBJ1", nt, A),   # converged: still two events
              RecombinationEvent("d1", "TRBV1", "TRBJ1", nt, A),   # same key: one event
              RecombinationEvent("d3", "TRBV1", "TRBJ1", "ACG" * len(FAR), FAR)]
    out = event_ratio([A], events, r=0)
    assert out["denominator"] == 3
    assert out["matched"] == 2


def test_event_ratio_radius_1_is_a_superset_of_exact_matching():
    nt = "TGT" * len(B)
    events = [RecombinationEvent("d1", "TRBV1", "TRBJ1", nt, B)]
    assert event_ratio([A], events, r=0)["matched"] == 0
    assert event_ratio([A], events, r=1)["matched"] == 1


def test_event_ratio_shares_one_denominator_across_epitopes():
    nt = "TGT" * len(A)
    events = [RecombinationEvent("d1", "TRBV1", "TRBJ1", nt, A)]
    out = event_ratio({"e1": [A], "e2": [FAR]}, events, r=0)
    assert out["epitopes"]["e1"]["matched"] == 1
    assert out["epitopes"]["e2"]["matched"] == 0
    assert out["denominator"] == 1


# --------------------------------------------------------------------------- unseen
def test_coverage_correction_degenerates_loudly_on_all_singletons(model):
    out = coverage_corrected_mass(model, [A, B, C], [1, 1, 1], n_units=3)
    assert out["degenerate"] and "singleton" in out["reason"]
    assert out["corrected"] == out["observed"]        # the bound, not an infinity


def test_coverage_correction_needs_two_capture_units(model):
    out = coverage_corrected_mass(model, [A, B], [1, 1], n_units=1)
    assert out["degenerate"] and "capture units" in out["reason"]


def test_coverage_correction_rejects_impossible_multiplicities(model):
    with pytest.raises(ValueError):
        coverage_corrected_mass(model, [A, B], [0, 1], n_units=2)
    with pytest.raises(ValueError):
        coverage_corrected_mass(model, [A, B], [5, 1], n_units=2)
    with pytest.raises(ValueError):
        coverage_corrected_mass(model, [A, B], [1], n_units=2)


def test_unseen_junctions_reports_richness_alongside_mass(model):
    seqs = [A, B, C, FAR, "CASSLAPGATNEKLWF", "CASSPAPGATNEKLFF"]
    out = unseen_junctions(model, seqs, [4, 3, 2, 1, 1, 1], n_units=5)
    if out["degenerate"]:
        pytest.skip(f"capture curve unidentified on this toy set: {out['reason']}")
    assert out["n_total"] >= out["n_observed"]
    assert out["n_unseen"] == pytest.approx(out["n_total"] - out["n_observed"], rel=1e-9)
    assert out["unseen_mass"] >= 0
    # the whole point: unseen members are the rare ones, so they are individually smaller
    assert out["rarity_ratio"] is None or out["rarity_ratio"] > 1.0


# --------------------------------------------------------------------------- the bridge
def test_precursor_frequency_is_linear_in_q_and_carries_its_calibration(model):
    raw = precursor_frequency(model, [A, B], r=1)
    cal = precursor_frequency(model, [A, B], r=1, q=9.41)
    assert raw["q"] == 1.0 and cal["q"] == 9.41
    assert cal["F"] == pytest.approx(9.41 * raw["F"], rel=1e-12)
    assert raw["F"] == pytest.approx(raw["retained"], rel=1e-12)


def test_precursor_frequency_omits_poisson_fields_without_n_eff(model):
    assert precursor_frequency(model, [A], r=1)["lambda"] is None
    out = precursor_frequency(model, [A], r=1, n_eff=1e8)
    assert out["lambda"] == pytest.approx(1e8 * out["F"], rel=1e-12)
    assert 0.0 <= out["p_ge_10"] <= out["p_ge_1"] <= 1.0


def test_occupancy_saturates_in_depth_and_splits_seen_from_unseen(model):
    from vdjmatch.precursor import occupancy

    seqs = [A, B, C, FAR]
    o = occupancy(model, seqs, r=1, n_eff=1e8)
    assert o["S"] > len(seqs)                       # the ball is bigger than the observed set
    assert 0 < o["n_seen"] < o["S"]
    assert o["n_unseen"] == pytest.approx(o["S"] - o["n_seen"], rel=1e-9)
    assert o["seen_fraction"] == pytest.approx(o["n_seen"] / o["S"], rel=1e-9)
    # deeper sampling can only reveal more, and it saturates at S rather than growing past it
    deep = occupancy(model, seqs, r=1, n_eff=1e14)
    assert deep["n_seen"] > o["n_seen"]
    assert deep["n_seen"] <= deep["S"] * (1 + 1e-9)
    assert occupancy(model, seqs, r=1)["n_seen"] is None    # no depth given, no count claimed


def test_occupancy_S_is_the_cognacy_weighted_ball_size(model):
    from seqtree.distance import union_size

    from vdjmatch.precursor import occupancy

    seqs = [A, FAR]                                  # disjoint balls, so shells are unambiguous
    o = occupancy(model, seqs, r=1, alpha=1.0)
    assert o["S"] == pytest.approx(union_size(seqs, r=1), rel=1e-12)
    # alpha down-weights the outer shell, so S falls between the observed count and the full ball
    o_a = occupancy(model, seqs, r=1, alpha=0.1)
    assert len(seqs) < o_a["S"] < o["S"]


def test_union_counts_its_members_and_skips_the_count_when_asked(model):
    u = union_mass(model, [A, B], r=1)
    assert u["n_union"] == ball_mass(model, [A, B], r=1)["n_union"]
    assert u["n_union"] > 2                       # the uncatalogued candidates are the rest
    assert union_mass(model, [A, B], r=1, count_members=False)["n_union"] is None


# --------------------------------------------------------------------------- the report / CLI
def test_report_row_carries_the_spread_and_the_finite_ball_census(model):
    from vdjmatch.precursor import summarise_group
    row = summarise_group(model, [A, B, C, FAR, "ASSLAPGATNEKL"], group="e", chain="TRB")
    assert row["n_records"] == 5 and row["n_junctions"] == 4 and row["n_dropped"] == 1
    assert row["pgen_log10_span"] > 0 and 0 < row["top1_share"] <= row["top10_share"] <= 1.0
    assert row["n_unseen_ball"] == row["n_ball"] - row["n_junctions"]
    assert row["unseen_ball_mass"] == pytest.approx(row["union"] - row["observed_mass"], rel=1e-9)
    assert row["unseen_status"] == "no capture units given"


def test_cli_precursor_writes_one_row_per_group(tmp_path, model):
    import polars as pl

    from vdjmatch.cli.__main__ import main

    src = tmp_path / "tcrs.tsv"
    pl.DataFrame({"cdr3": [A, B, FAR, C], "epitope": ["e1", "e1", "e2", "e2"],
                  "ref": ["s1", "s2", "s1", "s1"]}).write_csv(src, separator="\t")
    out = tmp_path / "out.txt"
    assert main(["precursor", str(src), "--group-by", "epitope", "--capture-col", "ref",
                 "--q", "9.41", "-o", str(out)]) == 0
    got = pl.read_csv(out, separator="\t")
    assert got.height == 2
    assert set(got["group"]) == {"e1", "e2"}
    assert (got["q"] == 9.41).all()
    assert (got["F"] > 0).all()
