"""The bench hook for a geometry arm: contract, puncturing and scoring direction.

`bench/geometry_detect.py` is the only script in `bench/` that exercises a vdjmatch
junction geometry. Its task-level run needs the supplied samples (restricted, not in
CI), so these tests cover the pure pieces and the argument contract, which is where a
regression would actually hide.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bench"))

geometry_detect = pytest.importorskip("geometry_detect")

from vdjmatch.match.geometry import load_v_loops  # noqa: E402

MODEL = load_v_loops("human")
REFS = ["CASSLGQAYEQYF", "CASSAGQAYEQYF", "CASSIVGGNEQFF"]
RV = ["TRBV19*01"] * 3


def test_resolvable_keeps_modelled_alleles_and_counts_the_rest():
    cdr3 = ["CASSLGQAYEQYF", "CASSAGQAYEQYF", "CASSIVGGNEQFF"]
    vs = ["TRBV19*01", "nonsense-gene", "TRBV19*01"]
    kept, alleles, dropped = geometry_detect.resolvable(cdr3, vs, MODEL, "TRB")
    assert dropped == 1
    assert len(kept) == len(alleles) == 2
    assert all(MODEL[a]["locus"] == "TRB" for a in alleles)


def test_resolvable_drops_short_junctions():
    kept, _, dropped = geometry_detect.resolvable(
        ["CASSLGQAYEQYF", "CASSF"], ["TRBV19*01"] * 2, MODEL, "TRB")
    assert len(kept) == 1 and dropped == 1


def test_resolvable_rejects_the_wrong_locus():
    kept, _, dropped = geometry_detect.resolvable(
        ["CASSLGQAYEQYF"], ["TRBV19*01"], MODEL, "TRA")
    assert kept == [] and dropped == 1


def test_resolvable_handles_empty_input():
    assert geometry_detect.resolvable([], [], MODEL, "TRB") == ([], [], 0)


def _score(queries, qv, **geom):
    return geometry_detect.score_task(
        queries, qv, REFS, RV, "human", 1, 90.0,
        dict(ctrim=2, v_loop_weight=1, alphabet="full",
             stem_depth=geom.get("stem_depth", 0), stem_weight=1.0))


def test_exact_identity_is_punctured_not_scored_as_a_perfect_hit():
    """A query that IS a reference must not score itself; that is the honest-benchmark rule."""
    nearest, _ = _score(["CASSLGQAYEQYF"], ["TRBV19*01"])
    assert nearest[0] > 0, "the query's own junction was counted as its nearest neighbour"


def test_nearest_distance_is_finite_when_a_different_reference_exists():
    nearest, kernel = _score(["CASSLGQAYEQYF"], ["TRBV19*01"])
    assert np.isfinite(nearest[0])
    assert kernel[0] > 0


def test_a_query_matching_every_reference_exactly_has_no_support():
    nearest, kernel = geometry_detect.score_task(
        ["CASSLGQAYEQYF"], ["TRBV19*01"], ["CASSLGQAYEQYF"], ["TRBV19*01"],
        "human", 1, 90.0,
        dict(ctrim=2, v_loop_weight=1, alphabet="full", stem_depth=0, stem_weight=1.0))
    assert not np.isfinite(nearest[0]) and kernel[0] == 0


def test_stem_never_decreases_the_nearest_distance():
    """Same non-negativity the dense geometry guarantees, seen through the producer."""
    plain, _ = _score(["CASSLGQAYEQYF"], ["TRBV19*01"], stem_depth=0)
    stem, _ = _score(["CASSLGQAYEQYF"], ["TRBV19*01"], stem_depth=4)
    assert stem[0] >= plain[0]


def test_batching_does_not_change_the_result(monkeypatch):
    """The query batch size is a memory bound, never a scoring parameter."""
    queries = ["CASSLGQAYEQYF", "CASSIVGGNEQFF", "CASSMRSTGELFF"]
    qv = ["TRBV19*01"] * 3
    whole = _score(queries, qv)
    monkeypatch.setattr(geometry_detect, "CELL_BUDGET", 4)
    split = _score(queries, qv)
    assert np.array_equal(whole[0], split[0])
    assert np.allclose(whole[1], split[1])


def test_cli_defaults_to_the_shipped_geometry():
    import subprocess

    out = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parents[2]
                             / "bench" / "geometry_detect.py"), "--help"],
        capture_output=True, text=True, timeout=120)
    assert out.returncode == 0
    assert "--stem-depth" in out.stdout and "--compare" in out.stdout


def test_missing_test_data_fails_with_the_fix_not_a_traceback(tmp_path, monkeypatch):
    import subprocess

    env = dict(**{k: v for k, v in __import__("os").environ.items()},
               VDJMATCH_TESTDATA=str(tmp_path))
    out = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parents[2]
                             / "bench" / "geometry_detect.py"), "--loci", "TRB"],
        capture_output=True, text=True, env=env, timeout=600)
    assert out.returncode != 0
    assert "VDJMATCH_TESTDATA" in (out.stderr + out.stdout)
    assert "Traceback" not in out.stderr
