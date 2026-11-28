import numpy as np
import polars as pl

from vdjmatch.cli.__main__ import main


def test_vx_fixed_target_excludes_exact_junction_hits_only(tmp_path):
    rectangle = tmp_path / "rectangle.npz"
    np.savez(rectangle, query_id=["q", "no-support"], reference_key=["r1", "r2"],
             native=[[0, 12], [np.inf, np.inf]],
             query_junction=[["CASSYF"], ["CASLYF"]], reference_junction=[["CASSYF"], ["CASSLF"]],
             query_v=[["TRBV19*01"], ["TRBV19*01"]], reference_v=[["TRBV19*01"], ["TRBV19*01"]],
             background=[np.inf, np.inf], source_log2=[0, 0])
    prefix = tmp_path / "result"
    assert main(["vx-fixed-target", str(rectangle), "--locus", "TRB", "--exclude-exact",
                 "--output-prefix", str(prefix)]) == 0
    rows = pl.read_csv(prefix.with_suffix(".scores.tsv"), separator="\t")
    assert rows.rows() == [("q", -12.0, -12.0), ("no-support", None, None)]


def test_mock_panel_shortlist_preserves_inclusive_cutoff_and_empty_panel():
    from vdjmatch.aggregate.vx_panels import _PanelScorer
    from vdjmatch.cli.paired_vx_base import _regional_corrections
    from vdjmatch.match.geometry import precompute_v_loop_distances

    context = precompute_v_loop_distances(species="human")
    geometry = _regional_corrections(context, "TRB")
    allele = geometry[0][0]
    panels = _PanelScorer(
        ["CASSTLGYF", "CASWF"], [allele] * 2,
        [np.array([0]), np.array([1]), np.array([], dtype=int)],
        [np.zeros(1), np.zeros(1), np.zeros(0)],
        (["CASRRQGYF"], [allele]), context, geometry, 0, 1,
    )
    calls = np.zeros((3, 3), dtype=bool)
    for indices, passed in panels.batches(["CASSTLGYF", "CASWF", "CASLYF"], [allele] * 3):
        calls[indices] = passed
    np.testing.assert_array_equal(calls, [
        [True, False, False], [False, True, False], [False, False, False],
    ])
