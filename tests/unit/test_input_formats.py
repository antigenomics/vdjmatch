"""Synthetic shared annotation input formats preserve pairing and source identities."""

import polars as pl
import pytest

from vdjmatch import io


def write(tmp_path, frame, name="input.tsv"):
    path = tmp_path / name
    frame.write_csv(path, separator="\t")
    return path


@pytest.mark.parametrize("identity", ["cloneid", "clonotype_id"])
def test_simple_identity_aliases_and_metadata(tmp_path, identity):
    path = write(tmp_path, pl.DataFrame({identity: ["007", "008"], "locus": ["TRB"] * 2,
        "v": ["TRBV1", "TRBV2"], "j": ["TRBJ1"] * 2,
        "junction_aa": ["CASSF", None], "note": ["a", "b"]}))
    raw, report = io.normalize_annotation_input(path)
    assert raw[identity].to_list() == ["007", "008"]
    assert raw["clone_id"].to_list() == ["007", "008"]
    assert raw["v_call"].to_list() == ["TRBV1", "TRBV2"]
    assert raw["note"].to_list() == ["a", "b"]
    assert raw.height == 2 and report["missing_sequences"] == 1
    normalized = io.read_rearrangement(path)
    assert normalized[identity].to_list() == ["007"]
    assert normalized["pair_id"].to_list() == ["007"]


def wide_frame():
    return pl.DataFrame({"query_id": ["q007", "q008", "q009"], "cloneid": ["007", "008", "009"],
        "sequence_id": ["0001", "0002", "0003"], "species": ["human"] * 3, "note": ["a", "b", "c"],
        "junction_aa_alpha": ["CAVF", None, None], "junction_aa_beta": ["CASSF", "CASAF", None],
        "v_alpha": ["TRAV1", None, None], "j_alpha": ["TRAJ1", None, None],
        "v_beta": ["TRBV1", "TRBV2", None], "j_beta": ["TRBJ1", "TRBJ2", None]})


def test_wide_expands_source_rows_without_losing_metadata_or_ids(tmp_path):
    raw, report = io.normalize_annotation_input(write(tmp_path, wide_frame()))
    assert raw.height == 6 and raw["query_id"].n_unique() == 6
    assert raw["clone_id"].to_list() == ["007", "007", "008", "008", "009", "009"]
    assert raw["locus"].to_list() == ["TRA", "TRB"] * 3
    assert raw["junction_aa"].to_list() == ["CAVF", "CASSF", None, "CASAF", None, None]
    assert raw["source_query_id"].to_list() == ["q007", "q007", "q008", "q008", "q009", "q009"]
    assert raw["source_sequence_id"].to_list() == ["0001", "0001", "0002", "0002", "0003", "0003"]
    assert raw["note"].to_list() == ["a", "a", "b", "b", "c", "c"]
    assert report["source_rows"] == 3 and report["input_rows"] == report["output_rows"] == 6
    assert report["missing_sequences"] == 3


def test_wide_readers_keep_incomplete_pairs(tmp_path):
    path = write(tmp_path, wide_frame())
    cells = io.read_cell(path, link="clone_id", source="airr")
    assert cells["pair_id"].to_list() == ["007", "008", "009"]
    assert cells["pair_status"].to_list() == ["paired", "missing_alpha", "missing_alpha"]
    assert cells["cdr3a"].to_list() == ["CAVF", None, None]
    assert cells["cdr3b"].to_list() == ["CASSF", "CASAF", None]
    assert cells["source_query_ida"].to_list() == ["q007", None, None]
    assert cells["source_query_idb"].to_list() == ["q007", "q008", None]
    assert io.read_rearrangement(path)["cdr3"].to_list() == ["CAVF", "CASSF", "CASAF"]


def test_wide_generated_links_and_empty_input(tmp_path):
    frame = wide_frame().drop("cloneid")
    raw, _ = io.normalize_annotation_input(write(tmp_path, frame))
    assert raw["clone_id"].n_unique() == 3
    cells = io.read_cell(write(tmp_path, frame))
    assert cells.height == 3
    empty, report = io.normalize_annotation_input(write(tmp_path, frame.head(0)))
    assert empty.height == report["output_rows"] == 0
    assert io.read_cell(write(tmp_path, frame.head(0))).height == 0


@pytest.mark.parametrize("field", ["cloneid", "query_id"])
def test_wide_duplicate_or_missing_identity_rejected(tmp_path, field):
    frame = wide_frame().with_columns(pl.lit("same").alias(field))
    with pytest.raises(ValueError, match="unique"):
        io.normalize_annotation_input(write(tmp_path, frame))
    frame = wide_frame().with_columns(pl.lit(None, dtype=pl.String).alias(field))
    with pytest.raises(ValueError, match="nonempty|non-null"):
        io.normalize_annotation_input(write(tmp_path, frame))


def test_no_anchor_fabrication_or_ambiguous_wide_selection(tmp_path):
    bare = pl.DataFrame({"cdr3_aa_alpha": ["AV"], "cdr3_aa_beta": ["ASS"], "cloneid": ["007"]})
    with pytest.raises(ValueError, match="junction"):
        io.normalize_annotation_input(write(tmp_path, bare))
    mixed = wide_frame().with_columns(pl.lit("CASSF").alias("junction_aa"))
    with pytest.raises(ValueError, match="mixed|ambiguous"):
        io.normalize_annotation_input(write(tmp_path, mixed))
    airr = pl.DataFrame({"sequence_id": ["001"], "v_call": ["TRBV1"], "cdr3_aa": ["ASS"]})
    with pytest.raises(ValueError, match="junction_aa"):
        io.normalize_annotation_input(write(tmp_path, airr))


def test_airr_long_cloneid_and_raw_invalid_rows(tmp_path):
    frame = pl.DataFrame({"cloneid": ["007", "007"], "sequence_id": ["001", "002"],
        "locus": ["TRA", "TRB"], "junction_aa": ["CAVF", "CA*F"],
        "v_call": ["TRAV1", "TRBV1"], "j_call": ["TRAJ1", "TRBJ1"], "note": ["a", "b"]})
    path = write(tmp_path, frame)
    raw, report = io.normalize_annotation_input(path)
    assert raw["sequence_id"].to_list() == ["001", "002"] and raw.height == 2
    assert report["invalid_sequences"] == 1
    cells = io.read_cell(path)
    assert cells["pair_id"].to_list() == ["007"] and cells["pair_status"].to_list() == ["missing_beta"]
    assert io.read_rearrangement(path)["cloneid"].to_list() == ["007"]
    assert io.read_cell(path, link="clone_id")["pair_id"].to_list() == ["007"]


def test_conflicting_link_aliases_rejected(tmp_path):
    frame = wide_frame().with_columns(pl.lit("other").alias("clone_id"))
    with pytest.raises(ValueError, match="conflicting"):
        io.normalize_annotation_input(write(tmp_path, frame))


def test_cell_and_clone_id_retain_distinct_counting_units(tmp_path):
    frame = pl.DataFrame({"cell_id": ["cell-a", "cell-b"], "clone_id": ["clone-a"] * 2,
                          "locus": ["TRA", "TRB"], "junction_aa": ["CAVF", "CASSF"]})
    raw, _ = io.normalize_annotation_input(write(tmp_path, frame))
    assert raw["cell_id"].to_list() == ["cell-a", "cell-b"]
    assert raw["clone_id"].to_list() == ["clone-a", "clone-a"]


def test_clonotype_id_and_clone_id_are_independent(tmp_path):
    frame = wide_frame().rename({"cloneid": "clone_id"}).with_columns(
        pl.Series("clonotype_id", ["clonotype-x", "clonotype-y", "clonotype-z"]))
    path = write(tmp_path, frame)
    raw, _ = io.normalize_annotation_input(path)
    assert raw["clone_id"].to_list() == ["007", "007", "008", "008", "009", "009"]
    assert raw["clonotype_id"].to_list() == ["clonotype-x"] * 2 + ["clonotype-y"] * 2 + ["clonotype-z"] * 2
    assert io.read_rearrangement(path)["clonotype_id"].to_list() == ["clonotype-x", "clonotype-x", "clonotype-y"]


def test_wide_source_controls_counts_and_lane_ids(tmp_path):
    frame = wide_frame().with_columns(pl.Series("sequence_id_alpha", ["a007", "a008", "a009"]),
                                    pl.Series("sequence_id_beta", ["b007", "b008", "b009"]))
    raw, _ = io.normalize_annotation_input(write(tmp_path, frame))
    assert raw["sequence_id"].to_list() == ["a007", "b007", "a008", "b008", "a009", "b009"]
    with pytest.raises(ValueError, match="source"):
        io.normalize_annotation_input(write(tmp_path, frame), source="unknown")
    for count in ["1.5", "-1", None]:
        with pytest.raises(ValueError, match="counts"):
            io.normalize_annotation_input(write(tmp_path, frame.with_columns(pl.lit(count).alias("count"))))


def test_wide_cells_can_share_clone_identity(tmp_path):
    frame = wide_frame().head(2).with_columns(pl.lit("clone-x").alias("cloneid"),
        pl.Series("cell_id", ["cell-a", "cell-b"]), pl.Series("pair_id", ["pair-a", "pair-b"]))
    raw, _ = io.normalize_annotation_input(write(tmp_path, frame))
    assert raw["clone_id"].to_list() == ["clone-x"] * 4
    assert raw["cell_id"].to_list() == ["cell-a", "cell-a", "cell-b", "cell-b"]
    assert raw["pair_id"].to_list() == ["pair-a", "pair-a", "pair-b", "pair-b"]
    cells = io.read_cell(write(tmp_path, frame))
    assert cells["pair_id"].to_list() == ["pair-a", "pair-b"]
    assert cells["clone_idb"].to_list() == ["clone-x"] * 2
