"""Synthetic ingestion contracts; no reference downloads or scientific runs."""

import gzip

import polars as pl
import pytest

from vdjmatch.io import columns, read_cell, read_rearrangement


def write_table(tmp_path, frame, name="query.tsv"):
    path = tmp_path / name
    frame.write_csv(path, separator="\t")
    return path


def test_airr_prefers_junction_and_keeps_identity(tmp_path):
    frame = pl.DataFrame(
        {
            "sequence_id": ["001", "002", "003"],
            "junction_aa": ["CASSF", "CASSF", "CASSF"],
            "cdr3_aa": ["ASS", "ASS", "ASS"],
            "v_call": ["TRBV1", "TRBV1", "TRBV2"],
            "j_call": ["TRBJ1"] * 3,
            "duplicate_count": ["2", "3", "4"],
        }
    )
    path = write_table(tmp_path, frame)
    out = read_rearrangement(path)
    assert out["cdr3"].to_list() == ["CASSF"] * 3
    assert out["query_id"].to_list() == [0, 1, 2]
    assert out["sequence_id"].to_list() == ["001", "002", "003"]
    assert out["count"].to_list() == [2, 3, 4]
    merged, report = read_rearrangement(path, dedup=True, return_report=True)
    assert merged["count"].to_list() == [5, 4]
    assert merged["query_ids"].to_list() == [[0, 1], [2]]
    assert merged["sequence_ids"].to_list() == [["001", "002"], ["003"]]
    assert report["retained_rows"] == 3 and report["deduplicated_rows"] == 1


def test_sequence_conventions_are_explicit(tmp_path):
    airr = pl.DataFrame({"cdr3_aa": ["ASS"], "v_call": ["TRBV1"], "j_call": ["TRBJ1"]})
    with pytest.raises(ValueError, match="junction_aa"):
        read_rearrangement(write_table(tmp_path, airr))
    custom = pl.DataFrame({"cdr3_aa": ["CASSF"], "v_gene": ["TRBV1"]})
    with pytest.raises(ValueError, match="sequence_convention"):
        columns.normalize_query(custom, source="custom")
    assert columns.normalize_query(custom)["cdr3"].to_list() == ["CASSF"]
    assert (
        columns.normalize_query(
            custom, source="custom", sequence_convention="junction"
        ).height
        == 1
    )
    with pytest.raises(ValueError, match="bare CDR3"):
        columns.normalize_query(custom, sequence_convention="cdr3")


@pytest.mark.parametrize(
    "count", ["bad", "-1", "1.5", None, "nan", "9223372036854775808"]
)
def test_bad_supplied_counts_fail(count):
    with pytest.raises(ValueError, match="counts"):
        columns.normalize_query(pl.DataFrame({"cdr3": ["CASSF"], "count": [count]}))


def test_report_filters_sequences_without_losing_original_offsets(tmp_path):
    frame = pl.DataFrame(
        {
            "cdr3": ["CASSF", None, "", "CA*F", "cassw"],
            "v": ["TRBV1"] * 5,
            "count": [1, 2, 3, 4, 5],
        }
    )
    out, report = read_rearrangement(write_table(tmp_path, frame), return_report=True)
    assert out["query_id"].to_list() == [0, 4]
    assert report == {
        "input_rows": 5,
        "missing_sequences": 2,
        "invalid_sequences": 1,
        "retained_rows": 2,
        "dropped_rows": 3,
        "output_rows": 2,
        "deduplicated_rows": 0,
    }


def test_custom_cell_link_and_missing_chain_types(tmp_path):
    frame = pl.DataFrame(
        {
            "junction_aa": ["CASSF", "CAVVF", "CASSW"],
            "v_call": ["TRBV1", "TRAV1", "TRBV2"],
            "j_call": ["TRBJ1", "TRAJ1", "TRBJ2"],
            "cell_id": ["ignored"] * 3,
            "selected_link": ["cell-2", "cell-2", "cell-1"],
            "sequence_id": ["b2", "a2", "b1"],
            "duplicate_count": [4, 5, 6],
        }
    )
    out = read_cell(write_table(tmp_path, frame), link="selected_link")
    assert out["pair_id"].to_list() == ["cell-2", "cell-1"]
    assert out["pair_status"].to_list() == ["paired", "missing_alpha"]
    assert out["sequence_idb"].to_list() == ["b2", "b1"]
    assert out["counta"].to_list() == [5, None]
    assert out.schema["cdr3a"] == pl.String and out.schema["counta"] == pl.Int64


def test_cell_multichain_is_explicit(tmp_path):
    frame = pl.DataFrame(
        {
            "cdr3": ["CAVVF", "CAAAF"],
            "v": ["TRAV1", "TRAV2"],
            "clone_id": ["same", "same"],
        }
    )
    path = write_table(tmp_path, frame)
    with pytest.raises(ValueError, match="multiple chains"):
        read_cell(path)
    out, report = read_cell(path, ambiguity="first", return_report=True)
    assert out["cdr3a"].to_list() == ["CAVVF"]
    assert out["pair_status"].to_list() == ["missing_beta"]
    assert report["ambiguous_cell_locus_groups"] == 1


def test_gzip_csv_separator_and_pair_dedup(tmp_path):
    path = tmp_path / "query.csv.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(
            "cdr3,v,j,cell_id,count\nCASSF,TRBV1,TRBJ1,one,2\nCASSF,TRBV1,TRBJ1,two,3\n"
        )
    out = read_rearrangement(path, dedup=True)
    assert out["count"].to_list() == [2, 3]
    assert out["pair_id"].to_list() == ["one", "two"]


def test_native_repertoire_uses_published_reader(tmp_path, monkeypatch):
    from vdjtools import io as upstream

    path = tmp_path / "native.txt"
    path.write_text(
        "count\tfreq\tcdr3nt\tcdr3aa\tv\td\tj\n2\t1\tTGTGCT\tCASSF\tTRBV1\tTRBD1\tTRBJ1\n"
    )
    calls = []
    original = upstream.read

    def read(*args, **kwargs):
        calls.append(kwargs["fmt"])
        return original(*args, **kwargs)

    monkeypatch.setattr(upstream, "read", read)
    out = read_rearrangement(path)
    assert calls == ["vdjtools"]
    assert out["cdr3"].to_list() == ["CASSF"] and out["count"].to_list() == [2]


def test_counts_remain_exact_above_float_precision():
    out = columns.normalize_query(
        pl.DataFrame({"cdr3": ["CASSF", "CASSW"], "count": ["9007199254740993", "2.0"]})
    )
    assert out["count"].to_list() == [9007199254740993, 2]


def test_airr_keeps_ambiguous_calls_and_declared_locus(tmp_path):
    frame = pl.DataFrame(
        {
            "junction_aa": ["CASSF"],
            "v_call": ["TRBV1,TRBV2"],
            "j_call": ["TRBJ1,TRBJ2"],
            "locus": ["TRB"],
            "sequence_id": ["source"],
        }
    )
    out = read_rearrangement(write_table(tmp_path, frame))
    assert out["v"].to_list() == ["TRBV1,TRBV2"]
    assert out["j"].to_list() == ["TRBJ1,TRBJ2"]
    assert out["locus"].to_list() == ["TRB"]


def test_tcrvdb_retains_original_row_identity_and_reports_drops(tmp_path):
    from vdjmatch.io import read_tcrvdb

    frame = pl.DataFrame(
        {
            "cdr3_alpha_aa": ["CAVVF", "BAD*"],
            "cdr3_beta_aa": ["CASSF", "CASSF"],
            "sequence_id": ["p1", "p2"],
        }
    )
    out, report = read_tcrvdb(write_table(tmp_path, frame), return_report=True)
    assert out["query_id"].to_list() == [0]
    assert out["sequence_id"].to_list() == ["p1"]
    assert report["dropped_rows"] == 1


def test_mixcr_reader_preserves_exact_counts_and_order(tmp_path):
    from vdjmatch.io import read_rearrangement

    p = tmp_path / "mixcr.tsv"
    p.write_text(
        "cloneCount\tallVHitsWithScore\tallJHitsWithScore\tnSeqCDR3\taaSeqCDR3\n"
        "9007199254740993\tTRBV1*01(10)\tTRBJ1-1*01(10)\tTGTGCTTTT\tCAF\n"
        "2\tTRBV2*01(10)\tTRBJ2-1*01(10)\tTGTGGGTTT\tCGF\n"
    )
    q = read_rearrangement(p)
    assert q["count"].to_list() == [9007199254740993, 2]
    assert q["cdr3"].to_list() == ["CAF", "CGF"] and q["query_id"].to_list() == [0, 1]


def test_airr_cell_optional_count_absent_is_one_observation(tmp_path):
    frame = pl.DataFrame(
        {
            "sequence_id": ["a", "b"],
            "cell_id": ["pair", "pair"],
            "junction_aa": ["CAVF", "CASSF"],
            "v_call": ["TRAV1", "TRBV1"],
            "j_call": ["TRAJ1", "TRBJ1"],
            "locus": ["TRA", "TRB"],
        }
    )
    out = read_cell(write_table(tmp_path, frame))
    assert out["counta"].to_list() == [1] and out["countb"].to_list() == [1]
