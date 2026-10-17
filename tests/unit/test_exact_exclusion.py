"""Puncture target/control zero-edit matches consistently before candidate reduction."""

import polars as pl
import pytest
from seqtree import Index

from vdjmatch import Annotator
from vdjmatch.match import search_params
from vdjmatch.runner.multisample import annotate_sample


def annotator():
    return Annotator.from_frame(
        pl.DataFrame(
            {
                "gene": ["TRB"] * 3,
                "cdr3": ["CASSF", "CASSF", "CATSF"],
                "epitope": ["PEP"] * 3,
                "record_id": ["exact1", "exact2", "neighbour"],
                "species": ["HomoSapiens"] * 3,
            }
        )
    )


def test_target_and_control_exact_hits_are_removed_before_reduction():
    ann = annotator()
    ctrl = Index.build(["CASSF", "CASRF"], "aa")
    default = ann.candidates(["CASSF"], scope="1", control=ctrl)
    assert default["n_target"].to_list() == [2]
    assert default["n_control"].to_list() == [2]
    excluded = ann.candidates(["CASSF"], scope="1", control=ctrl, exclude_exact=True)
    row = excluded.row(0, named=True)
    assert (
        row["n_target"],
        row["n_control"],
        row["nearest_edits"],
        row["n_records"],
    ) == (1, 1, 1, 1)
    assert row["n_reference"] == 2 and row["E"] == 1.0
    hits = ann.hits(["CASSF"], scope="1", exclude_exact=True, align=True)
    assert hits["db_cdr3"].to_list() == ["CATSF"]
    assert hits["cigar"].null_count() == 0
    out = ann.annotate(["CASSF"], scope="0", exclude_exact=True)
    assert out["vdjmatch_status"].to_list() == ["no_hit"]


def test_runner_forwards_exact_exclusion_and_preserves_typed_empty_hits(tmp_path):
    ann = annotator()
    query = tmp_path / "query.tsv"
    pl.DataFrame({"junction_aa": ["CASSF"], "locus": ["TRB"]}).write_csv(
        query, separator="\t"
    )
    ctrl = Index.build(["CASSF", "CASRF"], "aa")
    params = search_params("1")
    ordinary = ann._index.annotate(
        pl.DataFrame({"cdr3": ["CASSF"]}), params, gene="TRB"
    )
    empty = ann._index.annotate(
        pl.DataFrame({"cdr3": ["CASSF"]}),
        search_params("0"),
        gene="TRB",
        exclude_exact=True,
    )
    assert not empty.height and empty.schema == ordinary.schema
    result = annotate_sample(
        ann._index, query, scope="1", control=ctrl, exclude_exact=True
    )
    assert result["hits"]["db_cdr3"].to_list() == ["CATSF"]
    assert result["candidates"]["n_control"].to_list() == [1]
    assert result["calls"]["vdjmatch_status"].to_list() == ["matched"]


def test_paired_runner_forwards_joint_puncture(tmp_path):
    from vdjmatch.match import search_params
    ref = pl.DataFrame({"gene": ["TRA", "TRB", "TRA", "TRB"],
                        "cdr3": ["CAVVF", "CASSF", "CAVAF", "CASSF"],
                        "epitope": ["PEP"] * 4, "complex_id": [1,1,2,2]})
    ann = Annotator.from_frame(ref)
    query = tmp_path / "paired.tsv"
    pl.DataFrame({"junction_aa": ["CAVVF", "CASSF"], "locus": ["TRA", "TRB"],
                  "cell_id": ["cell", "cell"]}).write_csv(query, separator="\t")
    ctrl = {"TRA": Index.build(["CAVVF", "CAVAF"], "aa"),
            "TRB": Index.build(["CASSF", "CASRF"], "aa")}
    result = annotate_sample(ann._index, query, scope="1", paired=True,
                             exclude_exact=True, control=ctrl)
    assert result["hits"]["complex_id"].to_list() == ["2"]
    assert result["candidates"]["n_control_joint"].to_list() == [3]
    assert result["calls"]["vdjmatch_status"].to_list() == ["matched"]


def test_seqtrie_exact_filter_uses_identity_not_zero_edit_metadata():
    ann = annotator()
    ctrl = Index.build(["CASSF", "CASRF"], "aa")
    candidates = ann.candidates(["CASSF"], scope=search_params("1", engine="seqtrie"),
                                control=ctrl, exclude_exact=True)
    assert candidates["n_target"].to_list() == [1]
    assert candidates["n_control"].to_list() == [1]
