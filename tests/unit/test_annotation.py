"""M1/M2 contracts: identities, competing pMHC evidence and actual reference pairs."""

import polars as pl
from seqtree import Index
from vdjmatch import Annotator
from vdjmatch.match import VdjdbIndex, search_params


def reference():
    return pl.DataFrame(
        {
            "gene": ["TRB", "TRB", "TRB"],
            "cdr3": ["CASSF", "CASSF", "CASSL"],
            "v": ["TRBV1", "TRBV2", "TRBV1"],
            "j": ["TRBJ1"] * 3,
            "epitope": ["ONE", "TWO", "ONE"],
            "mhc_a": ["A", "B", "A"],
            "mhc_class": ["MHCI"] * 3,
            "species": ["HomoSapiens"] * 3,
            "record_id": ["r1", "r2", "r3"],
            "reference_id": ["s1", "s2", "s1"],
        }
    )


def test_same_junction_different_genes_preserves_queries():
    a = Annotator.from_frame(reference())
    q = pl.DataFrame(
        {
            "id": [5, 3, 8],
            "cdr3": ["CASSF", "CASSF", "CWWWF"],
            "v": ["TRBV1", "TRBV2", "TRBV1"],
        }
    )
    out = a.annotate(q, v="v", match_v=True, scope="0", threads=1)
    assert out["id"].to_list() == [5, 3, 8]
    assert out["vdjmatch_epitope"].to_list() == ["ONE", "TWO", None]
    assert out["vdjmatch_status"].to_list() == ["matched", "matched", "no_hit"]


def test_candidate_evidence_keeps_mhc_competitors_and_records():
    a = Annotator.from_frame(reference())
    q = pl.DataFrame({"cdr3": ["CASSF"], "v": ["TRBV1"]})
    ctrl = Index.build(["CWWWF", "CYYYYF"], "aa")
    c = a.candidates(q, v="v", control={"TRB": ctrl}, scope="1", threads=1)
    assert c.height == 2
    assert set(c["mhc_a"]) == {"A", "B"}
    assert c["n_competing_clonotypes"].min() > 0
    assert c["E"].min() > 0 and c["p_enrichment"].min() > 0
    assert c["n_records"].max() == 2
    assert c["n_clonotypes"].max() == 2
    assert c["estimator"].unique().to_list() == ["ned-v1"]


def test_duplicate_observations_do_not_inflate_ned():
    r = reference()
    q = pl.DataFrame({"cdr3": ["CASSF"]})
    c1 = Annotator.from_frame(r).candidates(q, scope="0", threads=1)
    c2 = Annotator.from_frame(
        pl.concat(
            [r, r.with_columns((pl.col("record_id") + "-extra").alias("record_id"))]
        )
    ).candidates(q, scope="0", threads=1)
    assert c1["ned_score"].to_list() == c2["ned_score"].to_list()
    assert c2["n_records"].to_list() == [2, 2]


def test_empty_nohit_invalid_status_schema():
    a = Annotator.from_frame(reference().head(0))
    q = pl.DataFrame({"cdr3": ["CASSF", None, "bad*"]})
    out = a.annotate(q, threads=1)
    assert out["vdjmatch_status"].to_list() == [
        "no_reference",
        "invalid_query",
        "invalid_query",
    ]
    assert a.annotate(q.head(0), threads=1).schema == out.schema


def test_engine_preserves_metadata_and_query_id():
    idx = VdjdbIndex.build(reference())
    q = pl.DataFrame(
        {
            "query_id": [42],
            "cdr3": ["CASSF"],
            "v": ["TRBV1"],
            "j": ["TRBJ1"],
            "count": [3],
            "locus": ["TRB"],
        }
    )
    h = idx.annotate(q, search_params("0"), gene="TRB", threads=1)
    assert h["query_id"].unique().to_list() == [42]
    assert set(h["record_id"]) == {"r1", "r2"}
    assert h["query_locus"].unique().to_list() == ["TRB"]


def test_paired_requires_same_complex_not_label_agreement():
    r = pl.DataFrame(
        {
            "gene": ["TRA", "TRB", "TRA", "TRB"],
            "cdr3": ["CAVVF", "CASSF", "CAWWF", "CASLF"],
            "epitope": ["ONE"] * 4,
            "complex_id": [1, 1, 2, 2],
        }
    )
    a = Annotator.from_frame(r)
    out = a.annotate_paired(
        pl.DataFrame({"cdr3_alpha_aa": ["CAVVF"], "cdr3_beta_aa": ["CASLF"]}),
        scope="0",
        threads=1,
    )
    assert out["vdjmatch_epitope"][0] is None
    assert out["vdjmatch_status"][0] == "no_joint_hit"


def test_supplied_id_order_and_empty_metadata():
    a = Annotator.from_frame(reference().with_columns(pl.lit("extra").alias("custom")))
    q = pl.DataFrame({"query_id": ["z", "a"], "cdr3": ["CASSF", "CASSL"]})
    c = a.candidates(q, scope="0", threads=1)
    assert set(c["query_id"]) == {"z", "a"}
    assert a.annotate(q, scope="0")["query_id"].to_list() == ["z", "a"]
    assert a.hits(["CWWWF"], scope="0").schema == a.hits(["CASSF"], scope="0").schema


def test_ties_abstain_and_control_predicate_is_explicit():
    import pytest

    a = Annotator.from_frame(reference())
    out = a.annotate(["CASSF"], scope="0", soft_v=False)
    assert out["vdjmatch_status"][0] == "ambiguous"
    assert out["vdjmatch_epitope"][0] is None
    with pytest.raises(ValueError, match="V/J-labelled controls"):
        a.candidates(
            pl.DataFrame({"cdr3": ["CASSF"], "v": ["TRBV1"]}),
            match_v=True,
            control=Index.build(["CWWWF"], "aa"),
        )
    for x in [float("nan"), float("inf"), 0.0, -1.0]:
        with pytest.raises(ValueError, match="finite"):
            a.candidates([], score_scale=x)
    with pytest.raises(ValueError, match="AIRR"):
        a.annotate(
            pl.DataFrame({"cdr3_aa": ["ASS"], "v_call": ["TRBV1"]}),
            cdr3="cdr3_aa",
            sequence_convention="junction",
        )


def test_reference_invalid_sequences_are_retained_and_reported():
    r = pl.concat(
        [
            reference(),
            reference()
            .head(1)
            .with_columns(pl.lit("r4").alias("record_id"), pl.lit("").alias("cdr3")),
        ]
    )
    a = Annotator.from_frame(r)
    assert a.reference_report == {
        "reference_rows": 4,
        "searchable_rows": 3,
        "unsearchable_rows": 1,
    }
    assert a.hits(["CASSF"], scope="0").height == 2


def test_paired_duplicates_do_not_vote_and_finite_controls():
    r = pl.DataFrame(
        {
            "gene": ["TRA", "TRB", "TRA", "TRB"],
            "cdr3": ["CAVVF", "CASSF", "CAVVF", "CASSF"],
            "epitope": ["ONE"] * 4,
            "complex_id": [1, 1, 2, 2],
        }
    )
    q = pl.DataFrame({"cdr3_alpha_aa": ["CAVVF"], "cdr3_beta_aa": ["CASSF"]})
    ctrl = {"TRA": Index.build(["CWWWF"], "aa"), "TRB": Index.build(["CWWWF"], "aa")}
    c = Annotator.from_frame(r).paired_candidates(q, scope="0", control=ctrl)
    assert c["n_hits"][0] == 1 and c["n_records"][0] == 2 and c["ned_score"][0] == 1.0
    assert c["E"][0] > 0 and c["p_enrichment"][0] > 0 and c["rule_of_three"][0]


def test_runner_api_default_scoring_and_summary_units(tmp_path):
    from vdjmatch.runner.multisample import annotate_sample
    from vdjmatch.match.scoring import load_vdjam

    a = Annotator.from_frame(reference())
    q = pl.DataFrame({"cdr3": ["CASSF", "CASSF", "CASAF"], "v": ["TRBV1"] * 3})
    path = tmp_path / "query.tsv"
    q.write_csv(path, separator="\t")
    result = annotate_sample(
        a._index, path, matrix=load_vdjam(), with_evalue=False, threads=1
    )
    api = a.candidates(q, threads=1).sort(["query_id", "rank"])
    runner = result["candidates"].sort(["query_id", "rank"])
    assert api.select("query_id", "epitope", "ned_score").equals(
        runner.select("query_id", "epitope", "ned_score")
    )
    one = result["summary"].filter(pl.col("epitope") == "ONE")
    assert (
        one["unique"][0] == 2 and one["n_query_rows"][0] == 3 and one["reads"][0] == 3
    )


def test_paired_rejects_cardinality_and_empty_controls_before_no_hit():
    import pytest

    r = reference().drop("record_id").with_columns(pl.lit(1).alias("complex_id"))
    q = pl.DataFrame({"cdr3_alpha_aa": ["CAVF"], "cdr3_beta_aa": ["CWWWF"]})
    with pytest.raises(ValueError, match="conflicting|multiple"):
        Annotator.from_frame(r).paired_candidates(q)
    valid = reference().with_columns(pl.lit(0).alias("complex_id"))
    empty = Index.build([], "aa")
    with pytest.raises(ValueError, match="empty"):
        Annotator.from_frame(valid).paired_candidates(
            q, control={"TRA": empty, "TRB": empty}
        )


def test_airr_junction_and_sequence_id_preserved_in_public_api():
    q = pl.DataFrame(
        {
            "sequence_id": ["query-a", "query-b"],
            "junction_aa": ["CASSF", "CWWWF"],
            "v_call": ["TRBV1", "TRBV1"],
            "locus": ["TRB", "TRB"],
        }
    )
    out = Annotator.from_frame(reference()).annotate(q, scope="0")
    assert out["sequence_id"].to_list() == ["query-a", "query-b"]
    assert out["vdjmatch_status"].to_list() == ["matched", "no_hit"]
