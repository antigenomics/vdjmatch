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
    assert c["E"][0] == 1.0 and c["p_enrichment"][0] > 0 and c["rule_of_three"][0]


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


def _paired_detail_reference():
    return pl.DataFrame(
        {
            "gene": ["TRA", "TRB", "TRA", "TRB"],
            "cdr3": ["CAVVF", "CASSF", "CAVVF", "CASSF"],
            "epitope": ["ONE"] * 4,
            "complex_id": ["complex-one", "complex-one", "complex-two", "complex-two"],
            "record_id": ["complex-one", "complex-one", "complex-two", "complex-two"],
            "reference_id": ["study-one", "study-one", "study-two", "study-two"],
            "source_note": ["a1", "b1", "a2", "b2"],
        }
    )


def test_paired_details_preserve_metadata_without_repeated_search(monkeypatch):
    a = Annotator.from_frame(_paired_detail_reference())
    q = pl.DataFrame(
        {
            "query_id": ["query-one"],
            "cdr3_alpha_aa": ["CAVVF"],
            "cdr3_beta_aa": ["CASSF"],
        }
    )
    calls = []
    original = a._index.annotate

    def counted(*args, **kwargs):
        calls.append((kwargs["gene"], kwargs["align"]))
        return original(*args, **kwargs)

    monkeypatch.setattr(a._index, "annotate", counted)
    hits, c = a.paired_candidates(q, scope="0", align=True, return_hits=True)
    assert calls == [("TRA", True), ("TRB", True)]
    assert hits["complex_id"].to_list() == ["complex-one", "complex-two"]
    assert hits["alpha_record_id"].to_list() == ["complex-one", "complex-two"]
    assert hits["beta_record_id"].to_list() == ["complex-one", "complex-two"]
    assert hits["alpha_reference_id"].to_list() == ["study-one", "study-two"]
    assert hits["beta_source_note"].to_list() == ["b1", "b2"]
    assert hits["query_id"].to_list() == ["query-one"] * 2
    for chain in ("alpha", "beta"):
        assert hits[chain + "_n_subs"].to_list() == [0, 0]
        assert hits[chain + "_score"].to_list() == [0, 0]
        assert hits[chain + "_cigar"].null_count() == 0
        assert hits[chain + "_match"].null_count() == 0
    assert c["n_hits"][0] == 1 and c["n_records"][0] == 2
    assert c.equals(a.paired_candidates(q, scope="0"))
    empty_q = q.with_columns(pl.lit("CWWWF").alias("cdr3_beta_aa"))
    empty_hits, empty_c = a.paired_candidates(
        empty_q, scope="0", align=True, return_hits=True
    )
    assert empty_hits.height == empty_c.height == 0
    assert empty_hits.schema == hits.schema
    assert empty_c.schema == c.schema


def test_paired_unit_cost_scale_and_annotate_forwarding():
    import math
    import pytest

    a = Annotator.from_frame(_paired_detail_reference())
    q = pl.DataFrame({"cdr3_alpha_aa": ["CAVAF"], "cdr3_beta_aa": ["CASSF"]})
    scope = search_params("1", gap_open=1, gap_extend=1)
    hits, c = a.paired_candidates(q, scope=scope, score_scale=1, return_hits=True)
    assert hits["alpha_score"].to_list() == [1, 1]
    assert c["ned_score"][0] == pytest.approx(math.exp(-1))
    out = a.annotate_paired(q, scope=scope, score_scale=1, align=True)
    assert out["vdjmatch_score"][0] == c["ned_score"][0]
    for invalid in (0, -1, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="score_scale"):
            a.paired_candidates(q, score_scale=invalid)


def test_paired_runner_details_scale_and_unsupported_flags(tmp_path):
    import math
    import pytest
    from vdjmatch.runner.multisample import annotate_sample

    a = Annotator.from_frame(_paired_detail_reference())
    path = tmp_path / "paired.tsv"
    pl.DataFrame(
        {
            "pair_id": ["cell-one", "cell-one"],
            "locus": ["TRA", "TRB"],
            "cdr3": ["CAVAF", "CASSF"],
        }
    ).write_csv(path, separator="\t")
    res = annotate_sample(
        a._index,
        path,
        scope="1",
        matrix=None,
        with_evalue=False,
        paired=True,
        align=True,
    )
    assert res["hits"].height == 2
    assert "alpha_record_id" in res["hits"].columns
    assert "beta_cigar" in res["hits"].columns
    assert res["candidates"]["ned_score"][0] == pytest.approx(math.exp(-1))
    for flag in ("match_v", "match_j"):
        with pytest.raises(ValueError, match="paired matching does not support"):
            annotate_sample(a._index, path, paired=True, **{flag: True})


def test_hf_asset_forwarded_and_native_progress(monkeypatch, capsys):
    from vdjmatch import db

    seen = {}

    def fetch(**kw):
        seen.update(kw)
        return "example.tsv"

    monkeypatch.setattr(db, "fetch_hf", fetch)
    monkeypatch.setattr(Annotator, "from_path", classmethod(lambda cls, *a, **kw: None))
    Annotator.latest(source="hf", asset="slim", pin="tag")
    assert seen["asset"] == "slim" and seen["tag"] == "tag"
    Annotator.from_frame(reference()).hits(["CASSF"], progress=True)
    assert "one native batch" in capsys.readouterr().err


def test_standard_calibration_resolves_reference_species(monkeypatch):
    import pytest
    from vdjmatch.evalue import control as controls

    requests = []
    monkeypatch.setattr(
        controls,
        "background",
        lambda locus, species: (
            requests.append((locus, species)) or Index.build(["CASSF", "CAVVF"], "aa")
        ),
    )
    mouse = reference().with_columns(pl.lit("MusMusculus").alias("species"))
    ann = Annotator.from_frame(mouse)
    ann.annotate(["CASSF"], scope="0", calibrate=True)
    assert requests == [("TRB", "mouse")]
    with pytest.raises(ValueError, match="does not match"):
        ann.candidates(["CASSF"], scope="0", calibrate=True, species="human")
    mixed = mouse.with_columns(
        pl.Series("species", ["MusMusculus", "HomoSapiens", "MusMusculus"])
    )
    with pytest.raises(ValueError, match="single reference species"):
        Annotator.from_frame(mixed).annotate(["CASSF"], calibrate=True, species="human")
    unknown = Annotator.from_frame(reference().drop("species"))
    with pytest.raises(ValueError, match="explicit control species"):
        unknown.annotate(["CASSF"], calibrate=True)
    unknown.annotate(["CASSF"], scope="0", calibrate=True, species="human")
    assert requests[-1] == ("TRB", "human")
    other_locus = (
        reference()
        .head(1)
        .with_columns(
            pl.lit("TRA").alias("gene"),
            pl.lit("CAVVF").alias("cdr3"),
            pl.lit("alpha").alias("record_id"),
        )
    )
    Annotator.from_frame(pl.concat([mouse, other_locus])).annotate(
        ["CASSF"],
        scope="0",
        calibrate=True,
    )
    assert requests[-1] == ("TRB", "mouse")
    # Supplied custom controls have no encoded organism and remain caller-owned.
    unknown.annotate(["CASSF"], scope="0", control=Index.build(["CASSF"], "aa"))

    paired = Annotator.from_frame(
        _paired_detail_reference().with_columns(pl.lit("MusMusculus").alias("species"))
    )
    q = pl.DataFrame({"cdr3_alpha_aa": ["CAVVF"], "cdr3_beta_aa": ["CASSF"]})
    requests.clear()
    paired.annotate_paired(q, scope="0", calibrate=True)
    assert requests == [("TRA", "mouse"), ("TRB", "mouse")]
    with pytest.raises(ValueError, match="does not match"):
        paired.paired_candidates(q, calibrate=True, species="human")
    with pytest.raises(ValueError, match="explicit control species"):
        Annotator.from_frame(_paired_detail_reference()).paired_candidates(
            q, calibrate=True
        )


def test_explicit_control_mappings_require_active_loci(monkeypatch):
    import pytest
    from vdjmatch.evalue import control as controls

    def unexpected(*args):
        raise AssertionError("incomplete explicit mapping must not fetch controls")

    monkeypatch.setattr(controls, "background", unexpected)
    ann = Annotator.from_frame(reference())
    ctrl = Index.build(["CASSF"], "aa")
    for calibrate in (False, True):
        for supplied in ({}, {"TRA": ctrl}):
            with pytest.raises(ValueError, match="every active locus"):
                ann.annotate(["CASSF"], control=supplied, calibrate=calibrate)
        paired = Annotator.from_frame(_paired_detail_reference())
        q = pl.DataFrame({"cdr3_alpha_aa": ["CAVVF"], "cdr3_beta_aa": ["CASSF"]})
        with pytest.raises(ValueError, match="both TRA and TRB"):
            paired.paired_candidates(q, control={"TRB": ctrl}, calibrate=calibrate)
