"""Unified route geometry and evidence against an independent scalar oracle."""
import math
import polars as pl
from seqtree import Index, SubstitutionMatrix, gapblock
from vdjmatch.api import Annotator


def test_unified_routes_counts_identity_and_parallel_equality():
    short = "CASSLGQAYEQYF"
    long = "CASSLGQAYEQYSSSSSSSSSSSF"
    refs = [short, short.replace("Q", "R"), long, long.replace("Q", "R")]
    ctrl = [short, long, "CYYYYYYYYYYYYYYYYYYYYYYF"]
    ann = Annotator.from_frame(pl.DataFrame({"gene": ["TRB"] * 4,
        "cdr3": refs, "epitope": ["E"] * 4}))
    q = pl.DataFrame({"query_id": ["a", "b", "repeat"], "cdr3": [short, long, long], "locus": ["TRB"] * 3})
    control = Index.build(ctrl)
    one = ann.unified_candidates(q, control=control, threads=1, exclude_exact=True, distance="gapblock")
    many = ann.unified_candidates(q, control=control, threads=2, exclude_exact=True, distance="gapblock")
    assert one.equals(many)
    assert set(one["query_id"]) == {"a", "b", "repeat"}
    row = one.filter(pl.col("query_id") == "b").row(0, named=True)
    matrix = SubstitutionMatrix.blosum62()
    kw = dict(matrix=matrix, gap_open=2*matrix.scale(), gap_extend=matrix.scale(),
              gap_prior=gapblock.positions_prior((3, 4, -4, -3)))
    distances = [(r, gapblock.gapblock_score(long, r, **kw)[0]) for r in refs if r != long]
    background = [gapblock.gapblock_score(long, r, **kw)[0] for r in ctrl if r != long]
    retained = [d for r,d in distances if d <= 5*matrix.scale()]
    expected = len(ctrl)/len(refs) * sum(1/(1+sum(b <= d for b in background)) for d in retained)
    assert math.isclose(row["ned_score"], expected)
    assert row["n_target"] == len(retained)
    assert row["n_control"] == sum(b <= 5*matrix.scale() for b in background)
    assert row["estimator"] == "background-mass-v1"
    assert row["nearest_edits"] is None


def test_duplicate_controls_rejected_and_metadata_stable():
    import pytest
    seqs = ["CASSF", "CASSLGQAYEQYSSSSSSSSSSSF"]
    ann = Annotator.from_frame(pl.DataFrame({"gene": ["TRB"]*2, "cdr3": seqs,
        "epitope": ["E"]*2, "sequence_id": ["r1", "r2"], "pair_id": ["rp1", "rp2"]}))
    q = pl.DataFrame({"cdr3": seqs, "locus": ["TRB"]*2,
        "sequence_id": ["q1", "q2"], "pair_id": ["qp1", "qp2"]})
    with pytest.raises(ValueError, match="unique junctions"):
        ann.unified_candidates(q.head(1), control=Index.build([seqs[0]]*2))
    hits, _ = ann.unified_candidates(q, control=Index.build(seqs), return_hits=True, distance="gapblock")
    exact = hits.filter(pl.col("query_cdr3") == pl.col("db_cdr3")).sort("query_id")
    assert exact["query_sequence_id"].to_list() == ["q1", "q2"]
    assert exact["query_pair_id"].to_list() == ["qp1", "qp2"]
    assert exact["sequence_id"].to_list() == ["r1", "r2"]
    assert exact["pair_id"].to_list() == ["rp1", "rp2"]
