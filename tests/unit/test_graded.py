"""A graded query searches once and preserves each fixed-radius result."""

import polars as pl
import pytest
from seqtree import Index

from vdjmatch.api import Annotator


@pytest.mark.parametrize("weighted", [False, True])
@pytest.mark.parametrize("match_v,match_j,expected", [(True, False, 2), (False, True, 2), (True, True, 1)])
def test_graded_gene_predicates_equal_fixed_and_reject_sequence_only_calibration(weighted, match_v, match_j, expected):
    from vdjmatch.match.scoring import load_vdjam
    from vdjmatch.match import search_params

    ref = pl.DataFrame({"gene": ["TRB"] * 3, "cdr3": ["CASSF", "CASRF", "CASTF"],
                        "v": ["TRBV19", "TRBV5-1", "TRBV19"],
                        "j": ["TRBJ2-7", "TRBJ2-7", "TRBJ1-1"], "epitope": ["E"] * 3})
    query = pl.DataFrame({"cdr3": ["CASSF"], "v": ["TRBV19*01"],
                          "j": ["TRBJ2-7*01"], "locus": ["TRB"]})
    ann = Annotator.from_frame(ref)
    matrix = load_vdjam() if weighted else None
    got = ann.graded_candidates(query, match_v=match_v, match_j=match_j, matrix=matrix)
    for radius in range(1, 6):
        gap = 100 if weighted else 1
        params = search_params(f"{radius},{min(radius,2)},{min(radius,2)},{radius}",
                               matrix=matrix or "", gap_open=gap, gap_extend=gap)
        fixed = ann.candidates(query, scope=params, match_v=match_v, match_j=match_j,
                               score_scale=400.0 if weighted else 1.0)
        assert got.filter(pl.col("radius") == radius).drop("radius").equals(fixed)
        assert fixed["n_clonotypes"].to_list() == [expected]
        assert fixed["n_reference"].to_list() == [expected]
    with pytest.raises(ValueError, match="V/J-labelled controls"):
        ann.graded_candidates(query, match_v=match_v, match_j=match_j, matrix=matrix,
                              control=Index.build(["CASSF"]))


def test_graded_balls_equal_five_fixed_searches(monkeypatch):
    ref = pl.DataFrame({"gene": ["TRB"] * 5,
                        "cdr3": ["CAAAAAAF", "CAAAAACF", "CAAAACCF", "CAACCCCF", "CACCCCCF"],
                        "epitope": ["E"] * 5})
    ann = Annotator.from_frame(ref)
    control = Index.build(["CAAAAAAF", "CAAAACCF", "CYYYYYYF"])
    expected = [ann.candidates(["CAAAAAAF"], scope=f"{r},{min(r,2)},{min(r,2)},{r}",
                               control=control, soft_v=False) for r in range(1, 6)]
    # The native batch is invoked once for target and once for control.
    original = Index.search_batch
    calls = []
    def batch(self, queries, params, threads=0):
        calls.append((len(queries), params.max_total_edits))
        return original(self, queries, params, threads)
    monkeypatch.setattr(Index, "search_batch", batch)
    got = ann.graded_candidates(["CAAAAAAF"], control=control)
    assert calls == [(1, 5), (1, 5)]
    for r, fixed in enumerate(expected, 1):
        row = got.filter(pl.col("radius") == r)
        assert row["n_target"].to_list() == fixed["n_target"].to_list()
        assert row["n_control"].to_list() == fixed["n_control"].to_list()
        assert row["p_enrichment"].to_list() == fixed["p_enrichment"].to_list()


def test_weighted_radius_keeps_substitution_hidden_by_cheaper_wide_gap():
    from vdjmatch.match.scoring import load_vdjam
    ann = Annotator.from_frame(pl.DataFrame({"gene": ["TRB"], "cdr3": ["AWA"], "epitope": ["E"]}))
    hits, c = ann.graded_candidates(["ACA"], matrix=load_vdjam(), return_hits=True)
    one = hits.filter(pl.col("radius") == 1)
    assert one.height == 1 and one["n_subs"][0] == 1
    assert c.filter(pl.col("radius") == 1)["n_target"][0] == 1
    assert hits.filter(pl.col("radius") == 5)["n_ins"][0] == 1
