"""A graded query searches once and preserves each fixed-radius result."""

import polars as pl
from seqtree import Index

from vdjmatch.api import Annotator


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
