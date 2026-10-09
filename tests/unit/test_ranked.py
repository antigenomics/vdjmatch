"""Small full-ranking contracts; no biological cohort computation."""
from types import SimpleNamespace

import numpy as np
import polars as pl
import pytest

from vdjmatch.api import Annotator
from vdjmatch.match.engine import VdjdbIndex
from vdjmatch.match import ranked


class Control:
    def __init__(self, refs):
        self.refs = refs

    def __len__(self):
        return len(self.refs)

    def ref_seqs(self):
        return self.refs.copy()


@pytest.fixture(autouse=True)
def native_contract(monkeypatch, request):
    if request.node.name == "test_actual_native_wrappers":
        return
    def scores(q, r, **kw):
        kw.pop("exclude_exact", None)
        return np.asarray(ranked.gapblock.score_matrix(q, r, **kw))

    def top(q, r, k, exclude_exact=False, **kw):
        mat = scores(q, r, **kw)
        return [[SimpleNamespace(ref_id=j, score=int(mat[i, j])) for j in sorted(range(len(r)), key=lambda j: (mat[i, j], j)) if not exclude_exact or q[i] != r[j]][:k] for i in range(len(q))]

    def counts(q, r, thresholds, exclude_exact=False, **kw):
        mat = scores(q, r, **kw)
        return [[sum(int(mat[i, j]) <= t and (not exclude_exact or q[i] != r[j]) for j in range(len(r))) for t in row] for i, row in enumerate(thresholds)]

    def pairs(qa, qb, ra, rb, k, exclude_exact=False, **kw):
        a, b = scores(qa, ra, **kw), scores(qb, rb, **kw)
        return [[(j, int(a[i, j]), int(b[i, j])) for j in sorted(range(len(ra)), key=lambda j: (max(a[i, j], b[i, j]), j)) if not exclude_exact or (qa[i], qb[i]) != (ra[j], rb[j])][:k] for i in range(len(qa))]

    monkeypatch.setattr(ranked.gapblock, "topk_batch", top, raising=False)
    monkeypatch.setattr(ranked.gapblock, "count_batch", counts, raising=False)
    monkeypatch.setattr(ranked.gapblock, "paired_topk_batch", pairs, raising=False)


def reference():
    return pl.DataFrame({"gene": ["TRB"]*3, "cdr3": ["CASSF", "CASSF", "CASAF"], "epitope": ["E", "E", "F"], "reference_id": ["S1", "S2", "S3"]})


def test_single_unique_keys_expand_records_and_separate_global():
    a = Annotator(VdjdbIndex.build(reference()))
    hits, c, g = ranked.ranked_evidence(a, pl.DataFrame({"query_id": ["Q"], "cdr3": ["CASSF"]}), k=2, control=Control(["CWWWF", "CYYYF"]))
    assert hits.height == 3
    assert g["n_reference"].to_list() == [2]
    assert g["requested_k"].to_list() == [2]
    assert 0 < g["p_global"][0] <= 1
    assert c.filter(pl.col("epitope") == "E")["n_records"][0] == 2
    assert c["p_enrichment"].null_count() == c.height
    assert c["nearest_edits"].null_count() == c.height
    assert c["distance_sum"].null_count() == c.height
    assert c.filter(pl.col("epitope") == "E")["n_exact"][0] == 1


def test_fewer_k_and_exact_exclusion():
    a = Annotator(VdjdbIndex.build(reference()))
    h, _, g = ranked.ranked_evidence(a, ["CASSF"], k=3, exclude_exact=True, control=Control(["CASSF", "CWWWF"]))
    assert h["db_cdr3"].to_list() == ["CASAF"]
    assert g["threshold"][0] is None and g["p_global"][0] == 1
    assert g["calibration"][0] == "finite-sample-rank-v1"
    assert g["control_size"][0] == 2 and g["n_control"][0] is None


def test_exact_rank_retains_original_exposures_and_excludes_both_samples():
    a = Annotator(VdjdbIndex.build(reference()))
    _, _, g = ranked.ranked_evidence(a, ["CASSF"], k=1, exclude_exact=True,
                                    control=Control(["CASSF", "CWWWF"]))
    row = g.row(0, named=True)
    assert row["n_reference"] == row["control_size"] == 2
    assert row["n_control"] == 0 and row["threshold"] == 1
    assert row["p_global"] == pytest.approx(.5)
    assert row["E"] is row["p_upper"] is row["finite_control_delta"] is None


def test_paired_observations_and_original_exposure():
    ref = pl.DataFrame({"complex_id": ["R1", "R1", "R2", "R2", "R3", "R3"], "gene": ["TRA", "TRB"]*3,
                        "cdr3": ["CAVF", "CASSF", "CAVF", "CASSF", "CAYF", "CASAF"], "epitope": ["E"]*4+["F"]*2})
    a = Annotator(VdjdbIndex.build(ref))
    q = pl.DataFrame({"query_id": ["P"], "cdr3a": ["CAVF"], "cdr3b": ["CASSF"]})
    ctrl = {"TRA": Control(["CAVF", "CWWF"]), "TRB": Control(["CASSF", "CWWWF"])}
    h, c, g = ranked.ranked_evidence(a, q, k=1, paired=True, control=ctrl, exclude_exact=True)
    assert h.height == 1 and h["complex_id"][0] == "R3"
    assert g["n_reference"][0] == 2
    assert 0 < g["p_global"][0] <= 1
    assert "alpha_record_id" in h.columns and "beta_record_id" in h.columns
    assert c["p_enrichment"][0] is None


def test_empty_and_invalid_control():
    a = Annotator(VdjdbIndex.build(reference()))
    with pytest.raises(ValueError, match="empty"):
        ranked.ranked_evidence(a, ["CASSF"], control=Control([]))
    with pytest.raises(ValueError, match="unique"):
        ranked.ranked_evidence(a, ["CASSF"], control=Control(["CASSF"]*2))
    h, c, g = ranked.ranked_evidence(a, [], k=1)
    assert h.height == c.height == g.height == 0


def test_predeclared_k_and_all_control_ties(monkeypatch):
    calls = []
    def statistic(k, n, nc, m):
        calls.append((k, n, nc, m))
        return {"p_global": 0.5}
    monkeypatch.setattr(ranked, "order_statistic", statistic)
    a = Annotator(VdjdbIndex.build(reference()))
    _, _, g = ranked.ranked_evidence(a, ["CASGF"], k=1, control=Control(["CASVF", "CASYF", "CASWF"]))
    assert calls == [(1, 2, 3, 3)]
    assert g["n_retained"][0] == 1
    assert g["threshold"][0] == 1


def test_paired_cardinality_and_missing_controls():
    a = Annotator(VdjdbIndex.build(reference()))
    with pytest.raises(ValueError, match="every active locus"):
        ranked.ranked_evidence(a, ["CASSF"], control={})
    q = pl.DataFrame({"cdr3a": ["CAVF"], "cdr3b": ["CASSF"]})
    with pytest.raises(ValueError, match="both TRA and TRB"):
        ranked.ranked_evidence(a, q, paired=True, control={"TRA": Control(["CAVF"])})


def test_forbidden_layout_is_not_a_finite_hit(monkeypatch):
    monkeypatch.setattr(ranked.gapblock, "topk_batch", lambda *a, **k: [[SimpleNamespace(ref_id=0, score=ranked.gapblock.UNREACHABLE)]])
    a = Annotator(VdjdbIndex.build(reference()))
    h, c, g = ranked.ranked_evidence(a, ["CASSF"], k=1)
    assert h.height == c.height == 0
    assert g["n_reference"][0] == 2
    assert g["threshold"][0] is None and g["p_global"][0] == 1


def test_paired_competitors_and_observation_hit_units():
    ref = pl.DataFrame({"complex_id": ["R1", "R1", "R2", "R2", "R3", "R3", "R4", "R4"],
        "gene": ["TRA", "TRB"]*4, "cdr3": ["CAVF", "CASSF", "CAVF", "CASSF", "CAYF", "CASAF", "CWWF", "CWWWF"],
        "epitope": ["E", "E", "E", "E", "F", "F", "G", "G"], "vdjdb_score": [1, 2, 1, 2, 3, 3, 0, 0]})
    q = pl.DataFrame({"cdr3a": ["CAVF"], "cdr3b": ["CASSF"]})
    _, candidates, _ = ranked.ranked_evidence(Annotator(VdjdbIndex.build(ref)), q, k=3, paired=True)
    rows = candidates.sort("rank").to_dicts()
    assert len(rows) == 3
    assert rows[0]["competing_score"] == rows[1]["ned_score"]
    assert rows[1]["competing_score"] == rows[0]["ned_score"]
    assert rows[2]["competing_score"] == rows[0]["ned_score"]
    assert rows[1]["score_margin"] < 0 and rows[2]["score_margin"] < 0
    assert [r["n_hits"] for r in rows] == [4, 4, 4]
    assert [r["n_clonotypes"] for r in rows] == [1, 1, 1]
    assert rows[0]["n_records"] == 2 and rows[0]["db_score"] == 2
    assert all(r["rule_of_three"] is False for r in rows)


def test_unsupported_locus_does_not_require_a_control():
    a = Annotator(VdjdbIndex.build(reference()))
    q = pl.DataFrame({"cdr3": ["CAVF"], "locus": ["TRA"]})
    h, c, g = ranked.ranked_evidence(a, q, control={})
    assert h.height == c.height == g.height == 0


def test_paired_metadata_ignores_unlinked_reference_rows():
    ref = pl.DataFrame({"complex_id": ["R1", "R1", "0", "0", None],
        "gene": ["TRA", "TRB", "TRA", "TRA", "TRB"],
        "cdr3": ["CAVF", "CASSF", "CAYF", "CWWF", "CASAF"],
        "epitope": ["E"] * 5})
    q = pl.DataFrame({"cdr3a": ["CAVF"], "cdr3b": ["CASSF"]})
    h, c, g = ranked.ranked_evidence(Annotator(VdjdbIndex.build(ref)), q, k=1, paired=True)
    assert h["complex_id"].to_list() == ["R1"]
    assert c["n_records"].to_list() == [1]
    assert g["n_reference"].to_list() == [1]


@pytest.mark.skipif(not hasattr(ranked.gapblock, "topk_batch"), reason="requires native gap-block top-K extension")
def test_actual_native_wrappers():
    from seqtree import Index
    a = Annotator(VdjdbIndex.build(reference()))
    h, c, g = ranked.ranked_evidence(a, ["CASSF"], k=2, control=Index.build(["CWWWF", "CYYYF"]))
    assert h.height == 3 and c.height == 2
    assert g["threshold"][0] == 1 and g["p_global"][0] > 0
    ref = pl.DataFrame({"complex_id": ["R1", "R1", "R2", "R2"], "gene": ["TRA", "TRB"]*2,
                        "cdr3": ["CAVF", "CASSF", "CAYF", "CASAF"], "epitope": ["E", "E", "F", "F"]})
    q = pl.DataFrame({"cdr3a": ["CAVF"], "cdr3b": ["CASSF"]})
    ctrl = {"TRA": Index.build(["CAVF", "CWWF"]), "TRB": Index.build(["CASSF", "CWWWF"])}
    h, c, g = ranked.ranked_evidence(Annotator(VdjdbIndex.build(ref)), q, k=1, paired=True, control=ctrl, exclude_exact=True)
    assert h["complex_id"].to_list() == ["R2"]
    assert g["n_reference"][0] == 2 and g["p_global"][0] > 0
