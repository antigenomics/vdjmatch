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
    # One target hit batch and one bounded control histogram batch.
    original = Index.search_batch
    calls = []
    def batch(self, queries, params, threads=0):
        calls.append((len(queries), params.max_total_edits))
        return original(self, queries, params, threads)
    monkeypatch.setattr(Index, "search_batch", batch)
    got = ann.graded_candidates(["CAAAAAAF"], control=control)
    assert calls == [(1, 5)]
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


@pytest.mark.parametrize("exclude_exact", [False, True])
def test_graded_histogram_preserves_all_fixed_control_fields(exclude_exact):
    from polars.testing import assert_frame_equal
    from vdjmatch.match.scope import search_params
    ref = pl.DataFrame({'gene': ['TRB']*4, 'cdr3': ['CASSF','CASRF','CASTF','CARRF'],
                        'epitope': ['A','A','B','B']})
    ann = Annotator.from_frame(ref)
    # Repeated model memberships and repeated query junctions retain their roles.
    control = Index.build(['CASSF','CASSF','CASRF','CYYYYF'])
    query = pl.DataFrame({'query_id':['a','b','c'], 'cdr3':['CASSF','CASSF','CARRF']})
    graded = ann.graded_candidates(query, control=control, exclude_exact=exclude_exact)
    for radius in range(1,6):
        fixed = ann.candidates(query, scope=search_params(f'{radius},{min(radius,2)},{min(radius,2)},{radius}', gap_open=1, gap_extend=1),
            control=control, exclude_exact=exclude_exact, score_scale=1.0)
        got = graded.filter(pl.col('radius')==radius).drop('radius')
        assert_frame_equal(got.sort('query_id','epitope'),fixed.sort('query_id','epitope'))


def test_bounded_observation_expansion_keeps_one_native_batch_and_all_candidates(monkeypatch):
    from polars.testing import assert_frame_equal
    from vdjmatch.api import _prepare
    from vdjmatch.match.scope import search_params
    amino='ACDEFGHIKLMNPQRSTVWY'
    sequences=['C'+amino[i%20]+amino[i//20]+'SSF' for i in range(72)]
    ref=pl.DataFrame({'gene':['TRB']*72,'cdr3':sequences,'epitope':['E'+str(i%3) for i in range(72)],'reference_id':['R'+str(i%5) for i in range(72)]})
    ann=Annotator.from_frame(ref)
    _,q=_prepare(pl.DataFrame({'query_id':['Q'+str(i) for i in range(73)],'cdr3':sequences+['CYYYYYYYYF']}))
    params=search_params('5,2,2,5',gap_open=1,gap_extend=1)
    kwargs=dict(threads=1,radii=range(1,6),exclude_exact=True,score_scale=1.0)
    _,expected=ann._evidence(q,params,**kwargs)
    original=Index.search_batch;calls=[]
    def batch(self,queries,params,threads=0):
        calls.append(len(queries));return original(self,queries,params,threads)
    monkeypatch.setattr(Index,'search_batch',batch)
    hits,got=ann._evidence(q,params,retain_hits=False,**kwargs)
    assert hits.is_empty() and calls==[73]
    assert_frame_equal(got,expected,check_exact=True)


def test_retained_hits_keep_the_same_epitope_exclusion():
    from vdjmatch.api import _prepare
    from vdjmatch.match.scope import search_params
    ref=pl.DataFrame({'gene':['TRB','TRB'],'cdr3':['CASSF','CASRF'],'epitope':['A','B']})
    ann=Annotator.from_frame(ref)
    _,q=_prepare(['CASSF'])
    excluded=pl.DataFrame({'query_id':q['query_id'],'epitope':['A']})
    hits,candidates=ann._evidence(q,search_params('1',gap_open=1,gap_extend=1),exclude_epitopes=excluded)
    assert hits['epitope'].to_list()==['B']
    assert candidates['epitope'].to_list()==['B']
