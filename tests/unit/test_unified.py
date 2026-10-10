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


def test_gapped_extension_preserves_original_edges_and_skips_control_search(monkeypatch):
    from vdjmatch.match.unified import gapped_extension
    q='CASSLGQAYEQYF'
    queries=pl.DataFrame({'query_id':['a'],'cdr3':[q],'v':['TRBV19']})
    refs=pl.DataFrame({'cdr3':[q,q.replace('Q','R'),q[:3]+'TTTTT'+q[8:]],'v':['TRBV19']*3})
    def forbidden(*args,**kwargs):
        raise AssertionError('no new edges must not trigger a control CDF search')
    monkeypatch.setattr(gapblock,'count_batch',forbidden)
    out=gapped_extension(queries,refs,Index.build([q]))
    assert out['gapped_density'].to_list()==[0.0]
    assert out['availability'].to_list()==[True]


def test_gapped_extension_equal_length_six_substitutions_and_zero_cdf_floor():
    from vdjmatch.match.unified import gapped_extension
    from vdjmatch.match.tcrdist import distance_matrix
    q='CASSLGQAYEQYF';r=q[:3]+'TTTTTT'+q[9:]
    assert len(q)==len(r) and sum(a!=b for a,b in zip(q,r))==6
    queries=pl.DataFrame({'query_id':['a'],'cdr3':[q],'v':['TRBV19']})
    refs=pl.DataFrame({'cdr3':[r],'v':['TRBV19']})
    total=distance_matrix([q],[r],['TRBV19'],['TRBV19'])[0,0]
    assert total<=90
    out=gapped_extension(queries,refs,Index.build([q]))
    assert math.isclose(out['gapped_density'][0],100*math.exp(-int(total)/12))
    missing=gapped_extension(queries.with_columns(pl.lit(None).alias('v')),refs,Index.build([q]))
    assert missing['gapped_density'].to_list()==[0.0] and missing['availability'].to_list()==[False]


def test_gapped_extension_full_control_puncture_and_serial_equality(monkeypatch):
    from vdjmatch.match.unified import gapped_extension
    from vdjmatch.match.tcrdist import distance_matrix
    q='CASSLGQAYEQYF';r=q+'F';same_trim='CAVSLGQAYEQWF'
    assert q[3:-2]==same_trim[3:-2] and q!=same_trim
    queries=pl.DataFrame({'query_id':['a','repeat','unknown'],'cdr3':[q]*3,
                          'v':['TRBV19','TRBV19','unknown']})
    refs=pl.DataFrame({'cdr3':[r,r],'v':['TRBV19','TRBV7-9']})
    control=Index.build([q,q,same_trim,'CAAF'])
    calls=[];native=gapblock.count_batch
    def observed(seqs,refs,thresholds,**kwargs):
        calls.append((seqs,refs,thresholds,kwargs['exclude_exact']))
        return native(seqs,refs,thresholds,**kwargs)
    monkeypatch.setattr(gapblock,'count_batch',observed)
    one=gapped_extension(queries,refs,control,threads=1)
    many=gapped_extension(queries,refs,control,threads=2)
    assert one.equals(many)
    assert one['M_gap'].to_list()==[2]*3
    assert one['availability'].to_list()==[True,True,False]
    assert one['n_reference_unavailable'].to_list()==[0]*3
    total,cdr3=distance_matrix([q],[r],['TRBV19'],['TRBV19'],return_cdr3=True)
    expected=math.exp(-int(total[0,0])/12)/.5
    assert all(math.isclose(s,expected) for s in one['gapped_density'][:2])
    assert one['gapped_density'][2]==0
    assert calls==[([q[3:-2]],[q[3:-2],q[3:-2]],[[2*int(cdr3[0,0])//3]],False)]*2
    other=gapped_extension(queries.head(1),refs,Index.build([same_trim]))
    assert math.isclose(other['gapped_density'][0],math.exp(-int(total[0,0])/12))


def test_gapped_extension_uses_v_loops_but_cdr3_only_background_and_first_v():
    from vdjmatch.match.unified import gapped_extension
    from vdjmatch.match.tcrdist import distance_matrix
    q='CASSLGQAYEQYF';r=q+'F'
    queries=pl.DataFrame({'query_id':['a'],'cdr3':[q],'v':['TRBV19']})
    refs=pl.DataFrame({'cdr3':[r],'v':['TRBV19*02']})
    control=Index.build([q,'CAVSLGQAYEQWF'])
    total,cdr3=distance_matrix([q],[r],['TRBV19'],['TRBV19*02'],return_cdr3=True)
    assert total[0,0]>cdr3[0,0]
    out=gapped_extension(queries,refs,control)
    assert math.isclose(out['gapped_density'][0],2*math.exp(-int(total[0,0])/12))
    unknown=pl.DataFrame({'cdr3':[r,r],'v':['unknown','TRBV19']})
    unavailable=gapped_extension(queries,unknown,control)
    assert unavailable['gapped_density'].to_list()==[0.0]
    assert unavailable['availability'].to_list()==[False]
    assert unavailable['n_reference_unavailable'].to_list()==[1]
