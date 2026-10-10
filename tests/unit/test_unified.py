"""Unified route geometry and evidence against an independent scalar oracle."""
import math
import polars as pl
from seqtree import Index, SubstitutionMatrix, gapblock
from vdjmatch.api import Annotator



def test_unified_gap_rejects_one_row_over_matrix_budget_before_native_allocation(monkeypatch):
    from types import SimpleNamespace
    import pytest
    from vdjmatch.match.unified import unified_evidence
    class OversizedSequences:
        def __len__(self): return 64*1024**2//4+1
    class VirtualReference:
        def select(self,*args): return self
        def unique(self): return self
        def sort(self,*args): return self
        def with_row_index(self,*args): return self
        def __getitem__(self,key): return self
        def to_list(self): return OversizedSequences()
    ann=SimpleNamespace(loci=['TRB'],_index=SimpleNamespace(records_for=lambda locus:VirtualReference()))
    def forbidden(*args,**kwargs):
        raise AssertionError('oversized one-query matrix must be rejected before native allocation')
    monkeypatch.setattr(gapblock,'score_matrix',forbidden)
    monkeypatch.setattr(gapblock,'count_batch',forbidden)
    q=pl.DataFrame({'query_id':['q'],'cdr3':['CASSLGQAYEQYF'],'locus':['TRB']})
    with pytest.raises(ValueError,match='64MiB matrix budget for one query'):
        unified_evidence(ann,q,control=Index.build(['CASSLGQAYEQYF']),distance='gapblock')


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


def test_gapped_extension_categorical_v_prior_and_flat_diagnostics(monkeypatch):
    from vdjmatch.match.unified import gapped_extension
    from vdjmatch.match import tcrdist
    q='CASSLGQAYEQYF'
    same,cross='TRBV24-1*01','TRBV24/OR9-2*01'
    model=tcrdist.load_v_loops()
    assert same.split('*')[0]!=cross.split('*')[0]
    assert all(model[same][k]==model[cross][k] for k in ('cdr1','cdr2','cdr25'))
    refs=pl.DataFrame({'cdr3':[q+'F',q+'W'],'v':[same,cross]})
    queries=pl.DataFrame({'query_id':['q','unknown'],'cdr3':[q,q],'v':[same,'unknown']})
    total,cdr3=tcrdist.distance_matrix([q],refs['cdr3'].to_list(),[same],[same,cross],return_cdr3=True)
    assert total[0,0]==total[0,1] and (total==cdr3).all()
    calls=[];native_distance=tcrdist.distance_matrix;native_count=gapblock.count_batch
    def distance(*args,**kwargs):
        calls.append('distance');return native_distance(*args,**kwargs)
    def count(*args,**kwargs):
        calls.append('count');return native_count(*args,**kwargs)
    monkeypatch.setattr(tcrdist,'distance_matrix',distance)
    monkeypatch.setattr(gapblock,'count_batch',count)
    control=Index.build([q]);tau=15
    out=gapped_extension(queries,refs,control,temperature=tau)
    assert calls==['distance','count']
    row=out.row(0,named=True);base=math.exp(-int(total[0,0])/tau)/.01
    assert math.isclose(row['gapped_density_same_v_unweighted'],base)
    assert math.isclose(row['gapped_density_cross_v_unweighted'],base)
    assert math.isclose(row['gapped_density'],1.25*base)
    assert math.isclose(row['gapped_density'],row['gapped_density_same_v_unweighted']+.25*row['gapped_density_cross_v_unweighted'])
    assert row['gapped_edges_same_v']==row['gapped_edges_cross_v']==1
    assert row['gapped_floor_density']==row['gapped_density']
    assert row['gapped_best_total_distance']==int(total[0,0])
    assert row['gapped_best_cdr3_distance']==int(cdr3[0,0])
    assert row['gapped_best_vloop_distance']==0
    empty=gapped_extension(queries.head(1),refs.head(0),control)
    no_edges=gapped_extension(queries.head(1),pl.DataFrame({'cdr3':[q],'v':[same]}),control)
    assert out.schema==empty.schema==no_edges.schema
    assert empty['gapped_status'].to_list()==['no_usable_reference']
    for frame in (out.tail(1),empty,no_edges):
        for name in out.columns:
            if name.startswith('gapped_best_'):
                assert frame[name].to_list()==[None]
            elif name.startswith(('gapped_density','gapped_edges_','gapped_floor_')):
                assert frame[name].to_list()==[0]
    assert calls==['distance','count','distance']
