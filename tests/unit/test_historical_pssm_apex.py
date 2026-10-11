"""Tiny independent fixed-frame geometry checks; no benchmark predictions."""
import numpy as np
import polars as pl
from seqtree import Index,SubstitutionMatrix,gapblock,SearchParams
from vdjmatch.match.regions import significance_weights,significance_pssm


def test_single_start_and_short_clamp():
    prior=gapblock.positions_prior((6,))
    for shorter in (0,1,4,5,6,7,13,23):
        for d in (1,2,4):
            assert [i for i in range(shorter+1) if prior(i,d,shorter+d)==0]==[min(6,shorter)]
        assert all(prior(i,0,shorter)==0 for i in range(shorter+1))


def test_apex_native_symmetry_cdf_and_original_equal_length():
    B=SubstitutionMatrix.blosum62();q='CASSLGQAYEQYF';r=q[:6]+'A'+q[6:]
    short='CAGF';short_r=short+'A';old=q[:-2]+'AF'
    qs=[q,r,short,short_r];refs=[q,q,r,old,short,short_r]
    weights={n:[max(1,round(100*x)) for x in significance_weights(n)] for n in {len(x) for x in qs+refs}}
    opt=dict(matrix=B,gap_open=2800,gap_extend=1400,gap_prior=gapblock.positions_prior((6,)),position_weights_by_length=weights)
    a=np.asarray(gapblock.score_matrix(qs,refs,threads=1,**opt))
    assert a[0,2]==a[1,0]==2800
    assert a[2,5]==a[3,4]==2800
    np.testing.assert_array_equal(a,np.asarray(gapblock.score_matrix(refs,qs,threads=1,**opt)).T)
    np.testing.assert_array_equal(a,np.asarray(gapblock.score_matrix(qs,refs,threads=4,**opt)))
    cutoffs=[[0,2800,7000]]*len(qs)
    counts=gapblock.count_batch(qs,refs,cutoffs,exclude_exact=True,threads=1,**opt)
    assert counts==[[sum(qs[i]!=s and a[i,j]<=c for j,s in enumerate(refs)) for c in cutoffs[i]] for i in range(len(qs))]
    assert counts==gapblock.count_batch(qs,refs,cutoffs,exclude_exact=True,threads=4,**opt)
    oldopt={**opt,'gap_prior':gapblock.positions_prior((3,4,-4,-3))}
    oldmatrix=np.asarray(gapblock.score_matrix([q],[q,old],threads=1,**oldopt))
    newmatrix=np.asarray(gapblock.score_matrix([q],[q,old],threads=1,**opt))
    np.testing.assert_array_equal(oldmatrix,newmatrix)
    p=SearchParams(max_subs=5,max_ins=0,max_dels=0,max_total_edits=5,engine='seqtm');p.pos_matrix=significance_pssm(len(q))
    hits=Index.build([q,old],'aa').search_batch([q],p,threads=1)[0]
    assert sorted((h.ref_id,h.score) for h in hits)==[(j,int(newmatrix[0,j])) for j in range(2)]


def test_exact_and_original_edge_exclusion_retained():
    from vdjmatch.match.unified import historical_pssm_extension
    q='CASSLGQAYEQYF';old=q[:-2]+'AF'
    queries=pl.DataFrame({'query_id':['q','dup'],'cdr3':[q,q],'v':['TRBV19']*2})
    refs=pl.DataFrame({'cdr3':[q,q,old],'v':['TRBV19']*3})
    out=historical_pssm_extension(queries,refs,Index.build([q,q,old]),threads=1)
    assert out['query_id'].to_list()==['q','dup']
    assert out['gapped_density'].to_list()==[0.,0.]
    assert out['gapped_edges_same_v'].to_list()==[0,0]
