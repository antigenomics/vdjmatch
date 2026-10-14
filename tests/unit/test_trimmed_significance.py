"""Recorded restricted significance profile and CDR3 flooring, synthetic only."""
import numpy as np
import polars as pl
import pytest
from seqtree import SubstitutionMatrix, gapblock
from vdjmatch.match import tcrdist as td
from vdjmatch.match.regions import significance_weights


def recorded_distance(q,r,qv,rv,matrix,context):
    costs=td._substitution_costs(matrix=matrix)
    weights={}
    for length in {len(s) for s in q+r}:
        profile=significance_weights(length)[3:-3]
        mean=sum(profile)/len(profile)
        weights[length-6]=[max(1,round(100*x/mean)) for x in profile]
    raw=np.asarray(gapblock.score_matrix([s[3:-3] for s in q],[s[3:-3] for s in r],
        matrix=SubstitutionMatrix.from_similarity((-costs).tolist()),
        gap_open=200*int(costs.max()),gap_extend=200*int(costs.max()),
        position_weights_by_length=weights,gap_prior=td._gap_prior_trim3,threads=1),dtype=np.int64)
    genes=td.resolve_v_alleles(qv+rv,context['model']).to_list()
    loops,ids=td._v_loop_groups(genes,context['model'],costs,context)
    return 3*(raw//200)+loops[ids[:len(q),None],ids[None,len(q):]]+td.stem_penalty(q,r,depth=4,costs=costs)


def frames(paired):
    q=pl.DataFrame({'seq0':['CASSLGQAYEQYF','CASSLGRAYEQYF'],'v0':['TRAV12-1','TRAV12-2'] if paired else ['TRBV19','TRBV7-9']})
    r=pl.DataFrame({'seq0':['CASSLGQAYEQYF','WASSLGQAYEQYW','CASSLWGQAAYEQYF'],
        'v0':['TRAV12-2','TRAV12-1','TRAV12-1'] if paired else ['TRBV7-9','TRBV19','TRBV19']})
    if paired:
        q=q.with_columns(pl.Series('seq1',['CASSIVGGNEQFF','CASSMRSTGELFF']),pl.Series('v1',['TRBV19','TRBV7-9']))
        r=r.with_columns(pl.Series('seq1',['CASSIVGGNEQFF','WASSIVGGNEQFW','CASSIVGAGNEQFF']),
            pl.Series('v1',['TRBV7-9','TRBV19','TRBV19']))
    return q,r


@pytest.mark.parametrize('matrix',['tcrdist','blosum62'])
@pytest.mark.parametrize('paired',[False,True])
def test_trimmed_profile_exact_native_boundaries_context_and_threads(matrix,paired):
    q,r=frames(paired);context=td.precompute_v_loop_distances(matrix=matrix,v_regions='cdr12')
    distance=np.zeros((q.height,r.height),dtype=np.int64)
    for i in range(2 if paired else 1):
        distance+=recorded_distance(q['seq'+str(i)].to_list(),r['seq'+str(i)].to_list(),
            q['v'+str(i)].to_list(),r['v'+str(i)].to_list(),matrix,context)
    cuts=sorted({0,*[200*int(x)+delta for x in distance.flat for delta in (-1,0,1)]})
    for threads in (1,4):
        for exclude in (False,True):
            iterator=td.grouped_total_distance_count_batch(q,r,cuts,paired=paired,matrix=matrix,
                ctrim=3,position_weighting='significance-trimmed',stem_depth=4,
                v_regions='cdr12',v_loop_distances=context,threads=threads,exclude_exact=exclude,
                reference_target_offsets=[0,1,2,3],reference_target_ids=[0,0,0],
                target_peptide_ids=[0],reference_exclusive_peptide_ids=[0,0,0])
            got=np.concatenate([np.asarray(chunk.counts)[:,0,:] for chunk in iterator])
            eligible=np.ones(distance.shape,dtype=bool)
            if exclude:
                exact=np.ones(distance.shape,dtype=bool)
                for i in range(2 if paired else 1):
                    exact &= np.asarray(q['seq'+str(i)])[:,None]==np.asarray(r['seq'+str(i)])[None,:]
                eligible=~exact
            expected=((200*distance[:,:,None]<=np.asarray(cuts))&eligible[:,:,None]).sum(axis=1)
            np.testing.assert_array_equal(got,expected)


def test_rounding_before_v_and_stem_offsets_on_a_genuine_residual():
    costs=td._substitution_costs();symbols=td.amino_acids();profile=significance_weights(13)[3:-3]
    mean=sum(profile)/len(profile);weights=[max(1,round(100*x/mean)) for x in profile]
    site,aa,raw=next((i,a,2*int(costs[symbols.index('A'),symbols.index(a)])*w)
        for i,w in enumerate(weights) for a in 'RNDCEQGHILKMFPSTWYV'
        if (2*int(costs[symbols.index('A'),symbols.index(a)])*w)%200)
    q='CAA'+'A'*7+'AAF';r=q[:site+3]+aa+q[site+4:]
    context=td.precompute_v_loop_distances()
    oracle=recorded_distance([q],[r],['TRBV19'],['TRBV7-9'],'tcrdist',context)[0,0]
    assert raw%200 and 3*(raw//200)<3*raw/200
    cuts=[int(200*oracle)-1,int(200*oracle)]
    iterator=td.grouped_total_distance_count_batch(pl.DataFrame({'seq0':[q],'v0':['TRBV19']}),
        pl.DataFrame({'seq0':[r],'v0':['TRBV7-9']}),cuts,ctrim=3,
        position_weighting='significance-trimmed',stem_depth=4,
        reference_target_offsets=[0,1],reference_target_ids=[0],target_peptide_ids=[0],
        reference_exclusive_peptide_ids=[0])
    assert np.asarray(list(iterator)[0].counts)[0,0,:].tolist()==[0,1]


def test_uniform_default_bitwise_and_other_count_routes_not_extended():
    q,r=frames(False);kw=dict(ctrim=3,stem_depth=4,reference_target_offsets=[0,1,2,3],
        reference_target_ids=[0,0,0],target_peptide_ids=[0],reference_exclusive_peptide_ids=[0,0,0])
    default=list(td.grouped_total_distance_count_batch(q,r,[0,5000,10000],**kw))
    explicit=list(td.grouped_total_distance_count_batch(q,r,[0,5000,10000],position_weighting='uniform',**kw))
    np.testing.assert_array_equal(default[0].counts,explicit[0].counts)
    with pytest.raises(ValueError,match='position_weighting'):
        td.total_distance_count_batch([],[],[],[],[],position_weighting='significance-trimmed')
    with pytest.raises(ValueError,match='trim3'):
        td.grouped_total_distance_count_batch(q,r,[0],position_weighting='significance-trimmed')


def test_old_native_capability_fails_explicitly(monkeypatch):
    monkeypatch.setattr(td.gapblock,'grouped_count_batch',lambda q,r,cuts,**kw:None)
    q,r=frames(False)
    with pytest.raises(RuntimeError,match='quant'):
        td.grouped_total_distance_count_batch(q,r,[0],ctrim=3,position_weighting='significance-trimmed')


def test_kernel_selection_forwards_profile_into_actual_native_batch(monkeypatch):
    from vdjmatch.evalue.kernel import grouped_kernel_selection
    q,r=frames(False)
    q=q.rename({'seq0':'cdr3','v0':'_allele'}).with_row_index('query_id')
    r=r.rename({'seq0':'cdr3','v0':'_allele'}).with_columns(pl.Series('epitope',['A','A','B']))
    calls=[];native=td.grouped_total_distance_count_batch
    def capture(*args,**kwargs):
        calls.append(kwargs['position_weighting'])
        return native(*args,**kwargs)
    monkeypatch.setattr(td,'grouped_total_distance_count_batch',capture)
    groups=[({'epitope':epitope},q,r.filter(pl.col('epitope')==epitope)) for epitope in ('A','B')]
    state,hit_ids=grouped_kernel_selection(q,r,groups,query_keys=['cdr3','_allele'],
        search_keys=['cdr3','_allele'],paired=False,cutoffs=[0,10,30,100,300],
        radius=100,kernel_scale=90,scorer='v-bounded',ctrim=3,stem_depth=4,
        position_weighting='significance-trimmed')
    assert calls==['significance-trimmed']
    assert sorted(state['query_id'].to_list())==[0,1]
    assert state['kernel_score'].is_finite().all()
    assert hit_ids=={0,1}
