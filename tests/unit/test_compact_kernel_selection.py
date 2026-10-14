"""Compact maximum reductions against the unchanged public table oracle."""
import hashlib
from itertools import product

import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal

import vdjmatch.evalue.kernel as kernel


def layout(names):
    order=np.asarray(sorted(range(len(names)),key=lambda i:(names[i] is None,names[i] or '')))
    ordered=[names[i] for i in order]
    starts=np.asarray([i for i,e in enumerate(ordered) if i==0 or e!=ordered[i-1]])
    return order,starts,[ordered[i] for i in starts]


def table_oracle(scores,exact,order,starts,peptides):
    C,T=scores.shape
    names=np.empty(T,dtype=object)
    for i,start in enumerate(starts):
        names[order[start:starts[i+1] if i+1<len(starts) else T]]=peptides[i]
    frame=pl.DataFrame({'query_id':np.repeat(np.arange(C),T),
        'epitope':pl.Series(np.tile(names,C),dtype=pl.String),
        'kernel_score':scores.reshape(-1),'exact_reference_count':exact.reshape(-1)})
    out=kernel.select_kernel(pl.DataFrame({'query_id':np.arange(C)}),frame)
    return {c:out[c].to_numpy() for c in ('selected_epitope','kernel_score',
        'kernel_competing_score','tied_peptides','exact_peptides')}


@pytest.mark.parametrize('names',[['B','A','B',None,'A'],['only'],[None],['z',None,'a']])
def test_maxima_ties_missing_lexical_and_exact_counts(names):
    rng=np.random.default_rng(19)
    scores=rng.choice([0.,-1.,1.,3.,np.nan,np.inf,-np.inf],size=(83,len(names)))
    scores[0]=0;scores[1]=np.nan;scores[2]=3;scores[3]=-1
    exact=rng.integers(0,4,size=scores.shape,dtype=np.uint32)
    args=layout(names)
    actual=kernel._compact_kernel_selection(scores,exact,*args)
    expected=table_oracle(scores,exact,*args)
    for c in actual:
        if c=='selected_epitope':assert actual[c].tolist()==expected[c].tolist()
        else:np.testing.assert_array_equal(actual[c],expected[c])
    assert actual['tied_peptides'].dtype==np.uint32
    assert actual['exact_peptides'].dtype==np.uint32


@pytest.fixture(scope='module')
def native_inputs():
    seqs=['CASSL'+a+b+'AYEQYF' for a,b in list(product('ACDEFGHIKLMNPQRSTVWY',repeat=2))[:71]]
    seqs+=['CWWWWWWWWWWWF',seqs[0],seqs[0]]
    queries=pl.DataFrame({'query_id':['q'+str(i) for i in range(len(seqs))],
        'cdr3':seqs,'_allele':['TRBV19*01' if i%2 else 'TRBV7-9*01' for i in range(len(seqs))],
        'j':['TRBJ1-1','TRBJ2-1']*(len(seqs)//2)})
    refs=pl.DataFrame({'cdr3':[seqs[0],seqs[0],seqs[2],seqs[3],seqs[0]],
        '_allele':['TRBV19*01','TRBV19*01','TRBV7-9*01','TRBV19*01','TRBV7-9*01'],
        'epitope':['B','A','B','C','A'],'restriction':['r1','r1','r2','r1','r2']})
    groups=[({'epitope':e},queries,part) for (e,_),part in refs.group_by('epitope','restriction',maintain_order=True)]
    from vdjmatch.match.geometry import precompute_v_loop_distances
    return queries,refs,groups,precompute_v_loop_distances(matrix='blosum62')


@pytest.mark.parametrize('scorer',kernel.KERNEL_SCORERS)
@pytest.mark.parametrize('exclude',[False,True])
def test_native_batch_old_new_exact_bytes_and_threads(native_inputs,monkeypatch,scorer,exclude):
    queries,refs,groups,context=native_inputs
    reducer=kernel._compact_kernel_selection
    outputs=[];hit_sets=[]
    for implementation,threads in [(table_oracle,1),(reducer,1),(reducer,4)]:
        monkeypatch.setattr(kernel,'_compact_kernel_selection',implementation)
        state,hits=kernel.grouped_kernel_selection(queries,refs,groups,
            query_keys=['cdr3','_allele'],search_keys=['cdr3','_allele'],paired=False,
            cutoffs=[0,4,8,13],radius=4,kernel_scale=22.5,scorer=scorer,
            threads=threads,matrix='blosum62',ctrim=3,stem_depth=4,
            exclude_exact=exclude,v_loop_distances=context)
        result=kernel.finish_kernel_selection(queries.select('query_id'),state)
        outputs.append(result);hit_sets.append(hits)
    assert_frame_equal(outputs[0],outputs[1],check_exact=True)
    assert_frame_equal(outputs[0],outputs[2],check_exact=True)
    digests=[hashlib.sha256(out.write_csv(separator='\t').encode()).hexdigest() for out in outputs]
    assert len(set(digests))==1
    assert hit_sets[0]==hit_sets[1]==hit_sets[2]
    assert outputs[1]['query_id'].dtype==queries['query_id'].dtype
    assert outputs[1].height==queries.height


@pytest.mark.parametrize('dtype',[pl.UInt32,pl.UInt64])
def test_numeric_source_id_dtype_retained(native_inputs,dtype):
    queries,refs,groups,context=native_inputs
    queries=queries.with_columns(pl.Series('query_id',range(queries.height),dtype=dtype))
    groups=[(task,queries,part) for task,_,part in groups]
    state,hits=kernel.grouped_kernel_selection(queries,refs,groups,
        query_keys=['cdr3','_allele'],search_keys=['cdr3','_allele'],paired=False,
        cutoffs=[0,4],radius=4,kernel_scale=22.5,scorer='v-bounded',
        matrix='blosum62',ctrim=3,v_loop_distances=context)
    assert state['query_id'].dtype==dtype
    assert set(state['query_id'])==set(range(queries.height))
    assert hits<=set(range(queries.height))


@pytest.mark.parametrize('assayed',[False,True])
def test_native_chunk_budget_preserves_scores_duplicates_and_source_order(native_inputs,monkeypatch,assayed):
    from vdjmatch.match import geometry
    queries,refs,groups,context=native_inputs
    queries=queries.reverse()
    groups=[(task,queries,part) for task,_,part in groups]
    native=geometry.grouped_total_distance_count_batch
    outputs=[];hits=[];budgets=[]
    for small in (True,False):
        def counted(*args,**kwargs):
            budgets.append(kwargs['chunk_size'])
            if small:kwargs['chunk_size']=64
            return native(*args,**kwargs)
        monkeypatch.setattr(geometry,'grouped_total_distance_count_batch',counted)
        state,supported=kernel.grouped_kernel_selection(queries,refs,groups,
            query_keys=['cdr3','_allele'],search_keys=['cdr3','_allele'],paired=False,
            cutoffs=[0,4,8,13],radius=4,kernel_scale=22.5,scorer='v-bounded',assayed=assayed,
            threads=4,matrix='blosum62',ctrim=3,stem_depth=4,exclude_exact=True,
            v_loop_distances=context)
        outputs.append(kernel.finish_kernel_selection(queries.select('query_id'),state));hits.append(supported)
    assert budgets==[65536,65536]
    assert_frame_equal(outputs[0],outputs[1],check_exact=True)
    assert outputs[1]['query_id'].to_list()==queries['query_id'].to_list()
    assert hits[0]==hits[1]
