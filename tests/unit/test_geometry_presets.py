"""Fixed historical geometry recipes against independent synthetic alignments."""
import numpy as np
import polars as pl
import pytest
from types import SimpleNamespace
from seqtree import SubstitutionMatrix,gapblock,pairwise
from vdjmatch.match import tcrdist as td


def inputs(paired):
    q=pl.DataFrame({'seq0':['CASSLGQAYEQYF','CASSGRYEQYF'],
        'v0':['TRAV12-1','TRAV12-2'] if paired else ['TRBV19','TRBV7-9'],
        'j0':['J1*01','J2*01']})
    r=pl.DataFrame({'seq0':['WASSLGQAYEQYW','CASSLWGQAAYEQYF','CASSLGQAYEQYF'],
        'v0':['TRAV12-2','TRAV12-1','TRAV12-1'] if paired else ['TRBV7-9','TRBV19','TRBV19'],
        'j0':['J2*02','J1*01','J2*01']})
    if paired:
        q=q.with_columns(pl.Series('seq1',['CASSIVGGNEQFF','CASSMRSTGELFF']),
            pl.Series('v1',['TRBV19','TRBV7-9']),pl.Series('j1',['J1','J2']))
        r=r.with_columns(pl.Series('seq1',['CASSIVGGNEQFF','CASSIVGAGNEQFF','CASSMRSTGELFF']),
            pl.Series('v1',['TRBV7-9','TRBV19','TRBV19']),pl.Series('j1',['J2','J1','J2']))
    return q,r


def test_recorded_frame_argument_order_and_native_units():
    q,r=inputs(False);context=td.precompute_v_loop_distances()
    qv=td.resolve_v_alleles(q['v0'],context['model']).to_list()
    rv=td.resolve_v_alleles(r['v0'],context['model']).to_list()
    options={}
    td._grouped_preset_options('framev',q,r,[(q['seq0'].to_list(),r['seq0'].to_list(),qv,rv)],options,context,False)
    # frame_prior(lam=1,c=Pmax=4), evaluated at trimmed column5.
    # Unweighted prototype raw is lifted by100*3, then floored to600.
    assert options['gap_prior'](8,1,14)==300
    assert options['score_quantum']==600
    assert options['gap_prior'](0,1,14)==gapblock.UNREACHABLE


def oracle(q,r,preset,matrix,context,paired):
    P=td._substitution_costs(matrix=matrix);pmax=int(P.max());aa=td.amino_acids()
    out=np.zeros((q.height,r.height),dtype=np.int64)
    for lane in range(2 if paired else 1):
        qs,rs=q['seq'+str(lane)].to_list(),r['seq'+str(lane)].to_list()
        qv,rv=q['v'+str(lane)].to_list(),r['v'+str(lane)].to_list()
        qt,rt=[s[3:-3] for s in qs],[s[3:-3] for s in rs]
        if preset in ('sw','nw'):
            similarity=SubstitutionMatrix.blosum62()
            sim=np.asarray(pairwise.score_matrix(qt,rt,similarity,
                mode='local' if preset=='sw' else 'global',gap_open=max(1,round(pmax/2)),
                gap_extend=max(1,round(pmax/2))))
            diag=np.asarray([similarity.similarity(a,a) for a in aa])
            selfq=np.asarray([sum(diag[aa.index(a)] for a in s) for s in qt])
            selfr=np.asarray([sum(diag[aa.index(a)] for a in s) for s in rt])
            d=3*np.maximum(0,selfq[:,None]+selfr[None,:]-2*sim)
        else:
            prior=td._gap_prior_trim3
            if preset=='freegap':prior=None
            if preset in ('framev','framemid','framej'):
                locus='TRA' if paired and lane==0 else 'TRB'
                column={'TRA':(0,1,3),'TRB':(1,3,6)}[locus][('framev','framemid','framej').index(preset)]
                prior=gapblock.frame_prior(column,pmax)
            raw=np.asarray(gapblock.score_matrix(qt,rt,matrix=SubstitutionMatrix.from_similarity((-P).tolist()),
                gap_open=2*pmax,gap_extend=2*pmax,gap_prior=prior,threads=1))
            d=3*(raw//2)
        genes=td.resolve_v_alleles(qv+rv,context['model']).to_list()
        loops,ids=td._v_loop_groups(genes,context['model'],P,context)
        d+=loops[ids[:len(qs),None],ids[None,len(qs):]]
        if preset=='flank':
            for pos in (0,1,2,-3,-2,-1):
                qi=np.asarray([aa.index(s[pos]) for s in qs]);ri=np.asarray([aa.index(s[pos]) for s in rs])
                d+=3*P[qi[:,None],ri[None,:]]
        if preset=='stem4j':d+=td.stem_penalty(qs,rs,depth=4,costs=P)
        if preset in td.J_GEOMETRY_PRESETS:
            qj=np.asarray([s.split('*')[0] for s in q['j'+str(lane)]]);rj=np.asarray([s.split('*')[0] for s in r['j'+str(lane)]])
            d+=pmax*(qj[:,None]!=rj[None,:])
        if preset in ('jlen','vjlen'):
            d+=pmax*np.abs(np.asarray([len(s) for s in qs])[:,None]-np.asarray([len(s) for s in rs])[None,:])
        if preset=='vjlen':
            d+=pmax*(np.asarray(qv)[:,None]!=np.asarray(rv)[None,:])
        out+=d
    return out


@pytest.mark.parametrize('preset',td.GEOMETRY_PRESETS)
@pytest.mark.parametrize('paired',[False,True])
def test_fixed_presets_exact_boundaries_full_identity_and_threads(preset,paired):
    matrix='tcrdist' if preset in ('freegap','framev','framemid','framej') else 'blosum62'
    context=td.precompute_v_loop_distances(matrix=matrix)
    q,r=inputs(paired);distance=oracle(q,r,preset,matrix,context,paired)
    cuts=sorted({0,*[int(200*x)+delta for x in distance.flat for delta in (-1,0,1)]})
    for threads in (1,4):
        for exclude in (False,True):
            chunks=td.grouped_total_distance_count_batch(q,r,cuts,paired=paired,matrix=matrix,
                ctrim=3,geometry_preset=preset,stem_depth=4 if preset=='stem4j' else 0,
                v_loop_distances=context,threads=threads,exclude_exact=exclude,
                reference_target_offsets=[0,1,2,3],reference_target_ids=[0,0,0],
                target_peptide_ids=[0],reference_exclusive_peptide_ids=[0,0,0])
            got=np.concatenate([np.asarray(chunk.counts)[:,0,:] for chunk in chunks])
            eligible=np.ones(distance.shape,dtype=bool)
            if exclude:
                exact=np.ones(distance.shape,dtype=bool)
                for lane in range(2 if paired else 1):
                    exact&=np.asarray(q['seq'+str(lane)])[:,None]==np.asarray(r['seq'+str(lane)])[None,:]
                eligible=~exact
            expected=((200*distance[:,:,None]<=np.asarray(cuts))&eligible[:,:,None]).sum(axis=1)
            np.testing.assert_array_equal(got,expected)


@pytest.mark.parametrize('empty',['neither','reference','control'])
@pytest.mark.parametrize('preset',['norm','jnorm'])
@pytest.mark.parametrize('paired',[False,True])
def test_normalization_pre_puncture_same_v_mass_and_empty_semantics(monkeypatch,preset,empty,paired):
    from vdjmatch.evalue.kernel import grouped_kernel_selection
    q=pl.DataFrame({'query_id':['stable'],'cdr3':['CASSLGQAYEQYF'],'_allele':['TRBV19*01'],'j':['J1']})
    r=pl.DataFrame({'cdr3':['CASSLGQAYEQYF','CASSWGQAYEQYF','CASSIGQAYEQYF'],
        '_allele':['TRBV19*01']*3,'j':['J1']*3,'epitope':['A','A','B']})
    if paired:
        q=q.rename({'cdr3':'cdr3a','_allele':'_allelea','j':'ja'}).with_columns(
            pl.lit('CASSIVGGNEQFF').alias('cdr3b'),pl.lit('TRBV7-9*01').alias('_alleleb'),pl.lit('J2').alias('jb'))
        r=r.rename({'cdr3':'alpha','_allele':'_allele_alpha','j':'j_alpha'}).with_columns(
            pl.lit('CASSIVGGNEQFF').alias('beta'),pl.lit('TRBV7-9*01').alias('_allele_beta'),pl.lit('J2').alias('j_beta'))
    target=0 if empty=='reference' else 1
    control=0 if empty=='control' else 1
    curves=np.asarray([[[target+control]*2,[target]*2,[target]*2,[target]*2]],dtype=np.uint32)
    sizes=np.asarray([[target+control,target,target,target]],dtype=np.uint32)
    def counts(*args,**kwargs):
        assert kwargs['geometry_preset']==preset
        if preset=='jnorm':
            assert 'j0' in args[0].columns
            if paired:assert 'j1' in args[0].columns
        return [SimpleNamespace(start=0,stop=1,counts=curves,sizes=sizes,exact_counts=np.zeros((1,1),dtype=np.uint32))]
    monkeypatch.setattr(td,'grouped_total_distance_count_batch',counts)
    keys=(['cdr3a','cdr3b','_allelea','_alleleb'] if paired else ['cdr3','_allele'])
    refkeys=(['alpha','beta','_allele_alpha','_allele_beta'] if paired else ['cdr3','_allele'])
    if preset=='jnorm':
        keys+=['ja','jb'] if paired else ['j']
        refkeys+=['j_alpha','j_beta'] if paired else ['j']
    state,hit_ids=grouped_kernel_selection(q,r,[({'epitope':'A'},q,r.filter(pl.col('epitope')=='A'))],
        query_keys=keys,search_keys=refkeys,paired=paired,cutoffs=[0,10],radius=0,
        kernel_scale=90,scorer='v-bounded',geometry_preset=preset,exclude_exact=True)
    assert state['query_id'].to_list()==['stable']
    if target:
        # Pre-puncture same-V count2 includes the exact query reference.
        assert state['kernel_score'][0]==target/(1+control)/(.25*target+.75*2)
        assert hit_ids=={'stable'}
    else:
        assert state['kernel_score'][0] is None
        assert hit_ids==set()
