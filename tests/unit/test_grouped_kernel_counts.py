"""Grouped count preparation keeps the established single/linked geometry."""
import polars as pl
import pytest

from vdjmatch.match import geometry


@pytest.mark.parametrize('paired', [False, True])
def test_grouped_distance_prepares_one_population_and_gene_codes(monkeypatch, paired):
    name='grouped_paired_sum_count_batch' if paired else 'grouped_count_batch'
    calls=[]
    def native(*args, **kwargs):
        calls.append((args,kwargs))
        return iter(())
    monkeypatch.setattr(geometry.gapblock,name,native,raising=False)
    q=pl.DataFrame({'seq0':['CASSLGQAYEQYF'], 'v0':['TRBV19*01']})
    r=pl.DataFrame({'seq0':['CASSLGQAYEQWF'], 'v0':['TRBV19*02']})
    if paired:
        q=q.rename({'seq0':'seq1','v0':'v1'}).with_columns(pl.lit('CAVRDSNYQLIW').alias('seq0'),pl.lit('TRAV12-2*01').alias('v0'))
        r=r.rename({'seq0':'seq1','v0':'v1'}).with_columns(pl.lit('CAVRDSNYQLVW').alias('seq0'),pl.lit('TRAV12-2*02').alias('v0'))
    list(geometry.grouped_total_distance_count_batch(q,r,[0,117000],paired=paired,
        reference_target_offsets=[0,1],reference_target_ids=[0],target_peptide_ids=[0],
        reference_exclusive_peptide_ids=[0],threads=1,ctrim=3,stem_depth=4,matrix='blosum62'))
    assert len(calls)==1
    args,options=calls[0]
    assert options['query_gene_ids']==options['reference_gene_ids']
    assert options['reference_target_offsets']==[0,1]
    assert options['threads']==1
    assert options['feature_offsets'+('_alpha' if paired else '')] is not None
    assert args[-1]==[0,117000]


def test_kernel_hit_exclusion_keeps_reference_and_control_models():
    from vdjmatch.evalue.kernel import grouped_kernel_enrichment
    q=pl.DataFrame({'query_id':['q'],'cdr3':['CASSLGQAYEQYF'],'_allele':['TRBV19*01']})
    r=pl.DataFrame({'cdr3':['CASSLGQAYEQYF','CASSLGQAYEQWF','CWWWWWWWWWWWF'],
        '_allele':['TRBV19*01']*3,'epitope':['A','A','B']})
    groups=[({'epitope':'A','mhc_a':'A','mhc_b':'B2M','mhc_class':'MHCI'},q,r.head(2))]
    def score(exclude):
        return grouped_kernel_enrichment(q,r,groups,query_keys=['cdr3','_allele'],search_keys=['cdr3','_allele'],
            paired=False,cutoffs=[0.,100.],radius=10000.,kernel_scale=10000.,scorer='all-bounded',
            assayed=True,species='human',threads=1,ctrim=3,exclude_exact=exclude)[0]
    full,masked=score(False),score(True)
    assert full['reference_population'].to_list()==masked['reference_population'].to_list()==[2]
    assert full['control_population'].to_list()==masked['control_population'].to_list()==[1]
    assert full['n_reference'][0]==masked['n_reference'][0]+1
    assert masked['exact_reference_count'][0]==0


def test_grouped_exact_exclusion_rejects_old_native_model_semantics(monkeypatch):
    from seqtree import _core
    monkeypatch.setattr(_core,'exact_exclusion_hits_only',False,raising=False)
    empty=pl.DataFrame(schema={'seq0':pl.String,'v0':pl.String})
    with pytest.raises(RuntimeError,match='hit-list-only exact exclusion'):
        geometry.grouped_total_distance_count_batch(empty,empty,[0],exclude_exact=True)


@pytest.mark.parametrize('paired', [False, True])
def test_grouped_same_v_counts_match_existing_conditional_counts(paired):
    import numpy as np
    q=pl.DataFrame({'seq0':['CASSLGQAYEQYF'],'v0':['TRBV19*01']})
    r=pl.DataFrame({'seq0':['CASSLGQAYEQYF','CASSLGRAYEQYF','CASSWGRAYEQYF'],
        'v0':['TRBV19*02','TRBV19*01','TRBV20-1*01']})
    if paired:
        q=q.rename({'seq0':'seq1','v0':'v1'}).with_columns(pl.lit('CAVRDSNYQLIW').alias('seq0'),pl.lit('TRAV12-2*01').alias('v0'))
        r=r.rename({'seq0':'seq1','v0':'v1'}).with_columns(pl.lit('CAVRDSNYQLIW').alias('seq0'),pl.lit('TRAV12-2*02').alias('v0'))
    cuts=[0,600,6000,48000]
    for exclude in (False,True):
        options=dict(paired=paired,ctrim=3,v_loop_weight=0,same_v_gene=True,exclude_exact=exclude,
            reference_target_offsets=[0,1,2,3],reference_target_ids=[0,0,0],
            target_peptide_ids=[0],reference_exclusive_peptide_ids=[0,0,0])
        chunks=list(geometry.grouped_total_distance_count_batch(q,r,cuts,**options))
        expected=[1,1,2,2] if not exclude else [0,0,1,1]
        np.testing.assert_array_equal(np.asarray(chunks[0].counts)[0,1],expected)
    with pytest.raises(ValueError,match='native group barrier'):
        geometry.grouped_total_distance_count_batch(q,r,[2**28],**options)
