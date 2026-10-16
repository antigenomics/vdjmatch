"""Numeric raw components preserve all source cases without a string case grid."""
import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal
from seqtree import Index
from vdjmatch.aggregate.candidates import PMHC


@pytest.mark.parametrize('paired',[False,True])
@pytest.mark.parametrize('exclude',[False,True])
@pytest.mark.parametrize('species',['human','mouse'])
def test_numeric_bank_matches_materialized_components_and_keeps_repeated_sources(monkeypatch,paired,exclude,species):
    from vdjmatch.runner import combined
    from vdjmatch.runner import component_bank
    q=pl.DataFrame({'query_id':['q','copy','missing','far'],
        'cdr3':['CASSMGQETQYF','CASSMGQETQYF','CASSLGQETQYF','W'*200],
        'v':['TRBV1','TRBV1',None,'TRBV1'],'j':['TRBJ1-1']*4,'locus':['TRB']*4})
    r=pl.DataFrame({'cdr3':['CASSLGQETQYF','CASSIGQETQYF','CATSLGQETQYF'],
        'v':['TRBV1']*3,'j':['TRBJ1-1']*3,'gene':['TRB']*3,
        'species':['HomoSapiens' if species=='human' else 'MusMusculus']*3,
        'epitope':['A','A','B'],'mhc_a':[None]*3,'mhc_b':['B2M']*3,'mhc_class':['MHCI']*3,
        'complex_id':['1','2','0']}).with_columns(pl.col('mhc_a').cast(pl.String))
    context={'density_index':Index.build(r['cdr3'].to_list()+['CASSAGQETQYF']),
        'prior_background':r.select('cdr3','v','j')}
    locus='TRB'
    if paired:
        locus='paired'
        av='TRAV1-2' if species=='human' else 'TRAV1'
        q=q.drop('locus').rename({'cdr3':'cdr3b','v':'vb','j':'jb'}).with_columns(
            pl.lit('CAVRDSNYQLIW').alias('cdr3a'),pl.lit(av).alias('va'),pl.lit('TRAJ33').alias('ja'))
        alpha=r.head(2).with_columns(pl.Series('cdr3',['CAVRDTNYQLIW','CAVRDMNYQLIW']),
            pl.lit(av).alias('v'),pl.lit('TRAJ33').alias('j'),pl.lit('TRA').alias('gene'))
        r=pl.concat([r,alpha])
        context['density_index']={'TRA':Index.build(alpha['cdr3'].to_list()),'TRB':context['density_index']}
    expected=combined.score_components(q,r,context,species=species,locus=locus,exclude_exact=exclude)
    def forbidden(*args,**kwargs):
        pytest.fail('numeric raw scoring must not construct a whole string case grid')
    monkeypatch.setattr(combined,'_cases',forbidden)
    bank=component_bank.score_component_bank(q,r,context,species=species,locus=locus,exclude_exact=exclude)
    assert bank.mapping['_query_row'].n_unique()==3
    assert all(not isinstance(value,np.ndarray) for value in bank.values.values())
    monkeypatch.setattr(component_bank,'_BATCH_ROWS',4)
    parts=list(bank.frames())
    assert all(part.height<=4 for part in parts)
    actual=pl.concat(parts).sort('query_id',*PMHC)
    assert_frame_equal(actual,expected.select(actual.columns).sort('query_id',*PMHC),check_exact=True)
    subset=pl.concat(bank.frames(q.filter(pl.col('query_id').is_in(['far','q']))))
    assert subset.height==4
    assert_frame_equal(subset.sort('query_id',*PMHC),actual.filter(pl.col('query_id').is_in(['far','q'])).sort(
        'query_id',*PMHC),check_exact=True)
    assert all(actual.schema[c] in (pl.Float64,pl.Boolean) for c in bank.values)


def test_component_array_budget_fails_before_scoring(monkeypatch):
    from vdjmatch.runner import component_bank,combined
    q=pl.DataFrame({'query_id':['q'],'cdr3':['CASSF'],'v':['TRBV1'],'j':['TRBJ1-1']})
    r=pl.DataFrame({'cdr3':['CASSF'],'v':['TRBV1'],'j':['TRBJ1-1'],'gene':['TRB'],'epitope':['A']})
    monkeypatch.setattr(component_bank,'_ARRAY_LIMIT_BYTES',1)
    def forbidden(*args,**kwargs):
        pytest.fail('array allocation budget must fail before native preparation')
    monkeypatch.setattr(combined,'_original',forbidden)
    with pytest.raises(ValueError,match='component array budget'):
        component_bank.score_component_bank(q,r,{'density_index':None,'prior_background':None},
            species='human',locus='TRB')


@pytest.mark.parametrize('include_length',[False,True])
def test_fixed_prior_reconstruction_preserves_unknowns_weights_and_500_activation(include_length):
    from itertools import product
    from vdjmatch.runner.component_bank import _original_model,_SparseScores
    from vdjmatch.match.historical import _prepare_prior_background,germline_prior,unpaired_score
    sequences=['CASS'+''.join(x)+'ETQYF' for x in list(product('ACDEFGHIKLMNPQRSTVWY',repeat=3))[:500]]
    ref=pl.DataFrame({'cdr3':sequences+sequences[:499],'v':['TRBV1']*999,'j':['TRBJ1-1']*999,
        'epitope':['A']*500+['B']*499,'mhc_a':['A']*999,'mhc_b':['B2M']*999,'mhc_class':['MHCI']*999})
    targets=ref.select(PMHC).unique(maintain_order=True).with_row_index('_bank_target')
    q=pl.DataFrame({'query_id':range(6),'cdr3':['CASSLGQETQYF','CASSMGQETQYF','CAVRDSNYQLIW','W'*80,'CA*F',None],
        'v':['TRBV1*01','TRBV1','TRBV2',None,'TRBV1','TRBV1'],'j':['TRBJ1-1']*6})
    background=pl.DataFrame({'cdr3':['CASSLGQETQYF','CAVRDSNYQLIW'],'v':['TRBV1','TRBV2'],
        'j':['TRBJ1-1','TRBJ2-1'],'event_weight':[2**32+13,17]})
    model=_original_model(q,ref,targets,_prepare_prior_background(background),include_length,lambda size:None)
    # One supported density at a missing-gene query must not fabricate a score.
    model.density=_SparseScores(np.array([0,1,6,7],dtype=np.uint64),np.ones(4),2,0.)
    order=np.array([4,1,0,5,2,3,1],dtype=np.uint32)
    density,prior,score,support=model.evaluate(order)
    expected_prior=np.full((6,2),np.nan);expected_score=np.full((6,2),np.nan)
    valid=q.filter(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$')&pl.col('v').is_not_null())
    all_density=model.density.take(np.arange(6,dtype=np.uint32)).reshape(6,2)
    for i,label in enumerate(['A','B']):
        r=ref.filter(pl.col('epitope')==label)
        p=germline_prior(valid,r,background,include_length=include_length)
        indices=p['query_id'].to_numpy()
        expected_prior[indices,i]=p['germline_lr'].to_numpy()
        expected_score[indices,i]=unpaired_score(all_density[indices,i],p['germline_lr'],r.height,log_scale=True).to_numpy()
    np.testing.assert_array_equal(prior,expected_prior[order].reshape(-1))
    np.testing.assert_array_equal(score,expected_score[order].reshape(-1))
    np.testing.assert_array_equal(support,density>0)
    assert support.reshape(-1,2)[5].all() and np.isnan(score.reshape(-1,2)[5]).all()
