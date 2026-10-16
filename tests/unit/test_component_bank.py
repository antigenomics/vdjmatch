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
    assert all(value.shape==(3,2) for value in bank.values.values())
    monkeypatch.setattr(component_bank,'_BATCH_ROWS',4)
    parts=list(bank.frames())
    assert all(part.height<=4 for part in parts)
    actual=pl.concat(parts).sort('query_id',*PMHC)
    assert_frame_equal(actual,expected.select(actual.columns).sort('query_id',*PMHC),check_exact=True)
    subset=pl.concat(bank.frames(q.filter(pl.col('query_id').is_in(['far','q']))))
    assert subset.height==4
    assert_frame_equal(subset.sort('query_id',*PMHC),actual.filter(pl.col('query_id').is_in(['far','q'])).sort(
        'query_id',*PMHC),check_exact=True)
    assert all(value.dtype in (np.dtype('float64'),np.dtype('bool')) for value in bank.values.values())


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
