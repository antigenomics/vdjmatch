"""Declared output is sparse; all labels still determine both randomized maxima."""
import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal


@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('species', ['human', 'mouse'])
def test_actual_compact_matches_full_declared_and_database_tails(monkeypatch, paired, species):
    from vdjmatch.evalue import combined
    from vdjmatch.aggregate.candidates import PMHC
    q=pl.DataFrame({'query_id':['near','far'],'cdr3':['CASSMGQETQYF','W'*80],
        'v':['TRBV1']*2,'j':['TRBJ1-1']*2,'locus':['TRB']*2})
    r=pl.DataFrame({'cdr3':['CASSLGQETQYF','CASSIGQETQYF','CATSLGQETQYF'],
        'v':['TRBV1']*3,'j':['TRBJ1-1']*3,'gene':['TRB']*3,
        'species':['HomoSapiens' if species == 'human' else 'MusMusculus']*3,
        'epitope':['A','A','B'],'mhc_a':['HLA-A*02:01']*3,'mhc_b':['B2M']*3,'mhc_class':['MHCI']*3})
    control=pl.DataFrame({'junction_aa':['CASSAGQETQYF','CASSVGQETQYF','CASSTGQETQYF'],
        'v_call':['TRBV1']*3,'j_call':['TRBJ1-1']*3})
    monkeypatch.setattr(combined,'compute_pgen',lambda ref,events,**kw:pl.concat([
        ref.select('junction_aa'),events.select('junction_aa')]).unique().with_columns(pl.lit(.01).alias('pgen')))
    if paired:
        av = 'TRAV1-2' if species == 'human' else 'TRAV1'
        q = q.drop('locus').rename({'cdr3':'cdr3b','v':'vb','j':'jb'}).with_columns(
            pl.Series('cdr3a',['CAVRDSNYQLIW','W'*80]),pl.lit(av).alias('va'),pl.lit('TRAJ33').alias('ja'))
        r = r.with_columns(pl.Series('complex_id',['1','2','0']))
        alpha = r.head(2).with_columns(pl.Series('cdr3',['CAVRDTNYQLIW','CAVRDMNYQLIW']),
            pl.lit(av).alias('v'),pl.lit('TRAJ33').alias('j'),pl.lit('TRA').alias('gene'))
        r = pl.concat([r,alpha])
        control = pl.concat([control,control.with_columns(
            pl.Series('junction_aa',['CAVRDSNYQLIW','CAVRDTNYQLIW','CAVRDMNYQLIW']),
            pl.lit(av).alias('v_call'),pl.lit('TRAJ33').alias('j_call'))])
    args=dict(species=species,locus='paired' if paired else 'TRB',replicates=2,seed=7,exclude_exact=True)
    full,full_tails=combined.calibrate_combined(q,r,control,**args)
    from vdjmatch.runner import component_bank
    monkeypatch.setattr(component_bank,'_BATCH_ROWS',6)  # test and observed rows share a batch
    requested=q.select('query_id').join(r.head(1).select(PMHC),how='cross')
    detail,tails=combined.calibrate_combined(q,r,control,detail_targets=requested,**args)
    keys=['query_id',*PMHC]
    expected=full.join(requested,on=keys,how='semi',nulls_equal=True)
    assert_frame_equal(detail.select(full.columns).sort(keys),expected.sort(keys),check_exact=True)
    expected_tails=full_tails.filter((pl.col('test_scope')=='database')|(pl.col('epitope')=='A'))
    assert_frame_equal(tails.select(full_tails.columns).sort('test_scope',*keys),
        expected_tails.sort('test_scope',*keys),check_exact=True)
    assert detail['declared_target'].to_list()==[True,True]
    assert not detail['reported_specificity'].any()
    assert detail['p_any'].min()>=1/3
    # An absent declared label remains an explicit no-reference diagnostic.
    missing=requested.head(1).with_columns(pl.lit('MISSING').alias('epitope'))
    detail,tails=combined.calibrate_combined(q,r,control,detail_targets=missing,**args)
    assert detail.height==1 and detail['combined_status'].item()=='no_reference_target'
    assert detail['p_any'].item()==1. and detail['declared_target'].item()
    assert tails.filter(pl.col('test_scope')=='epitope')['p_both'].item()==1.
    # Empty declarations still calibrate every label and retain complete database tails.
    detail,tails=combined.calibrate_combined(q,r,control,detail_targets=requested.head(0),**args)
    assert detail.height == 0 and 'reported_specificity' in detail.columns
    assert_frame_equal(tails.select(full_tails.columns).sort('query_id'),
        full_tails.filter(pl.col('test_scope')=='database').sort('query_id'),check_exact=True)


def test_complete_case_bitset_rejects_duplicates_unknown_and_missing():
    import pytest
    from vdjmatch.evalue.compact_combined import _mark,_complete
    seen=np.zeros(2,dtype=np.uint8)
    _mark(seen,np.repeat(np.arange(3),3),np.tile(np.arange(3),3),3,3)
    _complete(seen,9)
    with pytest.raises(RuntimeError,match='repeats'):
        _mark(seen,np.array([0]),np.array([0]),3,3)
    with pytest.raises(RuntimeError,match='unknown'):
        _mark(seen,np.array([3]),np.array([0]),3,3)
    with pytest.raises(RuntimeError,match='unknown'):
        _mark(seen,np.array([-1]),np.array([0]),3,3)
    with pytest.raises(RuntimeError,match='unknown'):
        _mark(seen,np.array([0]),np.array([-1]),3,3)
    with pytest.raises(RuntimeError,match='complete observed'):
        _complete(np.zeros(2,dtype=np.uint8),9)


@pytest.mark.parametrize('threshold', [.004, .005, .01])
def test_reported_ties_and_undeclared_rival_enter_both_any_call_tails(threshold, monkeypatch):
    from vdjmatch.evalue.compact_combined import calibrate_compact
    from vdjmatch.runner.component_bank import ComponentBank
    from vdjmatch.aggregate.candidates import PMHC
    r=pl.DataFrame({'epitope':['A','B','C'],'mhc_a':['A']*3,'mhc_b':['B2M']*3,'mhc_class':['MHCI']*3,
        'model':[0]*3})
    q=pl.DataFrame({'query_id':['near','far'],'pgen_bin':['bin']*2,'pgen':[.01]*2})
    normalizers=pl.DataFrame({'query_id':['norm'+str(i) for i in range(199)],'replicate':range(199)})
    tests=pl.DataFrame({'query_id':['test'+str(i) for i in range(199)],'replicate':range(199),'pgen_bin':['bin']*199})
    requested=pl.DataFrame({'query_id':['near','near','far'],'epitope':['A','C','A'],
        'mhc_a':['A']*3,'mhc_b':['B2M']*3,'mhc_class':['MHCI']*3})
    original_full = np.full
    def no_dense_nulls(shape, *args, **kwargs):
        assert shape != (199, 3), 'query-null storage must omit the dense draw/label matrix'
        return original_full(shape, *args, **kwargs)
    monkeypatch.setattr(np, 'full', no_dense_nulls)
    for rival in (False,True,None):
        def components(rows,model):
            values=np.zeros((rows.height,3))
            if model['model'][0]==0:
                for i,name in enumerate(rows['query_id']):
                    if name=='near':
                        values[i]=[10,10,np.nan]
                    elif name=='far':
                        values[i]=np.nan
                    elif rival:
                        values[i,1]=20
                    elif rival is None:
                        values[i]=np.nan
            elif rival:
                values[:,1]=20
            elif rival is None:
                values[:]=np.nan
            return ComponentBank(rows.select('query_id').with_row_index('_query_row'),model.select(PMHC),
                {'score':values.copy(),'call_score':values})
        def combine(raw,norm,model,rows,*,_consume,**kwargs):
            for frame in raw.frames():
                _consume(frame)
        detail,tails=calibrate_compact(q,normalizers,tests,r,components=components,combine=combine,
            models=(r.with_columns(pl.lit(i).alias('model')) for i in range(1,200)),
            detail_targets=requested,replicates=199,report_p=threshold,pgen_columns=['pgen'],combine_options={})
        near_a=detail.filter((pl.col('query_id')=='near')&(pl.col('epitope')=='A'))
        assert near_a['p_any'].item()==(1. if rival else .005)
        passes = not rival and threshold >= .005
        assert near_a['reported_specificity'].item()==passes
        fixed=tails.filter((pl.col('query_id')=='near')&(pl.col('epitope')=='A'))
        assert fixed['p_both'].item()==.005  # Fixed-label evidence alone would miss the rival.
        assert fixed['query_hit_draws'].item() == (0 if rival is None else 199)
        assert fixed['reference_hit_draws'].item() == (0 if rival is None else 199)
        reported=detail.filter(pl.col('reported_specificity'))
        assert sorted(reported['epitope'])==(['A','B'] if passes else [])
        assert detail.filter(pl.col('declared_target')).height==3
        assert detail.height==(4 if passes else 3)


def test_sparse_query_null_budget_keeps_nohits_and_stops_before_overallocation(monkeypatch):
    from vdjmatch.evalue import compact_combined as compact
    from vdjmatch.runner.component_bank import ComponentBank
    from vdjmatch.aggregate.candidates import PMHC
    r=pl.DataFrame({'epitope':['A'],'mhc_a':['A'],'mhc_b':['B2M'],'mhc_class':['MHCI']})
    q=pl.DataFrame({'query_id':['q'],'pgen_bin':['bin']})
    normalizers=pl.DataFrame({'query_id':['norm'],'replicate':[0]})
    tests=pl.DataFrame({'query_id':['test'],'replicate':[0],'pgen_bin':['bin']})
    requested=q.select('query_id').join(r,how='cross')
    calls=[];finite=False
    def components(rows,model):
        calls.append(rows.height)
        scores=np.full((rows.height,1),0. if finite else np.nan)
        return ComponentBank(rows.select('query_id').with_row_index('_query_row'),model.select(PMHC),
            {'score':scores,'call_score':scores})
    def combine(raw,norm,model,rows,*,_consume,**kwargs):
        for frame in raw.frames():_consume(frame)
    def run():
        return compact.calibrate_compact(q,normalizers,tests,r,components=components,combine=combine,
            models=iter([r]),detail_targets=requested,replicates=1,report_p=.01,pgen_columns=[],combine_options={})
    base=compact._query_null_bytes(1,1,1)
    monkeypatch.setattr(compact,'_QUERY_NULL_LIMIT_BYTES',base-1)
    with pytest.raises(ValueError,match='query-null buffer budget'):run()
    assert not calls
    monkeypatch.setattr(compact,'_QUERY_NULL_LIMIT_BYTES',base)
    detail,tails=run()
    assert detail['p_any'].item()==1. and tails['query_draws'].to_list()==[1,1]
    assert tails['query_hit_draws'].to_list()==[0,0]
    finite=True;calls.clear()
    with pytest.raises(ValueError,match='query-null buffer budget'):run()
    assert len(calls)==1  # fail on finite storage before any mock is built
    assert compact._query_null_bytes(89393,309843,2201)<2*1024**3
