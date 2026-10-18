"""Opt-in paired scopes preserve separate established linked and marginal units."""
import polars as pl
from polars.testing import assert_frame_equal
import pytest
from seqtree import Index
from vdjmatch import Annotator
from vdjmatch.match import search_params,load_vdjam
from vdjmatch.runner.multisample import _paired_scope_calls,annotate_sample


def reference():
    pairs=[('p1','CAVVF','CASSF','A','TRAV1','TRBV1'),
        ('p2','CAVAF','CASRF','A','TRAV2','TRBV2'),
        ('p3','CAVVF','CASSF','D','TRAV1','TRBV1')]
    rows=[dict(complex_id=cid,gene=g,cdr3=seq,epitope=ep,v=v,j=j,
        mhc_a='A*02:01',reference_id=cid,record_id=cid+g)
        for cid,a,b,ep,va,vb in pairs for g,seq,v,j in [('TRA',a,va,'TRAJ1'),('TRB',b,vb,'TRBJ1')]]
    rows.extend([dict(complex_id='0',gene='TRA',cdr3='CAVTF',epitope='B',v='TRAV1',j='TRAJ1',mhc_a='A*02:01',reference_id='R1',record_id='u1'),
        dict(complex_id='0',gene='TRB',cdr3='CASTF',epitope='C',v='TRBV1',j='TRBJ1',mhc_a='A*02:01',reference_id='R2',record_id='u2')])
    return pl.DataFrame(rows+[dict(rows[-2],record_id='u3',reference_id='independent')])


def queries():
    return pl.DataFrame({'query_id':['linked','alpha','beta','none','invalid','duplicate'],
        'cdr3a':['CAVVF','CAVTF',None,'CYYYYYYYYYYYYF','*','CAVVF'],
        'cdr3b':['CASSF',None,'CASTF','CYYYYYYYYYYYYF','CASSF','CASSF'],
        'va':['TRAV1',None,None,'TRAV1','TRAV1','TRAV2'],
        'vb':['TRBV1',None,'TRBV1','TRBV1','TRBV1','TRBV2'],
        'ja':['TRAJ1']*6,'jb':['TRBJ1']*6})


@pytest.mark.parametrize('weighted',[False,True])
def test_paired_three_edit_cap_preserves_inner_balls(monkeypatch,weighted):
    ann=Annotator.from_frame(reference());matrix=load_vdjam() if weighted else None
    expected,_=ann.paired_scope_candidates(queries(),ball=True,matrix=matrix)
    original=Index.search_paired_batch_consume;calls=[]
    def consume(self,other,qa,qb,params,callback,**options):
        calls.append(params.max_total_edits)
        return original(self,other,qa,qb,params,callback,**options)
    monkeypatch.setattr(Index,'search_paired_batch_consume',consume)
    actual,_=ann.paired_scope_candidates(queries(),ball=True,matrix=matrix,max_radius=3)
    for axis in actual:
        assert_frame_equal(actual[axis],expected[axis].filter(pl.col('radius')<=3),check_exact=True)
    assert calls==([1,2,3] if weighted else [3])


def controls(labelled):
    out={}
    for chain,sequences in [('TRA',['CAVVF','CAVAF','CAVTF']),('TRB',['CASSF','CASRF','CASTF'])]:
        frame=pl.DataFrame({'cdr3':sequences,'v':[chain+'V1',chain+'V2',chain+'V1'],
            'j':[chain+'J1']*3,'locus':[chain]*3,'species':['human']*3})
        out[chain]=frame if labelled else Index.build(sequences)
    return out


@pytest.mark.parametrize('match_v',[False,True])
@pytest.mark.parametrize('weighted',[False,True])
@pytest.mark.parametrize('exclude_exact',[False,True])
def test_all_families_match_existing_marginal_reducers_and_one_population(monkeypatch,match_v,weighted,exclude_exact):
    ann=Annotator.from_frame(reference());data=queries();matrix=load_vdjam() if weighted else None
    ctrl=controls(match_v)
    # Prepare expected source views once; independent single-chain APIs own their formulas.
    axes,prepared=ann.paired_scope_candidates(data,ball=True,matrix=matrix,match_v=match_v,
        control=ctrl,calibrate=True,species='human',exclude_exact=exclude_exact)
    for chain,q in zip(('TRA','TRB'),prepared):
        expected=ann.graded_candidates(q,locus=chain,matrix=matrix,control=ctrl,
            calibrate=True,species='human',match_v=match_v,exclude_exact=exclude_exact)
        expected=expected.sort('radius','query_id','rank',nulls_last=True)
        assert_frame_equal(axes[chain],expected,check_exact=True)
    calls=[];original=Index.search_paired_batch_consume
    def consume(self,other,qa,qb,params,callback,**kwargs):
        calls.append((len(qa),len(qb),kwargs['skip_empty']))
        return original(self,other,qa,qb,params,callback,**kwargs)
    monkeypatch.setattr(Index,'search_paired_batch_consume',consume)
    monkeypatch.setattr(Index,'search_batch_consume',lambda *a,**k:pytest.fail('extra single-chain reference search'))
    actual,_=ann.paired_scope_candidates(data,ball=True,matrix=matrix,match_v=match_v,
        control=ctrl,calibrate=True,species='human',exclude_exact=exclude_exact,threads=2)
    for axis in axes:assert_frame_equal(actual[axis],axes[axis],check_exact=True)
    assert calls==[(6,6,True)]*(5 if weighted else 1)
    if match_v:
        restricted=actual['linked'].filter(pl.col('query_id')=='duplicate')
        assert restricted['epitope'].unique().to_list()==['A']
        assert restricted['n_hits'].to_list()==[1]*5 and restricted['n_records'].to_list()==[1]*5
        assert actual['linked'].filter(pl.col('query_id')=='linked').height==(0 if exclude_exact else 10)
    assert actual['linked']['E'].null_count()==actual['linked'].height
    assert actual['linked']['p_enrichment'].null_count()==actual['linked'].height
    assert actual['linked']['calibration'].unique().to_list()==['uncalibrated']


def test_fixed_linked_parity_and_label_contract_with_incomplete_pairs():
    ann=Annotator.from_frame(reference());data=queries();sp=search_params('0',matrix='',gap_open=1,gap_extend=1)
    axes,_=ann.paired_scope_candidates(data,scope=sp)
    expected=ann.paired_candidates(data,cdr3a='cdr3a',cdr3b='cdr3b',scope=sp,score_scale=1.)
    assert_frame_equal(axes['linked'].drop('radius'),expected,check_exact=True)
    calls=_paired_scope_calls(data,axes,species='human')
    assert set(calls.filter(pl.col('query_id')=='linked')['paired_epitope'])=={'A','D'}
    assert calls.filter(pl.col('query_id')=='linked')['paired_label_source'].unique().to_list()==['linked']
    assert calls.filter(pl.col('query_id')=='alpha')['paired_epitope'].to_list()==['B']
    assert calls.filter(pl.col('query_id')=='beta')['paired_epitope'].to_list()==['C']
    assert calls.filter(pl.col('query_id')=='none')['paired_label_status'].to_list()==['no_support']
    assert calls['E'].null_count()==calls.height and not any(c.startswith('kernel_') for c in calls.columns)


def test_runner_preserves_flags_nohits_and_opt_in(tmp_path):
    ann=Annotator.from_frame(reference());source=tmp_path/'pairs.tsv'
    q=queries();q.select('query_id',pl.col('cdr3a').alias('junction_aa_alpha'),pl.col('cdr3b').alias('junction_aa_beta'),
        pl.col('va').alias('v_alpha'),pl.col('vb').alias('v_beta'),pl.col('ja').alias('j_alpha'),pl.col('jb').alias('j_beta')).write_csv(source,separator='\t')
    result=annotate_sample(ann._index,source,paired=True,paired_marginals=True,search_mode='ball',
        with_evalue=False,match_v=True,species='human',retain_hits=False)
    calls=result['calls']
    assert set(calls['query_id'])==set(range(6)) and set(calls['radius'])==set(range(1,6))
    assert calls.filter(pl.col('source_query_ida')=='alpha')['alpha_v_imputed'].all()
    assert calls.filter(pl.col('source_query_ida')=='alpha')['beta_query_valid'].not_().all()
    assert calls.filter(pl.col('source_query_ida')=='none')['paired_label_status'].unique().to_list()==['no_support']
    with pytest.raises(ValueError,match='paired graded balls'):
        annotate_sample(ann._index,source,paired=True,search_mode='ball',with_evalue=False)


def test_unlinked_union_and_wholly_unsupported_ball_rows():
    ann=Annotator.from_frame(reference())
    q=pl.DataFrame({'query_id':['union','nohit'],'cdr3a':['CAVTF','CYYYYYYYYYYYYF'],
        'cdr3b':['CASTF','CYYYYYYYYYYYYF']})
    axes,_=ann.paired_scope_candidates(q,scope=search_params('0'))
    assert axes['linked'].is_empty()
    calls=_paired_scope_calls(q,axes,species='human')
    assert set(calls.filter(pl.col('query_id')=='union')['paired_epitope'])=={'B','C'}
    assert calls.filter(pl.col('query_id')=='union')['paired_label_source'].unique().to_list()==['marginal_fallback']
    nohit=q.tail(1);axes,_=ann.paired_scope_candidates(nohit,ball=True)
    calls=_paired_scope_calls(nohit,axes,species='human',radii=range(1,6))
    assert calls.height==5 and calls['radius'].to_list()==list(range(1,6))
    assert calls['paired_label_status'].to_list()==['no_support']*5


def test_shipping_cli_exports_opt_in_axes_and_explicit_uncalibrated_manifest(tmp_path):
    import json
    from vdjmatch.cli.__main__ import main
    ref=tmp_path/'reference.tsv';source=tmp_path/'cells.tsv';prefix=tmp_path/'out'
    reference().rename({'cdr3':'junction_aa','gene':'locus','v':'v_call','j':'j_call'}).with_columns(
        pl.lit('HomoSapiens').alias('species')).write_csv(ref,separator='\t')
    queries().select('query_id',pl.col('cdr3a').alias('junction_aa_alpha'),pl.col('cdr3b').alias('junction_aa_beta'),
        pl.col('va').alias('v_alpha'),pl.col('vb').alias('v_beta')).write_csv(source,separator='\t')
    assert main(['match',str(source),'--vdjdb',str(ref),'--fresh-reference','--paired','--paired-marginals',
        '--search-mode','ball','--max-radius','3','--matrix','none','--match-v','--no-evalue','--no-hit-table',
        '--output-prefix',str(prefix)])==0
    for axis in ('candidates','TRA.candidates','TRB.candidates','calls'):
        assert (tmp_path/f'out.cells.{axis}.txt').is_file()
    manifest=json.loads((tmp_path/'out.cells.manifest.json').read_text())
    assert manifest['paired_marginals'] and manifest['paired_calibration_model'] is None
    assert manifest['paired_scope_model']['selection_calibration'].startswith('unavailable')
    calls=pl.read_csv(tmp_path/'out.cells.calls.txt',separator='\t')
    assert manifest['radii']==[1,2,3] and manifest['ball_model']['max_radius']==3
    assert set(calls['radius'])=={1,2,3} and calls['E'].null_count()==calls.height


def test_missing_reference_axis_does_not_veto_supported_unlinked_marginal():
    ann=Annotator.from_frame(reference().filter(pl.col('gene')=='TRA').with_columns(pl.lit('0').alias('complex_id')))
    q=pl.DataFrame({'query_id':['one'],'cdr3a':['CAVTF'],'cdr3b':['CASTF']})
    axes,_=ann.paired_scope_candidates(q,scope=search_params('0'))
    assert axes['linked'].is_empty() and axes['TRB'].is_empty()
    calls=_paired_scope_calls(q,axes,species='human')
    assert calls['paired_epitope'].to_list()==['B'] and calls['paired_label_source'].to_list()==['marginal_fallback']


def test_paired_scope_capability_fails_before_control_search(monkeypatch):
    ann=Annotator.from_frame(reference())
    monkeypatch.setattr(Index,'search_paired_batch_consume',None)
    monkeypatch.setattr(Index,'edit_histogram_batch',lambda *a,**k:pytest.fail('controls searched before capability check'))
    with pytest.raises(RuntimeError,match='skip_empty'):
        ann.paired_scope_candidates(queries(),ball=True,control=controls(False))
