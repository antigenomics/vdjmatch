"""Fresh all-label policy routing equals the frozen selector, including tail defaults."""
import json
import numpy as np
import polars as pl
from polars.testing import assert_frame_equal
import pytest
from vdjmatch.cli.score_policy import FAST_ROUTES,apply_reference_policy
from vdjmatch.evalue import reference_rule_bank as bank


@pytest.mark.parametrize('paired',[False,True])
def test_stream_policy_equals_selector_and_keeps_finite_fallbacks(monkeypatch,paired):
    lane='PAIRED' if paired else 'TRB'
    policy=pl.DataFrame({'species':['human']*6,'locus':[lane]*6,'epitope':list('ABCDFG'),
        'reference_class':['v_only','motif_only','v_and_motif','no_supported_axis','motif_only','v_and_motif']})
    queries=pl.DataFrame({'query_id':pl.Series([17,99],dtype=pl.UInt64),'seq0':['CASSAAAAF','CASSCCCCF'],'v0':['TRBV19*01']*2})
    reference=pl.DataFrame({'seq0':['CASSAAAAF']*7,'v0':['TRBV19*01']*7,'epitope':list('ABCDEFG')})
    grid=pl.DataFrame({'query_id':pl.Series([17]*7+[99]*7,dtype=pl.UInt64),'epitope':list('ABCDEFG')*2,
        'rule_v_score':[2.,2.,2.,2.,2.,0.,2.]+[None]*7,
        'rule_exact3_score':[6.,6.,6.,6.,6.,6.,float('inf')]+[None]*7,
        'kernel_v_bounded':[10.,10.,10.,10.,50.,10.,10.]+[0.]*7,
        '_native_supported':[True]*7+[False]*7})
    calls=[]
    def rules(q,r,*,consume,**kwargs):
        calls.append('rules')
        for name in ('rule_v','rule_exact3','rule_rank'):
            column=name+'_score' if name!='rule_rank' else 'kernel_v_bounded'
            consume(name,grid.select('query_id','epitope',pl.col(column).alias('score'),pl.col(column).is_not_null().alias('supported')))
        return {'native_batches':1}
    def kernel(q,r,*,consume,**kwargs):
        calls.append('kernel')
        consume('kernel_v_bounded',grid.select('query_id','epitope',pl.col('kernel_v_bounded').alias('score'),pl.col('_native_supported').alias('supported')))
        consume('kernel_other',grid.select('query_id','epitope',pl.lit(999.).alias('score'),pl.lit(True).alias('supported')))
        return {'native_batches':1}
    monkeypatch.setattr(bank,'stream_reference_rules',rules);monkeypatch.setattr(bank,'stream_reference_kernels',kernel)
    emitted=[]
    info=bank.stream_reference_policy(queries,reference,policy,locus=lane,paired=paired,
        consume=lambda name,frame:emitted.append((name,frame)))
    assert calls==(['kernel'] if paired else ['rules','kernel'])
    assert [name for name,_ in emitted]==['reference_class_policy']
    expected=apply_reference_policy(grid.with_columns(pl.lit('HomoSapiens').alias('species'),
        pl.lit(lane).alias('locus')),policy,original_column='kernel_v_bounded').select(
            'query_id','epitope',pl.col('policy_score').alias('score'),
            (pl.col('policy_score').is_finite() & pl.when(pl.col('selected_score_column')=='kernel_v_bounded')
                .then(pl.col('_native_supported')).otherwise(pl.lit(True))).fill_null(False).alias('supported'),
            'reference_class','selected_score_column','policy_status')
    assert_frame_equal(emitted[0][1],expected,check_exact=True)
    assert emitted[0][1].filter(pl.col('query_id')==99)['supported'].not_().all()
    assert info['query_target_rows']==14


def _policy(tmp_path,source,policy):
    from vdjmatch import db
    from vdjmatch.db.cache import sha256
    directory=tmp_path/'model';directory.mkdir()
    model=directory/'reference-policy.tsv';policy.write_csv(model,separator='\t')
    (directory/'manifest.json').write_text(json.dumps(dict(reference_model_asset=True,no_benchmark_inputs=True,
        routes=FAST_ROUTES,source_reference_identity={'human':db.provenance(source)},
        outputs_sha256={model.name:sha256(model)})))
    return model


def _run(source,prefix,model,*extra):
    from argparse import ArgumentParser
    from vdjmatch.cli.reference_rule_scores import register
    parser=ArgumentParser();register(parser.add_subparsers())
    args=parser.parse_args(['reference-rule-scores','--vdjdb',str(source),'--locus','TRB',
        '--class-policy',str(model),'--output-prefix',str(prefix),*extra])
    return args.func(args)


@pytest.mark.parametrize('tail_score',[2.,50.])
def test_cli_tail_winner_all_cases_and_frozen_policy(tmp_path,monkeypatch,tail_score):
    # Eligible observations define queries; the sparse tail remains a competing model.
    raw=pl.DataFrame({'cdr3':['CASSLGQAYEQYF']*30+['CASSLGRAYEQYF'],
        'v':['TRBV19']*31,'j':['TRBJ2-7']*31,'gene':['TRB']*31,
        'species':['HomoSapiens']*31,'epitope':['PRIMARY']*30+['TAIL']})
    source=tmp_path/'reference.tsv';raw.write_csv(source,separator='\t')
    model=_policy(tmp_path,source,pl.DataFrame({'species':['human'],'locus':['TRB'],
        'epitope':['PRIMARY'],'reference_class':['v_only']}))
    def rules(q,r,*,consume,**kwargs):
        assert set(r['epitope'])=={'PRIMARY','TAIL'}
        grid=q.select('query_id').join(r.select('epitope').unique(),how='cross')
        for scorer in ('rule_v','rule_exact3'):
            consume(scorer,grid.with_columns(pl.lit(2.).alias('score'),pl.lit(True).alias('supported')))
        return {'native_batches':1}
    def kernel(q,r,*,consume,**kwargs):
        assert set(r['epitope'])=={'PRIMARY','TAIL'}
        grid=q.select('query_id').join(r.select('epitope').unique(),how='cross')
        consume('kernel_v_bounded',grid.with_columns(pl.when(pl.col('epitope')=='TAIL').then(tail_score).otherwise(1.).alias('score'),
            pl.lit(True).alias('supported')))
        return {'native_batches':1}
    monkeypatch.setattr(bank,'stream_reference_rules',rules);monkeypatch.setattr(bank,'stream_reference_kernels',kernel)
    prefix=tmp_path/'out';assert _run(source,prefix,model)==0
    winners=pl.read_csv(tmp_path/'out.scores.tsv',separator='\t')
    continuous=pl.read_parquet(tmp_path/'out.continuous/*.parquet')
    audit=pl.read_csv(tmp_path/'out.policy-audit.tsv',separator='\t')
    manifest=json.loads((tmp_path/'out.manifest.json').read_text())
    assert set(winners['epitope'])==({'TAIL'} if tail_score>2 else {'PRIMARY','TAIL'})
    assert winners['scorer'].unique().to_list()==['reference_class_policy']
    assert continuous['epitope'].to_list()==['PRIMARY'] and continuous['score'].to_list()==[2.]
    assert audit['query_target_cases'].sum()==2 and set(audit['epitope'])=={'PRIMARY','TAIL'}
    assert audit.filter(pl.col('epitope')=='TAIL')['policy_status'].to_list()==['no_reference_profile']
    assert manifest['scorers']==['reference_class_policy'] and manifest['candidate_targets']==2
    assert manifest['continuous_rows']==manifest['expected_continuous_rows']==1
    assert manifest['class_policy']['inputs_frozen'] and manifest['class_policy']['calibration_status'].startswith('OPEN')
    assert not manifest['class_policy']['default_changed']
    with pytest.raises(ValueError,match='separate calculations'):_run(source,tmp_path/'exclusive',model,'--kernel-only')
    original=model.read_text()
    model.write_text(original+'\n')
    with pytest.raises(ValueError,match='frozen reference-only'):_run(source,tmp_path/'changed',model)
    model.write_text(original)
    source.write_text(source.read_text()+'\n')
    with pytest.raises(ValueError,match='different raw reference'):_run(source,tmp_path/'changed-source',model)


def test_actual_fresh_traversals_match_selector_with_unavailable_coordinates(tmp_path,monkeypatch):
    raw=pl.DataFrame({'cdr3':['CASSLGQAYEQYF']*2+['CASSLGRAYEQYF']*2+['CASSWGRAYEQYF']+['CASSF']*2,
        'v':['TRBV19']*4+['TRBV5-1']+['TRBV19']*2,'j':['TRBJ2-7']*7,'gene':['TRB']*7,
        'species':['HomoSapiens']*7,'epitope':['A']*2+['B']*2+['TAIL']+['SHORT']*2})
    source=tmp_path/'reference.tsv';raw.write_csv(source,separator='\t')
    policy=pl.DataFrame({'species':['human']*3,'locus':['TRB']*3,'epitope':['A','B','SHORT'],
        'reference_class':['v_only','motif_only','v_and_motif']})
    model=_policy(tmp_path,source,policy)
    rules,kernels=bank.stream_reference_rules,bank.stream_reference_kernels
    captured={};calls=[]
    def rule(q,r,*,consume,**kwargs):
        calls.append('rules')
        def emit(name,frame):
            if name in ('rule_v','rule_exact3'):captured.setdefault(name,[]).append(frame)
            consume(name,frame)
        return rules(q,r,consume=emit,**kwargs)
    def kernel(q,r,*,consume,**kwargs):
        calls.append('kernel')
        def emit(name,frame):
            if name=='kernel_v_bounded':captured.setdefault(name,[]).append(frame)
            consume(name,frame)
        return kernels(q,r,consume=emit,**kwargs)
    monkeypatch.setattr(bank,'stream_reference_rules',rule);monkeypatch.setattr(bank,'stream_reference_kernels',kernel)
    prefix=tmp_path/'out';assert _run(source,prefix,model,'--min-records','2')==0
    assert calls==['rules','kernel']
    wide=pl.concat(captured['kernel_v_bounded']).rename({'score':'kernel_v_bounded','supported':'_native_supported'})
    for name in ('rule_v','rule_exact3'):
        wide=wide.join(pl.concat(captured[name]).select('query_id','epitope',pl.col('score').alias(name+'_score')),
            on=['query_id','epitope'],validate='1:1')
    expected=apply_reference_policy(wide.with_columns(pl.lit('HomoSapiens').alias('species'),
        pl.lit('TRB').alias('locus')),policy,original_column='kernel_v_bounded').filter(pl.col('epitope').is_in(['A','B']))
    continuous=pl.read_parquet(tmp_path/'out.continuous/*.parquet')
    columns=['query_id','epitope','reference_class','selected_score_column','policy_status']
    actual=continuous.join(expected.select('query_id','epitope'),on=['query_id','epitope'],how='semi')
    assert_frame_equal(actual.select(*columns,'score').sort('query_id','epitope'),
        expected.select(*columns,pl.col('policy_score').alias('score')).sort('query_id','epitope'),check_exact=True)
    q=pl.read_csv(tmp_path/'out.queries.tsv',separator='\t')
    assert q.height==3 and q['model_available'].sum()==2 and continuous.height==9
    short=q.filter(~pl.col('model_available'))['query_id'].item()
    assert continuous.filter(pl.col('query_id')==short)['policy_status'].unique().to_list()==['no_finite_score']
    winners=pl.read_csv(tmp_path/'out.scores.tsv',separator='\t')
    assert winners.filter(pl.col('query_id')==short)['epitope'].null_count()==1
    audit=pl.read_csv(tmp_path/'out.policy-audit.tsv',separator='\t')
    assert audit['query_target_cases'].sum()==3*4


def test_paired_cli_retains_configured_kernel_default_without_rule_fusion(tmp_path,monkeypatch):
    rows=[]
    for cid,alpha,beta,epitope in [('a','CAVRDSNYQLIW','CASSLGQAYEQYF','A'),
        ('b','CAVRDANYQLIW','CASSLGRAYEQYF','B')]:
        for gene,seq,v,j in [('TRA',alpha,'TRAV1-2','TRAJ33'),('TRB',beta,'TRBV19','TRBJ2-7')]:
            rows.append(dict(complex_id=cid,cdr3=seq,v=v,j=j,gene=gene,species='HomoSapiens',epitope=epitope))
    source=tmp_path/'reference.tsv';pl.DataFrame(rows).write_csv(source,separator='\t')
    model=_policy(tmp_path,source,pl.DataFrame({'species':['human'],'locus':['PAIRED'],
        'epitope':['A'],'reference_class':['v_only']}))
    monkeypatch.setattr(bank,'stream_reference_rules',lambda *a,**k:pytest.fail('paired policy invented rule fusion'))
    def kernel(q,r,*,consume,paired,**kwargs):
        assert paired and {'seq0','seq1','v0','v1'} <= set(r.columns)
        frame=q.select('query_id').join(r.select('epitope').unique(),how='cross').with_columns(
            pl.when(pl.col('query_id')==0).then(0.).when(pl.col('epitope')=='A').then(3.).otherwise(5.).alias('score'),
            (pl.col('query_id')!=0).alias('supported'))
        consume('kernel_v_bounded',frame)
        return {'native_batches':1}
    monkeypatch.setattr(bank,'stream_reference_kernels',kernel)
    prefix=tmp_path/'out';assert _run(source,prefix,model,'--locus','paired','--min-records','1')==0
    continuous=pl.read_parquet(tmp_path/'out.continuous/*.parquet')
    assert continuous.filter(pl.col('epitope')=='A')['policy_status'].unique().to_list()==['paired_default']
    assert continuous.filter(pl.col('epitope')=='B')['policy_status'].unique().to_list()==['no_reference_profile']
    assert continuous['selected_score_column'].unique().to_list()==['kernel_v_bounded']
    assert continuous.filter(pl.col('query_id')==0)['supported'].not_().all()
    scores=pl.read_csv(tmp_path/'out.scores.tsv',separator='\t')
    assert scores.filter(pl.col('query_id')==0)['epitope'].null_count()==1
    assert scores.filter(pl.col('query_id')==1)['epitope'].to_list()==['B']
