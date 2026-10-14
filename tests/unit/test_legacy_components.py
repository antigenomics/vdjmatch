"""Synthetic oracle for the owning held-out components, without scientific reruns."""
import ast
import bisect
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from seqtree import Index, SearchParams
from vdjmatch.evalue import first_hit
from vdjmatch.match import vgene
from vdjmatch.match.historical import legacy_components
from vdjmatch.match.regions import gene_family, significance_pssm


def owner_function(path, name):
    # Extract functions only: importing analysis modules would access private data.
    source=Path(__file__).parents[2]/'bench'/path
    node=next(n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)
    env={'math':math,'np':np,'bisect':bisect,'Counter':Counter,'defaultdict':defaultdict,
         '_vg':vgene,'SOFTV_BETA':.25,'PSSM_SCALE':400.}
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),env)
    return env[name]


def inputs():
    alphabet='ACDEFGHIKLMNPQRSTVWY'
    seqs=['CASSL'+a+'GQAYEQYF' for a in alphabet]
    reference=pl.DataFrame({'cdr3':seqs+['CASSL'+a+'GQAYEQHF' for a in alphabet],
        'v':['TRBV19*01']*20+['TRBV7-9*01']*20,'j':['TRBJ1-1*01']*20+['TRBJ2-7*01']*20,
        'epitope':['A']*20+['B']*20})
    reference=pl.concat([reference,reference.head(1).with_columns(pl.lit('CONFLICT').alias('epitope'))])
    queries=pl.DataFrame({'query_id':['a','b','c','d'],'cdr3':[seqs[0],seqs[0], 'CASSLGQAYEQYF','CWWWWWWWWWWWWF'],
        'v':['TRBV19*02','TRBV7-9','UNKNOWN','TRBV19'],
        'j':['TRBJ1-1*02','TRBJ2-7','UNKNOWN','TRBJ1-1']})
    background=pl.DataFrame({'cdr3':seqs[:3]+['CASSLGQAYEQYF'],
        'v':['TRBV19','TRBV7-9','TRBV7-9','TRBV19'],
        'j':['TRBJ2-7','TRBJ2-7','TRBJ1-1','TRBJ1-1']})
    control=Index.build(background['cdr3'].to_list(),'aa')
    return queries,reference,control,background


def test_native_components_match_extracted_owning_scalar_formulas():
    q,r,ctrl,bg=inputs()
    actual=legacy_components(q,r,ctrl,bg)
    r=r.unique('cdr3',keep='first',maintain_order=True)
    assert actual['epitope'].unique().sort().to_list()==['A','B']
    rs=r['cdr3'].to_list(); seqs=q['cdr3'].unique(maintain_order=True).to_list()
    idx=Index.build(rs,'aa')
    hs={}
    for length in {len(s) for s in seqs}:
        batch=[s for s in seqs if len(s)==length]
        p=SearchParams(max_subs=5,max_ins=0,max_dels=0,max_total_edits=5,engine='seqtm')
        p.pos_matrix=significance_pssm(length)
        hs.update(zip(batch,idx.search_batch(batch,p,1)))
    th,cc=first_hit.scan(idx,r['epitope'].to_list(),ctrl,seqs,params=first_hit.scope(5,2,2),threads=1,exclude_exact=True)
    target,controls=dict(zip(seqs,th)),dict(zip(seqs,cc))
    ned=owner_function('holdout_features.py','ned_sim')
    odds=owner_function('hardcase.py','vj_logodds')
    lens=owner_function('holdout_features.py','_len_logodds')
    for query in q.iter_rows(named=True):
        seq=query['cdr3']
        rec={'v':query['v'],'ctrl':controls[seq], 'hits':[(h.score,h.n_subs,r['epitope'][h.ref_id],r['v'][h.ref_id])
            for h in hs[seq] if rs[h.ref_id]!=seq]}
        expected=ned(rec)
        for e in ['A','B']:
            row=actual.filter((pl.col('query_id')==query['query_id'])&(pl.col('epitope')==e)).row(0,named=True)
            ref=r.filter(pl.col('epitope')==e)
            assert row['ned_score']==expected.get(e,0.)
            for name,col in [('v_logodds','v'),('j_logodds','j')]:
                oracle=odds([gene_family(v) for v in ref[col]],[gene_family(v) for v in bg[col]])
                assert row[name]==oracle.get(gene_family(query[col]),0.)
            assert row['length_logodds']==lens([len(s) for s in ref['cdr3']],[len(s) for s in bg['cdr3']]).get(len(seq),0.)
            near=[h for h in target[seq] if 0<h[0]<=1]
            gate=first_hit.pvalue(near,controls[seq],r.height,len(ctrl))
            assert row['first_hit_p_enrichment']==gate['p_enrichment']
            assert row['first_hit_E']==gate['E']
            assert row['first_hit_gate']==(gate['p_enrichment']<.001)
            assert row['reference_union_size']==40


def test_missing_model_and_empty_queries_are_explicit():
    q,r,c,b=inputs()
    small=r.filter(pl.col('epitope')=='A').head(2)
    out=legacy_components(q,small,c,b)
    assert out['feature_model_available'].to_list()==[False]*4
    assert out['v_logodds'].null_count()==4
    empty=legacy_components(q.head(0),r,c,b)
    assert empty.height==0 and 'ned_score' in empty.columns


def test_cli_components_keep_invalid_identity_and_all_valid_policy(tmp_path):
    from vdjmatch.cli.__main__ import main
    import json
    q,r,_,bg=inputs()
    sample=tmp_path/'query.tsv'; reference=tmp_path/'reference.tsv';control=tmp_path/'control.tsv';prior=tmp_path/'prior.tsv'
    raw=q.rename({'cdr3':'junction_aa','v':'v_call','j':'j_call'}).with_columns(pl.lit('TRB').alias('locus'))
    raw=pl.concat([raw,raw.head(1).with_columns(pl.lit('bad').alias('query_id'),pl.lit('BAD?').alias('junction_aa'))])
    raw.write_csv(sample,separator='\t')
    r=r.with_columns(pl.lit('HomoSapiens').alias('species'),pl.lit('TRB').alias('gene'),
        pl.lit('HLA-A*02:01').alias('mhc_a'),pl.lit('B2M').alias('mhc_b'),pl.lit('MHCI').alias('mhc_class'))
    # Canonical rows excluded only by the shared current curation remain available
    # to the explicitly selected old all-valid raw reference policy.
    r=r.with_columns(pl.when(pl.col('epitope')=='A').then(pl.lit('KLGGALQAK')).otherwise(pl.col('epitope')).alias('epitope'))
    r.write_csv(reference,separator='\t');bg.write_csv(prior,separator='\t')
    bg.select(pl.col('cdr3').alias('junction_aa')).write_csv(control,separator='\t')
    prefix=tmp_path/'out'
    args=['historical-density',str(sample),'--vdjdb',str(reference),'--locus','TRB','--legacy-components',
          '--reference-population','all-valid','--germline-background',str(prior),'--control',str(control),'--output-prefix',str(prefix)]
    assert main(args)==0
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert out.height==10
    assert out.filter(pl.col('query_id')=='bad')['status'].to_list()==['invalid_query']*2
    assert out['epitope'].unique().sort().to_list()==['B','KLGGALQAK']
    manifest=json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert manifest['reference_population']=='all-valid'
    assert len(manifest['models'])==2
    assert manifest['estimator']=='legacy-heldout-components-v1'
    # Default remains the current curated benchmark reference policy.
    args[args.index('all-valid')]='curated-benchmark'
    assert main(args)==0
    assert pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')['epitope'].unique().sort().to_list()==['B','CONFLICT']

    # No valid junctions still produce all raw query dispositions on the union axis.
    raw.tail(1).write_csv(sample,separator='\t')
    assert main(args)==0
    empty=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert empty.height==2 and empty['status'].to_list()==['invalid_query']*2


def test_first_hit_rule_of_three_and_legacy_rejects_alternate_parameters(tmp_path):
    from vdjmatch.cli.__main__ import main
    q,r,_,bg=inputs()
    control=Index.build(['CWWWWWWWWWWWWF'],'aa')
    out=legacy_components(q.head(1),r,control,bg)
    assert out['first_hit_n_control'].to_list()==[0,0]
    assert out['first_hit_E'].to_list()==[120.,120.]
    args=['historical-density',str(tmp_path/'missing'),'--vdjdb',str(tmp_path/'missingref'),
          '--locus','TRB','--legacy-components','--germline-background',str(tmp_path/'bg'),
          '--puncture-populations','--output-prefix',str(tmp_path/'out')]
    with pytest.raises(SystemExit) as error:main(args)
    assert error.value.code==2
