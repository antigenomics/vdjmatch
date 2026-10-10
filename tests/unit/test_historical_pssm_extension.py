"""Opt-in historical positional extension against small exhaustive scalar oracles."""
import math
from functools import wraps
import json
from pathlib import Path
import polars as pl
import pytest
from seqtree import Index,SubstitutionMatrix,gapblock
from vdjmatch.match.regions import significance_weights,gene_family
from vdjmatch.match.vgene import vsim


def _penalty(q,r):
    base=SubstitutionMatrix.blosum62();scale=base.scale();d=abs(len(q)-len(r))
    weights=[max(1,round(100*x)) for x in significance_weights(max(len(q),len(r)))]
    prior=gapblock.positions_prior((3,4,-4,-3));values=[]
    for start in range(min(len(q),len(r))+1):
        score=(2*scale*100+(d-1)*scale*100 if d else 0)+prior(start,d,max(len(q),len(r)))
        for j in range(min(len(q),len(r))):
            qi=j+(d if len(q)>len(r) and j>=start else 0)
            ri=j+(d if len(r)>len(q) and j>=start else 0)
            score+=weights[max(qi,ri)]*base.penalty(q[qi],r[ri])
        values.append(score)
    return min(values)


def test_original_edges_excluded_and_weighted_cdf_matches_brute_oracle(monkeypatch):
    from vdjmatch.match.unified import historical_pssm_extension
    q='CASSLGQAYEQYF';old=q[:-2]+'AF';r=q[:4]+'A'+q[4:];cross=q[:3]+'G'+q[3:];six='ATTTIA'+q[6:]
    queries=pl.DataFrame({'query_id':['q','duplicate','unknown'],'cdr3':[q,q,q],'v':['TRBV19','TRBV19*01','unknown']})
    refs=pl.DataFrame({'cdr3':[q,old,r,cross,r,six],'v':['TRBV19','TRBV19','TRBV19','TRBV5-1','TRBV5-1','TRBV19']})
    control=Index.build([q,r,r,old]);calls=[];original=gapblock.count_batch
    @wraps(original)
    def count(qs,rs,thresholds,**kwargs):
        calls.append((qs,rs,thresholds,kwargs));return original(qs,rs,thresholds,**kwargs)
    monkeypatch.setattr(gapblock,'count_batch',count)
    out=historical_pssm_extension(queries,refs,control,threads=1)
    n_reference=len(set(refs['cdr3']));controls=list(dict.fromkeys(control.ref_seqs()))
    expected=0
    assert sum(a!=b for a,b in zip(q,six))==6
    for ref,v in [(r,'TRBV19'),(cross,'TRBV5-1'),(six,'TRBV19')]:
        score=_penalty(q,ref);assert score<=7000
        nc=sum(c!=q and _penalty(q,c)<=score for c in controls)
        weight=1 if gene_family('TRBV19')==gene_family(v) else .25*vsim('TRBV19',v)
        expected+=weight*math.exp(-score/400)/max(n_reference/len(controls)*nc,.01)
    assert out['query_id'].to_list()==queries['query_id'].to_list()
    assert out['M_gap'].to_list()==[3]*3
    assert out['gapped_density'][0]==pytest.approx(expected)
    assert out['gapped_density'][1]==pytest.approx(expected)
    assert out['gapped_density'][2]==0
    assert out['gapped_edges_same_v'][0]==2 and out['gapped_edges_cross_v'][0]==1
    assert out['gapped_best_reference_junction'][0] in (r,cross,six)
    assert out['gapped_best_reference_v'][0] in ('TRBV19','TRBV5-1')
    assert out['gapped_best_query_v'].to_list()[:2]==['TRBV19','TRBV19*01']
    assert out['gapped_best_block_position'][0] is None
    assert out['gapped_best_alignment_status'][0]=='not_returned_by_native'
    assert len(calls)==1
    qs,cs,thresholds,options=calls[0]
    assert cs==controls and options['exclude_exact'] is True
    assert options['gap_open']==2800 and options['gap_extend']==1400
    assert options['position_weights_by_length'][len(q)]==[max(1,round(100*x)) for x in significance_weights(len(q))]
    parallel=historical_pssm_extension(queries,refs,control,threads=4)
    assert out.equals(parallel)


def test_no_length_switch_unchanged_original_and_invalid_id_retention():
    from vdjmatch.match.unified import historical_pssm_extension
    from vdjmatch.match.historical import density
    q='CASSLGQAYEQYF';old=q[:-2]+'AF';long='CASSLGQAYEQYFQAYEQYFF'
    refs=pl.DataFrame({'cdr3':[old,q[:4]+'A'+q[4:],long[:4]+'A'+long[4:]],'v':['TRBV19']*3,
        'epitope':['E']*3,'mhc_a':['HLA-A*02:01']*3,'mhc_b':['B2M']*3,'mhc_class':['MHCI']*3,'species':['HomoSapiens']*3,'gene':['TRB']*3})
    queries=pl.DataFrame({'query_id':['short','long'],'cdr3':[q,long],'v':['TRBV19']*2})
    control=Index.build([q,old,long]);before=density(queries,refs,control,threads=1,pool_reference=True)
    extension=historical_pssm_extension(queries,refs,control,threads=1)
    after=density(queries,refs,control,threads=1,pool_reference=True)
    assert before.equals(after)
    assert (extension['gapped_density']>0).all()
    invalid=historical_pssm_extension(pl.DataFrame({'query_id':['invalid','null'],'cdr3':['C*F',None],'v':['TRBV19',None]}),refs,control)
    assert invalid['query_id'].to_list()==['invalid','null']
    assert invalid['gapped_density'].to_list()==[0.,0.]
    assert invalid['availability'].to_list()==[False,False]
    empty=historical_pssm_extension(queries,refs.head(0),control)
    assert empty['gapped_status'].to_list()==['no_usable_reference']*2


def test_cli_geometry_keeps_original_density_and_records_explicit_contract(tmp_path):
    from vdjmatch.cli.__main__ import main
    q='CASSLGQAYEQYF';old=q[:-2]+'AF';r=q[:4]+'A'+q[4:]
    sample,reference,control,prefix=[tmp_path/n for n in ('q.tsv','r.tsv','c.tsv','out')]
    pl.DataFrame({'query_id':['q'],'junction_aa':[q],'locus':['TRB'],'v_call':['TRBV19'],'j_call':['TRBJ1-1']}).write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':[old,r],'gene':['TRB']*2,'species':['HomoSapiens']*2,'epitope':['E']*2,'mhc_a':['HLA-A*02:01']*2,'v':['TRBV19']*2}).write_csv(reference,separator='\t')
    pl.DataFrame({'junction_aa':[q,old,r]}).write_csv(control,separator='\t')
    args=['historical-density',str(sample),'--vdjdb',str(reference),'--locus','TRB','--epitope','E','--mhc-a','HLA-A*02:01','--pool-reference','--control',str(control),'--output-prefix',str(prefix)]
    assert main(args)==0
    original=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert main(args+['--gapped-extension','--gap-geometry','historical-pssm'])==0
    scores=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert scores['historical_density'].to_list()==original['score'].to_list()
    assert scores['p_enrichment'].to_list()==original['p_enrichment'].to_list()
    assert scores['score'][0]==pytest.approx(original['score'][0]+scores['gapped_density'][0])
    info=json.loads(Path(str(prefix)+'.manifest.json').read_text())['gapped_extension']
    assert info['geometry']=='historical-pssm'
    assert info['gap_open']==2800 and info['gap_extend']==1400 and info['cutoff']==7000
    assert info['kernel_scale']==400 and info['position_frame']=='longer full junction'
    assert len(info['native_sha256'])==64
    assert info['significance']=='p_enrichment remains the original component test; no combined P-value'


def test_missing_positional_kernel_fails_before_control_loading(monkeypatch):
    from types import SimpleNamespace
    from vdjmatch.match.unified import historical_pssm_extension
    from vdjmatch.cli.__main__ import main
    def old_count(queries,refs,thresholds,threads=0):
        raise AssertionError('old native code must not be called')
    monkeypatch.setattr(gapblock,'count_batch',old_count)
    def controls():
        raise AssertionError('controls must not be loaded')
    control=SimpleNamespace(ref_seqs=controls)
    with pytest.raises(RuntimeError,match='requires positional seqtree'):
        historical_pssm_extension(pl.DataFrame(),pl.DataFrame(),control)
    with pytest.raises(SystemExit) as error:
        main(['historical-density','absent.tsv','--vdjdb','absent-reference.tsv','--locus','TRB',
              '--epitope','E','--mhc-a','HLA-A*02:01','--pool-reference','--gapped-extension',
              '--gap-geometry','historical-pssm','--output-prefix','unused'])
    assert error.value.code==2


def test_stale_native_binding_is_rejected(monkeypatch):
    from seqtree import _core
    from vdjmatch.match.unified import require_historical_pssm_kernel
    monkeypatch.setattr(_core,'gapblock_matrix',lambda:None)
    with pytest.raises(RuntimeError,match='requires positional seqtree'):
        require_historical_pssm_kernel()


def test_original_raw_v_alias_weight_is_preserved():
    from vdjmatch.match.unified import historical_pssm_extension
    q='CASSLGQAYEQYF';r=q[:4]+'A'+q[4:]
    query=pl.DataFrame({'query_id':['alias'],'cdr3':[q],'v':['TRAV14']})
    ref=pl.DataFrame({'cdr3':[r],'v':['TRAV14/DV4']})
    assert vsim('TRAV14','TRAV14/DV4')==0
    out=historical_pssm_extension(query,ref,Index.build([q,r]))
    assert out['gapped_density'].to_list()==[0.]
