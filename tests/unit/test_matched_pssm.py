"""One positional predicate for target and control, with original evidence preserved."""
import math
import json
from pathlib import Path
import polars as pl
import pytest
from seqtree import Index
from vdjmatch.match.unified import historical_pssm_extension
from test_historical_pssm_extension import _penalty


def test_matched_kernel_includes_original_edges_with_same_control_cdf():
    q='CASSLGQAYEQYF';old=q[:-2]+'AF';gapped=q[:6]+'A'+q[6:]
    queries=pl.DataFrame({'query_id':['q','dup'],'cdr3':[q,q],'v':['TRBV19']*2})
    refs=pl.DataFrame({'cdr3':[q,old,gapped],'v':['TRBV19']*3})
    control=Index.build([q,q,old,gapped])
    result=historical_pssm_extension(queries,refs,control,threads=1,matched_background=True)
    expected=0
    for r in [old,gapped]:
        penalty=_penalty(q,r)
        nc=sum(c!=q and _penalty(q,c)<=penalty for c in [q,old,gapped])
        expected+=math.exp(-penalty/400)/max(nc,.01)
    assert result['gapped_density'].to_list()==pytest.approx([expected,expected])
    assert result['gapped_edges_same_v'].to_list()==[2,2]
    original=historical_pssm_extension(queries,refs,control,threads=1)
    assert original['gapped_edges_same_v'].to_list()==[1,1]
    for threads in [4]:
        other=historical_pssm_extension(queries,refs,control,threads=threads,matched_background=True)
        assert result.equals(other)


def test_matched_cli_preserves_original_and_records_nonadditive_formula(tmp_path):
    from vdjmatch.cli.__main__ import main
    q='CASSLGQAYEQYF';old=q[:-2]+'AF';gap=q[:6]+'A'+q[6:]
    sample,reference,control,prefix=[tmp_path/n for n in ['q.tsv','r.tsv','c.tsv','out']]
    pl.DataFrame({'query_id':['q'],'junction_aa':[q],'locus':['TRB'],'v_call':['TRBV19'],'j_call':['TRBJ1-1']}).write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':[old,gap],'gene':['TRB']*2,'species':['HomoSapiens']*2,'epitope':['E']*2,'mhc_a':['HLA-A*02:01']*2,'v':['TRBV19']*2}).write_csv(reference,separator='\t')
    pl.DataFrame({'cdr3':[q,old,gap]}).write_csv(control,separator='\t')
    args=['historical-density',str(sample),'--locus','TRB','--vdjdb',str(reference),'--epitope','E','--mhc-a','HLA-A*02:01','--pool-reference','--control',str(control),'--gapped-extension','--gap-geometry','matched-pssm','--output-prefix',str(prefix)]
    main(args)
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert out['historical_density'][0]>0
    assert out['score'][0]==out['gapped_density'][0]
    assert out['estimator'][0]=='matched-positional-kernel-apex6-v1'
    manifest=json.loads(Path(str(prefix)+'.manifest.json').read_text())['gapped_extension']
    assert manifest['score_composition']=='matched_kernel_only'
    assert manifest['significance']=='p_enrichment remains the original component test; no combined P-value'
