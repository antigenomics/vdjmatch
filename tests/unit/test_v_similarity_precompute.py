"""Run-local legacy V priors retain directed ratios and historical scores."""
from difflib import SequenceMatcher

import polars as pl
from seqtree import Index

from vdjmatch.match import vgene
from vdjmatch.match.regions import gene_family


def test_directed_table_matches_scalar_semantics_and_missing_calls():
    queries=['A*01','B',None,'UNKNOWN']
    refs=['A*02','B','','OTHER','UNKNOWN*03']
    loops={'A':'tide','B':'diet'}
    table=vgene.precompute_v_similarities(queries,refs,loops=loops)
    assert table[('A','B')] != table[('B','A')]
    for q in queries:
        for r in refs:
            q,r=gene_family(q),gene_family(r)
            expected=1.0 if q==r else (SequenceMatcher(None,loops[q],loops[r]).ratio()
                if q in loops and r in loops else 0.0)
            assert table[(q,r)]==expected
    assert table[('','')]==1
    assert table[('UNKNOWN','UNKNOWN')]==1
    actual=vgene.precompute_v_similarities(['TRBV5-1','TRBV19','UNKNOWN',None],
                                          ['TRBV5-8','TRBV19*02','UNKNOWN',''])
    assert all(value==vgene.vsim(q,r) for (q,r),value in actual.items())


def test_density_and_extension_reuse_exact_scalar_table(monkeypatch):
    from vdjmatch.match import historical, unified
    q='CASSLGQAYEQYF'
    refs=pl.DataFrame({'cdr3':[q[:-2]+'AF',q[:6]+'A'+q[6:]],
        'v':['TRBV19','TRBV5-1'], 'epitope':['E','E'],
        'mhc_a':['HLA-A*02:01']*2,'mhc_b':['B2M']*2,'mhc_class':['MHCI']*2})
    queries=pl.DataFrame({'query_id':['q','duplicate','missing'],'cdr3':[q]*3,'v':['TRBV19','TRBV19*02',None]})
    control=Index.build([q,refs['cdr3'][0],refs['cdr3'][1]])
    expected_density=historical.density(queries,refs,control)
    expected_extension=unified.historical_pssm_extension(queries,refs,control)
    table={(gene_family(qv),gene_family(rv)):vgene.vsim(qv,rv)
           for qv in queries['v'] for rv in refs['v']}
    def forbidden(*args,**kwargs):
        raise AssertionError('supplied table must avoid V-ratio recomputation')
    monkeypatch.setattr(historical,'precompute_v_similarities',forbidden)
    monkeypatch.setattr(vgene,'precompute_v_similarities',forbidden)
    assert historical.density(queries,refs,control,v_similarities=table).equals(expected_density)
    assert unified.historical_pssm_extension(queries,refs,control,v_similarities=table).equals(expected_extension)


def test_historical_cli_precomputes_once_for_multiple_targets(tmp_path,monkeypatch):
    from vdjmatch.cli.__main__ import main
    q='CASSLGQAYEQYF'
    sample=tmp_path/'sample.tsv'
    reference=tmp_path/'reference.tsv'
    controls=tmp_path/'controls.tsv'
    pl.DataFrame({'junction_aa':[q],'v_call':['TRBV19*01'],'j_call':['TRBJ2-7*01'],'locus':['TRB']}).write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':[q[:-2]+'AF',q[:-3]+'AF'], 'v':['TRBV19','TRBV5-1'], 'j':['TRBJ2-7']*2,
        'gene':['TRB']*2,'species':['HomoSapiens']*2,'epitope':['E','F'],
        'mhc_a':['HLA-A*02:01']*2,'mhc_b':['B2M']*2,'mhc_class':['MHCI']*2}).write_csv(reference,separator='\t')
    pl.DataFrame({'junction_aa':[q,q[:-2]+'AF']}).write_csv(controls,separator='\t')
    calls=[]
    original=vgene.precompute_v_similarities
    def precompute(qgenes,rgenes=None,**kwargs):
        calls.append(1)
        return original(qgenes,rgenes,**kwargs)
    monkeypatch.setattr(vgene,'precompute_v_similarities',precompute)
    assert main(['historical-density',str(sample),'--vdjdb',str(reference),'--locus','TRB',
        '--epitope','E','--epitope','F','--mhc-a','HLA-A*02:01','--control',str(controls),
        '--output-prefix',str(tmp_path/'out')])==0
    assert len(calls)==1
    assert pl.read_csv(tmp_path/'out.scores.tsv',separator='\t').height==2
