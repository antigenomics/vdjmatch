"""PSSM density reproduces the original formula without collapsing query identity."""
import math
import polars as pl
from seqtree import Index
from vdjmatch.api import Annotator, _prepare
from vdjmatch.match.historical import density
from vdjmatch.match.regions import significance_pssm
from seqtree import SearchParams


def test_historical_density_formula_and_duplicate_queries():
    q='CASSLGQAYEQYF'; r='CASSLGRAYEQYF'
    ann=Annotator.from_frame(pl.DataFrame({'gene':['TRB','TRB'],'cdr3':[r,r],
        'v':['TRBV19','TRBV19'],'epitope':['E','E']}))
    _,queries=_prepare(pl.DataFrame({'query_id':['positive','negative'],'cdr3':[q,q],
                                   'v':['TRBV19','TRBV19'],'locus':['TRB','TRB']}))
    p=SearchParams(max_subs=5,max_ins=0,max_dels=0,max_total_edits=5,engine='seqtm')
    p.pos_matrix=significance_pssm(len(q))
    score=Index.build([r]).search(q,p)[0].score
    got=density(queries,ann._index.records_for('TRB'),Index.build([q,r]),threads=1)
    expected=math.exp(-score/400)/.5
    assert got['query_id'].to_list()==['positive','negative']
    assert all(math.isclose(x,expected) for x in got['score'])
    assert got['n_reference'].to_list()==[1,1]
    assert all(math.isclose(x,1-math.exp(-.5),rel_tol=1e-14) for x in got['p_enrichment'])


def test_cli_rejects_missing_target_before_control_loading(tmp_path, capsys):
    import pytest
    from vdjmatch.cli.__main__ import main
    q=tmp_path/'q.tsv';r=tmp_path/'r.tsv'
    pl.DataFrame({'junction_aa':['CASSLGQAYEQYF'],'locus':['TRB']}).write_csv(q,separator='\t')
    pl.DataFrame({'cdr3':['CASSLGRAYEQYF'],'gene':['TRB'],'species':['HomoSapiens'],
                  'epitope':['E'],'mhc_a':['HLA-A*02:01']}).write_csv(r,separator='\t')
    with pytest.raises(SystemExit) as error:
        main(['historical-density',str(q),'--vdjdb',str(r),'--locus','TRB',
              '--epitope','E','--epitope','missing','--mhc-a','HLA-A*02:01',
              '--control',str(tmp_path/'nonexistent'),'--output-prefix',str(tmp_path/'out')])
    assert error.value.code==2
    assert 'absent under the selected reference restriction: missing' in capsys.readouterr().err
