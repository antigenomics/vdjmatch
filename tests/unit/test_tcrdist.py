"""Independent tcrdist3==0.3 default-model oracle, generated2026-10-10."""
import numpy as np
import pytest
from vdjmatch.match.tcrdist import distance_matrix


def test_native_matches_tcrdist3_germline_and_cdr3_oracle():
    junctions=['CTSKGMKF','CPNKNINGEEPQLCF','CGNKELKNMHITGNARHAGFF','CSTNRTCCVAFAATHEDQYMWMNQTSYRWF']
    genes=['TRBV19*01','TRBV7-9*01','TRBV5-1*01','TRBV20-1*01']
    expected=[[0,153,239,360],[153,0,214,357],[239,214,0,331],[360,357,331,0]]
    a=distance_matrix(junctions,junctions,genes,genes,threads=1)
    b=distance_matrix(junctions,junctions,genes,genes,threads=2)
    assert np.array_equal(a,expected) and np.array_equal(a,b)
    with pytest.raises(ValueError,match='unknown human V'):
        distance_matrix(junctions[:1],junctions[:1],['unknown'],genes[:1])


def test_tcrdist_cli_preserves_unavailable_queries(tmp_path):
    import polars as pl
    from vdjmatch.cli.__main__ import main
    q=tmp_path/'query.tsv';r=tmp_path/'reference.tsv';prefix=tmp_path/'out'
    pl.DataFrame({'sequence_id':['a','b','c','d','e'],
                  'junction_aa':['CASSLGQAYEQYF',None,'','INVALID!', 'CASSLGQAYEQYF'],
                  'v_call':['TRBV19*01']*4+[None],'locus':['TRB']*5}).write_csv(q,separator='\t')
    pl.DataFrame({'gene':['TRB'],'cdr3':['CASSLGRAYEQYF'],'v':['TRBV19*01'],
                  'species':['HomoSapiens'],'epitope':['E']}).write_csv(r,separator='\t')
    assert main(['tcrdist-neighbours',str(q),'--vdjdb',str(r),'--locus','TRB',
                 '--output-prefix',str(prefix),'--exclude-exact'])==0
    scores=pl.read_csv(str(prefix)+'.candidates.tsv',separator='\t')
    status=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t')
    assert scores['distance'].to_list()==[9]
    assert status['query_id'].to_list()==list(range(5))
    assert status['sequence_id'].to_list()==list('abcde')
    assert status['available'].to_list()==[True,False,False,False,False]
