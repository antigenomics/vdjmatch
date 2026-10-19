"""Candidate-only linked evidence equals the legacy detailed paired batch."""
import polars as pl
from polars.testing import assert_frame_equal
import pytest
from seqtree import Index
from vdjmatch import Annotator
from vdjmatch.match import search_params


def reference():
    pairs=[('self','CAVVF','CASSF','A'),('aexact','CAVVF','CASRF','A'),
        ('bexact','CAVAF','CASSF','A'),('near','CAVAF','CASRF','A'),
        ('repeat','CAVAF','CASRF','A'),('conflict','CAVAF','CASRF','B'),
        ('unlinked','CAVVF',None,'C')]
    rows=[dict(complex_id=cid,gene=g,cdr3=seq,epitope=ep,
        mhc_a='A*02:01' if ep!='B' else 'A*01:01',record_id=f'{cid}-{g}',reference_id=cid,
        source_note='x'*10000) for cid,a,b,ep in pairs for g,seq in [('TRA',a),('TRB',b)] if seq]
    return pl.DataFrame(rows+[dict(rows[0],record_id='additional-source',reference_id='other')])

@pytest.mark.parametrize('exclude_exact',[False,True])
@pytest.mark.parametrize('weighted',[False,True])
def test_paired_stream_all_fields_and_single_pool(monkeypatch,exclude_exact,weighted):
    ann=Annotator.from_frame(reference())
    queries=pl.DataFrame({'query_id':[f'q{i:03}' for i in range(73)],
        'cdr3_alpha_aa':['CAVVF','CAVAF','CYYYYYYYYF',None,'CAVVF']*14+['CAVVF']*3,
        'cdr3_beta_aa':['CASSF','CASRF','CYYYYYYYYF','CASSF','invalid']*14+['CASSF']*3})
    sp=search_params('1',**({} if weighted else {'matrix':'','gap_open':1,'gap_extend':1}))
    ctrl={'TRA':Index.build(['CAVVF','CAVAF','CAVTF','CAVTF','CYYYYF']),
          'TRB':Index.build(['CASSF','CASRF','CASTF','CASTF','CYYYYF'])}
    class MaterializedControl:
        def __init__(self,index):self.index=index
        def __len__(self):return len(self.index)
        def search_batch(self,*args):return self.index.search_batch(*args)
        def ref_seq(self,*args):return self.index.ref_seq(*args)
    oracle_controls={gene:MaterializedControl(index) for gene,index in ctrl.items()}
    _,expected=ann.paired_candidates(queries,scope=sp,control=oracle_controls,exclude_exact=exclude_exact,return_hits=True)
    original=Index.search_paired_batch_consume;calls=[]
    def consume(self,other,qa,qb,params,callback,**kwargs):
        calls.append((len(qa),len(qb)))
        return original(self,other,qa,qb,params,callback,**kwargs)
    monkeypatch.setattr(Index,'search_paired_batch_consume',consume)
    monkeypatch.setattr(Index,'search_batch',lambda *a,**k:pytest.fail('full native hit lists retained'))
    actual=ann.paired_candidates(queries,scope=sp,control=ctrl,exclude_exact=exclude_exact,threads=2)
    assert_frame_equal(actual,expected,check_exact=True)
    assert calls==[(59,59)]


def test_paired_empty_and_capability(monkeypatch):
    ann=Annotator.from_frame(reference())
    queries=pl.DataFrame({'cdr3_alpha_aa':['CYYYYYYYYF'],'cdr3_beta_aa':['CYYYYYYYYF']})
    _,expected=ann.paired_candidates(queries,return_hits=True)
    assert_frame_equal(ann.paired_candidates(queries),expected)
    monkeypatch.setattr(Index,'search_paired_batch_consume',None)
    with pytest.raises(RuntimeError,match='search_paired_batch_consume'):
        ann.paired_candidates(queries)
    assert_frame_equal(ann.paired_candidates(queries,return_hits=True)[1],expected)
