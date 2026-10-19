"""Primitive native blocks preserve identities, observations and complete queries."""
import polars as pl
from polars.testing import assert_frame_equal
import pytest
from seqtree import Index
from vdjmatch.api import Annotator, _prepare
from vdjmatch.match import VdjdbIndex, search_params


@pytest.mark.parametrize('exclude_exact',[False,True])
def test_native_consumer_equals_full_hits_for_duplicate_junction_query_ids(monkeypatch,exclude_exact):
    ref=pl.DataFrame({'gene':['TRB']*4,'cdr3':['CASSF','CASRF','CASRF','CASTF'],
        'v':['TRBV1','TRBV1','TRBV2','TRBV1'],'j':['TRBJ1']*4,
        'epitope':['A','A','B','B'],'reference_id':['R1','R1','R2','R3'],
        'arbitrary_metadata':['note1','note2','note3','note4']})
    index=VdjdbIndex.build(ref)
    _,queries=_prepare(pl.DataFrame({'query_id':['one','two','nohit'],
        'cdr3':['CASSF','CASSF','CYYYYYYYYF'],'v':['TRBV1','TRBV2',None]}))
    params=search_params('2',gap_open=1,gap_extend=1)
    expected=index.annotate(queries,params,gene='TRB',exclude_exact=exclude_exact)
    original=Index.search_batch_consume;calls=[]
    def consume(self,keys,params,callback,**kwargs):
        calls.append(len(keys))
        return original(self,keys,params,callback,**kwargs)
    monkeypatch.setattr(Index,'search_batch_consume',consume)
    monkeypatch.setattr(Index,'search_batch',lambda *a,**k:pytest.fail('full native hit lists retained'))
    parts=[]
    index.annotate(queries,params,gene='TRB',exclude_exact=exclude_exact,
        consumer=parts.append,chunk=1,threads=2)
    assert all(part['query_id'].n_unique() <= 1 for part in parts)
    actual=pl.concat(parts).sort('query_id','score','db_cdr3','db_v')
    assert_frame_equal(actual,expected.sort('query_id','score','db_cdr3','db_v'),check_exact=True)
    assert calls==[2]


def test_candidate_projection_keeps_source_multiplicity_without_wide_metadata(monkeypatch):
    ref=pl.DataFrame({'gene':['TRB']*3,'cdr3':['CASSF','CASRF','CASRF'],
        'v':['TRBV1']*3,'epitope':['A']*3,'reference_id':['R1','R2','R3'],
        'wide_metadata':['x'*10000]*3})
    ann=Annotator.from_frame(ref)
    from vdjmatch import api
    reduce=api.reduce_candidates;seen=[]
    def checked(hits,*args,**kwargs):
        assert 'wide_metadata' not in hits.columns
        seen.append(hits.height)
        return reduce(hits,*args,**kwargs)
    monkeypatch.setattr(api,'reduce_candidates',checked)
    result=ann.candidates(['CASSF'],scope='1')
    assert result['n_records'].to_list()==[3] and result['n_studies'].to_list()==[3]
    assert seen==[3]
