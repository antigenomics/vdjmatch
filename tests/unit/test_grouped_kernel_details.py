"""Fresh old/new detailed CLI equality on synthetic full pMHC/source identities."""
import importlib.util
from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal


def inputs(tmp_path,paired):
    path=Path(__file__).with_name('test_grouped_kernel_selection.py')
    spec=importlib.util.spec_from_file_location('detail_fixture',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    sample,reference=module.fixture(tmp_path,paired)
    raw=pl.read_csv(sample,separator='\t')
    first=raw.head(2 if paired else 1)
    orphan=first.with_columns((pl.lit('orphan-')+pl.col('sequence_id')).alias('sequence_id'),
        pl.lit('orphan').alias('clone_id'),
        pl.when(pl.col('locus')=='TRB').then(pl.lit('C'+'W'*90+'F')).otherwise(pl.col('junction_aa')).alias('junction_aa'))
    invalid=first.with_columns((pl.lit('invalid-')+pl.col('sequence_id')).alias('sequence_id'),
        pl.lit('invalid').alias('clone_id'),pl.lit('UNKNOWN').alias('v_call'))
    pl.concat([raw,orphan,invalid]).write_csv(sample,separator='\t')
    return sample,reference


@pytest.mark.parametrize('paired',[False,True])
@pytest.mark.parametrize('assayed',[False,True])
@pytest.mark.parametrize('exclude',[False,True])
def test_detailed_old_new_all_fields_and_native_batch_budget(tmp_path,monkeypatch,paired,assayed,exclude):
    from vdjmatch.cli import neighbours as cli
    from vdjmatch.match import geometry as td
    from vdjmatch.cli.__main__ import main
    from vdjmatch.evalue import kernel
    sample,reference=inputs(tmp_path,paired)
    args=['annotate','--details',str(sample),'--vdjdb',str(reference),'--locus','paired' if paired else 'TRB']
    if assayed:args+=['--targets-from-sample']
    if exclude:args+=['--exclude-exact']
    eligibility=cli._grouped_detailed_eligible
    native=td.grouped_total_distance_count_batch;calls=[]
    scorer=kernel._kernel_scores;math_rows=[];batches=[];enrich=kernel.grouped_kernel_enrichment
    def math(*args,**kwargs):
        math_rows.append(args[0].shape[0])
        return scorer(*args,**kwargs)
    def detailed(queries,references,groups,**kwargs):
        start=len(math_rows)
        result=enrich(queries,references,groups,**kwargs)
        expected=(sum(part.select(kwargs['query_keys']).unique().height for _,part,_ in groups)
            if kwargs['assayed'] else queries.select(kwargs['query_keys']).unique().height*len(groups))
        batches.append((sum(math_rows[start:]),expected))
        return result
    monkeypatch.setattr(kernel,'_kernel_scores',math)
    monkeypatch.setattr(kernel,'grouped_kernel_enrichment',detailed)
    def counted(*args,**kwargs):
        calls.append(kwargs['threads'])
        return native(*args,**kwargs)
    monkeypatch.setattr(td,'grouped_total_distance_count_batch',counted)
    for threads in (1,4):
        old,new=tmp_path/('old'+str(threads)),tmp_path/('new'+str(threads))
        monkeypatch.setattr(cli,'_grouped_detailed_eligible',lambda *a:False)
        assert main(args+['--threads',str(threads),'--output-prefix',str(old)])==0
        calls.clear();batches.clear()
        monkeypatch.setattr(cli,'_grouped_detailed_eligible',eligibility)
        assert main(args+['--threads',str(threads),'--output-prefix',str(new)])==0
        assert calls==[threads]*(3 if paired else 1)
        assert batches and all(actual==expected for actual,expected in batches)
        old_files=sorted(tmp_path.glob(old.name+'*.tsv'))
        assert old_files
        for path in old_files:
            other=path.with_name(path.name.replace(old.name,new.name,1))
            expected=pl.read_csv(path,separator='\t',infer_schema_length=0)
            actual=pl.read_csv(other,separator='\t',infer_schema_length=0)
            assert_frame_equal(actual,expected,check_exact=True)
            if threads==4:
                serial=other.with_name(other.name.replace('new4','new1',1))
                assert_frame_equal(actual,pl.read_csv(serial,separator='\t',infer_schema_length=0),check_exact=True)


def test_sparse_union_keeps_different_source_targets_and_skips_unrequested_chunk(monkeypatch):
    from vdjmatch.evalue import kernel
    from types import SimpleNamespace
    import numpy as np
    from vdjmatch.match import geometry as td
    meta={'mhc_a':'A','mhc_b':'B2M','mhc_class':'MHCI'}
    q=pl.DataFrame({'query_id':range(72),'cdr3':['CASSLGQAYEQYF']*3+[str(i) for i in range(69)],
        '_allele':['TRBV19*01']*72})
    r=pl.DataFrame({'cdr3':['CASSLGQAYEQYF']*2,'_allele':['TRBV19*01']*2,'epitope':['A','B']})
    groups=[({'epitope':'A',**meta},q.head(2),r.head(1)),
        ({'epitope':'B',**meta},q.slice(2,1),r.tail(1))]
    def counts(*args,**kwargs):
        for start,stop in ((0,64),(64,70)):
            n=stop-start
            yield SimpleNamespace(start=start,stop=stop,counts=np.ones((n,7,2),dtype=np.uint32),
                sizes=np.ones((n,7),dtype=np.uint32),exact_counts=np.ones((n,2),dtype=np.uint32))
    monkeypatch.setattr(td,'grouped_total_distance_count_batch',counts)
    native=kernel._kernel_scores;seen=[]
    def score(*args,**kwargs):
        seen.append(args[0].shape)
        return native(*args,**kwargs)
    monkeypatch.setattr(kernel,'_kernel_scores',score)
    out=kernel.grouped_kernel_enrichment(q,r,groups,query_keys=['cdr3','_allele'],
        search_keys=['cdr3','_allele'],paired=False,cutoffs=[0,10],radius=0,kernel_scale=90,
        scorer='v-bounded',assayed=True)
    assert seen==[(2,2)]
    assert out.select('query_id','epitope').rows()==[(0,'A'),(1,'A'),(2,'B')]


@pytest.mark.parametrize('assayed',[False,True])
def test_score_consumer_preserves_every_case_without_retaining_batches(monkeypatch,assayed):
    import weakref
    from vdjmatch.evalue import kernel
    from vdjmatch.match import geometry
    from vdjmatch.aggregate.candidates import PMHC
    q=pl.DataFrame({'query_id':range(700),'cdr3':['CASSLGQAYEQYF']*700,
        '_allele':['TRBV19*01']*700})
    r=pl.DataFrame({'cdr3':['CASSLGQAYEQYF']*128,'_allele':['TRBV19*01']*128,
        'epitope':['E'+str(i) for i in range(128)]})
    groups=[({'epitope':e,'mhc_a':None,'mhc_b':'B2M','mhc_class':'MHCI'},q,
        r.filter(pl.col('epitope')==e)) for e in r['epitope']]
    options=dict(query_keys=['cdr3','_allele'],search_keys=['cdr3','_allele'],paired=False,
        cutoffs=[0,4],radius=4,kernel_scale=22.5,scorer='v-bounded',assayed=assayed,
        matrix='blosum62',ctrim=3,stem_depth=4,_scores_only=True)
    expected=kernel.grouped_kernel_enrichment(q,r,groups,**options)
    native=geometry.grouped_total_distance_count_batch;calls=[]
    def counted(*args,**kwargs):
        calls.append(1)
        return native(*args,**kwargs)
    monkeypatch.setattr(geometry,'grouped_total_distance_count_batch',counted)
    retained=[];parts=[]
    def consume(frame):
        assert all(ref() is None for ref in retained)
        retained.append(weakref.ref(frame))
        assert frame.height<=65536
        assert frame.columns==['query_id',*PMHC,'kernel_score','_target','_query_order']
        assert frame['_query_order'].max()<q.height and frame['_target'].max()<len(groups)
        parts.append(frame.select('query_id',*PMHC,'kernel_score'))
    result=kernel.grouped_kernel_enrichment(q,r,groups,_consume=consume,**options)
    assert result is None and calls==[1]
    assert all(ref() is None for ref in retained)
    assert [part.height for part in parts]==[65536,24064]
    assert_frame_equal(pl.concat(parts).sort('query_id',*PMHC),
        expected.sort('query_id',*PMHC),check_exact=True)
