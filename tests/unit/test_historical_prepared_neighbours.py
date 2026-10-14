"""Ordered transient reference edges preserve the original scalar reduction exactly."""
import importlib.util
import math
from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal
from seqtree import Index
from vdjmatch.match import historical


def fixture():
    spec=importlib.util.spec_from_file_location('historical_fixture',Path(__file__).with_name('test_historical.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize('puncture',[False,True])
def test_prepared_edges_single_exact_and_half_reference_searches(monkeypatch,puncture):
    q=pl.DataFrame(dict(query_id=['a','b','c','orphan','otherV'],cdr3=['CASSLGQAYEQYF','CASSLAGQAYEQYF','CASSLGQAYEQYF','CAWWWWWWWWWF','CASSLGQAYEQYF'],
        v=['TRBV19','TRBV19','TRBV19','TRBV19','TRBV7-9'],j=['TRBJ1-1']*5))
    r=pl.DataFrame([dict(cdr3=s,v=v,j='TRBJ1-1',epitope=e,mhc_a=mhc,mhc_b='B2M',mhc_class='MHCI',species='HomoSapiens',gene='TRB')
        for e,mhc in [('E','HLA-A*02:01'),('F','HLA-A*03:01')] for s,v in [('CASSLGQAYEQYF','TRBV19'),('CASSLGQAYEQFF','TRBV7-9'),('CASSLAGQAYEQFF','TRBV19')]])
    control=Index.build(q['cdr3'].to_list()+['CASSLGQAYEQFF','CASSLAGQAYEQFF'])
    calls=[];native=historical._density_neighbours
    def observed(seqs,rs,threads):calls.append(threads);return native(seqs,rs,threads)
    monkeypatch.setattr(historical,'_density_neighbours',observed)
    original=historical.shared_control_histograms
    def old(*args,**kwargs):
        kwargs.pop('prepared_batches',None);original(*args,**kwargs)
        return None
    monkeypatch.setattr(historical,'shared_control_histograms',old)
    expected=historical.density(q,r,control,threads=2,puncture_populations=puncture)
    old_calls=len(calls);calls.clear()
    monkeypatch.setattr(historical,'shared_control_histograms',original)
    actual=historical.density(q,r,control,threads=2,puncture_populations=puncture)
    assert_frame_equal(actual,expected,check_exact=True)
    assert old_calls==4 and len(calls)==2 and calls==[2,2]
    if not puncture:
        # Independent owning scalar formula, using the ordered native edges.
        from vdjmatch.match.regions import gene_family
        from vdjmatch.match.vgene import precompute_v_similarities
        seqs=q['cdr3'].unique(maintain_order=True).to_list()
        for key,part in r.group_by(historical.PMHC,maintain_order=True):
            ref=part.unique('cdr3',keep='first',maintain_order=True);rs=ref['cdr3'].to_list()
            hits,_,radii=native(seqs,rs,2)
            cdf=historical.control_histograms(control,seqs,threads=2,radii=radii)
            similarities=precompute_v_similarities(q['v'],ref['v'])
            for query in q.iter_rows(named=True):
                total=0.0;v=gene_family(query['v']);seq=query['cdr3']
                for hit in hits[seq]:
                    if rs[hit.ref_id]==seq:continue
                    rv=gene_family(ref['v'][hit.ref_id]);weight=1.0 if v==rv else .25*similarities[v,rv]
                    total+=weight*math.exp(-hit.score/400.0)/max(len(rs)/len(control)*cdf[seq][hit.n_subs],.01)
                value=actual.filter((pl.col('query_id')==query['query_id']) &
                    pl.all_horizontal(*[pl.col(c)==v for c,v in zip(historical.PMHC,key)]))['score'][0]
                assert value==total


def test_prepared_axes_reject_wrong_population():
    q=pl.DataFrame(dict(query_id=['q'],cdr3=['CASSLGQAYEQYF'],v=['TRBV19'],j=['TRBJ1-1']))
    r=pl.DataFrame(dict(cdr3=['CASSLGQAYEQFF'],v=['TRBV19'],j=['TRBJ1-1'],epitope=['E'],mhc_a=['A'],mhc_b=['B2M'],mhc_class=['MHCI']))
    c=Index.build(['CASSLGQAYEQFF']);prepared=[]
    counts=historical.shared_control_histograms(c,[(q,r)],threads=2,prepared_batches=prepared)
    with pytest.raises(ValueError,match='prepared.*axes'):
        historical.density(q.with_columns(pl.lit('CASSLGQAYEQFF').alias('cdr3')),r,c,threads=2,pool_reference=True,
            control_counts=counts,prepared_neighbours=prepared[0])


@pytest.mark.parametrize('locus',['TRA','paired'])
@pytest.mark.parametrize('all_targets',[False,True])
def test_cli_prepared_vs_old_exact_and_native_call_count(tmp_path,monkeypatch,locus,all_targets):
    from vdjmatch.cli.__main__ import main
    module=fixture();raw,sample,prefix,args=module._paired_cli_inputs(tmp_path)
    copied=raw.with_columns((pl.lit('other-')+pl.col('query_id')).alias('query_id'),
        (pl.lit('other-')+pl.col('clone_id')).alias('clone_id'),pl.lit('F').alias('epitope'))
    pl.concat([raw,copied]).write_csv(sample,separator='\t')
    reference=Path(args[args.index('--vdjdb')+1]);ref=pl.read_csv(reference,separator='\t')
    pl.concat([ref,ref.with_columns(pl.lit('F').alias('epitope'),pl.lit('HLA-A*03:01').alias('mhc_a'))]).write_csv(reference,separator='\t')
    # Query restriction remains A02, so use shared observations at that restriction too.
    ref=pl.read_csv(reference,separator='\t');pl.concat([ref,ref.filter(pl.col('epitope')=='F').with_columns(pl.lit('HLA-A*02:01').alias('mhc_a'))]).write_csv(reference,separator='\t')
    options=list(args)+['--threads','2']
    if all_targets:
        options[options.index('--targets-from-sample')]='--all-targets'
        options.append('--likelihood-score')
    if locus=='TRA':options[options.index('--locus')+1]='TRA';options[options.index('--control')+1]=options[options.index('--alpha-control')+1]
    calls=[];search=historical._density_neighbours;shared=historical.shared_control_histograms
    def observed(seqs,rs,threads):calls.append(threads);return search(seqs,rs,threads)
    def old(*a,**k):
        k.pop('prepared_batches',None);shared(*a,**k)
        return None
    monkeypatch.setattr(historical,'_density_neighbours',observed)
    monkeypatch.setattr(historical,'shared_control_histograms',old)
    assert main(options)==0
    expected=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t',infer_schema_length=0);old_n=len(calls);calls.clear()
    monkeypatch.setattr(historical,'shared_control_histograms',shared)
    assert main(options)==0
    actual=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t',infer_schema_length=0)
    assert_frame_equal(actual,expected,check_exact=True)
    assert old_n==2*len(calls) and len(calls)==((6 if all_targets else 4) if locus=='paired' else (3 if all_targets else 2))
    assert set(calls)=={2}


def test_empty_queries_do_not_consume_next_prepared_target():
    from collections import deque
    from vdjmatch.cli.historical import _take_prepared
    q=pl.DataFrame({'cdr3':['CASSF']});r=q.clone();prepared=deque(['next-population'])
    assert _take_prepared(prepared,q.head(0),r,{'counts':1}) is None
    assert _take_prepared(prepared,q,r.head(0),{'counts':1}) is None
    assert len(prepared)==1
    assert _take_prepared(prepared,q,r,{'counts':1})=='next-population'
    with pytest.raises(ValueError,match='missing prepared'):
        _take_prepared(prepared,q,r,{'counts':1})


def test_empty_internal_multi_target_queries():
    q=pl.DataFrame(schema={'query_id':pl.String,'cdr3':pl.String,'v':pl.String,'j':pl.String})
    r=pl.DataFrame(dict(cdr3=['CASSF','CASSF'],v=['TRBV19']*2,j=['TRBJ1-1']*2,
        epitope=['E','F'],mhc_a=['A']*2,mhc_b=['B2M']*2,mhc_class=['MHCI']*2))
    out=historical.density(q,r,Index.build(['CASSF']),threads=2)
    assert out.height==0


def test_prepared_csr_zero_edges_has_fixed_width_query_metadata(monkeypatch):
    from array import array
    n=1000
    q=pl.DataFrame({'cdr3':[f'query-{i}' for i in range(n)]})
    r=pl.DataFrame({'cdr3':['reference']})
    def empty(seqs,rs,threads):
        return {s:[] for s in seqs},[[0,0] for s in seqs],dict.fromkeys(seqs,0)
    monkeypatch.setattr(historical,'_density_neighbours',empty)
    batches=[]
    historical.shared_control_histograms(Index.build(['CASSF']),[(q,r)],prepared_batches=batches)
    p=batches[0]
    assert len(p)==8 and all(isinstance(a,array) for a in p[2:])
    assert len(p[2])==n+1 and not any(p[2])
    assert all(len(a)==0 for a in p[3:6])
    assert len(p[6])==len(p[7])==n
    assert historical._prepared_size(p[0],p[1],0)>n*145


def test_prepared_budget_rejects_before_search_for_zero_edge_metadata(monkeypatch):
    q=pl.DataFrame({'cdr3':[f'query-{i}' for i in range(1000)]})
    r=pl.DataFrame({'cdr3':['reference']})
    monkeypatch.setattr(historical,'_PREPARED_LIMIT_BYTES',1024)
    def forbidden(*a):pytest.fail('oversized metadata must fail before native search')
    monkeypatch.setattr(historical,'_density_neighbours',forbidden)
    with pytest.raises(ValueError,match='2 GiB.*restrict'):
        historical.shared_control_histograms(Index.build(['CASSF']),[(q,r)],prepared_batches=[])


def test_prepared_budget_includes_temporary_hits_and_fails_before_packing(monkeypatch):
    from types import SimpleNamespace
    q=pl.DataFrame({'cdr3':['CASSF']});r=pl.DataFrame({'cdr3':['CASRF']})
    minimum=historical._prepared_size(q['cdr3'].to_list(),r['cdr3'].to_list(),0)
    monkeypatch.setattr(historical,'_PREPARED_LIMIT_BYTES',minimum+20)
    def found(seqs,rs,threads):
        return {seqs[0]:[SimpleNamespace(ref_id=0,n_subs=1,score=100)]},[[0,1]],{seqs[0]:1}
    monkeypatch.setattr(historical,'_density_neighbours',found)
    def forbidden(*a):pytest.fail('packing must not allocate above budget')
    monkeypatch.setattr(historical,'_pack_prepared',forbidden)
    batches=[]
    with pytest.raises(ValueError,match='2 GiB.*restrict'):
        historical.shared_control_histograms(Index.build(['CASSF']),[(q,r)],prepared_batches=batches)
    assert batches==[]


def test_prepared_budget_shared_between_chain_collectors(monkeypatch):
    q=pl.DataFrame({'cdr3':['CASSF']});r=pl.DataFrame({'cdr3':['CASRF']})
    def empty(seqs,rs,threads):return {s:[] for s in seqs},[[0,0]],dict.fromkeys(seqs,0)
    monkeypatch.setattr(historical,'_density_neighbours',empty)
    size=historical._prepared_size(q['cdr3'].to_list(),r['cdr3'].to_list(),0)
    monkeypatch.setattr(historical,'_PREPARED_LIMIT_BYTES',2*size-1)
    budget=historical.PreparedBudget();alpha=[];beta=[]
    control=Index.build(['CASSF'])
    historical.shared_control_histograms(control,[(q,r)],prepared_batches=alpha,prepared_budget=budget)
    assert budget.retained==size
    with pytest.raises(ValueError,match='2 GiB.*restrict'):
        historical.shared_control_histograms(control,[(q,r)],prepared_batches=beta,prepared_budget=budget)
    assert len(alpha)==1 and beta==[]
