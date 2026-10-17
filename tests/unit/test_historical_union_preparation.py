"""Union preparation preserves target-specific scalar order and fixed populations."""
import polars as pl
import pytest
from polars.testing import assert_frame_equal
from seqtree import Index

from vdjmatch.aggregate.candidates import PMHC
from vdjmatch.match import historical


def _inputs():
    sequences=['CASSLGQAYEQYF','CASSLGQAYEQFF','CASSLAGQAYEQYF','CASSGQAYEQYF',
               'CWWWWWWWWWWWF','CASSLGQAYEQYF']
    query=pl.DataFrame({'query_id':range(len(sequences)),'cdr3':sequences,
        'v':['TRBV19']*len(sequences),'j':['TRBJ1-1']*len(sequences),'locus':['TRB']*len(sequences)})
    rows=[]
    for target,values in [('A',[('CASSLGQAYEQYF','TRBV19'),('CASSLGQAYEQFF','TRBV7-9'),
                               ('CASSLGQAYEQYF','TRBV7-9')]),
                          ('B',[('CASSLGQAYEQFF','TRBV19'),('CASSLGQAYEQYF','TRBV7-9'),
                                ('CASSLAGQAYEQFF',None)]),
                          ('C',[('CWWWWWWWWWWWF','TRBV19')])]:
        rows.extend(dict(cdr3=s,v=v,j='TRBJ1-1',gene='TRB',species='HomoSapiens',epitope=target,
            mhc_a='A',mhc_b='B2M',mhc_class='MHCI') for s,v in values)
    reference=pl.DataFrame(rows)
    control=Index.build(sequences+reference['cdr3'].to_list()+['CASSLGKAYEQYF'])
    return query,reference,control


def test_union_prep_exact_scalar_csr_and_one_reference_search(monkeypatch):
    query,reference,control=_inputs()
    batches=[(query,part) for part in reference.partition_by(PMHC,maintain_order=True)]
    expected=[]
    counts=historical.shared_control_histograms(control,batches,prepared_batches=expected)
    native=historical._density_neighbours
    calls=[]
    def observed(seqs,rs,threads,**kwargs):
        calls.append((list(seqs),list(rs),kwargs.get('unit_reference_ids') is not None))
        return native(seqs,rs,threads,**kwargs)
    monkeypatch.setattr(historical,'_density_neighbours',observed)
    query_reads=[]
    getitem=pl.DataFrame.__getitem__
    def observed_getitem(self,item):
        if self is query and isinstance(item,str) and item=='cdr3':
            query_reads.append(1)
        return getitem(self,item)
    monkeypatch.setattr(pl.DataFrame,'__getitem__',observed_getitem)
    actual=[]
    union_counts=historical.shared_control_histograms(control,batches,prepared_batches=actual,identical_queries=True)
    assert counts==union_counts and len(calls)==1 and calls[0][2]
    assert query_reads==[1]
    assert len(actual)==len(expected)==3
    for got,want,(_,part) in zip(actual,expected,batches):
        assert got==want
        scalar=historical.density(query,part,control,pool_reference=True,control_counts=counts)
        prepared=historical.density(query,part,control,pool_reference=True,control_counts=union_counts,
                                    prepared_neighbours=got)
        assert_frame_equal(prepared,scalar,check_exact=True)
    # Indel-only queries have radius-one significance support in target A.
    alpha=actual[0]
    for index in (2,3):
        position=list(alpha.query_indices).index(index)
        assert alpha.unit_counts[position]>0 and alpha.radii[position]>=1
    assert all(part.queries is actual[0].queries for part in actual)


@pytest.mark.parametrize('reorder',[True,False])
def test_union_validates_junction_axis_before_search(monkeypatch,reorder):
    query,reference,control=_inputs()
    def forbidden(*a,**k):pytest.fail('axis mismatch must fail before searching')
    monkeypatch.setattr(historical,'_density_neighbours',forbidden)
    changed=query.with_columns(pl.col('cdr3').reverse() if reorder else
        pl.col('cdr3').replace('CASSLGQAYEQFF','CASSLGQAYEQAF'))
    with pytest.raises(ValueError,match='identical ordered unique junction axes'):
        historical.shared_control_histograms(control,[(query,reference),(changed,reference)],
            prepared_batches=[],identical_queries=True)


def test_union_mutating_generator_revalidates_the_same_dataframe(monkeypatch):
    query,reference,control=_inputs()
    def batches():
        yield query,reference
        query.replace_column(query.get_column_index('cdr3'),query['cdr3'].reverse())
        yield query,reference
    def forbidden(*a,**k):pytest.fail('mutated axes must fail before searching')
    monkeypatch.setattr(historical,'_density_neighbours',forbidden)
    with pytest.raises(ValueError,match='identical ordered unique junction axes'):
        historical.shared_control_histograms(control,batches(),prepared_batches=[],identical_queries=True)


def test_union_empty_sparse_targets_keep_axis_and_zero_cdf():
    query,reference,control=_inputs()
    query=query.filter(pl.col('query_id')==4)
    reference=reference.filter(pl.col('epitope')!='C')
    batches=[(query,part) for part in reference.partition_by(PMHC,maintain_order=True)]
    prepared=[]
    counts=historical.shared_control_histograms(control,batches,prepared_batches=prepared,identical_queries=True)
    assert counts=={query['cdr3'].item():[0]*6}
    assert len(prepared)==2 and prepared[0].queries is prepared[1].queries
    assert all(not part.query_indices and not part.refs and list(part.offsets)==[0] for part in prepared)
    assert historical.shared_control_histograms(control,batches,identical_queries=True)==counts


def test_union_budget_rejects_before_search_and_before_join(monkeypatch):
    query,reference,control=_inputs()
    native=historical._density_neighbours
    gib=1024**3
    budget=historical.PreparedBudget()
    budget.retained=gib
    budget.check(gib,6*gib,peak_limit_bytes=8*gib)
    with pytest.raises(ValueError,match='2 GiB peak'):
        budget.check(gib,1)
    with pytest.raises(ValueError,match='2 GiB owned-state'):
        budget.check(gib+1,peak_limit_bytes=8*gib)
    with pytest.raises(ValueError,match='8 GiB peak.*2 GiB owned-state'):
        budget.check(gib,6*gib+1,peak_limit_bytes=8*gib)
    monkeypatch.setattr(historical,'_PREPARED_LIMIT_BYTES',1)
    def forbidden(*a,**k):pytest.fail('metadata budget must precede search')
    monkeypatch.setattr(historical,'_density_neighbours',forbidden)
    with pytest.raises(ValueError,match='2 GiB'):
        historical.shared_control_histograms(control,[(query,reference)],prepared_batches=[],identical_queries=True)
    calls=[]
    def observed(*args,**kwargs):
        calls.append(1)
        return native(*args,**kwargs)
    monkeypatch.setattr(historical,'_density_neighbours',observed)
    axis=tuple(query['cdr3'].unique(maintain_order=True))
    refs=reference['cdr3'].unique(maintain_order=True).to_list()
    base=historical._query_axis_size(axis)+20*historical._prepared_size(
        axis,refs,0,active_queries=0,include_query_axis=False)
    monkeypatch.setattr(historical,'_PREPARED_LIMIT_BYTES',base+100)
    prepared=[]
    with pytest.raises(ValueError,match='2 GiB'):
        historical.shared_control_histograms(control,[(query,reference)]*20,
            prepared_batches=prepared,identical_queries=True)
    assert prepared==[] and calls==[1]
    # The conversion guard must include both native edge lanes while they live.
    check=historical.PreparedBudget.check
    checks=[]
    def recorded(self,allocation,temporary=0,**kwargs):
        checks.append((self.retained,allocation,temporary,kwargs))
        return check(self,allocation,temporary,**kwargs)
    monkeypatch.setattr(historical.PreparedBudget,'check',recorded)
    monkeypatch.setattr(historical,'_PREPARED_LIMIT_BYTES',2*1024**3)
    historical.shared_control_histograms(control,[(query,reference)],prepared_batches=[],identical_queries=True)
    hit_rows,unit_rows,_=native(list(axis),refs,1)
    edge_count=sum(map(len,hit_rows.values()))+sum(row[1] for row in unit_rows)
    assert checks[3][2]-checks[2][2]==23*edge_count+2*((edge_count+7)//8)+4096
    assert checks[0][3]=={}  # Shared raw-axis ownership retains its default ceiling.
    assert all(row[3]=={'peak_limit_bytes':8*gib} for row in checks[1:])
    limit=sum(checks[3][:3])-1
    assert max(sum(row[:3]) for row in checks[:3])<limit
    checks.clear()
    monkeypatch.setattr(historical,'_UNION_PEAK_LIMIT_BYTES',limit)
    import numpy as np
    def forbidden_columns(*args,**kwargs):
        pytest.fail('typed edge columns must not allocate before their budget check')
    monkeypatch.setattr(np,'empty',forbidden_columns)
    with pytest.raises(ValueError,match='peak allocation budget'):
        historical.shared_control_histograms(control,[(query,reference)],prepared_batches=[],identical_queries=True)
    assert len(checks)==4


def test_streamed_all_six_scores_equal_old_per_target_preparation(monkeypatch):
    from vdjmatch.runner.retrieval_components import stream_components
    query,reference,control=_inputs()
    motif=pl.concat([reference.head(1)]*7).with_columns(pl.Series('cdr3',['CASSAAAA'+a+'F' for a in 'AILVDEN']))
    reference=pl.concat([reference,motif])
    query=query.with_columns(pl.when(pl.col('query_id')==4).then(pl.lit('CASSAAAAAF'))
        .otherwise(pl.col('cdr3')).alias('cdr3'))
    context={'density_index':control,'prior_background':reference.filter(pl.col('v').is_not_null()).select('cdr3','v','j')}
    owner=historical.shared_control_histograms
    def old(*args,**kwargs):
        kwargs.pop('identical_queries',None)
        return owner(*args,**kwargs)
    monkeypatch.setattr(historical,'shared_control_histograms',old)
    before={}
    stream_components(query,reference,context,species='human',locus='TRB',
        consume=lambda name,frame:before.setdefault(name,[]).append(frame))
    monkeypatch.setattr(historical,'shared_control_histograms',owner)
    after={}
    stream_components(query,reference,context,species='human',locus='TRB',
        consume=lambda name,frame:after.setdefault(name,[]).append(frame))
    assert set(after)==set(before)=={'accepted_original','current_original','current_original_log','gapblock','stem4','motif'}
    for name in before:
        # Motif's independent Polars floating sum can change reduction order
        # between calls; every field affected by union preparation stays exact.
        assert_frame_equal(pl.concat(after[name]).sort('query_id',*PMHC),
            pl.concat(before[name]).sort('query_id',*PMHC),check_exact=name!='motif',
            rel_tol=0,abs_tol=1e-16)


def test_union_dense_edges_use_exact_active_query_target_metadata(monkeypatch):
    from types import SimpleNamespace
    n=20
    query=pl.DataFrame({'cdr3':['CASS'+a+'YEQYF' for a in 'ACDEFGHIKLMNPQRSTVWY']})
    reference=query.rename({'cdr3':'cdr3'})
    axis=tuple(query['cdr3'])
    refs=reference['cdr3'].to_list()
    def found(seqs,rs,threads,*,unit_reference_ids):
        # Only one query is active in each target, despite many expanded edges.
        rows={s:[] for s in seqs}
        rows[seqs[0]]=[SimpleNamespace(ref_id=i,n_subs=1,score=i*7) for i in reversed(range(n))]
        unit_reference_ids.extend([[0] if i==0 else [] for i in range(n)])
        return rows,[(0,1 if i==0 else 0) for i in range(n)],dict.fromkeys(seqs,1)
    monkeypatch.setattr(historical,'_density_neighbours',found)
    monkeypatch.setattr(historical,'control_histograms',lambda *a,**k:{s:[0]*6 for s in axis})
    base=historical._query_axis_size(axis)+n*historical._prepared_size(
        axis,refs,0,active_queries=0,include_query_axis=False)
    expanded=n*(n+1)
    exact_active=n
    limit=base+expanded*13+exact_active*149
    assert limit < base+expanded*13+min(expanded,n*n)*149
    monkeypatch.setattr(historical,'_PREPARED_LIMIT_BYTES',limit)
    budget=historical.PreparedBudget()
    prepared=[]
    historical.shared_control_histograms(Index.build(['CASSF']),[(query,reference)]*n,
        prepared_batches=prepared,prepared_budget=budget,identical_queries=True)
    assert len(prepared)==n and budget.retained<=limit
    for part in prepared:
        assert list(part.query_indices)==[0]
        assert list(part.refs)==list(reversed(range(n)))
        assert list(part.penalties)==[i*7 for i in reversed(range(n))]
        assert list(part.unit_counts)==[1]
        assert list(part.offsets)==[0,n]


def test_union_target_windows_preserve_order_and_empty_targets(monkeypatch):
    query,reference,control=_inputs()
    query=query.filter(pl.col('query_id')!=4)
    populated=reference.filter(pl.col('epitope').is_in(['A','B']))
    empty=populated.head(1).with_columns(pl.lit('CFFFFF').alias('cdr3'))
    batches=[(query,empty if i in (0,7,19) else populated) for i in range(20)]
    calls=[]
    native=historical._density_neighbours
    def searched(*args,**kwargs):
        calls.append(1)
        return native(*args,**kwargs)
    monkeypatch.setattr(historical,'_density_neighbours',searched)
    checks=[]
    check=historical.PreparedBudget.check
    def recorded(self,allocation,temporary=0,**kwargs):
        checks.append(self.retained+allocation+temporary)
        return check(self,allocation,temporary,**kwargs)
    monkeypatch.setattr(historical.PreparedBudget,'check',recorded)
    expected=[]
    counts=historical.shared_control_histograms(control,batches,prepared_batches=expected,identical_queries=True)
    assert calls==[1]
    # The first six checks cover search, typed conversion and exact key metadata.
    # Force the wide expansion into several windows without changing any model.
    initial_peak=max(checks[:6])
    wide_peak=checks[6]
    assert wide_peak>initial_peak
    limit=(initial_peak+wide_peak)//2
    monkeypatch.setattr(historical,'_UNION_PEAK_LIMIT_BYTES',limit)
    joins=[]
    join=pl.DataFrame.join
    def observed(self,other,*args,**kwargs):
        if '_unit' in self.columns and '_local_ref' in other.columns:
            joins.append((other['_target'].min(),other['_target'].max()))
        return join(self,other,*args,**kwargs)
    monkeypatch.setattr(pl.DataFrame,'join',observed)
    actual=[]
    assert historical.shared_control_histograms(control,batches,prepared_batches=actual,identical_queries=True)==counts
    assert calls==[1,1] and len(joins)>1
    assert actual==expected
    assert joins[0][0]==0 and joins[-1][1]==19
    assert all(a[1]+1==b[0] for a,b in zip(joins,joins[1:]))
    assert all(not actual[i].refs and list(actual[i].offsets)==[0] for i in (0,7,19))
