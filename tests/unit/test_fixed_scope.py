"""Naive fixed-ball counts, exact policy and operation-budget checks."""
import polars as pl
import pytest
from vdjmatch.runner.fixed_scope import fixed_scope_scores

CONTEXT=dict(epitope='E',mhc_a='HLA-A*02',mhc_b='B2M',mhc_class='MHCI')

def test_naive_counts_pool_compatible_restrictions_and_preserve_nohits(monkeypatch):
    from vdjmatch.runner import fixed_scope
    queries=pl.DataFrame([dict(query_id=str(i),cdr3=s,v='TRBV1',j='TRBJ1-1',locus='TRB',**CONTEXT)
        for i,s in enumerate(['CASSF','CWWWF',None])])
    reference=pl.DataFrame([dict(cdr3=s,v=v,j='TRBJ1-1',gene='TRB',complex_id='0',species='HomoSapiens',
        **(CONTEXT|{'mhc_a':m})) for s,v,m in [('CASSF','TRBV1','HLA-A*02:01'),('CARSF','TRBV1','HLA-A*02:01'),
        ('CARSF','TRBV2','HLA-A*02:02'),('CASSFF','TRBV1','HLA-A*02:01')]])
    scores,report=fixed_scope_scores(queries,reference,exclude_exact=True)
    assert report['native_batches']==1
    hit=scores.filter(pl.col('query_id')=='0')
    assert hit['beta_neighbor_junctions'].item()==1
    assert hit['beta_neighbor_clonotypes'].item()==2
    assert hit['beta_nearest_score'].item()==-1
    assert scores.filter(pl.col('query_id').is_in(['1','2']))['beta_neighbor_junctions'].to_list()==[0,0]
    assert scores.filter(pl.col('query_id').is_in(['1','2']))['beta_nearest_edits'].null_count()==2
    ordinary,_=fixed_scope_scores(queries,reference,exclude_exact=False)
    assert ordinary.filter(pl.col('query_id')=='0')['beta_nearest_score'].item()==0
    insertion,_=fixed_scope_scores(queries,reference,scope='0,1,1,1',exclude_exact=True)
    assert insertion.filter(pl.col('query_id')=='0')['beta_neighbor_junctions'].item()==1


def test_paired_linked_and_unlinked_axes_have_distinct_count_units():
    q=pl.DataFrame([dict(query_id=l,pair_id='p',cdr3=s,v=v,j=j,locus=l,**CONTEXT)
        for l,s,v,j in [('TRA','CCAF','TRAV1','TRAJ1'),('TRB','CASSF','TRBV1','TRBJ1-1')]])
    r=pl.DataFrame([dict(cdr3=s,v=v,j=j,gene=l,complex_id=c,species='HomoSapiens',**CONTEXT)
        for l,s,v,j,c in [('TRA','CAAF','TRAV1','TRAJ1','1'),('TRB','CASSF','TRBV1','TRBJ1-1','1'),
        ('TRB','CARSF','TRBV1','TRBJ1-1','0')]])
    f,report=fixed_scope_scores(q,r,paired=True,exclude_exact=True)
    assert report['native_batches']==2 and f['locus'].item()=='paired'
    assert f['linked_neighbor_junction_pairs'].item()==1
    assert f['linked_nearest_bottleneck_edits'].item()==1
    assert f['alpha_neighbor_junctions'].item()==1
    assert f['beta_neighbor_junctions'].item()==1
    assert f['any_hit'].item()
    unlinked,_=fixed_scope_scores(q,r.with_columns(pl.lit('0').alias('complex_id')),paired=True,exclude_exact=True)
    assert unlinked['linked_neighbor_junction_pairs'].item()==0
    assert unlinked['alpha_neighbor_junctions'].item()==1 and unlinked['beta_neighbor_junctions'].item()==1
    assert unlinked['any_hit'].item()


def test_shipping_cli_restores_raw_source_and_reports_scope(tmp_path):
    from argparse import Namespace
    from vdjmatch.cli.fixed_scope import main
    q=tmp_path/'query.tsv';r=tmp_path/'reference.tsv';prefix=tmp_path/'naive'
    pl.DataFrame([dict(query_id='q',junction_aa='CASSF',v_call='TRBV1',j_call='TRBJ1-1',locus='TRB',**CONTEXT)]).write_csv(q,separator='\t')
    pl.DataFrame([dict(cdr3=s,v='TRBV1',j='TRBJ1-1',gene='TRB',species='HomoSapiens',complex_id='0',**CONTEXT)
        for s in ['CASSF','CARSF']]).write_csv(r,separator='\t')
    assert main(Namespace(sample=str(q),vdjdb=str(r),species='human',paired=False,targets_from_sample=True,
        scope='1,0,0,1',exclude_exact=True,reference_population='all-valid',threads=1,output_prefix=str(prefix)))==0
    frame=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert frame['naive_neighbor_junctions'].item()==1 and frame['naive_nearest_score'].item()==-1
    assert frame['beta_full_reference_records'].item()==2 and frame['beta_reference_records'].item()==2


def test_empty_reference_and_invalid_source_cases_stay_present(tmp_path):
    from argparse import Namespace
    from vdjmatch.cli.fixed_scope import main
    q=tmp_path/'q.tsv';r=tmp_path/'r.tsv';prefix=tmp_path/'empty'
    pl.DataFrame([dict(query_id=str(i),junction_aa=s,v_call='TRBV1',j_call='TRBJ1-1',locus='TRB',**CONTEXT)
        for i,s in enumerate([None,'C?F','CASSF'])]).write_csv(q,separator='\t')
    pl.DataFrame([dict(cdr3='CASSF',v='TRBV1',j='TRBJ1-1',gene='TRB',species='HomoSapiens',complex_id='0',**(CONTEXT|{'epitope':'other'}))]).write_csv(r,separator='\t')
    assert main(Namespace(sample=str(q),vdjdb=str(r),species='human',paired=False,scope='1,0,0,1',
        exclude_exact=True,reference_population='all-valid',threads=1,output_prefix=str(prefix)))==0
    f=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert f.height==3 and f['naive_neighbor_junctions'].to_list()==[0,0,0]
    assert f['naive_nearest_score'].null_count()==3 and not f['any_hit'].any()


@pytest.mark.parametrize('species', [None, ['human', 'mouse'], ['human', None]])
def test_trained_cli_reuses_source_boundary_and_writes_flat_provenance(tmp_path, monkeypatch, species):
    import sys
    import types
    import pytest
    from argparse import Namespace
    from vdjmatch.cli.fixed_scope import main

    q=tmp_path/'q.tsv';r=tmp_path/'r.tsv';prefix=tmp_path/'trained'
    pl.DataFrame([dict(query_id=str(i),junction_aa=s,v_call='TRBV1',j_call='TRBJ1-1',
        locus='TRB',positive=i==0,**CONTEXT) for i,s in enumerate(['CASSF',None])]).write_csv(q,separator='\t')
    if species is not None:
        pl.read_csv(q,separator='\t').with_columns(pl.Series('species',species)).write_csv(q,separator='\t')
    pl.DataFrame([dict(cdr3='CARSF',v='TRBV1',j='TRBJ1-1',gene='TRB',
        species='HomoSapiens',complex_id='0',**CONTEXT)]).write_csv(r,separator='\t')
    module=types.ModuleType('vdjmatch.runner.classical')
    def score(queries,reference,**options):
        assert queries.height==2 and reference.height==1
        if species is not None:
            assert queries['species'].to_list()==species
        assert options==dict(species='human',paired=False,scope='4,2,2,6',max_indels=2,
            exclude_exact=True,threads=2)
        return queries.with_columns(pl.Series('classical_score',[1.5,None])),{'populations':[{'locus':'TRB','native_batches':1}]}
    module.classical_scores=score
    monkeypatch.setitem(sys.modules,module.__name__,module)
    arguments=Namespace(sample=str(q),vdjdb=str(r),species='human',paired=False,
        score_model='trained',scope='4,2,2,6',max_indels=2,exclude_exact=True,
        reference_population='all-valid',threads=2,output_prefix=str(prefix))
    assert main(arguments)==0
    result=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert result.height==2 and result['positive'].to_list()==[True,False]
    assert result['classical_score'].to_list()==[1.5,None]
    provenance=pl.read_csv(str(prefix)+'.provenance.tsv',separator='\t')
    assert provenance['max_indels'].item()==2 and provenance['native_batches'].item()==1
    assert not list(tmp_path.glob('trained*.json'))
    pl.read_csv(q,separator='\t').with_columns(pl.lit(0).alias('classical_score')).write_csv(q,separator='\t')
    with pytest.raises(ValueError,match='reserved'):
        main(arguments)


def test_exact_hit_mask_keeps_single_and_linked_reference_metadata_fixed():
    q = pl.DataFrame([dict(query_id=l, pair_id='p', cdr3=s, v=v, j=j, locus=l, **CONTEXT)
        for l, s, v, j in [('TRA', 'CCAF', 'TRAV1', 'TRAJ1'), ('TRB', 'CASSF', 'TRBV1', 'TRBJ1-1')]])
    r = pl.DataFrame([dict(cdr3=s, v=v, j=j, gene=l, complex_id=c, species='HomoSapiens', **CONTEXT)
        for l, s, v, j, c in [('TRA', 'CCAF', 'TRAV1', 'TRAJ1', '1'),
                             ('TRB', 'CASSF', 'TRBV1', 'TRBJ1-1', '1'),
                             ('TRA', 'CAAF', 'TRAV1', 'TRAJ1', '2'),
                             ('TRB', 'CASSF', 'TRBV1', 'TRBJ1-1', '2')]])
    ordinary, _ = fixed_scope_scores(q, r, paired=True)
    excluded, _ = fixed_scope_scores(q, r, paired=True, exclude_exact=True)
    fields = [c for c in ordinary.columns if '_reference_' in c and '_excluded_' not in c]
    assert ordinary.select(fields).equals(excluded.select(fields))
    assert excluded['linked_reference_complexes'].item() == 2
    assert excluded['linked_reference_junction_pairs'].item() == 2
    assert excluded['linked_excluded_reference_complexes'].item() == 1
    assert excluded['linked_neighbor_junction_pairs'].item() == 1
    assert excluded['alpha_excluded_reference_records'].item() == 1
    assert excluded['beta_excluded_reference_records'].item() == 2
    assert excluded['alpha_neighbor_junctions'].item() == 1
    assert excluded['beta_neighbor_junctions'].item() == 0
