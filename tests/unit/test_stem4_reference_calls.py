"""Synthetic final-score calibration and independent physical fold oracles."""
import argparse
import json
import inspect

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from vdjmatch.cli import reference_calls
from vdjmatch.evalue.kernel import grouped_kernel_selection, finish_kernel_selection, ladder_cutoffs


def args(tmp_path, *options):
    parser=argparse.ArgumentParser();reference_calls.register(parser.add_subparsers())
    return parser.parse_args(['reference-calls','--vdjdb',str(tmp_path/'ref.tsv'),'--locus','TRB',
        '--output-prefix',str(tmp_path/'out'),*options])


def test_defaults_and_component_selector_rejection(tmp_path):
    a=args(tmp_path)
    assert a.scorer=='original' and a.reference_population=='curated-benchmark' and a.exact_policy=='exclude'
    with pytest.raises(ValueError,match='final kernel_score'):
        reference_calls.main(args(tmp_path,'--scorer','stem4','--calibration-score','E'))

def test_reference_call_radius_uses_context_units(tmp_path,monkeypatch):
    from vdjmatch.match import geometry
    from vdjmatch.evalue import kernel
    reference(tmp_path)
    original=geometry.precompute_v_loop_distances
    def context(**kwargs):
        result=original(**kwargs).copy()
        result['distance_units']=(52,2)
        return result
    class Captured(Exception):
        pass
    def capture(*a,**kwargs):
        assert kwargs['radius']==kwargs['kernel_scale']==22.5*52
        assert kwargs['cutoffs']==ladder_cutoffs(52)
        raise Captured
    monkeypatch.setattr(geometry,'precompute_v_loop_distances',context)
    monkeypatch.setattr(kernel,'grouped_kernel_selection',capture)
    with pytest.raises(Captured):
        reference_calls.main(args(tmp_path,'--scorer','stem4'))


def reference(tmp_path):
    rows=[]
    for seq,v,j,e,n in [('CASSLAAAYEQYF','TRBV19','TRBJ1-1','A',30),
        ('CASSLAAAYEQYF','TRBV19','TRBJ2-1','B',30),
        ('CASSLACAYEQYF','TRBV19','TRBJ1-1','B',1),
        ('CASSLADAYEQYF','TRBV19','TRBJ1-1','SPARSE',1),
        ('CASSLAEAYEQYF','UNKNOWN','TRBJ1-1','A',1)]:
        rows += [dict(cdr3=seq,v=v,j=j,epitope=e,gene='TRB',species='HomoSapiens',
            mhc_class='MHCI',mhc_a='HLA-A*02:01',mhc_b='B2M')]*n
    pl.DataFrame(rows).write_csv(tmp_path/'ref.tsv',separator='\t')


def test_stem4_eligible_source_records_and_identity_restoration(tmp_path):
    reference(tmp_path)
    assert reference_calls.main(args(tmp_path,'--scorer','stem4'))==0
    calls=pl.read_csv(tmp_path/'out.calls.tsv',separator='\t')
    assert calls.height==4
    shared=calls.filter(pl.col('cdr3')=='CASSLAAAYEQYF')
    assert shared['j'].n_unique()==2 and shared['kernel_score'].n_unique()==1
    assert calls.filter(~pl.col('available'))['status'].to_list()==['invalid_query_or_gene']
    m=json.loads((tmp_path/'out.manifest.json').read_text())
    assert m['eligible_benchmark_records']==62 and m['ineligible_benchmark_records']==1
    assert m['all_reference_peptides']==3 and m['searchable_reference_keys']==3
    assert m['geometry']['matrix_scale']==26
    assert m['geometry']['radius']==m['geometry']['kernel_scale']==585
    assert 'fixed reference/control model populations' in m['self_exclusion']
    assert 'vdjmatch_E' not in calls.columns


def require_folds():
    import seqtree.gapblock as gapblock
    if 'query_withheld_peptide_ids' not in inspect.signature(gapblock.grouped_count_batch).parameters:
        pytest.skip('new native fold capability is tested with the isolated candidate wheel')


def test_grouped_folds_equal_physically_filtered_models():
    require_folds()
    r=pl.DataFrame({'cdr3':['CASSLAAAYEQYF','CASSLAAAYEQYF','CASSLACAYEQYF',
        'CASSLADAYEQYF','CASSLAEAYEQYF'], '_allele':['TRBV19*01']*5,
        'epitope':['W','T','T','SPARSE','W']})
    q=pl.DataFrame({'query_id':[0,1,2], 'cdr3':['CASSLAFAYEQYF','CASSLAFAYEQYF','CASSLAAAYEQYF'],
        '_allele':['TRBV19*01']*3, '_fold':['W','SPARSE','T']})
    options=dict(query_keys=['cdr3','_allele','_fold'],search_keys=['cdr3','_allele'],
        paired=False,cutoffs=ladder_cutoffs(26),radius=585,kernel_scale=585,scorer='v-bounded',
        species='human',matrix='blosum62',ctrim=3,stem_depth=4,exclude_exact=True)
    def groups(ref):
        # SPARSE belongs to the whole model but is deliberately not an output label.
        return [({'epitope':e},q,ref.filter(pl.col('epitope')==e)) for e in ['W','T']]
    state,_=grouped_kernel_selection(q,r,groups(r),withheld_column='_fold',**options)
    actual=finish_kernel_selection(q.select('query_id'),state)
    expected=[]
    for fold in ['W','SPARSE','T']:
        part=q.filter(pl.col('_fold')==fold);ref=r.filter(pl.col('epitope')!=fold)
        state,_=grouped_kernel_selection(part,ref,groups(ref),**options)
        expected.append(finish_kernel_selection(part.select('query_id'),state))
    assert_frame_equal(actual,pl.concat(expected),check_exact=True)


def test_stem4_withholding_never_selects_withheld_peptide(tmp_path):
    require_folds();reference(tmp_path)
    reference_calls.main(args(tmp_path,'--scorer','stem4','--leave-epitope-out'))
    calls=pl.read_csv(tmp_path/'out.calls.tsv',separator='\t')
    assert calls.height==4
    assert not calls.filter(pl.col('selected_epitope')==pl.col('epitope')).height
    cases=pl.read_csv(tmp_path/'out.calibration-cases.tsv',separator='\t')
    assert cases['outcome'].unique().to_list()==['null']
    assert cases.height==calls.height


def test_original_default_reference_policy_is_preserved(tmp_path):
    reference(tmp_path)
    ref=pl.read_csv(tmp_path/'ref.tsv',separator='\t')
    extra=ref.head(1).with_columns(pl.lit('KLGGALQAK').alias('epitope'))
    pl.concat([ref,extra]).write_csv(tmp_path/'ref.tsv',separator='\t')
    control=tmp_path/'control.tsv'
    pl.DataFrame({'junction_aa':['CASSLAFAYEQYF','CASSLAGAYEQYF']}).write_csv(control,separator='\t')
    for policy in [None,'curated-benchmark','all-valid']:
        prefix=tmp_path/(policy or 'default')
        opts=['--min-target-records','1','--control',str(control),'--output-prefix',str(prefix)]
        if policy:opts+=['--reference-population',policy]
        reference_calls.main(args(tmp_path,*opts))
    for suffix in ['calls','candidates','calibration-cases','truth','target-support']:
        assert_frame_equal(pl.read_csv(str(tmp_path/'default')+'.'+suffix+'.tsv',separator='\t'),
            pl.read_csv(str(tmp_path/'curated-benchmark')+'.'+suffix+'.tsv',separator='\t'),check_exact=True,check_row_order=False)
    default=json.loads((tmp_path/'default.manifest.json').read_text())
    all_valid=json.loads((tmp_path/'all-valid.manifest.json').read_text())
    assert default['valid_reference_records']==63
    assert all_valid['valid_reference_records']==64


def test_empty_withheld_model_retains_no_reference(tmp_path):
    require_folds();reference(tmp_path)
    ref=pl.read_csv(tmp_path/'ref.tsv',separator='\t').filter(pl.col('epitope')=='A')
    ref.write_csv(tmp_path/'ref.tsv',separator='\t')
    reference_calls.main(args(tmp_path,'--scorer','stem4','--leave-epitope-out'))
    calls=pl.read_csv(tmp_path/'out.calls.tsv',separator='\t')
    assert calls.filter(pl.col('available'))['status'].unique().to_list()==['no_reference']
    assert calls['kernel_score'].null_count()==calls.height


@pytest.mark.parametrize('option,value,match',[('--control','unused.tsv','--control'),
    ('--scope','1,1,1,1','--scope')])
def test_stem4_rejects_unused_original_options_before_io(tmp_path,option,value,match):
    with pytest.raises(ValueError,match=match):
        reference_calls.main(args(tmp_path,'--scorer','stem4',option,value))


@pytest.mark.parametrize('alias,canonical',[('human','HomoSapiens'),('mouse','MusMusculus')])
def test_stem4_species_alias_normalized_before_reference_load(tmp_path,monkeypatch,alias,canonical):
    from vdjmatch import db
    class Loaded(Exception):pass
    def load(path,*,species,gene):
        assert species==canonical and gene=='TRB'
        raise Loaded
    monkeypatch.setattr(db,'load',load)
    with pytest.raises(Loaded):reference_calls.main(args(tmp_path,'--scorer','stem4','--species',alias))


@pytest.mark.parametrize('paired',[False,True])
def test_new_python_with_stale_native_rejected_before_geometry(monkeypatch,paired):
    from vdjmatch.match import geometry
    def new_wrapper(*args,query_withheld_peptide_ids=None,**kwargs):
        raise AssertionError('wrapper must not execute')
    def stale_native():pass
    def blocked_geometry(**kwargs):raise AssertionError('geometry must not execute')
    method='grouped_paired_sum_count_batch' if paired else 'grouped_count_batch'
    monkeypatch.setattr(geometry.gapblock,method,new_wrapper)
    monkeypatch.setattr(geometry.gapblock,'_gapblock_grouped_count_batch',stale_native)
    monkeypatch.setattr(geometry,'precompute_v_loop_distances',blocked_geometry)
    with pytest.raises(RuntimeError,match='updated native seqtree'):
        geometry.grouped_total_distance_count_batch(None,None,[],paired=paired,query_withheld_peptide_ids=[])


@pytest.mark.parametrize('policy,excluded',[('exclude',True),('retain',False)])
def test_stem4_explicit_exact_policy_reaches_native_counts(tmp_path,monkeypatch,policy,excluded):
    reference(tmp_path)
    from vdjmatch.evalue import kernel
    native=kernel.grouped_kernel_selection
    flags=[]
    def counted(*args,**kwargs):
        flags.append(kwargs['exclude_exact'])
        return native(*args,**kwargs)
    monkeypatch.setattr(kernel,'grouped_kernel_selection',counted)
    reference_calls.main(args(tmp_path,'--scorer','stem4','--exact-policy',policy))
    assert flags==[excluded]
    manifest=json.loads((tmp_path/'out.manifest.json').read_text())
    assert manifest['geometry']['exclude_exact'] is excluded
    assert manifest['options']['exact_policy']==policy
    assert ('retained' in manifest['self_exclusion']) is (not excluded)
    assert pl.read_csv(tmp_path/'out.calls.tsv',separator='\t').height==4


def test_original_explicit_retained_policy_reaches_evidence(tmp_path,monkeypatch):
    from vdjmatch.api import Annotator
    reference(tmp_path)
    control=tmp_path/'control.tsv'
    pl.DataFrame({'junction_aa':['CASSLAFAYEQYF']}).write_csv(control,separator='\t')
    evidence=Annotator._evidence
    flags=[]
    def counted(self,*args,**kwargs):
        flags.append(kwargs['exclude_exact'])
        return evidence(self,*args,**kwargs)
    monkeypatch.setattr(Annotator,'_evidence',counted)
    reference_calls.main(args(tmp_path,'--control',str(control),'--exact-policy','retain'))
    assert flags==[False]
    manifest=json.loads((tmp_path/'out.manifest.json').read_text())
    assert manifest['options']['exact_policy']=='retain'
    assert 'retained' in manifest['self_exclusion']
