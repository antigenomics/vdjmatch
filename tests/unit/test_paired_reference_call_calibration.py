"""Genuine paired nulls and unchanged independent-axis selection oracles."""
import argparse
import copy
import json

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from vdjmatch.cli import calibration,reference_calls,neighbours as geometry
from vdjmatch.evalue.kernel import grouped_kernel_selection,finish_kernel_selection,ladder_cutoffs
from vdjmatch.evalue.paired import build_paired_ref
from vdjmatch.match.geometry import precompute_v_loop_distances,impute_reference_v_alleles


def reference(path):
    meta=dict(species='HomoSapiens',mhc_class='MHCI',mhc_a='HLA-A*02:01',mhc_b='B2M')
    rows=[]
    for pair,e,a,b,av in [('p1','A','CAVRDSNYQLIW','CASSLAAAYEQYF','TRAV12-2'),
        ('p2','B','CAVRDSNYQLVW','CASSLACAYEQYF','TRAV12-2'),
        ('p3','A','CAVRDSNYQLIW','CASSLADAYEQYF','UNKNOWN'),
        ('p4','SPARSE','CAVRDSNYQLIW','CASSLAEAYEQYF','TRAV12-2')]:
        rows += [dict(meta,complex_id=pair,epitope=e,gene='TRA',cdr3=a,v=av,j='TRAJ33'),
            dict(meta,complex_id=pair,epitope=e,gene='TRB',cdr3=b,v='TRBV19',j='TRBJ1-1')]
    # Independent observations establish eligibility, not unique keys or pair counts.
    rows += [dict(meta,complex_id='0',epitope='A',gene='TRA',cdr3='CAVRDSNYQLIW',v='TRAV12-2',j='TRAJ33')]*30
    rows += [dict(meta,complex_id='0',epitope='B',gene='TRB',cdr3='CASSLACAYEQYF',v='TRBV19',j='TRBJ1-1')]*30
    rows += [dict(meta,complex_id='0',epitope='ONLY-UNLINKED',gene='TRB',cdr3='CASSLAAAYEQYF',v='TRBV19',j='TRBJ1-1')]
    pl.DataFrame(rows).write_csv(path,separator='\t')


def arguments(tmp_path,*extra):
    p=argparse.ArgumentParser();reference_calls.register(p.add_subparsers())
    return p.parse_args(['reference-calls','--vdjdb',str(tmp_path/'ref.tsv'),'--locus','paired',
        '--scorer','stem4','--leave-epitope-out','--output-prefix',str(tmp_path/'out'),*extra])


@pytest.mark.parametrize('exact',['retain','exclude'])
@pytest.mark.parametrize('policy',['all-valid','curated-benchmark'])
def test_pair_source_eligibility_truth_and_physical_fold_oracle(tmp_path,exact,policy):
    reference(tmp_path/'ref.tsv')
    a=arguments(tmp_path,'--exact-policy',exact,'--reference-population',policy)
    assert reference_calls.main(a)==0
    calls=pl.read_csv(tmp_path/'out.calls.tsv',separator='\t')
    manifest=json.loads((tmp_path/'out.manifest.json').read_text())
    assert calls.height==3 and set(calls['epitope'])=={'A','B'}
    truth=pl.read_csv(tmp_path/'out.truth.tsv',separator='\t')
    assert set(truth['complex_id'])=={'p1','p2','p3'}
    assert manifest['native_population_batches']==3
    assert calls.filter(pl.col('v_alpha')=='UNKNOWN')['linked_available'].to_list()==[True]
    assert calls.filter(pl.col('v_alpha')=='UNKNOWN')['v_geometry_imputed_alpha'].to_list()==[True]
    assert calls.filter(pl.col('v_alpha')=='UNKNOWN')['beta_available'].to_list()==[True]
    for axis in ['linked','alpha','beta']:
        assert not calls.filter(pl.col('selected_'+axis+'_epitope')==pl.col('epitope')).height
    joint=calls.select(pl.max_horizontal('kernel_score','kernel_alpha_score','kernel_beta_score').alias('kernel_any_axis_score'))
    assert_frame_equal(joint,calls.select('kernel_any_axis_score'),check_exact=True)
    from vdjmatch import db
    r=db.load(tmp_path/'ref.tsv',species='HomoSapiens').filter(pl.col('reference_valid'))
    context=precompute_v_loop_distances(matrix='blosum62')
    def resolve(frame,seq,v,locus,column):
        return frame.with_columns(impute_reference_v_alleles(frame[v],context['model'],locus=locus)[0].alias(column)).filter(
            pl.col(column).is_not_null() & pl.col(seq).str.len_chars().ge(8))
    # Physical removal is independent of native withholding memberships.
    for row in calls.iter_rows(named=True):
        surviving=r.filter(pl.col('epitope')!=row['epitope'])
        for axis in ['linked','alpha','beta']:
            if axis=='linked':
                ref=build_paired_ref(surviving)
                ref=resolve(resolve(ref,'alpha','v_alpha','TRA','_allele_alpha'),'beta','v_beta','TRB','_allele_beta')
                q=pl.DataFrame({'query_id':[row['query_id']],'cdr3a':[row['alpha']],'cdr3b':[row['beta']],
                    'va':[row['v_alpha']],'vb':[row['v_beta']]})
                q=resolve(resolve(q,'cdr3a','va','TRA','_allelea'),'cdr3b','vb','TRB','_alleleb')
                query_keys=['cdr3a','cdr3b','_allelea','_alleleb'];search_keys=['alpha','beta','_allele_alpha','_allele_beta']
            else:
                locus='TRA' if axis=='alpha' else 'TRB'
                ref=resolve(surviving.filter(pl.col('gene')==locus),'cdr3','v',locus,'_allele')
                q=resolve(pl.DataFrame({'query_id':[row['query_id']],'cdr3':[row[axis]],'v':[row['v_'+axis]]}),
                    'cdr3','v',locus,'_allele')
                query_keys=search_keys=['cdr3','_allele']
            ref=ref.unique([*search_keys,'epitope','mhc_a','mhc_b','mhc_class'],maintain_order=True)
            groups=[({'epitope':e},q,part) for (e,),part in ref.group_by('epitope',maintain_order=True)]
            state,_=grouped_kernel_selection(q,ref,groups,query_keys=query_keys,search_keys=search_keys,
                paired=axis=='linked',cutoffs=ladder_cutoffs(26),radius=585,kernel_scale=585,scorer='v-bounded',
                matrix='blosum62',ctrim=3,stem_depth=4,exclude_exact=exact=='exclude',v_loop_distances=context)
            expected=finish_kernel_selection(q.select('query_id'),state)
            score='kernel_score' if axis=='linked' else 'kernel_'+axis+'_score'
            if not expected.height:assert row[score] is None
            else:
                value=expected['kernel_score'][0]
                assert row[score]==value
                assert row['selected_'+axis+'_epitope']==expected['selected_epitope'][0]
                ties='tied_peptides' if axis=='linked' else axis+'_tied_peptides'
                assert row[ties]==expected['tied_peptides'][0]


def test_paired_route_requires_explicit_withheld_stem4(tmp_path):
    a=arguments(tmp_path);a.leave_epitope_out=False
    with pytest.raises(ValueError,match='paired reference-calls'):reference_calls.main(a)
    a=arguments(tmp_path);a.scorer='original'
    with pytest.raises(ValueError,match='paired reference-calls'):reference_calls.main(a)


def test_native_paired_annotation_binding_and_linked_score_rejection(tmp_path):
    reference(tmp_path/'ref.tsv')
    reference_calls.main(arguments(tmp_path,'--exact-policy','retain','--reference-population','all-valid'))
    calls=pl.read_csv(tmp_path/'out.calls.tsv',separator='\t').filter(pl.col('v_alpha')!='UNKNOWN').head(1)
    row=calls.row(0,named=True)
    sample=tmp_path/'q.tsv'
    pl.DataFrame([dict(query_id=i,clone_id='pair',junction_aa=row[axis],v_call=row['v_'+axis],
        j_call=row['j_'+axis],locus=locus) for i,(axis,locus) in enumerate([('alpha','TRA'),('beta','TRB')])]).write_csv(sample,separator='\t')
    p=argparse.ArgumentParser();geometry.register(p.add_subparsers(),name='annotate')
    args=p.parse_args(['annotate',str(sample),'--vdjdb',str(tmp_path/'ref.tsv'),'--locus','paired',
        '--distance-matrix','blosum62','--junction-ends','trim3','--stem-depth','4','--kernel-scorer','v-bounded',
        '--selection-only','--all-targets','--output-prefix',str(tmp_path/'generated')])
    assert geometry.main(args)==0
    output=pl.read_csv(tmp_path/'generated.selection.tsv',separator='\t')
    assert output['kernel_any_axis_score'][0]==max(output[c][0] or 0 for c in ['kernel_score','kernel_alpha_score','kernel_beta_score'])
    manifest=json.loads((tmp_path/'out.manifest.json').read_text())
    generated=json.loads((tmp_path/'generated.manifest.json').read_text())
    binding=calibration.pipeline_binding(manifest)
    assert calibration.pipeline_binding(generated)==binding
    recipe=dict(pipeline_binding=binding,comparison='>=',action='apply_cutoff',cutoff=0.,score_column='kernel_any_axis_score')
    assert calibration.apply_calibration(output,recipe,input_manifest=generated).height==1
    for column in ['kernel_score','kernel_alpha_score','kernel_beta_score']:
        bad=copy.deepcopy(recipe);bad['score_column']=column
        with pytest.raises(ValueError,match='final higher-is-better'):calibration.apply_calibration(output,bad,input_manifest=generated)
    flat=output.select('query_id',pl.col('kernel_any_axis_score').alias('score'))
    recipe['score_column']='score'
    assert calibration.apply_calibration(flat,recipe,input_manifest=generated).height==1
    with pytest.raises(ValueError,match='input-manifest'):calibration.apply_calibration(flat,recipe)



def test_same_aa_pair_different_j_keeps_truth_complex_mapping(tmp_path):
    reference(tmp_path/'ref.tsv')
    rows=pl.read_csv(tmp_path/'ref.tsv',separator='\t')
    extra=rows.filter(pl.col('complex_id')=='p1').with_columns(pl.lit('p5').alias('complex_id'),
        pl.when(pl.col('gene')=='TRB').then(pl.lit('TRBJ2-1')).otherwise(pl.col('j')).alias('j'))
    pl.concat([rows,extra]).write_csv(tmp_path/'ref.tsv',separator='\t')
    reference_calls.main(arguments(tmp_path))
    calls=pl.read_csv(tmp_path/'out.calls.tsv',separator='\t')
    shared=calls.filter((pl.col('alpha')=='CAVRDSNYQLIW')&(pl.col('beta')=='CASSLAAAYEQYF'))
    assert shared.height==2 and shared['j_beta'].n_unique()==2
    for column in ['kernel_score','kernel_alpha_score','kernel_beta_score','kernel_any_axis_score']:
        assert shared[column].n_unique()==1
    truth=pl.read_csv(tmp_path/'out.truth.tsv',separator='\t')
    assert set(truth['complex_id'])=={'p1','p2','p3','p5'}
    assert truth['query_id'].n_unique()==4


def test_no_eligible_pairs_retains_empty_schema_and_batch_accounting(tmp_path):
    reference(tmp_path/'ref.tsv')
    reference_calls.main(arguments(tmp_path,'--min-target-records','100'))
    calls=pl.read_csv(tmp_path/'out.calls.tsv',separator='\t')
    assert not calls.height and 'kernel_any_axis_score' in calls.columns
    manifest=json.loads((tmp_path/'out.manifest.json').read_text())
    assert manifest['query_pairs']==0 and manifest['native_population_batches']==0


def test_conflicting_pair_annotations_rejected_before_output(tmp_path):
    reference(tmp_path/'ref.tsv')
    rows=pl.read_csv(tmp_path/'ref.tsv',separator='\t').with_columns(
        pl.when((pl.col('complex_id')=='p1')&(pl.col('gene')=='TRB')).then(pl.lit('HLA-A*03:01')).otherwise(pl.col('mhc_a')).alias('mhc_a'))
    rows.write_csv(tmp_path/'ref.tsv',separator='\t')
    with pytest.raises(ValueError,match='conflicting peptide/MHC'):reference_calls.main(arguments(tmp_path))
    assert not (tmp_path/'out.calls.tsv').exists()
