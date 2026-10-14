"""Bound final-score cutoffs reject incompatible score producers before export."""
import argparse
import copy
import hashlib
import json
from importlib.resources import files
from pathlib import Path

import polars as pl
import pytest

from vdjmatch.cli import calibration
from vdjmatch.evalue.kernel import ladder_cutoffs
from vdjmatch.match.tcrdist import load_v_loops


def manifests():
    model_file = Path(str(files('vdjmatch.resources') / 'tcrdist/human_v_loops.tsv'))
    file_hash = hashlib.sha256(model_file.read_bytes()).hexdigest()
    semantic = hashlib.sha256(json.dumps(load_v_loops('human'), sort_keys=True, default=str).encode()).hexdigest()
    sources = {k: str(i) * 64 for i, k in enumerate(('kernel', 'tcrdist', 'seqtree_native', 'seqtree_python'), 1)}
    geometry = dict(species='human', matrix='blosum62', ctrim=3, v_loop_weight=1,
        alphabet='full', stem_depth=4, stem_weight=1., exclude_exact=True, radius=585,
        kernel_scale=585, cutoffs=ladder_cutoffs(26), cross_v=.25, floor=.01,
        scorer='v-bounded', matrix_scale=26, v_regions='cdr12plus', position_weighting='uniform')
    reference = dict(options=dict(species='HomoSapiens', locus='TRB', scorer='stem4'),
        reference=dict(input_sha256='a' * 64), reference_population='curated-benchmark',
        geometry=geometry, source_sha256=sources, model_sha256=semantic,
        selection='maximum restriction per peptide, then maximum peptide; exact score ties remain ambiguous',
        calibration_statistic='final selected kernel_score; no component P/E transfer')
    annotation = dict(parameters=dict(species='HomoSapiens', locus='TRB', exclude_exact=True,
        all_targets=True, background='vdjdb-other', neighbour_weighting='kernel', position_weighting='uniform'),
        reference=reference['reference'], reference_population='curated-benchmark',
        distance_scale=dict(matrix='blosum62', alphabet='full', radius=585, v_loop_weight=1,
            stem_depth=4, stem_weight=1., matrix_unit=26, v_regions='cdr12plus'),
        junction_end_trim=dict(n_terminal=3, c_terminal=3),
        kernel=dict(kernel_scale=585, cutoffs=ladder_cutoffs(26), cross_v=.25, scorer='v-bounded', score_column='kernel_score'),
        software=dict(kernel_source_sha256=sources['kernel'], scorer_source_sha256=sources['tcrdist'],
            seqtree_native_sha256=sources['seqtree_native'], seqtree_python_source_sha256=sources['seqtree_python']),
        v_loop_model=dict(sha256=file_hash))
    binding = dict(context=['HomoSapiens','TRB'], reference_population='curated-benchmark',
        reference_sha256='a' * 64, geometry=geometry, sources=sources, model_file_sha256=file_hash,
        model_semantic_sha256=semantic, statistic='final selected kernel_score',
        selection='max restriction then max peptide; ties retained',
        exact_policy='full-junction puncture across V/J; no ordinary exact-retained transfer')
    return reference, annotation, binding


def test_canonical_binding_both_shipped_manifests_and_single_cohort():
    reference, annotation, binding = manifests()
    assert calibration.pipeline_binding(reference) == binding
    assert calibration.pipeline_binding(annotation) == binding
    assert calibration.pipeline_binding({'core_manifests':[annotation]}) == binding
    with pytest.raises(ValueError, match='binding'):
        calibration.pipeline_binding({'core_manifests':[annotation, annotation]})
    reference['geometry']['exclude_exact']=False
    annotation['parameters']['exclude_exact']=False
    retained=calibration.pipeline_binding(reference)
    assert calibration.pipeline_binding(annotation)==retained
    assert retained['exact_policy']=='exact evidence retained across all labels; ordinary annotation'
    assert retained!=binding


@pytest.mark.parametrize('path,value', [
    (('parameters','species'),'MusMusculus'), (('parameters','locus'),'TRA'),
    (('reference','input_sha256'),'b'*64), (('reference_population',),'all-valid'),
    (('parameters','exclude_exact'),False), (('distance_scale','stem_depth'),3),
    (('kernel','scorer'),'all-bounded'), (('kernel','score_column'),'E'),
    (('v_loop_model','sha256'),'b'*64), (('software','seqtree_native_sha256'),'b'*64),
    (('software','kernel_source_sha256'),'b'*64), (('parameters','all_targets'),False),
    (('parameters','locus'),'paired'), (('parameters','withhold_epitope'),['TEST']),
])
def test_bound_cli_mismatch_fails_before_any_output(tmp_path, path, value):
    _, manifest, binding = manifests()
    node=manifest
    for key in path[:-1]: node=node[key]
    node[path[-1]]=value
    cases=tmp_path/'cases.tsv';pl.DataFrame({'query_id':['q'],'score':[2.]}).write_csv(cases,separator='\t')
    recipe=tmp_path/'recipe.json';recipe.write_text(json.dumps(dict(action='apply_cutoff',cutoff=1.,comparison='>=',pipeline_binding=binding)))
    evidence=tmp_path/'evidence.json';evidence.write_text(json.dumps(manifest))
    target=tmp_path/'new'/'calls.tsv'
    args=argparse.Namespace(cases=str(cases), calibration=str(recipe), input_manifest=str(evidence), null_cases=None, output=str(target))
    with pytest.raises(ValueError, match='binding'):
        calibration.apply_main(args)
    assert not target.parent.exists()


def test_bound_requires_manifest_and_all_fields_and_matching_application(tmp_path):
    reference, manifest, binding = manifests()
    rows=pl.DataFrame({'query_id':['q','missing'],'score':[2.,None]})
    recipe=dict(action='apply_cutoff',cutoff=1.,comparison='>=',pipeline_binding=binding)
    with pytest.raises(ValueError, match='input.manifest'):
        calibration.apply_calibration(rows, recipe)
    for key in ('reference_population','reference','v_loop_model','software','kernel'):
        incomplete=copy.deepcopy(manifest);incomplete.pop(key)
        with pytest.raises(ValueError,match='binding'):
            calibration.apply_calibration(rows,recipe,input_manifest=incomplete)
    assert calibration.apply_calibration(rows,recipe,input_manifest=manifest)['call_status'].to_list()==['accepted','no_supported_match']
    for incompatible in ({'score_column':'E'}, {'comparison':'<='}):
        with pytest.raises(ValueError,match='binding'):
            calibration.apply_calibration(rows,recipe | incompatible,input_manifest=manifest)
    cases=tmp_path/'cases.tsv';rows.write_csv(cases,separator='\t')
    frozen=tmp_path/'recipe.json';frozen.write_text(json.dumps(recipe))
    evidence=tmp_path/'manifest.json';evidence.write_text(json.dumps(reference))
    parser=argparse.ArgumentParser();calibration.register(parser.add_subparsers())
    output=tmp_path/'calls.tsv'
    args=parser.parse_args(['apply-call-calibration',str(cases),'--calibration',str(frozen),
        '--input-manifest',str(evidence),'--output',str(output)])
    assert calibration.apply_main(args)==0
    applied=json.loads(Path(str(output)+'.manifest.json').read_text())
    assert applied['pipeline_binding']==binding
    assert applied['input_manifest_sha256']==hashlib.sha256(evidence.read_bytes()).hexdigest()


def test_actual_cli_producers_have_matching_excluded_and_retained_bindings(tmp_path):
    import os
    import subprocess
    import sys

    reference=tmp_path/'reference.tsv'
    pl.DataFrame(dict(cdr3=['CASSLAAAYEQYF','CASSLACAYEQYF'],
        v=['TRBV19','TRBV19'],j=['TRBJ1-1','TRBJ1-1'],epitope=['A','B'],
        gene=['TRB','TRB'],species=['HomoSapiens','HomoSapiens'],
        mhc_class=['MHCI','MHCI'],mhc_a=['HLA-A*02:01','HLA-A*02:01'],
        mhc_b=['B2M','B2M'])).write_csv(reference,separator='\t')
    query=tmp_path/'query.tsv'
    pl.DataFrame(dict(sequence_id=['q'],junction_aa=['CASSLAAAYEQYF'],
        v_call=['TRBV19'],j_call=['TRBJ1-1'],locus=['TRB'])).write_csv(query,separator='\t')
    environment=dict(os.environ,POLARS_MAX_THREADS='1',OPENBLAS_NUM_THREADS='1',
        MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',OMP_THREAD_LIMIT='1')
    bindings={}
    for exact in ('exclude','retain'):
        produced=[]
        for route in ('reference-calls','annotate'):
            prefix=tmp_path/(route+'-'+exact)
            args=[sys.executable,'-m','vdjmatch.cli',route,'--vdjdb',str(reference),
                '--locus','TRB','--threads','1','--reference-population','all-valid',
                '--output-prefix',str(prefix)]
            if route=='reference-calls':
                args+=['--scorer','stem4','--min-target-records','1','--exact-policy',exact]
            else:
                args+=[str(query),'--all-targets']
                if exact=='exclude':args+=['--exclude-exact']
            result=subprocess.run(args,capture_output=True,text=True,env=environment,timeout=30)
            assert result.returncode==0,result.stderr
            manifest=json.loads(Path(str(prefix)+'.manifest.json').read_text())
            produced.append(calibration.pipeline_binding(manifest))
        cases=tmp_path/('calibration-'+exact+'.tsv')
        pl.DataFrame(dict(query_id=['positive','null','nohit'],score=[2.,.5,None],
            outcome=['correct','null','null'],null_source=[None,'synthetic','synthetic'])).write_csv(cases,separator='\t')
        frozen=tmp_path/('frozen-'+exact)
        producer=tmp_path/('reference-calls-'+exact+'.manifest.json')
        result=subprocess.run([sys.executable,'-m','vdjmatch.cli','calibrate-calls',str(cases),
            '--target-fpr','.01','--input-manifest',str(producer),'--output-prefix',str(frozen)],
            capture_output=True,text=True,env=environment,timeout=30)
        assert result.returncode==0,result.stderr
        recipe=json.loads(Path(str(frozen)+'.calibration.json').read_text())
        assert recipe['pipeline_binding']==produced[0]
        assert recipe['input_manifest_sha256']==hashlib.sha256(producer.read_bytes()).hexdigest()
        applied=tmp_path/('applied-'+exact+'.tsv')
        result=subprocess.run([sys.executable,'-m','vdjmatch.cli','apply-call-calibration',str(cases),
            '--calibration',str(frozen)+'.calibration.json',
            '--input-manifest',str(tmp_path/('annotate-'+exact+'.manifest.json')),'--output',str(applied)],
            capture_output=True,text=True,env=environment,timeout=30)
        assert result.returncode==0,result.stderr
        assert pl.read_csv(applied,separator='\t')['call_status'].to_list()==['accepted','no_supported_match','no_supported_match']
        assert produced[0]==produced[1]
        bindings[exact]=produced[0]
    assert bindings['exclude']!=bindings['retain']
    assert bindings['exclude']['geometry']['exclude_exact'] is True
    assert bindings['retain']['geometry']['exclude_exact'] is False


@pytest.mark.parametrize('failure',['component','reverse','missing_identity'])
def test_bound_calibration_freeze_fails_before_outputs(tmp_path,failure):
    reference,_,_=manifests()
    if failure=='missing_identity':reference.pop('source_sha256')
    producer=tmp_path/'producer.json';producer.write_text(json.dumps(reference))
    cases=tmp_path/'cases.tsv'
    pl.DataFrame(dict(query_id=['label','null'],score=[2.,.5],E=[.1,.2],
        outcome=['correct','null'],null_source=[None,'synthetic'])).write_csv(cases,separator='\t')
    parser=argparse.ArgumentParser();calibration.register(parser.add_subparsers())
    prefix=tmp_path/'new'/'calibration'
    options=['--score-column','E'] if failure=='component' else ['--lower-is-better'] if failure=='reverse' else []
    args=parser.parse_args(['calibrate-calls',str(cases),'--target-fpr','.01',
        '--input-manifest',str(producer),'--output-prefix',str(prefix),*options])
    with pytest.raises(ValueError,match='binding'):
        calibration.main(args)
    assert not prefix.parent.exists()
