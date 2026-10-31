import json
import subprocess
import sys
from pathlib import Path
import polars as pl


def test_reference_panel_cli_preserves_rows_reads_zero_panels_and_fdr_groups(tmp_path, monkeypatch):
    (tmp_path/'sample.tsv').write_text('junction_aa\tv\tj\tcount\nCASSF\tTRBV2\tTRBJ1-1\t7\nCASSF\tTRBV2\tTRBJ1-1\t3\nCASRF\tTRBV2\tTRBJ1-1\t5\n')
    (tmp_path/'meta.tsv').write_text('sample_id\tpath\tfamily\nS\tsample.tsv\tF\n')
    (tmp_path/'reference.tsv').write_text('cdr3\tgene\tv\tj\tspecies\tepitope\tmhc_a\tmhc_b\tmhc_class\nCASSF\tTRB\tTRBV2\tTRBJ1-1\tHomoSapiens\tLLWNGPMAV\tHLA-A*02:01\tB2M\tMHCI\n')
    (tmp_path/'probes.tsv').write_text('junction_aa\tv\tj\tjoint_pgen\nCASSF\tTRBV2\tTRBJ1-1\t0.0001\nCASRF\tTRBV2\tTRBJ1-1\t0.0001\nCATSF\tTRBV2\tTRBJ1-1\t0.0001\n')
    (tmp_path/'model.tsv').write_text('v\tj\trecords\nTRBV2\tTRBJ1-1\t100\n')
    prefix=tmp_path/'nested'/'out'
    command=[sys.executable,'-m','vdjmatch.cli','repertoire-exact','--legacy',str(tmp_path/'meta.tsv'),'--vdjdb',str(tmp_path/'reference.tsv'),'--control',str(tmp_path/'probes.tsv'),'--model-vj-counts',str(tmp_path/'model.tsv'),'--threads','1','--fdr-group','family','--output-prefix',str(prefix)]
    run=subprocess.run(command,capture_output=True,text=True,timeout=30)
    assert run.returncode==0,run.stderr+run.stdout
    counts=pl.read_csv(str(prefix)+'.counts.tsv',separator='\t')
    assert counts.height==201
    assert counts.filter(pl.col('replicate')==-1).select('matched_records','matched_reads').row(0)==(2,10)
    assert set(counts.filter(pl.col('replicate')>=0)['matched_records'])=={0,1}
    summary=pl.read_csv(str(prefix)+'.summary.tsv',separator='\t').row(0,named=True)
    assert summary['p_records']==summary['p_reads']==1/201
    assert summary['significant_bh_records'] and summary['significant_bh_reads']
    assert summary['original_effect_records']==float('inf')
    assert summary['original_significant_bh_records']
    manifest=json.loads((tmp_path/'nested/out.manifest.json').read_text())
    assert manifest['inputs_frozen'] and manifest['source_frozen']
    assert manifest['fdr_family']==['group_id','family']
    assert manifest['legacy'] and manifest['selection_q'] == 1.
    # Exercise shipped profile loading, native dispatch, selected V/J audit and outputs.
    # The native batch is stubbed here; its weighted probabilities have separate exact tests.
    import argparse
    from vdjmatch.cli.repertoire_exact import register
    from vdjtools.model import native
    monkeypatch.setattr(native,'pgen_aa_gene_batch',lambda model,sequences,*args,**kwargs:[.0001]*len(sequences))
    profile_probes=tmp_path/'profile-probes.tsv'
    pl.read_csv(tmp_path/'probes.tsv',separator='\t').drop('joint_pgen').write_csv(profile_probes,separator='\t')
    parser=argparse.ArgumentParser()
    register(parser.add_subparsers())
    profile_prefix=tmp_path/'profile'
    args=parser.parse_args(['repertoire-exact',str(tmp_path/'meta.tsv'),'--vdjdb',str(tmp_path/'reference.tsv'),
        '--control',str(profile_probes),'--selection-profile','vj-length','--mock-replicates','5',
        '--threads','1','--output-prefix',str(profile_prefix)])
    assert args.func(args)==0
    modern=json.loads(Path(str(profile_prefix)+'.manifest.json').read_text())
    assert modern['selection_profile_kind']=='vj-length' and modern['selection_q_origin']=='profile'
    assert modern['model_source']=='isalgo' and not modern['legacy']
    assert pl.read_csv(str(profile_prefix)+'.counts.tsv',separator='\t').height==6
    assert 'model_probability' in pl.read_csv(str(profile_prefix)+'.vj-audit.tsv',separator='\t').columns
    pl.DataFrame({'v':['TRBV2'],'j':['TRBJ1-1'],'records':[1.5]}).write_parquet(tmp_path/'fractional.parquet')
    invalid=command.copy();invalid[invalid.index('--model-vj-counts')+1]=str(tmp_path/'fractional.parquet')
    invalid[invalid.index('--output-prefix')+1]=str(tmp_path/'invalid')
    failed=subprocess.run(invalid,capture_output=True,text=True,timeout=30)
    assert failed.returncode!=0 and 'integer' in failed.stderr.lower()


def test_exact_source_restricts_selected_locus(tmp_path):
    from vdjmatch.cli.repertoire_exact import _source
    path=tmp_path/'mixed.tsv'
    path.write_text('junction_aa\tv\tj\tcount\nCASSF\tTRBV2\tTRBJ1-1\t7\nCAVRF\tTRAV1-1\tTRAJ1\t9\n')
    assert _source(path,locus='TRB')['count'].to_list()==[7]
    assert _source(path,locus='TRA')['count'].to_list()==[9]


def test_joint_pgen_uses_shared_gene_batch_and_preserves_rows(monkeypatch):
    from vdjtools.model import native
    from vdjmatch.cli.repertoire_exact import _joint_pgen
    probes=pl.DataFrame({'junction_aa':['CASSYVGSEQFF']*2,'v':['TRBV19','TRBV7-9'],
                         'j':['TRBJ2-1']*2})
    original=native.pgen_aa_gene_batch
    calls=[]
    def batch(model,sequences,v,j,**kwargs):
        calls.append((sequences,v,j,kwargs))
        return original(model,sequences,v,j,**kwargs)
    monkeypatch.setattr(native,'pgen_aa_gene_batch',batch)
    result=_joint_pgen(probes,species='human',locus='TRB',threads=1)
    assert result.select(probes.columns).equals(probes)
    assert result['pgen'].min()>0
    assert len(calls)==1 and calls[0][3]=={'mismatches':0,'threads':1}



def test_compiled_defaults_use_model_specific_q_and_legacy_is_explicit():
    from types import SimpleNamespace
    import pytest
    from vdjtools.model import selection_factor
    from vdjmatch.cli.repertoire_exact import _selection_options
    for species in ('human', 'mouse'):
        for locus in ('TRA', 'TRB'):
            args = SimpleNamespace(species=species, locus=locus, legacy=False,
                                   source=None, selection_q=None)
            source, q, origin = _selection_options(args)
            assert source == 'isalgo' and origin == 'calibrated'
            assert q == selection_factor(locus, organism=species)
            with pytest.raises(ValueError, match='explicit'):
                _selection_options(args, external_pgen=True)
    args = SimpleNamespace(species='human', locus='TRB', legacy=True,
                           source=None, selection_q=None)
    assert _selection_options(args) == ('olga-legacy', 1., 'legacy')
    args.selection_q = 2.
    assert _selection_options(args, external_pgen=True) == ('olga-legacy', 2., 'explicit')


def test_compiled_command_has_modern_defaults():
    import argparse
    from vdjmatch.cli.repertoire_exact import register
    parser = argparse.ArgumentParser()
    register(parser.add_subparsers())
    args = parser.parse_args(['repertoire-exact', 'samples.tsv', '--vdjdb', 'reference.tsv',
                              '--control', 'control.tsv', '--model-vj-counts', 'model.tsv',
                              '--output-prefix', 'out'])
    assert args.selection_q is None and not args.legacy and args.source is None


def test_profile_joint_pgen_weights_before_aa_max(monkeypatch):
    from vdjmatch.cli.repertoire_exact import _joint_pgen
    from vdjtools.model import native
    monkeypatch.setattr(native,'pgen_aa_gene_batch',lambda *args,**kwargs:[.1,.2])
    shape={'kind':'vj-length','length_min':5,'log_length':[0.], 'log_normalizer':0.,
        'log_vj':[{'v_call':'TRBV19','j_call':'TRBJ2-1','log_q':0.},
                  {'v_call':'TRBV7-9','j_call':'TRBJ2-1','log_q':-2.}]}
    probes=pl.DataFrame({'junction_aa':['CASSF']*2,'v':['TRBV19','TRBV7-9'],'j':['TRBJ2-1']*2})
    result=_joint_pgen(probes,species='human',locus='TRB',threads=1,profile=shape)
    import numpy as np
    np.testing.assert_allclose(result['pgen'],[.1,.2*np.exp(-2)])
    assert result['pgen'].max()==.1
