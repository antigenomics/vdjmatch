import json
import subprocess
import sys
import polars as pl


def test_reference_panel_cli_preserves_rows_reads_zero_panels_and_fdr_groups(tmp_path):
    (tmp_path/'sample.tsv').write_text('junction_aa\tv\tj\tcount\nCASSF\tTRBV2\tTRBJ1-1\t7\nCASSF\tTRBV2\tTRBJ1-1\t3\nCASRF\tTRBV2\tTRBJ1-1\t5\n')
    (tmp_path/'meta.tsv').write_text('sample_id\tpath\tfamily\nS\tsample.tsv\tF\n')
    (tmp_path/'reference.tsv').write_text('cdr3\tgene\tv\tj\tspecies\tepitope\tmhc_a\tmhc_b\tmhc_class\nCASSF\tTRB\tTRBV2\tTRBJ1-1\tHomoSapiens\tLLWNGPMAV\tHLA-A*02:01\tB2M\tMHCI\n')
    (tmp_path/'probes.tsv').write_text('junction_aa\tv\tj\tjoint_pgen\nCASSF\tTRBV2\tTRBJ1-1\t0.0001\nCASRF\tTRBV2\tTRBJ1-1\t0.0001\nCATSF\tTRBV2\tTRBJ1-1\t0.0001\n')
    (tmp_path/'model.tsv').write_text('v\tj\trecords\nTRBV2\tTRBJ1-1\t100\n')
    prefix=tmp_path/'nested'/'out'
    command=[sys.executable,'-m','vdjmatch.cli','repertoire-exact',str(tmp_path/'meta.tsv'),'--vdjdb',str(tmp_path/'reference.tsv'),'--control',str(tmp_path/'probes.tsv'),'--model-vj-counts',str(tmp_path/'model.tsv'),'--threads','1','--fdr-group','family','--output-prefix',str(prefix)]
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
