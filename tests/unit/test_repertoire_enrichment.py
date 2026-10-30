import numpy as np
import polars as pl
import pytest
from vdjmatch.aggregate.repertoire_enrichment import matched_enrichment


def test_exact_bin_count_matching_and_both_endpoints():
    observations=pl.DataFrame({"sample_id":["s"]*3,"query_id":[0,1,2],"count":[7,3,5],"pgen_bin":[-8,-8,-12]})
    control=pl.DataFrame({"query_id":[3,4,5],"event_weight":[1,1,1],"pgen_bin":[-8,-8,-12]})
    targets=pl.DataFrame({"epitope":["X","Y"]})
    # Duplicate reference observations must not multiply a matched receptor.
    hits=pl.DataFrame({"query_id":[0,0,1,2,3],"epitope":["X","X","X","X","X"]})
    summary,mocks,audit=matched_enrichment(observations,control,hits,targets,replicates=200,seed=7)
    x=summary.filter(pl.col("epitope")=="X").row(0,named=True)
    y=summary.filter(pl.col("epitope")=="Y").row(0,named=True)
    assert (x["matched_records"],x["matched_reads"])==(3,15)
    assert x["p_records"]==x["p_reads"]==1/201
    assert x["excess_records"]==2 and x["enrichment_records"]==3
    assert y["matched_records"]==y["matched_reads"]==0 and y["p_records"]==y["p_reads"]==1
    assert mocks.height==400 and set(mocks["mock_records"])=={0,1}
    assert set(mocks.filter(pl.col("epitope")=="X")["mock_reads"])=={3,7}
    assert audit["requested_records"].sum()==3 and audit["exact_histogram_and_counts"].all()
    assert mocks.equals(matched_enrichment(observations,control,hits,targets,replicates=200,seed=7)[1])
    with pytest.raises(ValueError,match="Pgen bin -12"):
        matched_enrichment(observations,control.filter(pl.col("pgen_bin")==-8),hits,targets,replicates=200)
    with pytest.raises(ValueError,match="Pgen bin -8"):
        matched_enrichment(observations,control.filter(pl.col("query_id")!=4),hits,targets,replicates=200)


@pytest.mark.parametrize("mismatches,missing_metadata",[(0,False),(1,False),(0,True),(0,"unknown")])
def test_repertoire_enrichment_cli_end_to_end(tmp_path,mismatches,missing_metadata):
    import os,subprocess,sys,json
    from pathlib import Path
    sample=tmp_path/"sample.tsv"
    sample.write_text("junction_aa\tv_call\tj_call\tlocus\tduplicate_count\nCASSF\tTRBV2\tTRBJ1-1\tTRB\t3\nCASSF\tTRBV2\tTRBJ1-1\tTRB\t4\n")
    if missing_metadata:
        replacement='UNKNOWN' if missing_metadata=='unknown' else ''
        sample.write_text(sample.read_text().replace('CASSF\tTRBV2\tTRBJ1-1\tTRB\t4',f'CASSF\t{replacement}\tTRBJ1-1\t\t4'))
    metadata=tmp_path/"samples.tsv";metadata.write_text("sample_id\tpath\nS0\tsample.tsv\n")
    control=tmp_path/"control.tsv";control.write_text("junction_aa\tv\tj\toccurrences\nCASSF\tTRBV2\tTRBJ1-1\t100\n")
    reference=tmp_path/"reference.tsv";reference.write_text("cdr3\tgene\tv\tj\tspecies\tepitope\tmhc_a\tmhc_b\tmhc_class\nCASSF\tTRB\tTRBV2\tTRBJ1-1\tHomoSapiens\tLLWNGPMAV\tHLA-A*02:01\tB2M\tMHCI\n")
    if missing_metadata=='unknown':
        reference.write_text(reference.read_text().replace('LLWNGPMAV','KLGGALQAK'))
    prefix=tmp_path/"fresh"
    command=[sys.executable,"-m","vdjmatch.cli","repertoire-enrichment",str(metadata),"--vdjdb",str(reference),"--control",str(control),"--threads","1","--pgen-mismatches",str(mismatches),"--output-prefix",str(prefix)]
    r=subprocess.run(command,text=True,capture_output=True,timeout=60)
    assert r.returncode==0,r.stdout+r.stderr
    summary=pl.read_csv(str(prefix)+".summary.tsv",separator="\t")
    assert summary["matched_records"].item()==2 and summary["matched_reads"].item()==7
    assert summary["p_records"].item()==summary["p_reads"].item()==1
    manifest=json.loads(Path(str(prefix)+".manifest.json").read_text())
    assert manifest["source_frozen"] and manifest["inputs_frozen"] and manifest["native_search_batches"]==201
    assert manifest["automatic_detection_threshold"] is None
    assert manifest["pgen_mismatches"]==mismatches
    ingestion=pl.read_csv(str(prefix)+".ingestion.tsv",separator="\t")
    assert ingestion["selected_records"].item()==2 and ingestion["imputed_v_records"].item()==int(bool(missing_metadata))


def test_tcr_source_counts_are_exact(tmp_path):
    from vdjmatch import io
    p=tmp_path/"source.tsv"
    p.write_text("Read.count\tCDR3.nucleotide.sequence\tCDR3.amino.acid.sequence\tbestVGene\tbestJGene\n9007199254740993\tTGTGCCAGCAGCTTT\tCASSF\tTRBV2\tTRBJ1-1\n1\tTGTGCCAGCAGCTTT\tCASSF\tTRBV2,TRBV3\t\n")
    frame=io.read_rearrangement(p,valid_aa=False)
    assert frame["count"].to_list()==[9007199254740993,1]
    assert frame["v"].to_list()==["TRBV2","TRBV2,TRBV3"]


def test_pgen_batch_respects_requested_mass(monkeypatch):
    from vdjmatch.evalue.mock_reference import compute_pgen
    from vdjtools.model import native
    observed=[]
    def batch(model,sequences,*,mismatches,threads):
        observed.append((sequences,mismatches,threads))
        return [0.001*(mismatches+1)]*len(sequences)
    monkeypatch.setattr(native,"pgen_aa_batch",batch)
    keys=pl.DataFrame({"junction_aa":["CASSF","CASSF"]})
    for mm in (0,1):
        result=compute_pgen(keys,keys.head(0),mismatches=mm,threads=2)
        assert result.height==1 and result["pgen"].item()==0.001*(mm+1)
    assert observed==[(["CASSF"],0,2),(["CASSF"],1,2)]


def test_v_ambiguity_policy():
    from vdjmatch.cli.repertoire_enrichment import _single_v
    frame=pl.DataFrame({'v':['TRBV2*01,TRBV2*02','TRBV2,TRBV3',None,''], 'count':[1,2,3,4]})
    result=_single_v(frame,'v')
    assert result['v'].to_list()==['TRBV2']*4
    assert result['original_v'].to_list()==frame['v'].to_list()
    assert result['count'].to_list()==[1,2,3,4]
    assert result['v_imputed'].to_list()==[False,True,True,True]
    all_missing=_single_v(pl.DataFrame({'v':[None,'']}),'v')
    assert all_missing['v'].str.starts_with('TRBV').all() and all_missing['v_imputed'].all()


def test_streamed_mock_measurements_equal_precomputed_memberships():
    observations=pl.DataFrame(dict(sample_id=['s']*3,query_id=[0,1,2],count=[7,3,5],pgen_bin=[-8,-8,-12]))
    control=pl.DataFrame(dict(query_id=[3,4,5],event_weight=[1,1,1],pgen_bin=[-8,-8,-12]))
    targets=pl.DataFrame(dict(epitope=['X','Y']))
    hits=pl.DataFrame(dict(query_id=[0,0,1,2,3],epitope=['X']*5))
    sizes=[]
    def measure(rows,counts):
        sizes.append(len(rows))
        matched=np.isin(rows,[0,1,2,3])
        return np.array([matched.sum(),0]),np.array([counts[matched].sum(),0])
    expected=matched_enrichment(observations,control,hits,targets,seed=7)
    actual=matched_enrichment(observations,control,None,targets,seed=7,measure=measure)
    assert all(a.equals(b) for a,b in zip(expected,actual))
    assert sizes==[3]*201


def test_floored_pgen_binning_includes_zero_and_extreme_tail():
    from vdjmatch.cli.repertoire_enrichment import _repertoire_strata
    x=pl.DataFrame(dict(junction_aa=['CASSF','CATSF','CARTF','CASRF'],pgen=[0.,2.**-112,2.**-61,2.**-20]))
    result=_repertoire_strata(x)
    assert dict(zip(result['junction_aa'],result['pgen_bin']))==dict(CASSF=-60,CATSF=-60,CARTF=-59,CASRF=-20)


def test_control_allele_pooling_preserves_complete_baseline_and_missing_keys():
    from vdjmatch.cli.repertoire_enrichment import _single_v
    from vdjmatch.evalue.mock_reference import control_events
    frame=pl.DataFrame({'junction_aa':['CASSF']*2,'junction':['TGTGCCAGCAGCTTT']*2,
        'v':['TRBV19*01','TRBV19*02'],'j':['TRBJ1-1']*2})
    complete=control_events(_single_v(frame,'v'))
    assert complete['event_weight'].sum()==1
    missing=control_events(_single_v(frame.with_columns(pl.Series('v',['TRBV19',None])),'v'))
    assert missing['event_weight'].sum()==2 and missing['v_imputed'].all()


def test_mixed_control_population_selected_before_imputation():
    from vdjmatch.cli.repertoire_enrichment import _single_v
    from vdjmatch.evalue.mock_reference import control_events
    raw=pl.DataFrame({'junction_aa':['CASSF']*3,'v':['TRAV1-1*01','TRBV19*01',None],
        'j':['TRAJ1','TRBJ1-1','TRBJ1-1'],'occurrences':[10,2,3]})
    pooled=_single_v(raw,'v',impute=False)
    events=_single_v(control_events(pooled,locus='TRB'),'v',impute=False)
    assert events['event_weight'].sum()==5 and events['v'].to_list()==['TRBV19']
    assert events['v_imputed'].all()


def test_reference_mock_enrichment_counts_receptors_once_and_retains_no_hits():
    from vdjmatch.aggregate import repertoire_enrichment as module
    assert hasattr(module,'reference_enrichment'), 'M3 reference enrichment is missing'
    sample=pl.DataFrame(dict(receptor_id=['r1','r1','r2'],locus=['TRA','TRB','TRB'],
        junction_aa=['CAVSF','CASSF','CATSF'],v=['TRAV1','TRBV2','TRBV2'],j=['TRAJ1','TRBJ1','TRBJ1'],count=[3,3,2]))
    reference=pl.DataFrame(dict(group_id=['X','X','X','Y'],locus=['TRA','TRB','TRB','TRB'],
        junction_aa=['CAVSF','CASSF','CASSF','CWWEF'],v=['TRAV1','TRBV2','TRBV2','TRBV3'],j=['TRAJ1','TRBJ1','TRBJ1','TRBJ2'],pgen_bin=[-10]*4))
    controls={locus:pl.DataFrame(dict(junction_aa=['CQQQF','CTTTF','CPPPF','CNNNF'],
        v=[v]*4,j=[j]*4,pgen_bin=[-10]*4)) for locus,v,j in [('TRA','TRAV2','TRAJ2'),('TRB','TRBV4','TRBJ2')]}
    summary,mocks,hits=module.reference_enrichment(sample,reference,controls,replicates=3,threads=1)
    x=summary.filter(pl.col('group_id')=='X').row(0,named=True)
    y=summary.filter(pl.col('group_id')=='Y').row(0,named=True)
    assert (x['matched_receptors'],x['matched_count'],x['p_receptors'])==(2,5,.25)
    assert (y['matched_receptors'],y['p_receptors'])==(0,1.)
    assert mocks.height==6 and mocks['mock_receptors'].sum()==0
    assert hits.select('receptor_id','group_id').n_unique()==2
    assert x['n_receptors']==2 and x['total_count']==5
    parallel=module.reference_enrichment(sample,reference,controls,replicates=3,threads=2)
    assert all(a.equals(b) for a,b in zip((summary,mocks,hits),parallel))
    with pytest.raises(ValueError,match='count'):
        module.reference_enrichment(sample.with_columns(pl.Series('count',[3,4,2])),reference,controls,replicates=3)


@pytest.mark.parametrize("supported", [True, False])
def test_reference_mode_cli_and_significance_options(tmp_path,supported):
    import subprocess,sys,json
    from pathlib import Path
    sample=tmp_path/'sample.tsv'
    sample.write_text('junction_aa\tv_call\tj_call\tlocus\tduplicate_count\nCASSF\tTRBV2\tTRBJ1-1\tTRB\t7\n')
    (tmp_path/'metadata.tsv').write_text('sample_id\tpath\nS\tsample.tsv\n')
    (tmp_path/'reference.tsv').write_text('cdr3\tgene\tv\tj\tspecies\tepitope\tmhc_a\tmhc_b\tmhc_class\nCASSF\tTRB\tTRBV2\tTRBJ1-1\tHomoSapiens\tSYNTHETIC\tHLA-A*02:01\tB2M\tMHCI\n')
    (tmp_path/'control.tsv').write_text('junction_aa\tv\tj\toccurrences\nCASSF\t'+('TRBV3' if supported else 'TRBV2')+'\tTRBJ1-1\t100\n')
    command=[sys.executable,'-m','vdjmatch.cli','repertoire-enrichment',str(tmp_path/'metadata.tsv'),
        '--null','reference','--vdjdb',str(tmp_path/'reference.tsv'),'--control',str(tmp_path/'control.tsv'),
        '--pgen-mismatches','0','--threads','1','--significance','both','--output-prefix',str(tmp_path/'result')]
    result=subprocess.run(command,capture_output=True,text=True,timeout=30)
    assert result.returncode==0,result.stdout+result.stderr
    summary=pl.read_csv(tmp_path/'result.summary.tsv',separator='\t')
    assert summary['matched_receptors'].item()==1
    if supported:
        assert summary['p_receptors'].item()==1/201
        assert summary['significant_nominal'].item() and summary['significant_bh'].item()
    else:
        assert summary['status'].item()=='untestable_null'
        assert summary['p_receptors'].item() is None and summary['q_receptors'].item() is None
        assert summary['significant_nominal'].item() is None and summary['significant_bh'].item() is None
        assert pl.read_parquet(tmp_path/'result.mocks.parquet')['mock_receptors'].null_count()==200
    manifest=json.loads((tmp_path/'result.manifest.json').read_text())
    assert manifest['null']=='reference' and manifest['inputs_frozen'] and manifest['source_frozen']
    assert manifest['mock_policy']=='pgen' and manifest['counting_unit']=='unique receptor identity'
    assert pl.read_csv(tmp_path/'result.hits.tsv',separator='\t').height==1


def test_reference_cli_rejects_unattainable_tail_resolution(tmp_path):
    from argparse import Namespace
    from vdjmatch.cli.repertoire_reference import main
    with pytest.raises(ValueError,match='resolution'):
        main(Namespace(mock_replicates=200,seed=0,threads=1,nominal_alpha=.001,fdr_alpha=.1,significance='nominal'))


def test_projected_control_read_preserves_population_and_keys(tmp_path):
    from vdjmatch.io.airr import _read_table
    from vdjmatch.evalue.mock_reference import control_events
    import polars.selectors as cs
    p=tmp_path/'control.tsv'
    p.write_text('cdr3aa\tv\tj\toccurrences\tspecies\tlocus\tunused\nCASSF\tTRBV2\tTRBJ1-1\t3\thuman\tTRB\textra\nCATSF\tTRBV3\tTRBJ1-1\t0\thuman\tTRB\textra\n')
    projected=_read_table(p,columns=cs.exclude('unused'))
    assert projected.columns==['cdr3aa','v','j','occurrences','species','locus']
    assert control_events(projected).equals(control_events(_read_table(p)))
