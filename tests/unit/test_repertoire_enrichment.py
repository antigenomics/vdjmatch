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
    for mode,expected in [('tra',1),('trb',2),('intersection',1)]:
        separate=module.reference_enrichment(sample,reference,controls,replicates=3,threads=1,chain_mode=mode)[0]
        assert separate.filter(pl.col('group_id')=='X')['matched_receptors'].item()==expected
    combined=module.reference_enrichment(sample,reference,controls,replicates=3,chain_mode='all')
    for mode in ['tra','trb','union','intersection']:
        separate=module.reference_enrichment(sample,reference,controls,replicates=3,chain_mode=mode)
        for actual,expected in zip(combined,separate):
            selected=actual.filter(pl.col('chain_mode')==mode).drop('chain_mode').select(expected.columns)
            assert selected.sort(expected.columns).equals(expected.sort(expected.columns))
    parallel=module.reference_enrichment(sample,reference,controls,replicates=3,threads=2)
    assert all(a.equals(b) for a,b in zip((summary,mocks,hits),parallel))
    with pytest.raises(ValueError,match='count'):
        module.reference_enrichment(sample.with_columns(pl.Series('count',[3,4,2])),reference,controls,replicates=3)


@pytest.mark.parametrize("supported", [True, False])
@pytest.mark.parametrize("chain_mode", ["union", "all"])
def test_reference_mode_cli_and_significance_options(tmp_path,supported,chain_mode):
    import subprocess,sys,json
    from pathlib import Path
    sample=tmp_path/'sample.tsv'
    sample.write_text('junction_aa\tv_call\tj_call\tlocus\tduplicate_count\nCASSF\tTRBV2\tTRBJ1-1\tTRB\t7\n')
    (tmp_path/'metadata.tsv').write_text('sample_id\tpath\nS\tsample.tsv\n')
    (tmp_path/'reference.tsv').write_text('cdr3\tgene\tv\tj\tspecies\tepitope\tmhc_a\tmhc_b\tmhc_class\nCASSF\tTRB\tTRBV2\tTRBJ1-1\tHomoSapiens\tSYNTHETIC\tHLA-A*02:01\tB2M\tMHCI\n')
    (tmp_path/'control.tsv').write_text('junction_aa\tv\tj\toccurrences\nCASSF\t'+('TRBV3' if supported else 'TRBV2')+'\tTRBJ1-1\t100\n')
    command=[sys.executable,'-m','vdjmatch.cli','repertoire-enrichment',str(tmp_path/'metadata.tsv'),
        '--null','reference','--vdjdb',str(tmp_path/'reference.tsv'),'--control',str(tmp_path/'control.tsv'),
        '--pgen-mismatches','0','--threads','1','--significance','both','--chain-mode',chain_mode,'--output-prefix',str(tmp_path/'result')]
    if chain_mode=='all':command.extend(['--model-source','isalgo'])
    result=subprocess.run(command,capture_output=True,text=True,timeout=30)
    assert result.returncode==0,result.stdout+result.stderr
    summary=pl.read_csv(tmp_path/'result.summary.tsv',separator='\t')
    if chain_mode=='all':
        assert summary.height==4
        assert summary.filter(pl.col('chain_mode').is_in(['tra','intersection']))['matched_receptors'].sum()==0
        summary=summary.filter(pl.col('chain_mode')=='union')
    assert summary['matched_receptors'].item()==1
    if supported:
        assert summary['p_receptors'].item()==1/201
        assert summary['significant_nominal'].item() and summary['significant_bh'].item()
    else:
        assert summary['status'].item()=='untestable_null'
        assert summary['p_receptors'].item() is None and summary['q_receptors'].item() is None
        assert summary['significant_nominal'].item() is None and summary['significant_bh'].item() is None
        assert pl.read_parquet(tmp_path/'result.mocks.parquet')['mock_receptors'].null_count()==(800 if chain_mode=='all' else 200)
    manifest=json.loads((tmp_path/'result.manifest.json').read_text())
    assert manifest['null']=='reference' and manifest['inputs_frozen'] and manifest['source_frozen']
    if chain_mode=='all':assert manifest['models']['TRB']['source']=='isalgo'
    assert manifest['mock_policy']=='pgen' and manifest['counting_unit']=='unique receptor identity'
    assert pl.read_csv(tmp_path/'result.hits.tsv',separator='\t').height==(2 if chain_mode=='all' else 1)


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


def test_exact_reference_panel_counts_preserve_source_rows_and_reads():
    from vdjmatch.aggregate.repertoire_enrichment import exact_reference_counts
    obs=pl.DataFrame(dict(junction_aa=['CASSF','CASSF','CASSF','CASRF'],
        v=['V1','V1','V1','V1'],j=['J1','J1','J2','J1'],count=[7,3,20,4]))
    panels=pl.DataFrame(dict(group_id=['X','X','X','Y'],replicate=[-1,0,1,-1],
        junction_aa=['CASSF','CASRF','CATSF','CATSF'],v=['V1']*4,j=['J1']*4))
    out=exact_reference_counts(obs,panels)
    assert out.select('group_id','replicate','matched_records','matched_reads').rows()==[
        ('X',-1,2,10),('X',0,1,4),('X',1,0,0),('Y',-1,0,0)]
    with pytest.raises(ValueError,match='unique'):
        exact_reference_counts(obs,pl.concat([panels,panels.head(1)]))
    with pytest.raises(ValueError,match='nonnegative'):
        exact_reference_counts(obs.with_columns(pl.lit(-1).alias('count')),panels)


@pytest.mark.parametrize('mismatches,expected',[(0,['exact']),(1,['exact','one'])])
def test_reference_lookup_zero_and_one_mismatch_memberships(mismatches,expected):
    from vdjmatch.aggregate.repertoire_enrichment import reference_enrichment
    sample=pl.DataFrame(dict(receptor_id=['exact','one','two','wrong_v'],locus=['TRB']*4,
        junction_aa=['CASSF','CATSF','CATRF','CASSF'],v=['TRBV2']*3+['TRBV3'],j=['TRBJ1']*4,count=[2,3,5,7]))
    reference=pl.DataFrame(dict(group_id=['X'],locus=['TRB'],junction_aa=['CASSF'],v=['TRBV2'],j=['TRBJ1'],pgen_bin=[-10]))
    controls={'TRB':pl.DataFrame(dict(junction_aa=['CQQQF','CWWWF','CNNNF'],v=['TRBV2']*3,j=['TRBJ1']*3,pgen_bin=[-10]*3))}
    summary,mocks,hits=reference_enrichment(sample,reference,controls,mismatches=mismatches,replicates=3,threads=1)
    assert hits['receptor_id'].sort().to_list()==expected
    assert summary['matched_count'].item()==(2 if mismatches==0 else 5)
    assert mocks['mock_receptors'].sum()==0


def test_strict_beta_with_fuzzy_alpha_preserves_pairing():
    from vdjmatch.aggregate.repertoire_enrichment import reference_enrichment
    sample=pl.DataFrame(dict(receptor_id=['pair','pair','beta_only','alpha_only'],locus=['TRA','TRB','TRB','TRA'],
        junction_aa=['CAVRDSNFQLIW','CASSLGRAYEQYF','CASSLGRAYEQYF','CAVRDSNFQLIW'],
        v=['TRAV12-2*01','TRBV19*01','TRBV19*01','TRAV12-2*01'],j=['TRAJ33*01','TRBJ2-7*01','TRBJ2-7*01','TRAJ33*01'],count=[5,5,3,2]))
    reference=pl.DataFrame(dict(group_id=['X','X'],locus=['TRA','TRB'],junction_aa=['CAVRDSNYQLIW','CASSLGQAYEQYF'],
        v=['TRAV12-2*01','TRBV19*01'],j=['TRAJ33*01','TRBJ2-7*01'],pgen_bin=[-10,-10]))
    controls={locus:pl.DataFrame(dict(junction_aa=['CAVRRNTGGFKTIF','CAVRDDNTGGFKTIF'],v=[v]*2,j=[j]*2,pgen_bin=[-10]*2))
        for locus,v,j in [('TRA','TRAV1-2*01','TRAJ9*01'),('TRB','TRBV4-1*01','TRBJ1-1*01')]}
    summary,mocks,hits=reference_enrichment(sample,reference,controls,replicates=2,threads=1,
        chain_mode='intersection',tra_matcher='stem4',mismatches=1)
    assert summary['matched_receptors'].item()==1 and summary['matched_count'].item()==5
    assert set(hits['receptor_id'])=={'pair'} and set(hits['locus'])=={'TRA','TRB'}
    assert hits.filter(pl.col('locus')=='TRA')['chain_score'].min()>0
    assert mocks['mock_receptors'].sum()==0
    parallel=reference_enrichment(sample,reference,controls,replicates=2,threads=2,chain_mode='intersection',tra_matcher='stem4',mismatches=1)
    assert all(a.sort(a.columns).equals(b.sort(b.columns)) for a,b in zip((summary,mocks,hits),parallel))
    combined=reference_enrichment(sample,reference,controls,replicates=2,threads=1,chain_mode='all',tra_matcher='stem4')
    for actual,expected in zip(combined,(summary,mocks,hits)):
        selected=actual.filter(pl.col('chain_mode')=='intersection').drop('chain_mode').select(expected.columns)
        assert selected.sort(expected.columns).equals(expected.sort(expected.columns))
    exact=reference_enrichment(sample,reference,controls,replicates=2,threads=1,
        chain_mode='intersection',tra_matcher='stem4',mismatches=0)[0]
    assert exact['matched_receptors'].item()==0
