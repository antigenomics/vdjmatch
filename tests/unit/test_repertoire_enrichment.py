import numpy as np
import polars as pl
import pytest
from vdjmatch.aggregate.repertoire_enrichment import matched_enrichment


def test_requested_targets_preserve_full_fuzzy_model_and_mock_draws():
    from vdjmatch.aggregate.repertoire_enrichment import reference_enrichment
    sample=pl.DataFrame(dict(receptor_id=['q','outer_ladder'],locus=['TRB']*2,
        junction_aa=['CASSLGRAYEQYF','C'+'W'*18+'F'],
        v=['TRBV19*01']*2,j=['TRBJ2-7*01']*2,count=[3,5]))
    reference=pl.DataFrame(dict(group_id=['X','Y'],locus=['TRB']*2,
        junction_aa=['CASSLGQAYEQYF','CASSLGRAYEQYF'],v=['TRBV19*01']*2,j=['TRBJ2-7*01']*2,pgen_bin=[-10]*2))
    controls={'TRB':reference.drop('group_id')}
    options=dict(replicates=2,trb_matcher='stem4',chain_mode='trb')
    full=reference_enrichment(sample,reference,controls,**options)
    selected=reference_enrichment(sample,reference,controls,target_groups=['X'],**options)
    for all_rows,selected_rows in zip(full,selected):
        assert all_rows.filter(pl.col('group_id')=='X').equals(selected_rows)
    truncated=reference_enrichment(sample,reference.head(1),controls,**options)
    assert selected[2].filter(pl.col('receptor_id')=='q')['chain_score'].item()<truncated[2].filter(pl.col('receptor_id')=='q')['chain_score'].item()
    # Outer ladder support extends beyond the nominal radius and must survive.
    assert selected[2].filter(pl.col('receptor_id')=='outer_ladder')['chain_score'].item()>0


def test_single_chain_source_rows_preserve_clonotype_endpoint():
    from vdjmatch.cli.repertoire_reference import _observations
    from vdjmatch.aggregate.repertoire_enrichment import reference_enrichment
    frame=pl.DataFrame(dict(locus=['TRB']*3,cdr3=['CASSLGQAYEQYF']*3,
        v=['TRBV19*01']*3,j=['TRBJ2-7*01']*3,count=[7,3,2]))
    _,collapsed=_observations(frame,locus='TRB',species='human')
    _,rows=_observations(frame,locus='TRB',species='human',single_chain_unit='source-row')
    assert collapsed.height==1 and rows.height==3
    reference=rows.head(1).select('locus','junction_aa','v','j').with_columns(pl.lit('X').alias('group_id'),pl.lit(-10).alias('pgen_bin'))
    controls={'TRB':reference.drop('group_id').with_columns(pl.lit('CASSWWWWYEQYF').alias('junction_aa'))}
    for observations,n in [(collapsed,1),(rows,3)]:
        summary,_,_=reference_enrichment(observations,reference,controls,replicates=2,mismatches=0)
        assert summary['matched_receptors'].item()==n
        assert summary['n_receptors'].item()==n
        assert summary['matched_count'].item()==12
    with pytest.raises(ValueError,match='single-chain'):
        _observations(frame,locus='paired',species='human',single_chain_unit='source-row')


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
@pytest.mark.parametrize("batch", [False, True])
def test_reference_mode_cli_and_significance_options(tmp_path,supported,chain_mode,batch):
    import subprocess,sys,json
    from pathlib import Path
    sample=tmp_path/'sample.tsv'
    sample.write_text('junction_aa\tv_call\tj_call\tlocus\tduplicate_count\nCASSIIDAYEQYF\tTRBV2\tTRBJ1-1\tTRB\t7\n')
    (tmp_path/'metadata.tsv').write_text('sample_id\tpath\nS\tsample.tsv\n')
    (tmp_path/'reference.tsv').write_text('cdr3\tgene\tv\tj\tspecies\tepitope\tmhc_a\tmhc_b\tmhc_class\nCASSIIDAYEQYF\tTRB\tTRBV2\tTRBJ1-1\tHomoSapiens\tSYNTHETIC\tHLA-A*02:01\tB2M\tMHCI\n')
    (tmp_path/'control.tsv').write_text('junction_aa\tv\tj\toccurrences\nCASSIIDAYEQYF\t'+('TRBV3' if supported else 'TRBV2')+'\tTRBJ1-1\t100\n')
    command=[sys.executable,'-m','vdjmatch.cli','repertoire-enrichment',str(tmp_path/'metadata.tsv'),
        '--null','reference','--vdjdb',str(tmp_path/'reference.tsv'),'--control',str(tmp_path/'control.tsv'),
        '--pgen-mismatches','0','--threads','1','--significance','both','--chain-mode',chain_mode,'--output-prefix',str(tmp_path/'result')]
    if chain_mode=='all':command.extend(['--model-source','isalgo'])
    if batch:command.append('--batch-samples')
    if batch and chain_mode=='all':command.extend(['--trb-matcher','stem4','--trb-min-score','0.50'])
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
        assert summary['q_count'].item()==1/201
        assert summary['significant_nominal_count'].item() and summary['significant_bh_count'].item()
    else:
        assert summary['status'].item()=='untestable_null'
        assert summary['p_receptors'].item() is None and summary['q_receptors'].item() is None
        assert summary['significant_nominal'].item() is None and summary['significant_bh'].item() is None
        assert summary['q_count'].item() is None and summary['significant_bh_count'].item() is None
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


@pytest.mark.parametrize("kernel", ["stem4", "gapblock"])
@pytest.mark.parametrize("fuzzy_locus", ["TRA", "TRB"])
def test_strict_beta_with_fuzzy_alpha_preserves_pairing(monkeypatch, fuzzy_locus, kernel):
    from vdjmatch.runner import combined
    original=combined._kernel;scopes=[]
    def capture(queries,reference,cases,**options):
        scopes.append((queries.height,cases is not None))
        return original(queries,reference,cases,**options)
    monkeypatch.setattr(combined,"_kernel",capture)
    from vdjmatch.aggregate.repertoire_enrichment import reference_enrichment
    sample=pl.DataFrame(dict(receptor_id=['pair','pair','beta_only','alpha_only'],locus=['TRA','TRB','TRB','TRA'],
        junction_aa=['CAVRDSNFQLIW','CASSLGRAYEQYF','CASSLGRAYEQYF','CAVRDSNFQLIW'],
        v=['TRAV12-2*01','TRBV19*01','TRBV19*01','TRAV12-2*01'],j=['TRAJ33*01','TRBJ2-7*01','TRBJ2-7*01','TRAJ33*01'],count=[5,5,3,2]))
    reference=pl.DataFrame(dict(group_id=['X','X'],locus=['TRA','TRB'],junction_aa=['CAVRDSNYQLIW','CASSLGQAYEQYF'],
        v=['TRAV12-2*01','TRBV19*01'],j=['TRAJ33*01','TRBJ2-7*01'],pgen_bin=[-10,-10]))
    controls={locus:pl.DataFrame(dict(junction_aa=['CAVRRNTGGFKTIF','CAVRDDNTGGFKTIF'],v=[v]*2,j=[j]*2,pgen_bin=[-10]*2))
        for locus,v,j in [('TRA','TRAV1-2*01','TRAJ9*01'),('TRB','TRBV4-1*01','TRBJ1-1*01')]}
    options={fuzzy_locus.lower()+'_matcher':kernel}
    if fuzzy_locus=='TRB':
        def swap(frame):
            replacements={'TRAV12-2*01':'TRBV19*01','TRBV19*01':'TRAV12-2*01',
                'TRAJ33*01':'TRBJ2-7*01','TRBJ2-7*01':'TRAJ33*01',
                'TRAV1-2*01':'TRBV4-1*01','TRBV4-1*01':'TRAV1-2*01',
                'TRAJ9*01':'TRBJ1-1*01','TRBJ1-1*01':'TRAJ9*01'}
            frame=frame.with_columns(pl.col('v','j').replace(replacements))
            return frame.with_columns(pl.col('locus').replace({'TRA':'TRB','TRB':'TRA'})) if 'locus' in frame.columns else frame
        sample,reference=swap(sample),swap(reference)
        controls={('TRB' if locus=='TRA' else 'TRA'):swap(frame) for locus,frame in controls.items()}
    summary,mocks,hits=reference_enrichment(sample,reference,controls,replicates=2,threads=1,
        chain_mode='intersection',**options,mismatches=1)
    assert summary['matched_receptors'].item()==1 and summary['matched_count'].item()==5
    assert set(hits['receptor_id'])=={'pair'} and set(hits['locus'])=={'TRA','TRB'}
    assert hits.filter(pl.col('locus')==fuzzy_locus)['chain_score'].min()>0
    assert mocks['mock_receptors'].sum()==0
    assert scopes[0]==(1,True)
    parallel=reference_enrichment(sample,reference,controls,replicates=2,threads=2,chain_mode='intersection',**options,mismatches=1)
    assert all(a.sort(a.columns).equals(b.sort(b.columns)) for a,b in zip((summary,mocks,hits),parallel))
    combined=reference_enrichment(sample,reference,controls,replicates=2,threads=1,chain_mode='all',**options)
    for actual,expected in zip(combined,(summary,mocks,hits)):
        selected=actual.filter(pl.col('chain_mode')=='intersection').drop('chain_mode').select(expected.columns)
        assert selected.sort(expected.columns).equals(expected.sort(expected.columns))
    from polars.testing import assert_frame_equal
    inputs={'A':sample,'B':sample.with_columns(pl.col('count')+7)}
    cohort=pl.concat([part.with_columns(pl.lit(name).alias('sample_id')) for name,part in inputs.items()])
    batched=reference_enrichment(cohort,reference,controls,replicates=2,chain_mode='all',**options)
    for name,part in inputs.items():
        separate=reference_enrichment(part,reference,controls,replicates=2,chain_mode='all',**options)
        for actual,expected in zip(batched,separate):
            selected=actual.filter(pl.col('sample_id')==name).drop('sample_id').select(expected.columns)
            assert_frame_equal(selected.sort(expected.columns),expected.sort(expected.columns),check_dtypes=False)
    exact=reference_enrichment(sample,reference,controls,replicates=2,threads=1,
        chain_mode='intersection',**options,mismatches=0)[0]
    assert exact['matched_receptors'].item()==0
    cutoff=float(hits.filter(pl.col('locus')==fuzzy_locus)['chain_score'].max())
    rejected=reference_enrichment(sample,reference,controls,replicates=2,chain_mode='intersection',**options,**{fuzzy_locus.lower()+'_min_score':cutoff})
    assert rejected[0]['matched_receptors'].item()==0
    assert rejected[2].is_empty()


def test_abundance_enrichment_can_differ_from_receptor_enrichment():
    from vdjmatch.aggregate.repertoire_enrichment import reference_enrichment
    sample=pl.DataFrame(dict(receptor_id=['expanded','small'],locus=['TRB']*2,
        junction_aa=['CASSF','CATSF'],v=['TRBV2']*2,j=['TRBJ1']*2,count=[7,1]))
    reference=pl.DataFrame(dict(group_id=['X'],locus=['TRB'],junction_aa=['CASSF'],
        v=['TRBV2'],j=['TRBJ1'],pgen_bin=[-10]))
    controls={'TRB':pl.DataFrame(dict(junction_aa=['CATSF'],v=['TRBV2'],j=['TRBJ1'],pgen_bin=[-10]))}
    summary,_,_=reference_enrichment(sample,reference,controls,replicates=200,mismatches=0)
    assert summary['matched_receptors'].item()==1 and summary['mock_mean_receptors'].item()==1
    assert summary['p_receptors'].item()==1
    assert summary['matched_count'].item()==7 and summary['mock_mean_count'].item()==1
    assert summary['p_count'].item()==1/201

@pytest.mark.parametrize('cutoff', [-1., float('nan'), float('inf')])
def test_fuzzy_cutoff_rejected_before_cli_input_loading(cutoff):
    from types import SimpleNamespace
    from vdjmatch.cli.repertoire_enrichment import main
    with pytest.raises(ValueError, match='finite nonnegative'):
        main(SimpleNamespace(tra_min_score=cutoff, tra_matcher='stem4'))


@pytest.mark.parametrize('matcher,bound', [('hamming',None)])
def test_bulk_panel_matchers_count_each_source_row_once(matcher,bound):
    from vdjmatch.aggregate.repertoire_enrichment import exact_reference_counts
    obs=pl.DataFrame(dict(junction_aa=['CASSIIDAYEQYF','CASSIIDAYEQYF','CASSMMDAYEQYF','CASSIIRAYEQYF','CASSIIDAYEQYF'],
        v=['V']*5,j=['J']*4+['wrong'],count=[2,3,7,11,100]))
    panels=pl.DataFrame(dict(group_id=['X','X','X','none'],replicate=[-1,-1,0,0],
        junction_aa=['CASSIIDAYEQYF','CASSMIDAYEQYF','CASSIIDAYEQYF','CWWWWWWWWWWWF'],v=['V']*4,j=['J']*4))
    result=exact_reference_counts(obs,panels,matcher=matcher)
    assert result.select('matched_records','matched_reads').rows()==(
        [(4,23),(3,16),(0,0)] if matcher=='hamming' else [(3,12),(3,12),(0,0)])
    parallel=exact_reference_counts(obs,panels,matcher=matcher,threads=2)
    assert result.equals(parallel)
    empty=exact_reference_counts(obs.head(0),panels,matcher=matcher)
    assert empty['matched_records'].sum()==empty['matched_reads'].sum()==0


@pytest.mark.parametrize('kernel',['stem4','gapblock'])
@pytest.mark.parametrize('targets',[None,['X']])
def test_fuzzy_same_v_reducer_preserves_model_and_excludes_cross_v_support(kernel,targets):
    from vdjmatch.aggregate.repertoire_enrichment import reference_enrichment
    obs=pl.DataFrame(dict(receptor_id=['same','different'],locus=['TRB']*2,
        junction_aa=['CASSIIDAYEQYF']*2,v=['TRBV19','TRBV7-9'],j=['TRBJ2-7']*2,count=[3,7]))
    ref=obs.head(1).select('locus','junction_aa','v','j').with_columns(pl.lit('X').alias('group_id'),pl.lit(-10).alias('pgen_bin'))
    control=ref.drop('group_id').with_columns(pl.lit('CASSMMDAYEQYF').alias('junction_aa'))
    options=dict(replicates=2,trb_matcher=kernel,chain_mode='trb',target_groups=targets)
    soft=reference_enrichment(obs,ref,{'TRB':control},**options)
    strict=reference_enrichment(obs,ref,{'TRB':control},kernel_scorer='vonly-bounded',**options)
    assert set(soft[2]['receptor_id'])=={'same','different'}
    assert strict[2]['receptor_id'].to_list()==['same']
    assert strict[0]['matched_count'].item()==3
    assert strict[0]['n_receptors'].item()==2


def test_fixed_kernel_panel_counts_keep_cohort_rows_and_other_peptide_context(monkeypatch):
    from vdjmatch.aggregate.repertoire_enrichment import kernel_reference_counts
    obs=pl.DataFrame(dict(sample_id=['A','A','A','B','B'],junction_aa=['CASSIIDAYEQYF','CASSIIDAYEQYF','CASSMMDAYEQYF','CASSMMDAYEQYF','CWWWWWWWWWWWF'],v=['TRBV19']*5,j=['TRBJ2-7']*5,count=[2,3,7,11,13]))
    panels=pl.DataFrame(dict(group_id=['X','X'],replicate=[-1,0],junction_aa=['CASSIIDAYEQYF','CASSMMDAYEQYF'],v=['TRBV19']*2,j=['TRBJ2-7']*2))
    inventory=pl.DataFrame(dict(group_id=['X'],epitope=['LLWNGPMAV'],mhc_a=['HLA-A*02:01'],mhc_b=['B2M'],mhc_class=['MHCI']))
    full=pl.DataFrame(dict(cdr3=['CASSIIDAYEQYF','CWWWWWWWWWWWF'],v=['TRBV19']*2,j=['TRBJ2-7']*2,locus=['TRB']*2,epitope=['LLWNGPMAV','OTHER'],mhc_a=['HLA-A*02:01']*2,mhc_b=['B2M']*2,mhc_class=['MHCI']*2))
    far='C'+'W'*70+'F'
    obs=pl.concat([obs,obs.head(1).with_columns(pl.lit('B').alias('sample_id'),pl.lit(far).alias('junction_aa'),pl.lit(17,dtype=pl.Int64).alias('count'))])
    import vdjmatch.runner.combined as combined
    original_kernel=combined._kernel;seen=[]
    def capture(queries,*args,**kwargs):
        seen.extend(queries['cdr3'].to_list())
        return original_kernel(queries,*args,**kwargs)
    monkeypatch.setattr(combined,'_kernel',capture)
    original=full.clone()
    result=kernel_reference_counts(obs,panels,inventory,full,locus='TRB',species='human',min_score=.95)
    assert full.equals(original)
    assert result.select('sample_id','replicate','matched_records','matched_reads').rows()==[('A',-1,2,5),('A',0,1,7),('B',-1,0,0),('B',0,1,11)]
    for sample in ['A','B']:
        separate=kernel_reference_counts(obs.filter(pl.col('sample_id')==sample),panels,inventory,full,locus='TRB',species='human',min_score=.95)
        assert separate.equals(result.filter(pl.col('sample_id')==sample))
    assert seen and far not in seen


@pytest.mark.parametrize('anchored',[False,True])
@pytest.mark.parametrize('target_groups',[None,['X']])
def test_fuzzy_preserves_peptide_background_across_mhc_restrictions(anchored,target_groups):
    from vdjmatch.aggregate.repertoire_enrichment import reference_enrichment
    from vdjmatch.runner.combined import _kernel
    from vdjmatch.match.geometry import precompute_v_loop_distances
    from vdjmatch.aggregate.candidates import PMHC
    sample=pl.DataFrame(dict(receptor_id=['q'],locus=['TRB'],junction_aa=['CASSLGRAYEQYF'],v=['TRBV19*01'],j=['TRBJ2-7*01'],count=[3]))
    reference=pl.DataFrame(dict(group_id=['X','Y','Z'],locus=['TRB']*3,
        junction_aa=['CASSLGQAYEQYF','CASSLGRAYEQYF','CASSWWWWYEQYF'],v=['TRBV19*01']*3,j=['TRBJ2-7*01']*3,pgen_bin=[-10]*3,
        epitope=['SYNTHETIC','SYNTHETIC','OTHER'],mhc_a=['HLA-A*02:01','HLA-A*03:01','HLA-A*02:01'],mhc_b=['B2M']*3,mhc_class=['MHCI']*3))
    reference=reference.with_columns(pl.lit(None).alias('mhc_b'))
    control=reference.head(1).select('locus','junction_aa','v','j','pgen_bin').with_columns(pl.lit('CASSYYYYYEQYF').alias('junction_aa'))
    query=sample.rename({'junction_aa':'cdr3'}).with_row_index('query_id')
    cases=query.select('query_id').join(reference.head(1).select(PMHC),how='cross')
    expected=_kernel(query,reference.rename({'junction_aa':'cdr3'}),cases,species='human',threads=1,paired=False,
        exclude_exact=False,contexts={'blosum62':precompute_v_loop_distances(species='human',alphabet='full',matrix='blosum62',v_regions='cdr12plus')},kernels=('stem4',))['stem4_score'].item()
    controls={'TRB':control}
    if anchored:
        alpha=lambda frame:frame.with_columns(pl.lit('TRA').alias('locus'),pl.lit('CAVRDSNYQLIW').alias('junction_aa'),pl.lit('TRAV1-1*01').alias('v'),pl.lit('TRAJ33*01').alias('j'))
        sample=pl.concat([sample,alpha(sample)])
        reference=pl.concat([reference,alpha(reference.head(1))])
        controls['TRA']=alpha(control).with_columns(pl.lit('CAVRDSNYQLIF').alias('junction_aa'))
    _,_,hits=reference_enrichment(sample,reference,controls,replicates=2,trb_matcher='stem4',chain_mode='intersection' if anchored else 'trb',target_groups=target_groups)
    scored=hits.filter((pl.col('group_id')=='X')&(pl.col('locus')=='TRB'))
    assert scored.height==1
    assert scored['chain_score'].item()==pytest.approx(expected,rel=0,abs=1e-15)


def test_fuzzy_keeps_reference_rows_requiring_production_v_geometry_imputation(monkeypatch):
    from vdjmatch.aggregate.repertoire_enrichment import reference_enrichment
    import vdjmatch.evalue.kernel as kernel
    original=kernel.grouped_kernel_enrichment
    observed=[]
    def capture(queries,references,*args,**kwargs):
        observed.append(references.height)
        return original(queries,references,*args,**kwargs)
    monkeypatch.setattr(kernel,'grouped_kernel_enrichment',capture)
    sample=pl.DataFrame(dict(receptor_id=['q'],locus=['TRA'],junction_aa=['CAVRDSNYQLIW'],v=['TRAV1-1'],j=['TRAJ33'],count=[1]))
    reference=pl.DataFrame(dict(group_id=['X','Y'],locus=['TRA']*2,junction_aa=['CAVRDSNYQLIW','CAVRDSNFQLIW'],v=['TRAV1-1','TRAV14-1'],j=['TRAJ33']*2,pgen_bin=[-10]*2))
    control=reference.drop('group_id').with_columns(pl.Series('junction_aa',['CAVRDSNQQLIW','CAVRDSNHQLIW']))
    reference_enrichment(sample,reference,{'TRA':control},replicates=2,tra_matcher='stem4',chain_mode='tra',target_groups=['X'])
    assert observed[0]==reference.height


@pytest.mark.parametrize('matcher,factor',[('blosum62',.55),('positional-blosum62',.77)])
def test_dimensionless_blosum_radius_recovers_low_cost_multiple_substitutions(matcher,factor):
    from vdjmatch.aggregate.repertoire_enrichment import exact_reference_counts
    obs=pl.DataFrame(dict(junction_aa=['CASSIIDAYEQYF','CASSMMDAYEQYF','CASSIIRAYEQYF'],v=['TRBV19']*3,j=['TRBJ2-7']*3,count=[2,3,5]))
    panels=obs.head(1).select('junction_aa','v','j').with_columns(pl.lit('X').alias('group_id'),pl.lit(-1).alias('replicate'))
    result=exact_reference_counts(obs,panels,matcher=matcher,radius_factor=factor)
    assert result.select('matched_records','matched_reads').row(0)==(2,5)
    strict=exact_reference_counts(obs,panels,matcher='hamming')
    assert strict.select('matched_records','matched_reads').row(0)==(2,7)
    exact=exact_reference_counts(obs,panels,matcher=matcher,radius_factor=.0001)
    assert exact.select('matched_records','matched_reads').row(0)==(1,2)


@pytest.mark.parametrize('matcher,factor',[('blosum62',.55),('positional-blosum62',.77)])
@pytest.mark.parametrize('locus,v,j',[('TRA','TRAV12-2','TRAJ33'),('TRB','TRBV19','TRBJ2-7')])
def test_dimensionless_blosum_reference_matching_is_symmetric(matcher,factor,locus,v,j):
    from vdjmatch.aggregate.repertoire_enrichment import reference_enrichment
    sample=pl.DataFrame(dict(receptor_id=['exact','two_low','one_high'],locus=[locus]*3,junction_aa=['CASSIIDAYEQYF','CASSMMDAYEQYF','CASSIIRAYEQYF'],v=[v]*3,j=[j]*3,count=[2,3,5]))
    reference=sample.head(1).select('locus','junction_aa','v','j').with_columns(pl.lit('X').alias('group_id'),pl.lit(-10).alias('pgen_bin'))
    controls={locus:reference.drop('group_id').with_columns(pl.lit('CASSMMDAYEQYF').alias('junction_aa'))}
    options={locus.lower()+'_matcher':matcher,locus.lower()+'_radius_factor':factor}
    summary,_,hits=reference_enrichment(sample,reference,controls,replicates=2,**options)
    assert sorted(hits['receptor_id'])==['exact','two_low']
    assert summary['matched_count'].item()==5


@pytest.mark.parametrize('matcher',['blosum62','positional-blosum62'])
@pytest.mark.parametrize('factor',[None,0.,-1.,float('nan'),float('inf'),True])
def test_dimensionless_blosum_rejects_invalid_factors(matcher,factor):
    from vdjmatch.match.scope import blosum_radius_params
    with pytest.raises(ValueError,match='dimensionless radius factor'):
        blosum_radius_params(matcher,factor,12)


@pytest.mark.parametrize('matcher',['blosum62','positional-blosum62'])
def test_dimensionless_blosum_length_batches_preserve_identity(matcher):
    from vdjmatch.aggregate.repertoire_enrichment import exact_reference_counts
    obs=pl.DataFrame(dict(junction_aa=['CASSIIDAYEQYF','CASSIIDAYEQYFF'],v=['TRBV19']*2,j=['TRBJ2-7']*2,count=[2,3]))
    panels=obs.select('junction_aa','v','j').with_columns(pl.Series('group_id',['short','long']),pl.lit(-1).alias('replicate'))
    result=exact_reference_counts(obs,panels,matcher=matcher,radius_factor=.10)
    assert result.select('group_id','matched_records','matched_reads').rows()==[('long',1,3),('short',1,2)]
