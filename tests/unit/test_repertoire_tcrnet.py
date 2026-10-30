"""TCRnet filtering retains original receptor denominators and distinct match evidence."""
import polars as pl


def test_filtered_flexible_matches_keep_pairing_orphans_and_zero_targets():
    from vdjmatch.aggregate.repertoire_tcrnet import TcrnetMatcher
    reference=pl.DataFrame(dict(cdr3=['CYSSF','CYSRF','CWWWF'],gene=['TRB']*3,
        v=['TRBV3']*3,j=['TRBJ1-1']*3,species=['HomoSapiens']*3,
        epitope=['X','X','Y'],mhc_a=['HLA-A*02:01']*3,mhc_b=['B2M']*3,mhc_class=['MHCI']*3))
    controls={'TRB':pl.DataFrame(dict(junction_aa=['CVVVF','CVVWF','CFVVF'])),
              'TRA':pl.DataFrame(dict(junction_aa=['CVVVF','CVVWF','CFVVF']))}
    observations=pl.DataFrame(dict(receptor_id=['p','p','o','n'],locus=['TRA','TRB','TRB','TRB'],
        junction_aa=['CASSF','CASSF','CATSF','CWWWF'],v=['TRAV1','TRBV2','TRBV2','TRBV2'],
        j=['TRAJ1','TRBJ1-1','TRBJ1-1','TRBJ1-1'],count=[5,5,2,9]))
    matcher=TcrnetMatcher(reference,controls,scope='3,0,0,3',threads=1)
    summary,hits,filtering=matcher.match(observations)
    x=summary.filter(pl.col('epitope')=='X').row(0,named=True)
    y=summary.filter(pl.col('epitope')=='Y').row(0,named=True)
    assert (x['n_receptors'],x['total_count'],x['retained_receptors'])==(3,16,2)
    assert (x['matched_receptors'],x['matched_count'])==(2,7)
    assert y['matched_receptors']==0 and y['reference_retained_keys']==0
    assert hits['receptor_id'].unique().sort().to_list()==['o','p']
    assert not any(c.startswith(('p_','q_','significant')) for c in summary.columns)
    assert filtering.filter(pl.col('locus')=='TRA')['retained'].sum()==0
    assert matcher.reference_filtering['retained'].sum()==2
    parallel=TcrnetMatcher(reference,controls,scope='3,0,0,3',threads=2).match(observations)
    assert all(a.sort(a.columns).equals(b.sort(b.columns)) for a,b in zip((summary,hits,filtering),parallel))

    empty=TcrnetMatcher(reference.tail(1),controls,threads=1)
    empty_summary,empty_hits,_=empty.match(observations)
    assert empty_hits.is_empty() and empty_summary['matched_receptors'].item()==0
    assert empty_summary['n_receptors'].item()==3


def test_filtered_mode_cli_exports_counts_and_filter_evidence(tmp_path):
    import subprocess,sys,json
    (tmp_path/'samples.tsv').write_text('sample_id\tpath\nS\tsample.tsv\n')
    (tmp_path/'sample.tsv').write_text('junction_aa\tv_call\tj_call\tlocus\tduplicate_count\nCASSF\tTRBV2\tTRBJ1-1\tTRB\t5\nCATSF\tTRBV2\tTRBJ1-1\tTRB\t2\nCWWWF\tTRBV2\tTRBJ1-1\tTRB\t9\n')
    (tmp_path/'reference.tsv').write_text('cdr3\tgene\tv\tj\tspecies\tepitope\tmhc_a\tmhc_b\tmhc_class\nCYSSF\tTRB\tTRBV3\tTRBJ1-1\tHomoSapiens\tX\tHLA-A*02:01\tB2M\tMHCI\nCYSRF\tTRB\tTRBV3\tTRBJ1-1\tHomoSapiens\tX\tHLA-A*02:01\tB2M\tMHCI\n')
    (tmp_path/'control.tsv').write_text('junction_aa\tv\tj\nCVVVF\tTRBV3\tTRBJ1-1\nCVVWF\tTRBV3\tTRBJ1-1\nCFVVF\tTRBV3\tTRBJ1-1\n')
    command=[sys.executable,'-m','vdjmatch.cli','repertoire-tcrnet',str(tmp_path/'samples.tsv'),
        '--vdjdb',str(tmp_path/'reference.tsv'),'--control',str(tmp_path/'control.tsv'),
        '--threads','1','--scope','3,0,0,3','--output-prefix',str(tmp_path/'result')]
    result=subprocess.run(command,capture_output=True,text=True,timeout=30)
    assert result.returncode==0,result.stdout+result.stderr
    summary=pl.read_csv(tmp_path/'result.summary.tsv',separator='\t')
    assert summary['matched_receptors'].item()==2 and summary['n_receptors'].item()==3
    assert summary['matched_count'].item()==7 and summary['total_count'].item()==16
    assert not any(c.startswith(('p_','q_','significant')) for c in summary.columns)
    assert pl.read_csv(tmp_path/'result.sample-filter.tsv',separator='\t').height==3
    manifest=json.loads((tmp_path/'result.manifest.json').read_text())
    assert manifest['inputs_frozen'] and manifest['source_frozen'] and not manifest['repertoire_p_values']
