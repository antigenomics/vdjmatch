"""Same-predicate count enrichment and gene-control provenance contracts."""
import json
import polars as pl
import pytest
from vdjmatch.cli.__main__ import main
from vdjmatch.evalue.control import raw_gene_background


def test_gene_controls_keep_distinct_v_keys_and_report_unavailable(tmp_path):
    p=tmp_path/'b.tsv'
    pl.DataFrame({'junction_aa':['CASSLGQAYEQYF']*5+['INVALID!'],
        'v_call':['TRBV19','TRBV19*01','TRBV7-9','TRBV19,TRBV7-9',None,'TRBV19']}).write_csv(p,separator='\t')
    keys,report=raw_gene_background(p)
    assert keys.height==2 and report['eligible_rows']==3
    assert report['excluded_gene_rows']==2 and report['excluded_sequence_rows']==1
    assert report['source_rows']==6 and report['unique_keys']==2
    pl.DataFrame({'cdr3_aa':['CASSLGQAYEQYF'],'v_call':['TRBV19']}).write_csv(p,separator='\t')
    with pytest.raises(ValueError,match='junction_aa'):
        raw_gene_background(p)


def test_gene_controls_filter_declared_metadata_before_shared_gene_resolution(tmp_path):
    p=tmp_path/'mixed.tsv'
    pl.DataFrame({'junction_aa':['CASSLGQAYEQYF']*5,'v_call':['TRBV19']*5,
                  'species':['HomoSapiens','mouse','human',None,'unknown'],
                  'locus':['TRB','TRB','TRA','TRB','TRB']}).write_csv(p,separator='\t')
    keys,report=raw_gene_background(p)
    assert keys.height==report['eligible_rows']==1
    assert report['excluded_metadata_rows']==4


@pytest.mark.parametrize('background',['real','generative','vdjdb-other'])
def test_cli_counts_and_jeffreys_ratio_use_same_full_distance(tmp_path,background):
    q='CASSLGQAYEQYF';r=q[:-3]+'W'+q[-2:]
    sample=tmp_path/'q.tsv';reference=tmp_path/'r.tsv';control=tmp_path/'b.tsv';prefix=tmp_path/'out'
    pl.DataFrame({'junction_aa':[q,q,q],'v_call':['TRBV19']*3,'locus':['TRB']*3,
        'species':['human']*3,'epitope':['E','absent','exact-only'],
        'mhc_a':['HLA-A*02']*3,'mhc_b':['B2M']*3,'mhc_class':['MHCI']*3}).write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':[q,r,r,q],'v':['TRBV19']*4,'gene':['TRB']*4,
        'species':['HomoSapiens']*4,'epitope':['E','E','other','exact-only'],
        'mhc_a':['HLA-A*02:01']*4,'mhc_b':['B2M']*4,'mhc_class':['MHCI']*4}).write_csv(reference,separator='\t')
    pl.DataFrame({'junction_aa':[q,r],'v_call':['TRBV19']*2}).write_csv(control,separator='\t')
    args=['tcrdist-neighbours',str(sample),'--vdjdb',str(reference),'--locus','TRB',
          '--targets-from-sample','--junction-ends','trim3','--radius','0','--exclude-exact',
          '--background',background,'--output-prefix',str(prefix)]
    if background!='vdjdb-other':args+=['--control',str(control)]
    assert main(args)==0
    scores=pl.read_csv(str(prefix)+'.enrichment.tsv',separator='\t')
    e=scores.filter(pl.col('epitope')=='E').row(0,named=True)
    assert e['n_reference']==e['reference_population']==1
    assert e['n_control']==e['control_population']==1
    assert e['expected_count']==.75 and e['enrichment']==pytest.approx(4/3)
    assert not e['zero_control_hits']
    unavailable=scores.filter(pl.col('epitope')!='E')
    assert unavailable['enrichment'].null_count()==2
    assert unavailable['reference_population'].to_list()==[0,0]
    manifest=json.loads((tmp_path/'out.manifest.json').read_text())
    assert manifest['calibration']=='none'
    assert manifest['background_enrichment']['geometry']=='trim3/3'


def test_background_options_reject_undefined_joint_and_ignored_controls():
    base=['tcrdist-neighbours','missing','--vdjdb','missing','--output-prefix','out']
    for flags in (['--locus','paired','--targets-from-sample','--background','real','--control','b'],
                  ['--locus','TRB','--control','b'],
                  ['--locus','TRB','--background','real','--targets-from-sample']):
        with pytest.raises(SystemExit) as error:main(base+flags)
        assert error.value.code==2
