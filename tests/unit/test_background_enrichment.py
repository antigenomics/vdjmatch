"""Same-predicate count enrichment and gene-control provenance contracts."""
import json
import polars as pl
import pytest
from vdjmatch.cli.__main__ import main
from vdjmatch.evalue.control import raw_gene_background


def test_gene_control_population_ignores_j_conflicts_but_preserves_consensus_inputs(tmp_path):
    path=tmp_path/'control.tsv'
    pl.DataFrame({'junction_aa':['CASSLGQAYEQYF']*3,'v_call':['TRBV19']*3,
                  'j_call':['TRBJ1-1','TRBJ2-1','TRBJ1-1']}).write_csv(path,separator='\t')
    keys,report=raw_gene_background(path,include_j=True)
    assert sorted(keys['j'].to_list())==['TRBJ1-1','TRBJ2-1']
    assert report['unique_keys']==report['population_size']==1
    draws,report=raw_gene_background(path,include_j=True,deduplicate=False)
    assert draws.height==report['population_size']==3 and report['unique_keys']==1


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
@pytest.mark.parametrize('weighting',['ball','local-rank'])
def test_cli_counts_and_jeffreys_ratio_use_same_full_distance(tmp_path,background,weighting):
    q='CASSLGQAYEQYF';r=q[:-3]+'W'+q[-2:]
    sample=tmp_path/'q.tsv';reference=tmp_path/'r.tsv';control=tmp_path/'b.tsv';prefix=tmp_path/'out'
    pl.DataFrame({'junction_aa':[q,q,q],'v_call':['TRBV19']*3,'locus':['TRB']*3,
        'species':['human']*3,'epitope':['E','absent','exact-only'],
        'mhc_a':['HLA-A*02']*3,'mhc_b':['B2M']*3,'mhc_class':['MHCI']*3}).write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':[q,r,r,q],'v':['TRBV19']*4,'gene':['TRB']*4,
        'species':['HomoSapiens']*4,'epitope':['E','E','other','exact-only'],
        'mhc_a':['HLA-A*02:01']*4,'mhc_b':['B2M']*4,'mhc_class':['MHCI']*4}).write_csv(reference,separator='\t')
    control_sequences=[q,r,r] if background=='generative' else [q,r]
    pl.DataFrame({'junction_aa':control_sequences,'v_call':['TRBV19']*len(control_sequences)}).write_csv(control,separator='\t')
    args=['neighbours',str(sample),'--vdjdb',str(reference),'--locus','TRB',
          '--targets-from-sample','--junction-ends','trim3','--radius','0','--exclude-exact',
          '--background',background,'--output-prefix',str(prefix)]
    args+=['--neighbour-weighting',weighting,'--rank-signal-prior','.3']
    if background!='vdjdb-other':args+=['--control',str(control)]
    assert main(args)==0
    scores=pl.read_csv(str(prefix)+'.enrichment.tsv',separator='\t')
    e=scores.filter(pl.col('epitope')=='E').row(0,named=True)
    # E has two model keys (q,r), but only r remains in its hit list.
    N = 2
    M = 3 if background == 'generative' else 2
    m = 2 if background == 'generative' else 1
    assert e['n_reference'] == 1 and e['reference_population'] == N
    assert e['n_control'] == m and e['control_population'] == M
    expected = N * (m + .5) / (M + 1)
    assert e['expected_count'] == expected
    assert e['enrichment'] == pytest.approx(1 / expected)
    assert not e['zero_control_hits']
    absent = scores.filter(pl.col('epitope') == 'absent').row(0, named=True)
    exact_only = scores.filter(pl.col('epitope') == 'exact-only').row(0, named=True)
    assert absent['enrichment'] is None and absent['reference_population'] == 0
    assert exact_only['reference_population'] == 1 and exact_only['n_reference'] == 0
    assert exact_only['enrichment'] == 0
    manifest=json.loads((tmp_path/'out.manifest.json').read_text())
    assert manifest['calibration']=='none'
    assert manifest['background_enrichment']['geometry']=='trim3/3'
    if weighting=='local-rank':
        # One near target averages inverse ranks1..m+1; the other model
        # target stays censored, averaging inverse ranks m+1..M+1.
        H = sum(1 / k for k in range(1, M + 2))
        near = sum(1 / k for k in range(1, m + 2)) / (m + 1)
        far = sum(1 / k for k in range(m + 1, M + 2)) / (M - m + 1)
        bf = (M + 1) * (near + (N - 1) * far) / (N * H)
        assert e['rank_bayes_factor'] == pytest.approx(bf)
        assert e['p_rank_bound'] == pytest.approx(min(1, 1 / bf))
        assert e['posterior_rank_signal'] == pytest.approx(.3 * bf / (.7 + .3 * bf))
        assert absent['rank_bayes_factor'] is None
        assert absent['posterior_rank_signal'] == pytest.approx(.3)
        exact_bf = (M + 1) * far / H
        assert exact_only['rank_bayes_factor'] == pytest.approx(exact_bf)
        assert exact_only['posterior_rank_signal'] == pytest.approx(.3 * exact_bf / (.7 + .3 * exact_bf))
        assert manifest['local_rank_model']['prior']==.3


def test_background_options_reject_undefined_joint_and_ignored_controls():
    base=['neighbours','missing','--vdjdb','missing','--output-prefix','out']
    for flags in (['--locus','paired','--targets-from-sample','--background','real','--control','b'],
                  ['--locus','TRB','--control','b'],
                  ['--locus','TRB','--background','real','--targets-from-sample']):
        with pytest.raises(SystemExit) as error:main(base+flags)
        assert error.value.code==2


@pytest.mark.parametrize('weighting',['linear','local-rank'])
@pytest.mark.parametrize('geometry',['uniform','ndn'])
def test_linear_paired_score_keeps_linkage_and_masks_hits_with_fixed_populations(tmp_path,weighting,geometry):
    if geometry=='ndn':
        pytest.importorskip('arda.cdr3fix')
        from seqtree import gapblock
        if 'query_position_weights' not in gapblock._gapblock_count_batch.__doc__:pytest.skip('native per-input geometry unavailable')
    a='CASSLGQAYEQYF';b='CASSLGRAYEQYF'
    sample=tmp_path/'q.tsv';reference=tmp_path/'r.tsv';prefix=tmp_path/'out'
    pl.DataFrame({'query_id':['a','b'],'clone_id':['001']*2,'junction_aa':[a,a],
        'v_call':['TRAV1-1','TRBV19'],'locus':['TRA','TRB'],'species':['human']*2,
        'epitope':['E']*2,'mhc_a':['HLA-A*02']*2,'mhc_b':['B2M']*2,'mhc_class':['MHCI']*2}).write_csv(sample,separator='\t')
    pl.DataFrame({'complex_id':['1','1','2','2','3','3','4','4'],'gene':['TRA','TRB']*4,
        'cdr3':[a,a,a,b,b,a,b,b],'v':['TRAV1-1','TRBV19']*4,'species':['HomoSapiens']*8,
        'epitope':['E']*4+['other']*4,'mhc_a':['HLA-A*02']*8,'mhc_b':['B2M']*8,'mhc_class':['MHCI']*8}).write_csv(reference,separator='\t')
    args=['neighbours',str(sample),'--vdjdb',str(reference),'--locus','paired',
        '--targets-from-sample','--background','vdjdb-other','--neighbour-weighting',weighting,
        '--exclude-exact','--radius','18','--position-weighting',geometry,'--output-prefix',str(prefix)]
    assert main(args)==0
    row=pl.read_csv(str(prefix)+'.enrichment.tsv',separator='\t').row(0,named=True)
    # Exact/exact hit is removed; both full linked model pairs remain.
    assert row['reference_population'] == 2 and row['n_reference'] == 1
    assert row['control_population'] == row['n_control'] == 2
    if weighting=='linear':
        assert row['reference_weight']==pytest.approx(.5)
        assert row['control_weight']==pytest.approx(.5)  # linked distances9 and18
        # N=2, M=2, target weight1/2, control weight1/2.
        assert row['enrichment'] == pytest.approx(.5 / (2 * (.5 + .5) / 3))
    else:
        # Near target inverse-rank mean3/4; censored second target1/3.
        bf = 3 * (3/4 + 1/3) / (2 * (1 + 1/2 + 1/3))
        assert row['rank_bayes_factor'] == pytest.approx(bf)
        assert row['p_rank_bound'] == pytest.approx(min(1, 1 / bf))
        assert row['p_nearest_rank'] == pytest.approx(5/6)
        assert row['posterior_rank_signal'] == pytest.approx(bf / (1 + bf))
    status=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t',infer_schema_length=0)
    assert status['pair_id'].to_list()==['001']


@pytest.mark.parametrize('flags',[['--radius','0','--background','vdjdb-other','--targets-from-sample'],[]])
def test_linear_rejects_zero_radius_or_missing_background(flags):
    with pytest.raises(SystemExit):
        main(['neighbours','missing','--vdjdb','missing','--locus','TRB',
              '--neighbour-weighting','linear','--output-prefix','out',*flags])


@pytest.mark.parametrize('factor',[5,5.1])
def test_matrix_relative_radius_and_bins_preserve_evidence(monkeypatch,tmp_path,factor):
    from vdjmatch.match import geometry
    sample=tmp_path/'sample.tsv';reference=tmp_path/'reference.tsv'
    q='CASSLGQAYEQYF';near='CASSLGQAWAYEQYF';far='CASSLWWWWGQAYEQYF'
    pl.DataFrame({'junction_aa':[q],'v_call':['TRBV19'],'locus':['TRB'],
        'species':['human'],'epitope':['E'],'mhc_a':['HLA-A*02'],
        'mhc_b':['B2M'],'mhc_class':['MHCI']}).write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':[near,far],'v':['TRBV19','TRBV7-9'],'gene':['TRB']*2,
        'species':['HomoSapiens']*2,'epitope':['E','other'],'mhc_a':['HLA-A*02:01']*2,
        'mhc_b':['B2M']*2,'mhc_class':['MHCI']*2}).write_csv(reference,separator='\t')
    args=['neighbours',str(sample),'--vdjdb',str(reference),'--locus','TRB',
        '--targets-from-sample','--background','vdjdb-other','--neighbour-weighting','local-rank',
        '--unbounded-nearest','--nearest-rank','--junction-ends','trim3',
        '--radius-factor',str(factor)]
    assert main(args+['--output-prefix',str(tmp_path/'base')])==0
    matrix=geometry._substitution_costs()
    monkeypatch.setattr(geometry,'_substitution_costs',lambda alphabet='full',matrix_name='clipped':2*matrix)
    assert main(args+['--output-prefix',str(tmp_path/'scaled')])==0
    base=pl.read_csv(tmp_path/'base.enrichment.tsv',separator='\t')
    scaled=pl.read_csv(tmp_path/'scaled.enrichment.tsv',separator='\t')
    fields=['n_reference','n_control','reference_population','control_population',
        'rank_bayes_factor','conditional_rank_bayes_factor']
    assert base.select(fields).equals(scaled.select(fields))
    a=json.loads((tmp_path/'base.manifest.json').read_text())
    b=json.loads((tmp_path/'scaled.manifest.json').read_text())
    assert b['distance_scale']['radius']==2*a['distance_scale']['radius']
    assert b['distance_scale']['bin_width']==2*a['distance_scale']['bin_width']
    assert b['distance_scale']['gap_cost']==2*a['distance_scale']['gap_cost']
