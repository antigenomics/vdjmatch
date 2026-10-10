"""PSSM density reproduces the original formula without collapsing query identity."""
import math
from pathlib import Path
import polars as pl
from seqtree import Index
from vdjmatch.api import Annotator, _prepare
from vdjmatch.match.historical import density, germline_prior, unpaired_score, paired_score
from vdjmatch.match.regions import significance_pssm
from seqtree import SearchParams


def test_historical_density_formula_and_duplicate_queries():
    q='CASSLGQAYEQYF'; r='CASSLGRAYEQYF'
    ann=Annotator.from_frame(pl.DataFrame({'gene':['TRB','TRB'],'cdr3':[r,r],
        'v':['TRBV19','TRBV19'],'epitope':['E','E']}))
    _,queries=_prepare(pl.DataFrame({'query_id':['positive','negative'],'cdr3':[q,q],
                                   'v':['TRBV19','TRBV19'],'locus':['TRB','TRB']}))
    p=SearchParams(max_subs=5,max_ins=0,max_dels=0,max_total_edits=5,engine='seqtm')
    p.pos_matrix=significance_pssm(len(q))
    score=Index.build([r]).search(q,p)[0].score
    got=density(queries,ann._index.records_for('TRB'),Index.build([q,r]),threads=1)
    expected=math.exp(-score/400)/.5
    assert got['query_id'].to_list()==['positive','negative']
    assert all(math.isclose(x,expected) for x in got['score'])
    assert got['n_reference'].to_list()==[1,1]
    assert all(math.isclose(x,1-math.exp(-.5),rel_tol=1e-14) for x in got['p_enrichment'])


def test_target_first_controls_equal_full_cdf_and_skip_unsupported_queries():
    from vdjmatch.match.historical import control_histograms
    q='CASSLGQAYEQYF';r='CASSLGRAYEQYF';unsupported='CWWWWWWWWWWWF'
    reference=Annotator.from_frame(pl.DataFrame({'gene':['TRB'],'cdr3':[r],
        'v':['TRBV19'],'epitope':['E']}))._index.records_for('TRB')
    _,queries=_prepare(pl.DataFrame({'cdr3':[q,unsupported],'v':['TRBV19']*2,'locus':['TRB']*2}))
    index=Index.build([q,r,unsupported])
    full=control_histograms(index,[q,unsupported])
    class ObservedControl:
        calls=[]
        def __len__(self):return len(index)
        def edit_histogram_batch(self,sequences,params,threads,exclude_exact):
            self.calls.append((list(sequences),params.max_total_edits))
            return index.edit_histogram_batch(sequences,params,threads,exclude_exact)
    control=ObservedControl()
    narrow=density(queries,reference,control,pool_reference=True)
    expected=density(queries,reference,index,pool_reference=True,control_counts=full)
    assert narrow.equals(expected)
    assert control.calls==[([q],1)]


def test_cli_rejects_missing_target_before_control_loading(tmp_path, capsys):
    import pytest
    from vdjmatch.cli.__main__ import main
    q=tmp_path/'q.tsv';r=tmp_path/'r.tsv'
    pl.DataFrame({'junction_aa':['CASSLGQAYEQYF'],'locus':['TRB']}).write_csv(q,separator='\t')
    pl.DataFrame({'cdr3':['CASSLGRAYEQYF'],'gene':['TRB'],'species':['HomoSapiens'],
                  'epitope':['E'],'mhc_a':['HLA-A*02:01']}).write_csv(r,separator='\t')
    with pytest.raises(SystemExit) as error:
        main(['historical-density',str(q),'--vdjdb',str(r),'--locus','TRB',
              '--epitope','E','--epitope','missing','--mhc-a','HLA-A*02:01',
              '--control',str(tmp_path/'nonexistent'),'--output-prefix',str(tmp_path/'out')])
    assert error.value.code==2
    assert 'absent under the selected reference restriction: missing' in capsys.readouterr().err


def test_reference_pooling_is_explicit_and_preserves_default():
    q, r = 'CASSLGQAYEQYF', 'CASSLGRAYEQYF'
    ann = Annotator.from_frame(pl.DataFrame({'gene': ['TRB', 'TRB'], 'cdr3': [r, r],
        'v': ['TRBV19', 'TRBV19'], 'epitope': ['E', 'E'],
        'mhc_a': ['HLA-A*02', 'HLA-A*02:01']}))
    _, queries = _prepare(pl.DataFrame({'cdr3': [q], 'v': ['TRBV19'], 'locus': ['TRB']}))
    reference = ann._index.records_for('TRB')
    control = Index.build([q, r])
    separate = density(queries, reference, control)
    pooled = density(queries, reference, control, pool_reference=True)
    assert separate.height == 2 and pooled.height == 1
    assert pooled['mhc_a'].to_list() == [None]
    assert pooled['n_reference'].to_list() == [1]
    assert all(math.isclose(s, pooled['score'][0]) for s in separate['score'])
    import pytest
    with pytest.raises(ValueError, match='one epitope'):
        density(queries, reference.with_columns(pl.Series('epitope', ['E', 'F'])),
                control, pool_reference=True)


def test_original_prior_counts_raw_background_and_reference_representatives():
    queries = pl.DataFrame({'query_id': ['positive', 'negative'], 'cdr3': ['CAAF', 'CAAF'],
                            'v': ['TRBV1*02', 'TRBV1'], 'j': ['TRBJ1', 'TRBJ1']})
    reference = pl.DataFrame({'cdr3': ['CAGF', 'CAGF'], 'v': ['TRBV1', 'TRBV2'],
                              'j': ['TRBJ1', 'TRBJ2']})
    background = pl.DataFrame({'cdr3': ['CATF', 'CATF', 'CAAF'],
                               'v': ['TRBV2', 'TRBV2', 'TRBV1'],
                               'j': ['TRBJ2', 'TRBJ2', 'TRBJ1']})
    got = germline_prior(queries, reference, background)
    expected = 2 * math.log((1.5 / 31) / (1.5 / 33)) + math.log((1.5 / 31) / (3.5 / 33))
    assert got['query_id'].to_list() == ['positive', 'negative']
    assert got['n_reference'].to_list() == [1, 1]
    assert got['background_size'].to_list() == [3, 3]
    assert all(math.isclose(x, expected) for x in got['germline_lr'])
    # Alpha's original denominator is deliberately the same beta raw background.
    alpha_q = queries.with_columns(pl.lit('TRAV1').alias('v'), pl.lit('TRAJ1').alias('j'))
    alpha_r = reference.with_columns(pl.lit('TRAV1').alias('v'), pl.lit('TRAJ1').alias('j'))
    alpha = germline_prior(alpha_q, alpha_r, background, include_length=False)
    assert all(math.isclose(x, 2 * math.log((1.5 / 31) / (.5 / 33))) for x in alpha['germline_lr'])


def test_original_sparse_boundary_and_paired_average_rank_fusion():
    sparse = unpaired_score([0, 2], [1, -1], 499)
    assert math.isclose(sparse[0], math.log(1e-6) + 1)
    assert math.isclose(sparse[1], math.log(2 + 1e-6) - 1)
    assert unpaired_score([0, 2], [1, -1], 500).to_list() == [0, 2]
    assert paired_score([0, 0, 2], [2, 0, 0], [0, 1, 2]).to_list() == [2.5, 2, 4.5]
    import pytest
    with pytest.raises(ValueError, match='finite'):
        paired_score([0], [float('nan')], [0])
    with pytest.raises(ValueError, match='aligned'):
        paired_score([0], [0, 1], [0])


def test_cli_assayed_targets_preserve_rows_restrictions_and_statuses(tmp_path):
    from vdjmatch.cli.__main__ import main
    q, r = 'CASSLGQAYEQYF', 'CASSLGRAYEQYF'
    sample, reference, control, prefix = [tmp_path / n for n in ('q.tsv', 'r.tsv', 'c.tsv', 'out')]
    pl.DataFrame({'query_id': ['pos', 'neg', 'missing_ref', 'bad', 'empty'],
        'junction_aa': [q, q, q, 'CAXF', None], 'locus': ['TRB']*5,
        'v_call': ['TRBV19']*5, 'j_call': ['TRBJ1-1']*5,
        'species': ['human']*5, 'epitope': ['E', 'E', 'absent', 'E', 'E'],
        'mhc_a': ['HLA-A*02']*5, 'mhc_b': ['B2M']*5, 'mhc_class': ['MHCI']*5,
        'binder': [1, 0, 1, 0, 0]}).write_csv(sample, separator='\t')
    pl.DataFrame({'cdr3': [r, r, 'CASSLGQSYEQYF'], 'gene': ['TRB']*3,
        'v': ['TRBV19']*3, 'j': ['TRBJ1-1']*3, 'species': ['HomoSapiens']*3,
        'epitope': ['E']*3, 'mhc_a': ['HLA-A*02', 'HLA-A*02:01', 'HLA-A*02:01'],
        'mhc_b': ['B2M', 'B2M', 'other'], 'mhc_class': ['MHCI']*3}).write_csv(reference, separator='\t')
    pl.DataFrame({'junction_aa': [q, r]}).write_csv(control, separator='\t')
    assert main(['historical-density', str(sample), '--vdjdb', str(reference),
                 '--locus', 'TRB', '--species', 'human', '--targets-from-sample',
                 '--control', str(control), '--output-prefix', str(prefix)]) == 0
    out = pl.read_csv(str(prefix)+'.scores.tsv', separator='\t')
    assert out['query_id'].to_list() == ['pos', 'neg', 'missing_ref', 'bad', 'empty']
    assert out['binder'].to_list() == [1, 0, 1, 0, 0]
    assert out['status'].to_list() == ['scored', 'scored', 'no_reference', 'invalid_query', 'invalid_query']
    assert out['n_reference'].to_list()[:2] == [1, 1]
    assert out['mhc_a'].to_list() == ['HLA-A*02']*5
    assert math.isclose(out['score'][0], out['score'][1])


def test_cli_sample_targets_reject_explicit_selector(tmp_path, capsys):
    import pytest
    from vdjmatch.cli.__main__ import main
    with pytest.raises(SystemExit) as error:
        main(['historical-density', str(tmp_path/'absent.tsv'), '--vdjdb', str(tmp_path/'r.tsv'),
              '--locus', 'TRB', '--targets-from-sample', '--epitope', 'E',
              '--output-prefix', str(tmp_path/'out')])
    assert error.value.code == 2
    assert 'mutually exclusive' in capsys.readouterr().err


def test_cli_sample_targets_all_invalid_keep_score_schema(tmp_path):
    from vdjmatch.cli.__main__ import main
    sample, reference, control, prefix = [tmp_path / n for n in ('q.tsv', 'r.tsv', 'c.tsv', 'out')]
    pl.DataFrame({'query_id': ['bad'], 'junction_aa': ['CAXF'], 'locus': [None],
        'v_call': ['TRBV19'], 'j_call': ['TRBJ1-1'], 'species': ['human'],
        'epitope': ['E'], 'mhc_a': ['HLA-A*02'], 'mhc_b': ['B2M'],
        'mhc_class': ['MHCI']}).write_csv(sample, separator='\t')
    pl.DataFrame({'cdr3': ['CASSLGRAYEQYF'], 'gene': ['TRB'], 'v': ['TRBV19'],
        'species': ['HomoSapiens'], 'epitope': ['E'], 'mhc_a': ['HLA-A*02'],
        'mhc_b': ['B2M'], 'mhc_class': ['MHCI']}).write_csv(reference, separator='\t')
    pl.DataFrame({'junction_aa': ['CASSLGRAYEQYF']}).write_csv(control, separator='\t')
    assert main(['historical-density', str(sample), '--vdjdb', str(reference),
                 '--locus', 'TRB', '--targets-from-sample', '--control', str(control),
                 '--output-prefix', str(prefix)]) == 0
    out = pl.read_csv(str(prefix)+'.scores.tsv', separator='\t')
    assert out['query_id'].to_list() == ['bad']
    assert out['status'].to_list() == ['invalid_query']
    assert out['score'].to_list() == [None]


def _paired_cli_inputs(tmp_path):
    import polars as pl
    alpha, alpha_mut = 'CAVRDSNYQLIW', 'CAVRDTNYQLIW'
    beta, beta_mut = 'CASSLGQAYEQYF', 'CASSLGRAYEQYF'
    sample, reference, control, alpha_control, prior, prefix = [tmp_path/n for n in
        ('pairs.tsv','r.tsv','c.tsv','ca.tsv','prior.tsv','out')]
    raw=pl.DataFrame({'query_id':['a1','b1','a2','b2','a3','a4','b4'],
        'clone_id':['one','one','two','two','missing','invalid','invalid'],
        'junction_aa':[alpha,beta_mut,alpha_mut,beta,alpha,'CAXF',beta_mut],
        'locus':['TRA','TRB','TRA','TRB','TRA','TRA','TRB'],
        'v_call':['TRAV1','TRBV19','TRAV1','TRBV19','TRAV1','TRAV1','TRBV19'],
        'j_call':['TRAJ1','TRBJ1-1','TRAJ1','TRBJ1-1','TRAJ1','TRAJ1','TRBJ1-1'],
        'species':['human']*7,'epitope':['E']*7,'mhc_a':['HLA-A*02']*7,
        'mhc_b':['B2M']*7,'mhc_class':['MHCI']*7,'binder':[1,1,0,0,1,1,1]})
    raw.write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':[alpha,beta], 'gene':['TRA','TRB'], 'v':['TRAV1','TRBV19'],
        'j':['TRAJ1','TRBJ1-1'],'species':['HomoSapiens']*2,'epitope':['E']*2,
        'mhc_a':['HLA-A*02:01']*2,'mhc_b':['B2M']*2,'mhc_class':['MHCI']*2}).write_csv(reference,separator='\t')
    pl.DataFrame({'junction_aa':[beta,beta_mut]}).write_csv(control,separator='\t')
    pl.DataFrame({'junction_aa':[alpha,alpha_mut]}).write_csv(alpha_control,separator='\t')
    pl.DataFrame({'cdr3':[beta,beta_mut,beta_mut],'v':['TRBV19']*3,'j':['TRBJ1-1']*3}).write_csv(prior,separator='\t')
    args=['historical-density',str(sample),'--vdjdb',str(reference),'--locus','paired',
          '--targets-from-sample','--control',str(control),'--alpha-control',str(alpha_control),
          '--germline-background',str(prior),'--output-prefix',str(prefix)]
    return raw,sample,prefix,args


def test_cli_original_paired_rank_ties_independent_exact_exclusion_and_missing(tmp_path):
    from vdjmatch.cli.__main__ import main
    _,_,prefix,args=_paired_cli_inputs(tmp_path)
    assert main(args)==0
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert out['query_id'].to_list()==['one','two','missing','invalid']
    assert out['pair_id'].to_list()==['one','two','missing','invalid']
    assert out['binder'].to_list()==[1,0,1,1]
    assert out['status'].to_list()==['scored','scored','missing_beta','invalid_alpha']
    assert out['query_id_alpha'].to_list()==['a1','a2','a3','a4']
    assert out['query_id_beta'].to_list()==['b1','b2',None,'b4']
    assert out['density_score_alpha'][0]==0 and out['density_score_beta'][0]>0
    assert out['density_score_beta'][1]==0 and out['density_score_alpha'][1]>0
    assert out['p_enrichment_alpha'][0]==1 and out['p_enrichment_beta'][1]==1
    assert out['score'].to_list()==[1.5,1.5,None,None]
    assert out['rank_germline_prior'].to_list()[:2]==[.5,.5]
    assert 'p_enrichment' not in out.columns
    import json
    manifest=json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert 'beta raw background' in manifest['paired_contract']['alpha_prior_background']


def test_cli_original_paired_rejects_linked_metadata_conflict(tmp_path,capsys):
    import pytest
    from vdjmatch.cli.__main__ import main
    raw,sample,_,args=_paired_cli_inputs(tmp_path)
    raw.with_columns(pl.when(pl.col('query_id')=='b1').then(pl.lit('different'))
                     .otherwise(pl.col('epitope')).alias('epitope')).write_csv(sample,separator='\t')
    with pytest.raises(SystemExit) as error:
        main(args)
    assert error.value.code==2
    assert 'conflicting' in capsys.readouterr().err


def test_cli_extension_retains_original_score_and_source_identities(tmp_path):
    import json
    from vdjmatch.cli.__main__ import main
    q='CASSLGQAYEQYF';r=q+'F'
    sample,reference,control,prefix=[tmp_path/n for n in ('q.tsv','r.tsv','c.tsv','out')]
    pl.DataFrame({'query_id':['a'],'junction_aa':[q],'locus':['TRB'],
                  'v_call':['TRBV19'],'j_call':['TRBJ1-1']}).write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':[r],'gene':['TRB'],'species':['HomoSapiens'],
                  'epitope':['E'],'mhc_a':['HLA-A*02:01'],'v':['TRBV19']}).write_csv(reference,separator='\t')
    pl.DataFrame({'junction_aa':[q]}).write_csv(control,separator='\t')
    assert main(['historical-density',str(sample),'--vdjdb',str(reference),'--locus','TRB',
                 '--epitope','E','--mhc-a','HLA-A*02:01','--pool-reference',
                 '--control',str(control),'--gapped-extension','--output-prefix',str(prefix)])==0
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert out['historical_density'].to_list()==[0.0]
    assert out['gapped_density'][0]>0 and out['score'][0]==out['gapped_density'][0]
    info=json.loads(Path(str(prefix)+'.manifest.json').read_text())['gapped_extension']
    assert all(len(info[k])==64 for k in ['source_sha256','distance_source_sha256','model_sha256'])



def test_cli_sample_targets_retains_unselected_species_rows(tmp_path):
    from vdjmatch.cli.__main__ import main
    sample, reference, control, prefix = [tmp_path/n for n in ('q.tsv','r.tsv','c.tsv','out')]
    q, r = 'CASSLGQAYEQYF', 'CASSLGRAYEQYF'
    pl.DataFrame({'query_id':['selected','unselected'], 'junction_aa':[q,q], 'locus':['TRB','TRB'],
        'v_call':['TRBV19','TRBV19'], 'j_call':['TRBJ1-1','TRBJ1-1'], 'species':['human','mouse'],
        'epitope':['E','E'], 'mhc_a':['HLA-A*02','H-2-Kb'], 'mhc_b':['B2M','B2M'],
        'mhc_class':['MHCI','MHCI'], 'binder':[1,0]}).write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':[r], 'gene':['TRB'], 'v':['TRBV19'], 'species':['HomoSapiens'],
        'epitope':['E'], 'mhc_a':['HLA-A*02'], 'mhc_b':['B2M'], 'mhc_class':['MHCI']}).write_csv(reference,separator='\t')
    pl.DataFrame({'junction_aa':[q,r]}).write_csv(control,separator='\t')
    assert main(['historical-density',str(sample),'--vdjdb',str(reference),'--locus','TRB',
        '--targets-from-sample','--control',str(control),'--output-prefix',str(prefix)])==0
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert out['query_id'].to_list()==['selected','unselected']
    assert out['binder'].to_list()==[1,0]
    assert out['status'].to_list()==['scored','unselected_species']
    assert out['score'][1] is None
    assert out['mhc_a'].to_list()==['HLA-A*02','H-2-Kb']


def test_cli_original_paired_retains_unselected_species_and_chain_ids(tmp_path):
    from vdjmatch.cli.__main__ import main
    raw,sample,prefix,args=_paired_cli_inputs(tmp_path)
    mouse=raw.head(2).with_columns(pl.Series('query_id',['ma','mb']),pl.lit('mouse-pair').alias('clone_id'),
                                 pl.lit('mouse').alias('species'),pl.lit('H-2-Kb').alias('mhc_a'))
    pl.concat([raw,mouse]).write_csv(sample,separator='\t')
    assert main(args)==0
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert out['query_id'].to_list()==['one','two','missing','invalid','mouse-pair']
    assert out['status'][-1]=='unselected_species' and out['score'][-1] is None
    assert out['query_id_alpha'][-1]=='ma' and out['query_id_beta'][-1]=='mb'
    assert out['species'][-1]=='mouse' and out['mhc_a'][-1]=='H-2-Kb'


def test_cli_all_invalid_extension_prior_has_complete_null_schema(tmp_path):
    from vdjmatch.cli.__main__ import main
    sample,reference,control,prior,prefix=[tmp_path/n for n in ('q.tsv','r.tsv','c.tsv','prior.tsv','out')]
    q='CASSLGQAYEQYF'
    pl.DataFrame({'query_id':['invalid'], 'junction_aa':['CAXF'], 'locus':['TRB'],
        'v_call':['TRBV19'], 'j_call':['TRBJ1-1'], 'species':['human'], 'epitope':['E'],
        'mhc_a':['HLA-A*02'], 'mhc_b':['B2M'], 'mhc_class':['MHCI']}).write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':[q], 'gene':['TRB'], 'v':['TRBV19'], 'j':['TRBJ1-1'],
        'species':['HomoSapiens'], 'epitope':['E'], 'mhc_a':['HLA-A*02'],
        'mhc_b':['B2M'], 'mhc_class':['MHCI']}).write_csv(reference,separator='\t')
    pl.DataFrame({'junction_aa':[q]}).write_csv(control,separator='\t')
    pl.DataFrame({'cdr3':[q], 'v':['TRBV19'], 'j':['TRBJ1-1']}).write_csv(prior,separator='\t')
    assert main(['historical-density',str(sample),'--vdjdb',str(reference),'--locus','TRB',
        '--targets-from-sample','--control',str(control),'--germline-background',str(prior),
        '--gapped-extension','--output-prefix',str(prefix)])==0
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert out['status'].to_list()==['invalid_query']
    required={'historical_score','historical_density','gapped_density','density_score','germline_lr',
              'availability','gapped_status','M_gap','n_reference_unavailable','score','p_enrichment'}
    assert required<=set(out.columns)
    assert all(out[name].to_list()==[None] for name in required)
