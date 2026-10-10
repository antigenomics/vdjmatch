"""PSSM density reproduces the original formula without collapsing query identity."""
import math
import pytest
from pathlib import Path
import polars as pl
from seqtree import Index
from vdjmatch.api import Annotator, _prepare
from vdjmatch.match.historical import density, germline_prior, unpaired_score, paired_score
from vdjmatch.match.regions import significance_pssm
from seqtree import SearchParams





@pytest.mark.parametrize('name',['score','cdr3a','status','cohort_size','pair_id','species',
    'v_call','junction_aa','binder','germline_lr_alpha','unpaired_score_beta','status_alpha',
    'rank_germline_prior','n_reference_alpha','control_size_beta','estimator_alpha',
    'query_ida','countb','_locus'])
def test_cli_rejects_reserved_cohort_column_before_pair_parsing(tmp_path,monkeypatch,capsys,name):
    from vdjmatch.cli.__main__ import main
    from vdjmatch import io,db
    from vdjmatch.evalue import control
    raw,sample,_,args=_paired_cli_inputs(tmp_path)
    if name not in raw.columns:
        raw=raw.with_columns(pl.lit('A').alias(name))
    raw.write_csv(sample,separator='\t')
    def forbidden(*args,**kwargs):
        raise AssertionError('reserved cohort name must fail before pairing or scoring')
    monkeypatch.setattr(io,'read_cell',forbidden)
    monkeypatch.setattr(db,'load',forbidden)
    monkeypatch.setattr(control,'raw_background',forbidden)
    args[args.index('--vdjdb')+1]=str(tmp_path/'absent-reference.zip')
    with pytest.raises(SystemExit) as error:
        main(args+['--cohort-column',name])
    assert error.value.code==2
    assert '--cohort-column collides with a canonical input, parsed or score column' in capsys.readouterr().err


def test_cli_cohort_sharing_preserves_independent_ranks_and_unpaired_scores(tmp_path,monkeypatch):
    import pytest
    from vdjmatch.cli.__main__ import main
    from vdjmatch.match import historical
    raw,sample,prefix,args=_paired_cli_inputs(tmp_path)
    assert main(args)==0
    default=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert 'unpaired_score_alpha' not in default.columns
    raw=raw.with_columns(pl.lit('A').alias('benchmark_cohort'))
    solo=raw.head(2).with_columns(pl.Series('query_id',['solo-a','solo-b']),
        pl.lit('solo').alias('clone_id'),pl.lit('B').alias('benchmark_cohort'))
    isolated=raw.head(1).with_columns(pl.lit('isolated-a').alias('query_id'),
        pl.lit('isolated').alias('clone_id'),pl.lit('chain-view').alias('benchmark_cohort'))
    combined=pl.concat([raw,solo,isolated]);combined.write_csv(sample,separator='\t')
    calls=[];native=historical.density
    def density(queries,reference,control,**kwargs):
        calls.append(queries.height)
        return native(queries,reference,control,**kwargs)
    monkeypatch.setattr(historical,'density',density)
    assert main(args+['--cohort-column','benchmark_cohort'])==0
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert len(calls)==2 # One shared target union calculation per chain, not per view.
    assert out['query_id'].to_list()==['one','two','missing','invalid','solo','isolated']
    assert out['score'].to_list()==default['score'].to_list()+[0,None]
    assert out['cohort_size'].to_list()==[2,2,2,2,1,0]
    assert out['benchmark_cohort'].to_list()==['A']*4+['B','chain-view']
    assert out['status'][-1]=='missing_beta'
    assert out['status_alpha'][-1]=='scored' and out['status_beta'][-1]=='missing_beta'
    assert out['unpaired_score_alpha'][-1] is not None and out['unpaired_score_beta'][-1] is None
    # The chain component equals the shipping single-chain CLI calculation exactly.
    isolated.write_csv(sample,separator='\t')
    alpha_control=args[args.index('--alpha-control')+1]
    single=['historical-density',str(sample),'--vdjdb',args[args.index('--vdjdb')+1],
        '--locus','TRA','--targets-from-sample','--control',alpha_control,
        '--germline-background',args[args.index('--germline-background')+1],
        '--output-prefix',str(prefix)]
    assert main(single)==0
    alpha=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert alpha['score'][0]==out['unpaired_score_alpha'][-1]
    assert alpha['p_enrichment'][0]==out['p_enrichment_alpha'][-1]
    solo.tail(1).write_csv(sample,separator='\t')
    beta_args=single.copy();beta_args[beta_args.index('--locus')+1]='TRB'
    beta_args[beta_args.index('--control')+1]=args[args.index('--control')+1]
    assert main(beta_args)==0
    beta=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert beta['score'][0]==out['unpaired_score_beta'][-2]
    assert beta['p_enrichment'][0]==out['p_enrichment_beta'][-2]
    # Cohort is linked metadata: conflicts or missing values are errors.
    conflict=combined.with_columns(pl.when(pl.col('query_id')=='b1').then(pl.lit('B'))
                                   .otherwise(pl.col('benchmark_cohort')).alias('benchmark_cohort'))
    conflict.write_csv(sample,separator='\t')
    with pytest.raises(SystemExit):main(args+['--cohort-column','benchmark_cohort'])
    combined.with_columns(pl.when(pl.col('query_id')=='b1').then(pl.lit(None,dtype=pl.String))
                          .otherwise(pl.col('benchmark_cohort')).alias('benchmark_cohort')).write_csv(sample,separator='\t')
    with pytest.raises(SystemExit):main(args+['--cohort-column','benchmark_cohort'])


def test_cli_single_prior_missing_genes_keep_ids_and_density_only_contract(tmp_path):
    from vdjmatch.cli.__main__ import main
    q,r='CASSLGQAYEQYF','CASSLGRAYEQYF'
    sample,reference,control,prior,prefix=[tmp_path/n for n in ('q.tsv','r.tsv','c.tsv','prior.tsv','out')]
    raw=pl.DataFrame({'query_id':['null-v','blank-j','null-j','unknown','valid','invalid'],
        'junction_aa':[q]*5+['CAXF'],'locus':['TRB']*6,
        'v_call':[None,'TRBV19','TRBV19','TRBVunknown','TRBV19',None],
        'j_call':['TRBJ1-1','   ',None,'TRBJunknown','TRBJ1-1','TRBJ1-1'],
        'species':['human']*6,'epitope':['E']*6,'mhc_a':['HLA-A*02']*6,
        'mhc_b':['B2M']*6,'mhc_class':['MHCI']*6,'binder':[0,1,0,1,0,1]})
    raw.write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':[r],'gene':['TRB'],'v':['TRBV19'],'j':['TRBJ1-1'],
        'species':['HomoSapiens'],'epitope':['E'],'mhc_a':['HLA-A*02'],'mhc_b':['B2M'],
        'mhc_class':['MHCI']}).write_csv(reference,separator='\t')
    pl.DataFrame({'junction_aa':[q,r]}).write_csv(control,separator='\t')
    pl.DataFrame({'cdr3':[q],'v':['TRBV19'],'j':['TRBJ1-1']}).write_csv(prior,separator='\t')
    args=['historical-density',str(sample),'--vdjdb',str(reference),'--locus','TRB',
          '--targets-from-sample','--control',str(control),'--output-prefix',str(prefix)]
    assert main(args+['--germline-background',str(prior)])==0
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert out['query_id'].to_list()==raw['query_id'].to_list()
    assert out['binder'].to_list()==raw['binder'].to_list()
    assert out['status'].to_list()==['missing_v_or_j']*3+['scored','scored','invalid_query']
    assert out['germline_lr'][:3].to_list()==out['score'][:3].to_list()==[None]*3
    assert out['density_score'][:5].null_count()==0
    assert out['germline_lr'][3] is not None and out['score'][3] is not None
    assert main(args)==0
    density_only=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert density_only['status'].to_list()==['scored']*5+['invalid_query']
    assert density_only['score'][:5].to_list()==out['density_score'][:5].to_list()
    # Explicit pooled selection has the same missing-gene full-score disposition.
    explicit=[a for a in args if a!='--targets-from-sample']+['--pool-reference','--epitope','E','--mhc-a','HLA-A*02',
                                                              '--germline-background',str(prior)]
    assert main(explicit)==0
    pooled=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert pooled['query_id'].to_list()==raw['query_id'][:5].to_list()
    assert pooled['status'].to_list()==['missing_v_or_j']*3+['scored']*2
    assert pooled['score'][:3].to_list()==[None]*3
    raw.head(3).write_csv(sample,separator='\t')
    assert main(args+['--germline-background',str(prior),'--gapped-extension'])==0
    missing=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert missing['query_id'].to_list()==['null-v','blank-j','null-j']
    assert missing['status'].to_list()==['missing_v_or_j']*3
    for name in ('score','germline_lr','historical_score'):
        assert missing[name].to_list()==[None]*3


def test_cli_paired_missing_genes_excluded_from_prior_and_fusion(tmp_path):
    from vdjmatch.cli.__main__ import main
    raw,sample,prefix,args=_paired_cli_inputs(tmp_path)
    alpha=raw.head(2).with_columns(pl.Series('query_id',['a-null','b-null']),pl.lit('null-alpha').alias('clone_id'),
        pl.when(pl.col('locus')=='TRA').then(pl.lit(None,dtype=pl.String)).otherwise(pl.col('v_call')).alias('v_call'))
    beta=raw.head(2).with_columns(pl.Series('query_id',['a-blank','b-blank']),pl.lit('blank-beta').alias('clone_id'),
        pl.when(pl.col('locus')=='TRB').then(pl.lit('  ')).otherwise(pl.col('j_call')).alias('j_call'))
    raw=raw.with_columns(pl.when(pl.col('query_id').is_in(['a3','a4'])).then(pl.lit(None,dtype=pl.String))
                        .otherwise(pl.col('v_call')).alias('v_call'))
    pl.concat([raw,alpha,beta]).write_csv(sample,separator='\t')
    assert main(args)==0
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert out['query_id'].to_list()==['one','two','missing','invalid','null-alpha','blank-beta']
    assert out['status'].to_list()==['scored','scored','missing_beta','invalid_alpha',
                                    'missing_v_or_j_alpha','missing_v_or_j_beta']
    assert out['score'].to_list()==[1.5,1.5,None,None,None,None]
    assert out['cohort_size'].to_list()==[2]*6
    assert out['query_id_alpha'][-2:].to_list()==['a-null','a-blank']
    assert out['query_id_beta'][-2:].to_list()==['b-null','b-blank']
    assert out['germline_lr_alpha'][-2] is None and out['germline_lr_beta'][-1] is None
    assert out['germline_lr'][-2:].to_list()==[None,None]
    pl.concat([alpha,beta]).write_csv(sample,separator='\t')
    assert main(args)==0
    missing=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert missing.columns==out.columns
    assert missing['status'].to_list()==['missing_v_or_j_alpha','missing_v_or_j_beta']
    assert missing['cohort_size'].to_list()==[0,0]
    assert missing['score'].to_list()==missing['rank_germline_prior'].to_list()==[None,None]


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
    manifest=json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert len(manifest['software']['seqtree_native_sha256'])==64



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
              'availability','gapped_status','M_gap','n_reference_unavailable','score','p_enrichment',
              'gapped_density_same_v_unweighted','gapped_density_cross_v_unweighted',
              'gapped_edges_same_v','gapped_edges_cross_v','gapped_floor_density',
              'gapped_best_total_distance','gapped_best_cdr3_distance','gapped_best_vloop_distance'}
    assert required<=set(out.columns)
    assert all(out[name].to_list()==[None] for name in required)

    # A valid query with no compatible target reference keeps the same null schema.
    pl.DataFrame({'query_id':['no-reference'], 'junction_aa':[q], 'locus':['TRB'],
        'v_call':['TRBV19'], 'j_call':['TRBJ1-1'], 'species':['human'], 'epitope':['other'],
        'mhc_a':['HLA-A*02'], 'mhc_b':['B2M'], 'mhc_class':['MHCI']}).write_csv(sample,separator='\t')
    assert main(['historical-density',str(sample),'--vdjdb',str(reference),'--locus','TRB',
        '--targets-from-sample','--control',str(control),'--germline-background',str(prior),
        '--gapped-extension','--output-prefix',str(prefix)])==0
    empty=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert empty['status'].to_list()==['no_reference']
    assert required<=set(empty.columns)
    assert all(empty[name].to_list()==[None] for name in required)
