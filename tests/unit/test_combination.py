"""Global evidence formula convergence, monotonicity and missing dispositions."""
import numpy as np
import polars as pl
import pytest
from vdjmatch.aggregate.combination import fit_evidence_combination,combine_evidence,fixed_evidence_combination


def test_combination_fit_and_frozen_prediction():
    rng=np.random.default_rng(14);n=300
    l=rng.normal(size=n);r=rng.exponential(size=n);m=rng.uniform(size=n)
    q=pl.DataFrame({'likelihood_score':l,'conditional_rank_score':r,'motif_score':m,
        'reference_coverage':rng.uniform(size=n),'knee_shape':rng.uniform(size=n),
        'knee_radius_ratio':rng.uniform(size=n),'original_support':rng.uniform(size=n),'paired':np.arange(n)%2==0})
    y=(l+np.log1p(r)+m+rng.normal(scale=.7,size=n)>1).astype(int)
    model=fit_evidence_combination(q,y)
    score=combine_evidence(q,model)
    assert np.isfinite(score).all()
    assert score[y==1].mean()>score[y==0].mean()
    for c in ('likelihood_score','conditional_rank_score'):
        stronger=q.with_columns((pl.col(c)+1).alias(c))
        assert np.all(combine_evidence(stronger,model)>=score-1e-12)
    unsupported=q.head(1).with_columns(pl.lit(None,dtype=pl.Float64).alias(c)
        for c in ('likelihood_score','conditional_rank_score'))
    assert np.isnan(combine_evidence(unsupported,model)[0])
    bad=q.with_columns(pl.lit(2).alias('knee_shape'))
    with pytest.raises(ValueError,match=r'\[0,1\]'):combine_evidence(bad,model)


def test_cli_frozen_groups_and_unsupported_rows(tmp_path):
    from vdjmatch.cli.__main__ import main
    train=pl.DataFrame({'likelihood_score':[-2.,-1.,1.,2.],
        'conditional_rank_score':[.1,.2,1.,2.], 'motif_score':[0.,0.,1.,1.],
        'reference_coverage':[.5]*4,'knee_shape':[.2]*4,'knee_radius_ratio':[.2]*4,
        'original_support':[1.]*4,'paired':[False]*4,
        'epitope':['TRAIN']*4,'label':[0,0,1,1]})
    source=tmp_path/'train.tsv';model=tmp_path/'model.json'
    train.write_csv(source,separator='\t')
    original=source.read_bytes()
    with pytest.raises(SystemExit):
        main(['fit-evidence',str(source),'--model',str(source)])
    assert source.read_bytes()==original
    assert main(['fit-evidence',str(source),'--model',str(model)])==0
    with pytest.raises(SystemExit) as rejected:
        main(['combine-evidence',str(source),'--model',str(model),'--require-heldout','--output',str(tmp_path/'bad.tsv')])
    assert rejected.value.code==2
    test=train.with_columns(pl.lit('HELDOUT').alias('epitope')).with_columns(
        pl.when(pl.int_range(pl.len())==0).then(None).otherwise(pl.col(c)).alias(c)
        for c in ('likelihood_score','conditional_rank_score'))
    path=tmp_path/'test.tsv';output=tmp_path/'scores.tsv';test.write_csv(path,separator='\t')
    assert main(['combine-evidence',str(path),'--model',str(model),'--require-heldout','--output',str(output)])==0
    scored=pl.read_csv(output,separator='\t')
    assert scored.height==4 and scored['combined_status'][0]=='unsupported'
    assert scored['combined_score'].null_count()==1
    fixed=tmp_path/'fixed.tsv'
    assert main(['combine-evidence',str(path),'--fixed','--output',str(fixed)])==0
    assert pl.read_csv(fixed,separator='\t')['combined_score'].equals(scored['fixed_score'])


def test_fixed_formula_dense_sparse_and_missing():
    frame=pl.DataFrame({'likelihood_score':[2.,2.,None,2.,None],
        'conditional_rank_score':[4.,4.,4.,None,None],
        'reference_coverage':[1.]*5,'knee_shape':[1.]*5,'knee_radius_ratio':[0.]*5,
        'original_support':[1.,0.,1.,0.,0.],'paired':[True]*5})
    actual=fixed_evidence_combination(frame)
    np.testing.assert_allclose(actual[:4],[2.,2.+np.log(4),np.log(4),2.])
    assert np.isnan(actual[-1])
    # At fixed reference structure the original term always retains unit weight.
    stronger=frame.with_columns((pl.col('likelihood_score')+3).alias('likelihood_score'))
    np.testing.assert_allclose(fixed_evidence_combination(stronger)[[0,1,3]]-actual[[0,1,3]],3.)


def test_null_normalized_components_use_same_target_ties_sources_and_nohits(tmp_path):
    from vdjmatch.aggregate.combination import null_normalized_evidence
    from vdjmatch.cli.__main__ import main
    context=dict(species='human',locus='TRB',epitope='E',mhc_a='HLA-A*02:01',mhc_b='B2M',mhc_class='MHCI')
    q=pl.DataFrame([dict(query_id='strong',likelihood_score=100.,conditional_rank_score=1.1,**context),
                    dict(query_id='empty',likelihood_score=None,conditional_rank_score=None,**context)])
    nulls=pl.DataFrame([dict(query_id=str(i),likelihood_score=l,conditional_rank_score=r,
        outcome='null',null_source=source,**context) for i,(l,r,source) in enumerate([
            (50.,1.0,'withheld'),(100.,1.1,'withheld'),(None,None,'withheld'),(None,None,'withheld'),
            (200.,1.2,'generated'),(None,None,'generated')])])
    result=null_normalized_evidence(q,nulls)
    assert result['p_component_likelihood'][0]==pytest.approx(2/3)
    assert result['p_component_rank'][0]==pytest.approx(2/3)
    assert result['normalized_score'][0]==pytest.approx(-np.log(2/3))
    assert result['normalized_score'][1] is None
    assert result['p_component_likelihood'][1]==result['p_component_rank'][1]==1
    # Scaling either component by a strictly increasing affine map is neutral.
    for column in ['likelihood_score','conditional_rank_score']:
        transformed_q=q.with_columns((pl.col(column)*100+3).alias(column))
        transformed_nulls=nulls.with_columns((pl.col(column)*100+3).alias(column))
        assert null_normalized_evidence(transformed_q,transformed_nulls)['normalized_score'].equals(result['normalized_score'])
    with pytest.raises(ValueError,match='calibration group'):
        null_normalized_evidence(q,nulls.with_columns(pl.lit('other').alias('epitope')))
    source=tmp_path/'evidence.tsv';control=tmp_path/'nulls.tsv';output=tmp_path/'scores.tsv'
    q.write_csv(source,separator='\t');nulls.write_csv(control,separator='\t')
    assert main(['combine-evidence',str(source),'--null-evidence',str(control),'--output',str(output)])==0
    saved=pl.read_csv(output,separator='\t')
    assert saved['combined_score'].equals(result['normalized_score'])
    assert saved['combined_status'].to_list()==['scored','unsupported']
    assert saved.select(q.columns).equals(q)

    # A query with no reference label must survive in rejection denominators.
    absent=q.tail(1).with_columns(pl.lit(None,dtype=pl.String).alias(c)
        for c in ('epitope','mhc_a','mhc_b','mhc_class'))
    absent.write_csv(source,separator='\t')
    assert main(['combine-evidence',str(source),'--null-evidence',str(control),'--output',str(output)])==0
    saved=pl.read_csv(output,separator='\t')
    assert saved['query_id'].to_list()==['empty'] and saved['combined_status'].to_list()==['unsupported']

    for invalid in [nulls.with_columns(pl.lit('correct').alias('outcome')),
                    nulls.with_columns(pl.lit(0.).alias('conditional_rank_score')),
                    nulls.with_columns(pl.lit(float('nan')).alias('likelihood_score'))]:
        with pytest.raises(ValueError):null_normalized_evidence(q,invalid)


def test_structure_conditioned_null_selection_preserves_original_and_recovers_floor(tmp_path):
    from vdjmatch.aggregate.combination import null_normalized_evidence
    from vdjmatch.cli.__main__ import main
    context=dict(species='human',locus='paired',epitope='E',mhc_a='HLA-A*02:01',mhc_b='B2M',mhc_class='MHCI')
    frame=pl.DataFrame([dict(query_id=str(i),likelihood_score=l,conditional_rank_score=r,
        original_support=u,reference_coverage=0. if i==0 else 1.,knee_shape=1.,knee_radius_ratio=0.,paired=True,**context)
        for i,(l,r,u) in enumerate([(-30.,10.,0.),(-30.,10.,1.),(10.,1.,0.),(None,10.,1.),(None,None,0.)])])
    nulls=pl.DataFrame([dict(query_id=str(i),likelihood_score=float(i),conditional_rank_score=float(i+1),
        outcome='null',null_source='reference',**context) for i in range(9)])
    scored=null_normalized_evidence(frame,nulls,structure_conditioned=True)
    with pytest.raises(ValueError,match='missing columns: knee_shape'):
        null_normalized_evidence(frame.drop('knee_shape'),nulls,structure_conditioned=True)
    np.testing.assert_allclose(scored['normalized_score'].to_list()[:4],[np.log(10),0.,np.log(10),np.log(10)])
    assert scored['normalized_score'][4] is None
    assert np.all(scored['normalized_score'].to_numpy()[:4]>=-np.log(scored['p_component_likelihood'].to_numpy()[:4]))
    for column in ('likelihood_score','conditional_rank_score'):
        transformed=lambda f:f.with_columns((pl.col(column)*100+3).alias(column))
        assert null_normalized_evidence(transformed(frame),transformed(nulls),structure_conditioned=True)['normalized_score'].equals(scored['normalized_score'])
    source=tmp_path/'evidence.tsv';control=tmp_path/'null.tsv';output=tmp_path/'scores.tsv'
    frame.write_csv(source,separator='\t');nulls.write_csv(control,separator='\t')
    args=['combine-evidence',str(source),'--structure-conditioned','--output',str(output)]
    with pytest.raises(SystemExit):main(args+['--fixed'])
    assert main(args+['--null-evidence',str(control)])==0
    assert pl.read_csv(output,separator='\t')['combined_score'].equals(scored['normalized_score'])


def test_structure_conditioned_selection_has_no_query_support_switch():
    """The same component evidence and repertoire must produce the same score."""
    from vdjmatch.aggregate.combination import null_normalized_evidence
    context=dict(species='human',locus='TRA',epitope='E',mhc_a='A*02:01',mhc_b='B2M',mhc_class='MHCI')
    q=pl.DataFrame([dict(query_id=str(i),likelihood_score=-30.,conditional_rank_score=10.,
        original_support=float(i),reference_coverage=1.,knee_shape=1.,knee_radius_ratio=0.,**context)
        for i in range(2)])
    nulls=pl.DataFrame([dict(query_id=str(i),likelihood_score=float(i),conditional_rank_score=float(i+1),
        outcome='null',null_source='reference',**context) for i in range(9)])
    dense=null_normalized_evidence(q,nulls,structure_conditioned=True)
    np.testing.assert_allclose(dense['normalized_score'].to_numpy(),[0.,0.])
    assert null_normalized_evidence(q.drop('original_support'),nulls,structure_conditioned=True)['normalized_score'].equals(dense['normalized_score'])
    sparse=null_normalized_evidence(q.with_columns(pl.lit(0.).alias('reference_coverage')),
                                    nulls,structure_conditioned=True)
    np.testing.assert_allclose(sparse['normalized_score'].to_numpy(),np.log(10))


def test_three_channel_sum_keeps_scales_weights_ties_and_missing_evidence(tmp_path):
    from vdjmatch.aggregate.combination import null_normalized_evidence, THREE_COMPONENTS, THREE_WEIGHTS
    from vdjmatch.cli.__main__ import main
    context=dict(species='human',locus='TRB',epitope='E',mhc_a='HLA-A*02',mhc_b='B2M',mhc_class='MHCI')
    q=pl.DataFrame([dict(query_id=f'q:{i}',original_score=o,gapblock_score=g,stem4_score=s,
        original_weight=.5,gapblock_weight=.25,stem4_weight=.25,**context)
        for i,(o,g,s) in enumerate([(2.,20.,200.),(None,20.,None),(None,None,None)])])
    nulls=pl.DataFrame([dict(query_id=str(i),original_score=o,gapblock_score=g,stem4_score=s,
        outcome='null',null_source=source,**context) for i,(o,g,s,source) in enumerate([
            (1.,10.,100.,'withheld'),(2.,20.,200.,'withheld'),(None,None,None,'withheld'),
            (3.,30.,300.,'generated'),(None,None,None,'generated')])])
    result=null_normalized_evidence(q,nulls,three_channel=True)
    np.testing.assert_allclose(result['normalized_score'].to_list()[:2],[-np.log(2/3),-.25*np.log(2/3)])
    assert result['normalized_score'][2] is None
    contributions=result.select(pl.sum_horizontal(pl.col('^contribution_.*$'))).to_series()
    np.testing.assert_allclose(contributions.head(2),result['normalized_score'].head(2))
    for column in THREE_COMPONENTS:
        transform=lambda f:f.with_columns((pl.col(column)*17-13).alias(column))
        assert null_normalized_evidence(transform(q),transform(nulls),three_channel=True)['normalized_score'].equals(result['normalized_score'])
    zero=q.with_columns(pl.lit(0.).alias('gapblock_weight'),pl.lit(.5).alias('stem4_weight'))
    assert null_normalized_evidence(zero,nulls,three_channel=True)['normalized_score'][1] is None
    for weights in ([1.,0.,.1],[-1.,1.,1.],[float('nan'),0.,1.],[None,0.,1.]):
        invalid=q.with_columns(pl.lit(w).alias(c) for c,w in zip(THREE_WEIGHTS,weights))
        with pytest.raises(ValueError,match='weights'):
            null_normalized_evidence(invalid,nulls,three_channel=True)
    source=tmp_path/'evidence.tsv';control=tmp_path/'null.tsv';output=tmp_path/'scores.tsv'
    q.write_csv(source,separator='\t');nulls.write_csv(control,separator='\t')
    args=['combine-evidence',str(source),'--three-channel','--output',str(output)]
    with pytest.raises(SystemExit):main(args+['--fixed'])
    assert main(args+['--null-evidence',str(control)])==0
    saved=pl.read_csv(output,separator='\t')
    assert saved['combined_score'].equals(result['normalized_score'])
    assert saved.select(q.columns).equals(q)
