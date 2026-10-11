"""Global evidence formula convergence, monotonicity and missing dispositions."""
import numpy as np
import polars as pl
import pytest
from vdjmatch.aggregate.combination import fit_evidence_combination,combine_evidence


def test_combination_fit_and_frozen_prediction():
    rng=np.random.default_rng(14);n=300
    l=rng.normal(size=n);r=rng.exponential(size=n);m=rng.uniform(size=n)
    q=pl.DataFrame({'likelihood_score':l,'conditional_rank_score':r,'motif_score':m,
        'reference_coverage':rng.uniform(size=n),'knee_shape':rng.uniform(size=n),'paired':np.arange(n)%2==0})
    y=(l+np.log1p(r)+m+rng.normal(scale=.7,size=n)>1).astype(int)
    model=fit_evidence_combination(q,y)
    score=combine_evidence(q,model)
    assert np.isfinite(score).all()
    assert score[y==1].mean()>score[y==0].mean()
    for c in ('likelihood_score','conditional_rank_score','motif_score'):
        stronger=q.with_columns((pl.col(c)+1).alias(c))
        assert np.all(combine_evidence(stronger,model)>=score-1e-12)
    unsupported=q.head(1).with_columns(pl.lit(None,dtype=pl.Float64).alias(c)
        for c in ('likelihood_score','conditional_rank_score','motif_score'))
    assert np.isnan(combine_evidence(unsupported,model)[0])
    bad=q.with_columns(pl.lit(2).alias('knee_shape'))
    with pytest.raises(ValueError,match=r'\[0,1\]'):combine_evidence(bad,model)


def test_cli_frozen_groups_and_unsupported_rows(tmp_path):
    from vdjmatch.cli.__main__ import main
    train=pl.DataFrame({'likelihood_score':[-2.,-1.,1.,2.],
        'conditional_rank_score':[0.,0.,1.,2.], 'motif_score':[0.,0.,1.,1.],
        'reference_coverage':[.5]*4,'knee_shape':[.2]*4,'paired':[False]*4,
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
        for c in ('likelihood_score','conditional_rank_score','motif_score'))
    path=tmp_path/'test.tsv';output=tmp_path/'scores.tsv';test.write_csv(path,separator='\t')
    assert main(['combine-evidence',str(path),'--model',str(model),'--require-heldout','--output',str(output)])==0
    scored=pl.read_csv(output,separator='\t')
    assert scored.height==4 and scored['combined_status'][0]=='unsupported'
    assert scored['combined_score'].null_count()==1
