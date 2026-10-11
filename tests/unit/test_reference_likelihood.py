"""Analytic Gaussian evidence, exact model puncturing and missing-feature routes."""
import numpy as np
import polars as pl
import pytest
from scipy.special import ndtr

from vdjmatch.match.reference_likelihood import reference_likelihood


def test_joint_reference_likelihood_matches_gaussian_and_physical_deletion():
    def frame(values, junctions):
        return pl.DataFrame({'query_id': [str(i) for i in range(len(values))],
            'junction_alpha': junctions, 'junction_beta': ['CASSF']*len(values),
            'original_alpha_density_score': [None]*len(values),
            'original_beta_density_score': [None]*len(values),
            'original_germline_lr': values, 'conditional_score': [None]*len(values)})
    target = frame([1.,2.,3.], ['CAVF','CAAF','CAGF'])
    background = frame([-1.,0.,1.], ['CATF','CAWF','CAYF'])
    query = frame([2.], ['CASF'])
    got = reference_likelihood(query, target, background, paired=True)
    assert got['log_likelihood_ratio'].item() == pytest.approx(1.)
    assert got['p_gaussian_model'].item() == pytest.approx(ndtr(-np.sqrt(2.)))
    assert got['model_target_count'].item()==got['model_background_count'].item()==3
    shifted = [t.with_columns((pl.col('original_germline_lr')+1e9).alias('original_germline_lr'))
               for t in [query,target,background]]
    assert reference_likelihood(*shifted,paired=True)['score'].item()==pytest.approx(got['score'].item())
    # Two identities have conflicting feature records. Remove every matching junction,
    # not just a row or matching gene value, from both classes.
    target = pl.concat([target, target.head(1).with_columns(pl.lit('duplicate').alias('query_id'),
        pl.lit(4.).alias('original_germline_lr'))])
    background = pl.concat([background, target.head(1).with_columns(pl.lit('conflict').alias('query_id'))])
    query = pl.concat([query, target.head(1).with_columns(pl.lit('matched').alias('query_id'))])
    got = reference_likelihood(query, target, background, paired=True)
    for i, row in enumerate(query.iter_rows(named=True)):
        keep = ~((pl.col('junction_alpha')==row['junction_alpha']) &
                 (pl.col('junction_beta')==row['junction_beta']))
        expected = reference_likelihood(query.slice(i,1), target.filter(keep), background.filter(keep), paired=True)
        assert got.row(i)==expected.row(0)
    assert got['model_target_count'].to_list()==[4,2]
    assert got['model_background_count'].to_list()==[4,3]
    # Query outcomes do not influence fitting, score, or the conditional model tail.
    labelled = reference_likelihood(query.with_columns(pl.Series('binder',[0,1])),
        target.with_columns(pl.lit(0).alias('binder')),background.with_columns(pl.lit(1).alias('binder')),paired=True)
    assert labelled.drop('binder').equals(got)
    unavailable = reference_likelihood(query.with_columns(pl.lit(None,dtype=pl.Float64).alias('original_germline_lr')),
        target, background, paired=True)
    assert unavailable['score'].to_list()==[None,None]
    assert unavailable['p_gaussian_model'].to_list()==[1.,1.]
    depleted = reference_likelihood(target.head(1), target.head(1), background, paired=True)
    assert depleted['score'].item() is None and depleted['likelihood_status'].item()=='no_target_model'
    with pytest.raises(ValueError,match='finite or null'):
        reference_likelihood(query.with_columns(pl.lit(float('inf')).alias('original_germline_lr')),target,background,paired=True)


def test_joint_reference_likelihood_handles_correlation_and_missing_conditional():
    def frame(values, keys):
        return pl.DataFrame({'query_id':keys,'junction_aa':keys,
            'original_score':values,'conditional_score':np.expm1(np.arcsinh(values))})
    target=frame(np.array([1.,2.,3.]),['CAVF','CAAF','CAGF'])
    background=frame(np.array([0.,.1,.2]),['CATF','CAWF','CAYF'])
    query=frame(np.array([1.5,2.5]),['CASF','CALF'])
    full=reference_likelihood(query,target,background)
    missing=reference_likelihood(query.with_columns(pl.lit(None,dtype=pl.Float64).alias('conditional_score')),target,background)
    # Exactly duplicated transformed information must not be multiplied as independent evidence.
    np.testing.assert_allclose(full['score'],missing['score'],rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(full['p_gaussian_model'],missing['p_gaussian_model'],rtol=1e-12,atol=1e-12)
    assert missing['likelihood_status'].to_list()==['available','available']


def test_reference_likelihood_shipping_cli_and_input_protection(tmp_path):
    from vdjmatch.cli.__main__ import main
    import json
    def frame(values, keys):
        return pl.DataFrame({'query_id':keys,'junction_aa':keys,'original_score':np.sinh(values),
            'conditional_score':[None]*len(keys),'binder':[0]*len(keys)})
    cases,target,background,out=[tmp_path/name for name in ['q.tsv','t.tsv','b.tsv','out.tsv']]
    frame([2.],['CASF']).write_csv(cases,separator='\t')
    frame([1.,2.,3.],['CAVF','CAAF','CAGF']).write_csv(target,separator='\t')
    frame([-1.,0.,1.],['CATF','CAWF','CAYF']).write_csv(background,separator='\t')
    args=['composite-evidence',str(cases),'--background',str(background),
        '--target-reference',str(target),'--output',str(out)]
    assert main(args)==0
    result=pl.read_csv(out,separator='\t')
    assert result['score'].item()==pytest.approx(1.)
    assert result['p_gaussian_model'].item()==pytest.approx(ndtr(-np.sqrt(2.)))
    assert result['query_id'].item()=='CASF' and result['binder'].item()==0
    manifest=json.loads((tmp_path/'out.tsv.manifest.json').read_text())
    assert manifest['target_reference_queries']==3 and manifest['query_labels_used'] is False
    before=target.read_bytes()
    with pytest.raises(SystemExit):
        main(args[:-1]+[str(target)])
    assert target.read_bytes()==before
