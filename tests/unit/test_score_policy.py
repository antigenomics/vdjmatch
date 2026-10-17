"""A single runnable contract for blind selection, fallbacks and input validation."""
import polars as pl
import pytest
from vdjmatch.cli import score_policy as policy_module

def test_policy_preserves_rows_and_fallbacks():
    policy=pl.DataFrame(dict(species=['HomoSapiens']*4,locus=['TRA','TRB','paired','TRA'],
        epitope=['A','A','A','B'],reference_class=['V_dominant','junction_projection_dominant','general_fallback','V_dominant']))
    scores=pl.DataFrame(dict(species=['human']*6,locus=['TRB','TRA','paired','TRA','TRA','TRA'],
        epitope=['A','A','A','B','missing','B'],case_id=[5,4,3,2,1,0],label=[1,0,1,0,1,0],
        original_score=[-1.,-2.,-3.,-4.,-5.,None],rule_v_score=[2.,3.,4.,float('nan'),1.,None],rule_rank_score=[8.,7.,6.,5.,4.,None]))
    result=policy_module.apply_reference_policy(scores,policy)
    assert result['case_id'].to_list()==[5,4,3,2,1,0]
    assert result['policy_score'].to_list()==[8.,3.,-3.,None,-5.,None]
    assert result['policy_status'].to_list()==['selected_component','selected_component','reference_fallback','component_unavailable_no_score','no_reference_profile','no_finite_score']
    flipped=policy_module.apply_reference_policy(scores.with_columns(1-pl.col('label')),policy)
    assert result['policy_score'].equals(flipped['policy_score'])
    with pytest.raises(ValueError,match='unique'):
        policy_module.apply_reference_policy(scores,pl.concat([policy,policy.head(1)]))
    with pytest.raises(ValueError,match='unknown'):
        policy_module.apply_reference_policy(scores,policy.with_columns(pl.lit('invented').alias('reference_class')))


def test_fast_classes_factorization_fallback_and_pairs():
    classes=['v_only','motif_only','v_and_motif','no_supported_axis',
             'motif_only','motif_only','v_and_motif','v_only','motif_only','v_only']
    policy=pl.DataFrame(dict(species=['HomoSapiens']*10,locus=['TRA']*7+['paired','TRA+TRB','TRA'],
        epitope=list('ABCDEFGHIJ'),reference_class=classes))
    scores=pl.DataFrame(dict(case_id=list(range(11)),species=['human']*11,
        locus=['TRA']*7+['paired','TRA+TRB','TRA','TRA'],epitope=list('ABCDEFGHIJK'),
        label=[0]*11,stem4_score=[10.]*9+[None,11.],
        rule_v_score=[2.,2.,2.,2.,0.,float('nan'),2.,2.,2.,None,2.],
        rule_exact3_score=[6.,6.,6.,6.,6.,6.,float('inf'),6.,6.,6.,6.]))
    result=policy_module.apply_reference_policy(scores,policy,original_column='stem4_score')
    assert result['case_id'].equals(scores['case_id'])
    assert result['policy_score'].to_list()==[2.,3.,6.,10.,10.,10.,10.,10.,10.,None,11.]
    assert result['selected_score_column'].to_list()==[
        'rule_v_score','rule_exact3_score / rule_v_score','rule_exact3_score',
        'stem4_score','stem4_score','stem4_score','stem4_score','stem4_score',
        'stem4_score','rule_v_score','stem4_score']
    assert result['policy_status'].to_list()==[
        'selected_component','selected_component','selected_component','reference_fallback',
        'component_unavailable_default_fallback','component_unavailable_default_fallback',
        'component_unavailable_default_fallback','paired_default','paired_default',
        'no_finite_score','no_reference_profile']
    blind=policy_module.apply_reference_policy(scores.drop('label'),policy,original_column='stem4_score')
    assert result.drop('label').equals(blind)
    with pytest.raises(ValueError,match='reserved'):
        policy_module.apply_reference_policy(result,policy,original_column='stem4_score')
    with pytest.raises(ValueError,match='declared'):
        policy_module.apply_reference_policy(scores.drop('rule_exact3_score'),policy,original_column='stem4_score')
    v_policy=policy.filter(pl.col('reference_class').is_in(['v_only','no_supported_axis']))
    policy_module.apply_reference_policy(scores.drop('rule_exact3_score'),v_policy,original_column='stem4_score')


def test_fast_and_legacy_unavailable_components_remain_distinct():
    policy=pl.DataFrame(dict(species=['human']*2,locus=['TRA']*2,epitope=['A','B'],
        reference_class=['V_dominant','v_only']))
    scores=pl.DataFrame(dict(species=['human']*2,locus=['TRA']*2,epitope=['A','B'],
        original_score=[5.,5.],rule_v_score=[None,None]))
    result=policy_module.apply_reference_policy(scores,policy)
    assert result['policy_score'].to_list()==[None,5.]
    assert result['policy_status'].to_list()==['component_unavailable_no_score','component_unavailable_default_fallback']


@pytest.mark.parametrize('denominator',[0.,-1.,float('nan'),float('inf'),None])
def test_conditional_motif_requires_finite_positive_denominator(denominator):
    policy=pl.DataFrame(dict(species=['human'],locus=['TRA'],epitope=['A'],reference_class=['motif_only']))
    scores=pl.DataFrame(dict(species=['human'],locus=['TRA'],epitope=['A'],
        original_score=[7.],rule_v_score=[denominator],rule_exact3_score=[6.]))
    result=policy_module.apply_reference_policy(scores,policy)
    assert result['policy_score'].item()==7.
    assert result['selected_score_column'].item()=='original_score'
    assert result['policy_status'].item()=='component_unavailable_default_fallback'


def test_fast_policy_cli_manifest(tmp_path):
    import json
    from argparse import Namespace
    from vdjmatch.db.cache import sha256
    model=tmp_path/'reference-policy.tsv'
    pl.DataFrame(dict(species=['human'],locus=['TRA'],epitope=['A'],
        reference_class=['no_supported_axis'])).write_csv(model,separator='\t')
    (tmp_path/'manifest.json').write_text(json.dumps(dict(reference_model_asset=True,
        no_benchmark_inputs=True,outputs_sha256={model.name:sha256(model)})))
    scores=tmp_path/'scores.tsv'
    pl.DataFrame(dict(species=['human'],locus=['TRA'],epitope=['A'],
        default_score=[7.])).write_csv(scores,separator='\t')
    prefix=tmp_path/'selected'
    assert policy_module.main(Namespace(scores=str(scores),policy=str(model),
        original_column='default_score',output_prefix=str(prefix)))==0
    report=json.loads((tmp_path/'selected.manifest.json').read_text())
    assert report['labels_used'] is report['default_changed'] is False
    assert 'core-exact3' in report['policy']
    assert 'do not transfer' in report['calibration_status']
    assert pl.read_csv(tmp_path/'selected.scores.tsv',separator='\t')['policy_score'].item()==7.
