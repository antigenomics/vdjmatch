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
