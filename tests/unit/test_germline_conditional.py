"""Fixed conditional J frequency models and linked-pair categories."""
import polars as pl
import pytest
from vdjmatch.evalue.germline import conditional_j_usage

@pytest.mark.parametrize('paired',[False,True])
def test_conditional_j_counts_stay_fixed_under_hit_exclusion(paired):
    def frame(rows,query=False):
        f=pl.DataFrame(rows,schema=['seq0','v0','j0'],orient='row')
        if paired:f=f.with_columns(pl.col('seq0').alias('seq1'),pl.lit('TRAV1*01').alias('v1'),pl.lit('TRAJ1*01').alias('j1'))
        return f.with_row_index('query_id') if query else f
    q=frame([('CAAAAF','TRBV1*01','TRBJ1-1*01'),('CAAAGF','TRBV1*01',None)],True)
    r=frame([('CAAAAF','TRBV1*01','TRBJ1-1*01'),('CAAACF','TRBV1*01','TRBJ1-1*01'),('CAAADF','TRBV2*01','TRBJ2-1*01')])
    b=frame([('CAAAAF','TRBV1*01','TRBJ1-1*01'),('CAAAEF','TRBV1*01','TRBJ2-1*01'),('CAAAGF','TRBV2*01','TRBJ2-1*01')])
    out=conditional_j_usage(q,r,b,paired=paired,exclude_exact=True)
    assert out['j_reference_population'][0]==out['j_control_population'][0]==3
    assert out['j_reference_stratum'][0]==out['j_control_stratum'][0]==2
    assert out['j_reference_joint'][0]==2 and out['j_control_joint'][0]==1
    # Two observed J categories plus one unseen category; dense hand oracle.
    assert out['conditional_j_usage_ratio'][0]==pytest.approx((2+2.5/4.5)/(1+1.5/4.5))
    assert out.equals(conditional_j_usage(q,r,b,paired=paired,exclude_exact=False))
    assert out['conditional_j_usage_ratio'][1]==1 and out['conditional_j_status'][1]=='missing_query_j_neutral'
    neutral=conditional_j_usage(q,r,r,paired=paired)
    assert neutral['conditional_j_usage_ratio'].to_list()==[1,1]
    # Conflicting source J calls stay one ambiguous search key, not a first label.
    conflict=pl.concat([r,r.head(1).with_columns(pl.lit('TRBJ2-1*01').alias('j0'))])
    c=conditional_j_usage(q,conflict,b,paired=paired)
    assert c['j_reference_population'][0]==3 and c['j_reference_joint'][0]==1


def test_j_predictive_is_independent_of_other_query_categories():
    r=pl.DataFrame({'seq0':['CAAAAF'],'v0':['TRBV1*01'],'j0':['TRBJ1-1*01']})
    b=r.with_columns(pl.lit('TRBJ2-1*01').alias('j0'))
    q=r.with_row_index('query_id')
    extra=q.with_columns(pl.lit(1).cast(pl.UInt32).alias('query_id'),pl.lit('TRBJ2-7*01').alias('j0'))
    alone=conditional_j_usage(q,r,b)
    together=conditional_j_usage(pl.concat([q,extra]),r,b)
    assert alone.equals(together.head(1))


def test_linked_j_categories_are_joint_not_independent_chain_products():
    r=pl.DataFrame({'seq0':['CAAAAF','CAAACF'],'v0':['TRAV1*01']*2,'j0':['TRAJ1*01','TRAJ2*01'],'seq1':['CASSAF','CASSCF'],'v1':['TRBV1*01']*2,'j1':['TRBJ1-1*01','TRBJ2-1*01']})
    b=r.with_columns(pl.Series('j1',['TRBJ2-1*01','TRBJ1-1*01']))
    q=r.head(1).with_row_index('query_id')
    out=conditional_j_usage(q,r,b,paired=True)
    # Chain marginals coincide; only their same-pair J co-occurrence differs.
    assert out['conditional_j_usage_ratio'][0]==pytest.approx(12)
