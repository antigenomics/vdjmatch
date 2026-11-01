"""Cohort batching preserves separate seeded calls, including shared receptor IDs."""
import pytest
import polars as pl
from polars.testing import assert_frame_equal
from vdjmatch.aggregate.repertoire_enrichment import reference_enrichment


def test_cohort_matches_separate_samples():
    sample=pl.DataFrame(dict(receptor_id=['r','r','orphan'],locus=['TRA','TRB','TRB'],
        junction_aa=['CAVSF','CASSF','CQQQF'],v=['TRAV1','TRBV2','TRBV2'],
        j=['TRAJ1','TRBJ1','TRBJ1'],count=[3,3,2]))
    reference=pl.DataFrame(dict(group_id=['X','X','Y'],locus=['TRA','TRB','TRB'],
        junction_aa=['CAVSF','CASSF','CWWEF'],v=['TRAV1','TRBV2','TRBV3'],
        j=['TRAJ1','TRBJ1','TRBJ2'],pgen_bin=[-10]*3))
    controls={locus:pl.DataFrame(dict(junction_aa=['CQQQF','CTTTF','CPPPF','CNNNF'],
        v=[v]*4,j=[j]*4,pgen_bin=[-10]*4))
        for locus,v,j in [('TRA','TRAV1','TRAJ1'),('TRB','TRBV2','TRBJ1')]}
    other=sample.with_columns(pl.col('count')+4).filter(pl.col('locus')=='TRB')
    inputs={'first':sample,'second':other}
    cohort=pl.concat([frame.with_columns(pl.lit(name).alias('sample_id')) for name,frame in inputs.items()])
    actual=reference_enrichment(cohort,reference,controls,replicates=3,chain_mode='all')
    for name,frame in inputs.items():
        expected=reference_enrichment(frame,reference,controls,replicates=3,chain_mode='all')
        for batch,single in zip(actual,expected):
            selected=batch.filter(pl.col('sample_id')==name).drop('sample_id').select(single.columns)
            assert_frame_equal(selected.sort(single.columns),single.sort(single.columns),check_dtypes=False)

    invalid=cohort.with_columns(pl.when(pl.col('locus')=='TRA').then(None).otherwise(pl.col('receptor_id')).alias('receptor_id'))
    with pytest.raises(ValueError,match='non-null'):
        reference_enrichment(invalid,reference,controls,replicates=3)
