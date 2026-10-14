"""Grouped count preparation keeps the established single/linked geometry."""
import polars as pl
import pytest

from vdjmatch.match import tcrdist


@pytest.mark.parametrize('paired', [False, True])
def test_grouped_distance_prepares_one_population_and_gene_codes(monkeypatch, paired):
    name='grouped_paired_sum_count_batch' if paired else 'grouped_count_batch'
    calls=[]
    def native(*args, **kwargs):
        calls.append((args,kwargs))
        return iter(())
    monkeypatch.setattr(tcrdist.gapblock,name,native,raising=False)
    q=pl.DataFrame({'seq0':['CASSLGQAYEQYF'], 'v0':['TRBV19*01']})
    r=pl.DataFrame({'seq0':['CASSLGQAYEQWF'], 'v0':['TRBV19*02']})
    if paired:
        q=q.rename({'seq0':'seq1','v0':'v1'}).with_columns(pl.lit('CAVRDSNYQLIW').alias('seq0'),pl.lit('TRAV12-2*01').alias('v0'))
        r=r.rename({'seq0':'seq1','v0':'v1'}).with_columns(pl.lit('CAVRDSNYQLVW').alias('seq0'),pl.lit('TRAV12-2*02').alias('v0'))
    list(tcrdist.grouped_total_distance_count_batch(q,r,[0,117000],paired=paired,
        reference_target_offsets=[0,1],reference_target_ids=[0],target_peptide_ids=[0],
        reference_exclusive_peptide_ids=[0],threads=1,ctrim=3,stem_depth=4,matrix='blosum62'))
    assert len(calls)==1
    args,options=calls[0]
    assert options['query_gene_ids']==options['reference_gene_ids']
    assert options['reference_target_offsets']==[0,1]
    assert options['threads']==1
    assert options['feature_offsets'+('_alpha' if paired else '')] is not None
    assert args[-1]==[0,117000]
