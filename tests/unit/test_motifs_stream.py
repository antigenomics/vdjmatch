"""The MF×IKF consumer retains fixed model gates and the materialized score law."""
import polars as pl
import pytest
from polars.testing import assert_frame_equal

from vdjmatch.aggregate.candidates import PMHC
from vdjmatch.match import motifs


@pytest.mark.parametrize('exclude', [False, True])
@pytest.mark.parametrize('targeted', [False, True])
def test_streamed_motif_blocks_match_materialized_model(monkeypatch, exclude, targeted):
    rows=[{'cdr3':'CASSAAAA'+a+'F','v':'TRBV1*01','epitope':'A'} for a in 'AILVDEN']
    rows += [{'cdr3':'CASSGGGG'+a+'F','v':'TRBV1*01','epitope':'B'} for a in 'AILVDEN']
    reference=pl.DataFrame(rows).with_columns(pl.lit('A').alias('mhc_a'),
        pl.lit('B').alias('mhc_b'),pl.lit('I').alias('mhc_class'))
    query=pl.DataFrame({'query_id':range(9),'cdr3':[rows[0]['cdr3']]*9,'v':['TRBV1*02']*9})
    target=None
    if targeted:
        query=query.with_columns(*[pl.lit(value).alias(name) for name,value in
            zip(PMHC,['A','A','B','I'])])
        target=reference.filter(pl.col('epitope')=='A')
    expected=motifs.masked_motif_scores(query,reference,exclude_exact=exclude,target_reference=target)
    monkeypatch.setattr(motifs,'_STREAM_TOKEN_ROWS',64)
    monkeypatch.setattr(motifs,'_STREAM_OUTPUT_ROWS',1)
    parts=[]
    assert motifs.masked_motif_scores(query,reference,exclude_exact=exclude,
        target_reference=target,_consume=parts.append) is None
    assert len(parts)>=9 and all(frame.height<=1 for frame in parts)
    assert_frame_equal(pl.concat(parts).sort('query_id',*PMHC),
                       expected.sort('query_id',*PMHC),check_exact=True)


@pytest.mark.parametrize('budget,message',[(1,'token budget'),(64,'expanded-hit budget')])
def test_streamed_motif_budget_fails_before_oversized_query_join(monkeypatch,budget,message):
    reference=pl.DataFrame({'cdr3':['CASSAAAA'+a+'F' for a in 'AILVDEN'],
        'v':['TRBV1']*7,'epitope':['A']*7,'mhc_a':['A']*7,'mhc_b':['B']*7,'mhc_class':['I']*7})
    query=reference.head(1).select('cdr3','v').with_row_index('query_id')
    reference=pl.concat([reference.with_columns(pl.lit(str(i)).alias('epitope')) for i in range(10)])
    monkeypatch.setattr(motifs,'_STREAM_TOKEN_ROWS',budget)
    with pytest.raises(ValueError,match=message):
        motifs.masked_motif_scores(query,reference,_consume=lambda frame:pytest.fail('over budget'))


@pytest.mark.parametrize('streamed',[False,True])
@pytest.mark.parametrize('targeted',[False,True])
def test_fixed_motif_model_excludes_missing_reference_v_only(streamed,targeted):
    reference=pl.DataFrame({'cdr3':['CASSAAAA'+a+'F' for a in 'AILVDEN'],
        'v':['TRBV1']*7,'epitope':['A']*7,'mhc_a':['A']*7,'mhc_b':['B']*7,'mhc_class':['I']*7})
    query=reference.head(1).with_row_index('query_id')
    missing=reference.head(2).with_columns(pl.Series('v',[None,'  ']))
    mixed=pl.concat([reference,missing])
    expected=motifs.masked_motif_scores(query,reference,exclude_exact=True,
        target_reference=reference if targeted else None)
    parts=[]
    actual=motifs.masked_motif_scores(query,mixed,exclude_exact=True,
        target_reference=mixed if targeted else None,_consume=parts.append if streamed else None)
    if streamed:
        assert actual is None
        actual=pl.concat(parts)
    assert_frame_equal(actual.sort('query_id',*PMHC),expected.sort('query_id',*PMHC),check_exact=True)
    assert motifs.masked_motif_scores(query,missing).is_empty()
    with pytest.raises(ValueError,match='valid complete junction'):
        motifs.masked_motif_scores(query,missing.with_columns(pl.lit('BAD?').alias('cdr3')))
    with pytest.raises(ValueError,match='complete junction and V'):
        motifs.masked_motif_scores(query.with_columns(pl.lit(None).alias('v')),reference)
