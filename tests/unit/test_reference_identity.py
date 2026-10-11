"""Declared pair links survive source normalization and index identity generation."""
import polars as pl
import pytest
from vdjmatch import Annotator
from vdjmatch.db import schema
from vdjmatch.match import VdjdbIndex


def pairs(**extra):
    return pl.DataFrame({'gene':['TRA','TRB'],'cdr3':['CAVVF','CASSF'],
                         'epitope':['PEP','PEP'],**extra})


def test_reference_curation_preserves_rows_and_natural_same_peptide():
    raw = pl.DataFrame({
        'gene': ['TRB'] * 5, 'cdr3': ['CASSF'] * 5,
        'epitope': ['KLGGALQAK', 'SLLMWITQV', 'SLLMWITQV', 'SLLMWITQV', 'PEP'],
        'reference.id': [None, 'PMID:40498839', 'PMID:32102902',
                         '["https://pubmed.ncbi.nlm.nih.gov/40498839/"]', 'PMID:140498839'],
    })
    normalized = schema.normalize(raw)
    assert normalized.height == raw.height
    assert normalized['reference_valid'].to_list() == [False, False, True, False, True]
    assert normalized['reference_exclusion_reason'].to_list() == [
        'KLGGALQAK_artifact', 'engineered_phage_library_PMID40498839', None,
        'engineered_phage_library_PMID40498839', None]
    assert schema.normalize(normalized).equals(normalized)
    idx = VdjdbIndex.build(raw)
    assert idx.records_for('TRB')['epitope'].to_list() == ['SLLMWITQV', 'PEP']


def test_paired_builder_applies_curation_to_targets_and_complement():
    from vdjmatch.evalue.paired import build_paired_ref
    raw=pl.concat([
        pairs(complex_id=['artifact']*2).with_columns(pl.lit('KLGGALQAK').alias('epitope'),pl.lit(None,dtype=pl.String).alias('reference.id')),
        pairs(complex_id=['engineered']*2).with_columns(pl.lit('SLLMWITQV').alias('epitope'),pl.Series('reference.id',['PMID:40498839','PMID:32102902'])),
        pairs(complex_id=['natural']*2).with_columns(pl.lit('SLLMWITQV').alias('epitope'),pl.lit('PMID:32102902').alias('reference.id')),
    ])
    ref=schema.normalize(raw)
    assert ref.height==6  # Keep all source observations and exclusion reasons.
    out=build_paired_ref(ref)
    assert out['complex_id'].to_list()==['natural']
    assert out['epitope'].to_list()==['SLLMWITQV']


def test_string_complex_without_record_ids_survives_index_and_joint_query():
    raw=pairs(complex_id=['self','self'])
    normalized=schema.normalize(raw)
    assert normalized['complex_id'].to_list()==['self','self']
    idx=VdjdbIndex.build(raw)
    assert idx.records_for('TRA')['complex_id'].to_list()==['self']
    assert idx.records_for('TRB')['complex_id'].to_list()==['self']
    generated=pl.concat([idx.records_for('TRA'),idx.records_for('TRB')])
    assert generated['record_id'].n_unique()==2
    assert schema.normalize(generated)['complex_id'].to_list()==['self','self']
    assert schema.normalize(normalized).equals(normalized)
    ann=Annotator.from_frame(raw)
    q=pl.DataFrame({'cdr3_alpha_aa':['CAVVF'],'cdr3_beta_aa':['CASSF']})
    assert ann.paired_candidates(q,scope='0')['n_hits'].to_list()==[1]


def test_explicit_pair_links_override_distinct_chain_record_ids():
    raw=pairs(complex_id=['complex','complex'],record_id=['alpha-record','beta-record'])
    out=schema.normalize(raw)
    assert out['complex_id'].to_list()==['complex','complex']
    assert out['record_id'].to_list()==['alpha-record','beta-record']
    assert schema.normalize(out).equals(out)


@pytest.mark.parametrize('identity',['record_id','cell_id'])
def test_missing_links_infer_observation_and_inherit_single_declared_link(identity):
    raw=pairs(**{identity:['observation','observation']})
    assert schema.normalize(raw)['complex_id'].to_list()==['observation','observation']
    raw=raw.with_columns(pl.Series('complex_id',['declared',None]))
    assert schema.normalize(raw)['complex_id'].to_list()==['declared','declared']
    raw=raw.with_columns(pl.Series('complex_id',['0','0']))
    assert schema.normalize(raw)['complex_id'].to_list()==['observation','observation']


def test_conflicting_declared_links_within_shared_observation_fail():
    with pytest.raises(ValueError,match='conflicting.*complex'):
        schema.normalize(pairs(record_id=['observation','observation'],complex_id=['one','two']))


def test_inferred_link_cannot_collide_with_unrelated_declared_link():
    raw=pl.concat([pairs(record_id=['collision','collision'],complex_id=['0','0']),
                   pairs(record_id=['other','other'],complex_id=['collision','collision'])])
    with pytest.raises(ValueError,match='collides'):
        schema.normalize(raw)


def test_legacy_integer_links_preserved_without_lossy_string_coercion():
    raw=pairs(**{'complex.id':[7,7]})
    out=schema.normalize(raw)
    assert out.schema['complex_id']==pl.Int64
    assert out['complex_id'].to_list()==[7,7]
    assert schema.normalize(out).equals(out)
    raw=pairs(complex_id=['01','01'])
    assert schema.normalize(raw)['complex_id'].to_list()==['01','01']
    empty=pairs(complex_id=[None,None])
    assert schema.normalize(empty)['complex_id'].to_list()==[0,0]


@pytest.mark.parametrize('link_column',['complex_id','complex.id'])
def test_wide_tables_preserve_declared_pair_links(link_column):
    raw=pl.DataFrame({'cdr3.alpha':['CAVVF','CAVAF'], 'cdr3.beta':['CASSF','CASRF'],
                      'epitope':['PEP','PEP'],link_column:['source-pair','0']})
    out=schema.normalize(raw)
    assert out['complex_id'].to_list()==['source-pair','__wide_row:1','source-pair','__wide_row:1']
    assert schema.normalize(out).equals(out)


@pytest.mark.parametrize('identity',['record_id','cell_id'])
@pytest.mark.parametrize('value',[7,'7','01'])
def test_numeric_inferred_identity_is_idempotent(identity,value):
    raw=pairs(**{identity:[value,value]})
    first=schema.normalize(raw)
    second=schema.normalize(first)
    assert second.equals(first)
    assert first[identity].to_list()==[value,value]
    if value=='01':
        assert first['complex_id'].to_list()==['01','01']
    else:
        assert first.schema['complex_id']==pl.Int64
        assert first['complex_id'].to_list()==[7,7]


@pytest.mark.parametrize('identity',['record_id','cell_id'])
@pytest.mark.parametrize('value',[0,'0'])
def test_zero_observation_identity_is_retained_and_paired(identity,value):
    raw=pairs(**{identity:[value,value]})
    first=schema.normalize(raw)
    assert first['complex_id'].to_list()==[f'__{identity}:0']*2
    assert first[identity].to_list()==[value,value]
    assert schema.normalize(first).equals(first)
    ann=Annotator.from_frame(raw)
    q=pl.DataFrame({'cdr3_alpha_aa':['CAVVF'],'cdr3_beta_aa':['CASSF']})
    assert ann.paired_candidates(q,scope='0')['n_hits'].to_list()==[1]


def test_zero_inference_namespace_collision_fails():
    raw=pl.concat([pairs(record_id=['0','0'],complex_id=['0','0']),
                   pairs(record_id=['other','other'],complex_id=['__record_id:0']*2)])
    with pytest.raises(ValueError,match='collides'):
        schema.normalize(raw)


@pytest.mark.parametrize('identity',['record_id','cell_id'])
def test_zero_namespace_cannot_merge_another_inferred_observation(identity):
    raw=pl.concat([pairs(**{identity:['0','0']}),
                   pairs(**{identity:[f'__{identity}:0']*2})])
    with pytest.raises(ValueError,match='collide between observations'):
        schema.normalize(raw)
