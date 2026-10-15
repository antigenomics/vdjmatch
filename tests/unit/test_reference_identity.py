"""Declared pair links survive source normalization and index identity generation."""
import polars as pl
import pytest
from vdjmatch import Annotator
from vdjmatch.db import schema
from vdjmatch.match import VdjdbIndex


def test_replicated_support_preserves_pmhc_species_and_curation():
    from vdjmatch import db
    raw = pl.DataFrame({
        'gene':['TRB']*8, 'cdr3':['CASSF']*8,
        'v':['TRBV5-1']*8, 'j':['TRBJ1-1']*8,
        'epitope':['PEP']*8, 'species':['HomoSapiens']*6+['MusMusculus']*2,
        'mhc_a':['A*02:01']*3+['A*01:01']*3+['A*02:01']*2,
        'mhc_b':['B2M']*8, 'mhc_class':['MHCI']*8,
        'reference_id':[' PMID:1 ','PMID:1','PMID:2','PMID:3','PMID:4','   ','PMID:5','PMID:6'],
        'reference_valid':[True,True,True,True,False,True,True,True],
    })
    out = db.replicated(raw)
    assert out.select('species','mhc_a','mhc_b','n_refs').sort('species').rows() == [
        ('HomoSapiens','A*02:01','B2M',2), ('MusMusculus','A*02:01','B2M',2)]
    assert db.replicated(raw, min_refs=1).height == 3


@pytest.mark.parametrize('value',[0,-1,True,1.5])
def test_replicated_rejects_invalid_reference_count(value):
    from vdjmatch import db
    with pytest.raises(ValueError, match='positive integer'):
        db.replicated(pl.DataFrame(), min_refs=value)


def test_reference_support_cli(tmp_path):
    import subprocess
    import sys
    raw = pairs(v=['TRAV1-2','TRBV5-1'],j=['TRAJ33','TRBJ1-1']).with_columns(
        pl.lit('HomoSapiens').alias('species'),pl.lit('A*02:01').alias('mhc_a'),
        pl.lit('B2M').alias('mhc_b'),pl.lit('MHCI').alias('mhc_class'))
    raw = pl.concat([raw.with_columns(pl.lit('PMID:1').alias('reference_id')),
                     raw.with_columns(pl.lit('PMID:2').alias('reference_id'))])
    source, output = tmp_path/'reference.tsv', tmp_path/'support.tsv'
    raw.write_csv(source, separator='\t')
    result = subprocess.run([sys.executable,'-m','vdjmatch.cli','reference-support',
        '--vdjdb',str(source),'--output',str(output)],capture_output=True,text=True,timeout=20)
    assert result.returncode == 0, result.stderr
    out = pl.read_csv(output, separator='\t')
    assert out['n_refs'].to_list() == [2,2]
    assert {'species','mhc_a','mhc_b','mhc_class'} <= set(out.columns)
    before=source.read_bytes()
    result = subprocess.run([sys.executable,'-m','vdjmatch.cli','reference-support',
        '--vdjdb',str(source),'--output',str(source)],capture_output=True,text=True,timeout=20)
    assert result.returncode!=0 and source.read_bytes()==before


@pytest.mark.parametrize('gene',['TRBV5-1','TRBV21/OR9-2'])
def test_replicated_gene_resolution_is_explicit(gene):
    from vdjmatch import db
    raw = pl.DataFrame({'gene':['TRB']*2,'cdr3':['CASSF']*2,
        'v':[gene+'*01',gene+'*02'],'j':['TRBJ1-1*01','TRBJ1-1*02'],
        'epitope':['PEP']*2,'reference_id':['PMID:1','PMID:2']})
    assert db.replicated(raw).is_empty()
    out = db.replicated(raw, gene_level=True)
    assert out.select('v','j','n_refs').rows() == [(gene,'TRBJ1-1',2)]


def test_replicated_gene_resolution_cannot_turn_ambiguous_calls_into_support():
    from vdjmatch import db
    raw = pl.DataFrame({'gene':['TRB']*4,'cdr3':['CASSF']*4,
        'v':['TRBV5-1*01,TRBV7-9*01']*2+['TRAV1-2*01']*2,
        'j':['TRBJ1-1*01']*4,'epitope':['PEP']*4,
        'reference_id':['PMID:1','PMID:2']*2})
    assert db.replicated(raw, gene_level=True).is_empty()


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


@pytest.mark.parametrize('link_column', ['pair_id', 'clone_id'])
def test_custom_long_junction_reference_preserves_declared_pair(link_column):
    raw = pl.DataFrame({'junction_aa': ['CAVVF', 'CASSF'], 'locus': ['TRA', 'TRB'],
        'epitope': ['PEP'] * 2, link_column: ['custom-pair'] * 2})
    out = schema.normalize(raw)
    assert out['complex_id'].to_list() == ['custom-pair'] * 2
    assert out[link_column].to_list() == ['custom-pair'] * 2
    assert schema.normalize(out).equals(out)


def test_custom_wide_junction_reference_loader_retains_sources_and_marginal(tmp_path):
    from vdjmatch import db
    raw = pl.DataFrame({'junction_aa_alpha': ['CAVVF', None], 'junction_aa_beta': ['CASSF', 'CASRF'],
        'v_alpha': ['TRAV1-2', None], 'j_alpha': ['TRAJ33', None],
        'v_beta': ['TRBV5-1'] * 2, 'j_beta': ['TRBJ1-1'] * 2,
        'pair_id': ['custom-pair', '0'], 'epitope': ['PEP'] * 2, 'source_note': ['paired', 'marginal']})
    path = tmp_path / 'custom.tsv'; raw.write_csv(path, separator='\t')
    out = db.load(path)
    assert out.select('gene', 'cdr3', 'v', 'j', 'complex_id', 'source_note').rows() == [
        ('TRA', 'CAVVF', 'TRAV1-2', 'TRAJ33', 'custom-pair', 'paired'),
        ('TRB', 'CASSF', 'TRBV5-1', 'TRBJ1-1', 'custom-pair', 'paired'),
        ('TRB', 'CASRF', 'TRBV5-1', 'TRBJ1-1', '0', 'marginal')]
    assert set(raw.columns) <= set(out.columns)
    assert schema.normalize(out).equals(out)


@pytest.mark.parametrize('dropped,expected', [([], 'canonical'), (['complex_id'], 'dot'),
    (['complex_id', 'complex.id'], 'pair'), (['complex_id', 'complex.id', 'pair_id'], 'clone')])
def test_custom_pair_link_field_precedence(dropped, expected):
    raw = pairs(complex_id=['canonical'] * 2, **{'complex.id': ['dot'] * 2},
        pair_id=['pair'] * 2, clone_id=['clone'] * 2).drop(dropped)
    assert schema.normalize(raw)['complex_id'].to_list() == [expected] * 2


def test_custom_beta_only_wide_is_unlinked_and_bare_cdr3_is_rejected():
    out = schema.normalize(pl.DataFrame({'junction_aa_beta': ['CASSF'], 'epitope': ['PEP']}))
    assert out.select('gene', 'cdr3', 'complex_id').rows() == [('TRB', 'CASSF', 0)]
    with pytest.raises(ValueError, match='anchors'):
        schema.normalize(pl.DataFrame({'cdr3_aa_alpha': ['AVV'], 'cdr3_aa_beta': ['ASS'], 'epitope': ['PEP']}))


@pytest.mark.parametrize('layout', ['wide', 'long'])
def test_custom_paired_reference_one_control_cli_preserves_mock_roles(tmp_path, monkeypatch, layout):
    from vdjmatch.cli.__main__ import main
    from vdjmatch.evalue import combined
    from vdjmatch.evalue.paired import build_paired_ref
    reference = tmp_path / 'custom-reference.tsv'
    pl.DataFrame({'junction_aa_alpha': ['CAVRDSNYQLIW', None],
        'junction_aa_beta': ['CASSLGQETQYF', 'CASSIGQETQYF'],
        'v_alpha': ['TRAV1-2', None], 'j_alpha': ['TRAJ33', None],
        'v_beta': ['TRBV1'] * 2, 'j_beta': ['TRBJ1-1'] * 2,
        'pair_id': ['reference-pair', '0'], 'epitope': ['E'] * 2}).write_csv(reference, separator='\t')
    if layout == 'long':
        pl.DataFrame({'junction_aa': ['CAVRDSNYQLIW', 'CASSLGQETQYF', 'CASSIGQETQYF'],
            'v_call': ['TRAV1-2', 'TRBV1', 'TRBV1'], 'j_call': ['TRAJ33', 'TRBJ1-1', 'TRBJ1-1'],
            'locus': ['TRA', 'TRB', 'TRB'], 'pair_id': ['reference-pair', 'reference-pair', '0'],
            'epitope': ['E'] * 3}).write_csv(reference, separator='\t')
    sample = tmp_path / 'sample.tsv'
    pl.DataFrame({'query_id': ['complete', 'partial'],
        'junction_aa_alpha': ['CAVRDMNYQLIW', None],
        'junction_aa_beta': ['CASSMGQETQYF', 'CASSVGQETQYF'],
        'v_alpha': ['TRAV1-2', None], 'j_alpha': ['TRAJ33', None],
        'v_beta': ['TRBV1'] * 2, 'j_beta': ['TRBJ1-1'] * 2}).write_csv(sample, separator='\t')
    control = tmp_path / 'control.tsv'
    pl.DataFrame({'junction_aa': ['CAVRDTNYQLIW', 'CAVRDMNYQLIW', 'CASSAGQETQYF', 'CASSVGQETQYF'],
        'locus': ['TRA', 'TRA', 'TRB', 'TRB'],
        'v_call': ['TRAV1-2'] * 2 + ['TRBV1'] * 2,
        'j_call': ['TRAJ33'] * 2 + ['TRBJ1-1'] * 2,
        'occurrences': [4] * 4}).write_csv(control, separator='\t')
    monkeypatch.setattr(combined, 'compute_pgen', lambda ref, events, **kw:
        pl.concat([ref.select('junction_aa'), events.select('junction_aa')]).unique().with_columns(pl.lit(.01).alias('pgen')))
    observed = []
    mock_models = combined.iter_mock_models
    def verify_models(model, *args, **kwargs):
        expected = model.select('gene', 'complex_id')
        paired_population = 'TRA' in args[0]
        roles = [('TRB', 'reference-pair'), ('TRB', '0')]
        assert expected.rows() == ([('TRA', 'reference-pair')] + roles if paired_population else roles)
        assert build_paired_ref(model).height == int(paired_population)
        for randomized in mock_models(model, *args, **kwargs):
            assert randomized.select('gene', 'complex_id').equals(expected)
            assert build_paired_ref(randomized).height == int(paired_population)
            observed.append(randomized.height)
            yield randomized
    monkeypatch.setattr(combined, 'iter_mock_models', verify_models)
    prefix = tmp_path / 'result'
    assert main(['annotate', str(sample), '--vdjdb', str(reference), '--control', str(control),
        '--mock-replicates', '2', '--seed', '7', '-o', str(prefix)]) == 0
    queries = pl.read_csv(str(prefix) + '.queries.tsv', separator='\t')
    assert queries.height == 4 and queries['calibration_status'].to_list().count('missing_sequence') == 1
    assert queries['calibration_pair_status'].to_list().count('complete_pair') == 2
    assert observed == [2, 2, 3, 3]  # Beta marginal models, then complete paired models.
    assert not list(tmp_path.glob('result*.json'))


@pytest.mark.parametrize('identity', ['record_id', 'cell_id'])
def test_clone_metadata_does_not_merge_declared_observation_pairs(identity):
    raw = pl.concat([pairs(**{identity: [name] * 2, 'clone_id': ['shared-clone'] * 2})
                     for name in ('observation-one', 'observation-two')])
    out = schema.normalize(raw)
    assert out['complex_id'].to_list() == ['observation-one'] * 2 + ['observation-two'] * 2
    assert out['clone_id'].to_list() == ['shared-clone'] * 4
    assert schema.normalize(out).equals(out)
    explicit = raw.with_columns(pl.Series('pair_id', ['explicit-one'] * 2 + ['explicit-two'] * 2))
    assert schema.normalize(explicit)['complex_id'].to_list() == ['explicit-one'] * 2 + ['explicit-two'] * 2
