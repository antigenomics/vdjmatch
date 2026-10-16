"""Synthetic identity/coverage checks for explicitly supplied durable Pgen features."""
import hashlib
import json

import polars as pl
import pytest

from vdjmatch.db.cache import sha256
from vdjmatch.evalue.mock_reference import load_pgen_features


def _bank(tmp_path, monkeypatch):
    import vdjmatch.evalue.mock_reference as mocks
    sequences = ['CAAAAAAAF', 'CAAAAACAF', 'CAGGGGGGF']
    frame = pl.DataFrame({'junction_aa': sequences, 'pgen': [.25, .125, 0.],
        'species': ['human']*3, 'locus': ['TRB']*3})
    path = tmp_path/'bank.features.parquet'
    frame.write_parquet(path)
    files = {'genes_v.tsv': 'a'*64, 'manifest.json': 'b'*64, 'genes_j.tsv': 'c'*64, 'v_choice.tsv': 'd'*64}
    metadata = dict(organism='human', locus='TRB', chain_type='VJ',
        events={'v_choice': {'kind': 'choice', 'given': []}})
    monkeypatch.setattr(mocks, '_pgen_model_identity', lambda *args: {'files': files, 'metadata': metadata})
    model_hash = hashlib.sha256(''.join(f'{name}\t{files[name]}\n' for name in sorted(files)).encode()).hexdigest()
    manifest = dict(features=dict(path=path.name, sha256=sha256(path)),
        options=dict(species='human', locus='TRB', mismatches=1),
        model=dict(files=files, sha256=model_hash, collapse_alleles=False,
            metadata=metadata))
    manifest_path = tmp_path/'bank.manifest.json'
    manifest_path.write_text(json.dumps(manifest))
    return path, manifest_path, frame, manifest


def test_feature_bank_checks_identity_and_complete_union(tmp_path, monkeypatch):
    path, manifest_path, frame, manifest = _bank(tmp_path, monkeypatch)
    reference, control = frame.head(1), frame.tail(2)
    result, identity = load_pgen_features(path, reference, control, species='HomoSapiens')
    assert result.equals(frame.select('junction_aa', 'pgen'))
    assert identity['sha256'] == sha256(path)
    assert identity['manifest_sha256'] == sha256(manifest_path)
    assert identity['model'] == manifest['model']
    missing = pl.DataFrame({'junction_aa': ['CREF']})
    with pytest.raises(ValueError, match='cover 1 required'):
        load_pgen_features(path, missing, control)
    manifest['features']['sha256'] = '0'*64
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='table SHA256'):
        load_pgen_features(path, reference, control)


@pytest.mark.parametrize('section,key,value', [
    ('options', 'mismatches', 0), ('options', 'mismatches', True),
    ('options', 'species', 'mouse'), ('options', 'locus', 'TRA'),
    ('model', 'collapse_alleles', True), ('model', 'sha256', '0'*64),
    ('model', 'files', {'unknown.tsv': 'a'*64}),
    ('model', 'metadata', {'organism': 'mouse', 'locus': 'TRB'}),
    ('features', 'path', 'wrong.features.parquet'),
    ('model', 'files', {'manifest.json': 'b'*64}),
])
def test_feature_bank_rejects_wrong_manifest(tmp_path, monkeypatch, section, key, value):
    path, manifest_path, frame, manifest = _bank(tmp_path, monkeypatch)
    manifest[section][key] = value
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        load_pgen_features(path, frame, frame)


@pytest.mark.parametrize('kind', ['duplicate', 'unsorted', 'invalid', 'null', 'nan', 'negative', 'above_one', 'species', 'missing'])
def test_feature_bank_rejects_malformed_rows(tmp_path, monkeypatch, kind):
    path, manifest_path, frame, manifest = _bank(tmp_path, monkeypatch)
    if kind == 'duplicate':
        frame = pl.concat([frame, frame.head(1)])
    elif kind == 'unsorted':
        frame = frame.reverse()
    elif kind == 'invalid':
        frame = frame.with_columns(pl.lit('C*F').alias('junction_aa'))
    elif kind == 'null':
        frame = frame.with_columns(pl.lit(None, dtype=pl.Float64).alias('pgen'))
    elif kind in ('nan', 'negative', 'above_one'):
        frame = frame.with_columns(pl.lit({'nan': float('nan'), 'negative': -.1, 'above_one': 1.1}[kind]).alias('pgen'))
    elif kind == 'species':
        frame = frame.with_columns(pl.lit('mouse').alias('species'))
    else:
        frame = frame.drop('locus')
    frame.write_parquet(path)
    manifest['features']['sha256'] = sha256(path)
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        load_pgen_features(path, frame.select('junction_aa'), frame.select('junction_aa'))


@pytest.mark.parametrize('command', ['reference-structure', 'reference-coverage', 'reference-diagnostics'])
def test_diagnostics_use_features_without_native_pgen(tmp_path, monkeypatch, command):
    if command != 'reference-structure':
        pytest.importorskip('leidenalg')
    import vdjmatch.evalue.mock_reference as mocks
    from vdjmatch.cli.__main__ import main
    path, _, frame, _ = _bank(tmp_path, monkeypatch)
    reference, control, prefix = (tmp_path/name for name in ('reference.tsv', 'control.tsv', 'out'))
    frame.select(pl.col('junction_aa').alias('cdr3')).with_columns(
        pl.lit('TRBV1').alias('v'), pl.lit('TRBJ1-1').alias('j'),
        pl.lit('HomoSapiens').alias('species'), pl.lit('TRB').alias('gene'),
        pl.lit('DEMO').alias('epitope')).write_csv(reference, separator='\t')
    frame.select('junction_aa').with_columns(pl.lit('TRBV1').alias('v_call'),
        pl.lit('TRBJ1-1').alias('j_call'), pl.lit(10).alias('occurrences')).write_csv(control, separator='\t')
    def forbidden(*args, **kwargs):
        pytest.fail('supplied validated features must not recompute Pgen')
    monkeypatch.setattr(mocks, 'compute_pgen', forbidden)
    args = [command, '--vdjdb', str(reference), '--control', str(control), '--locus', 'TRB',
        '--pgen-features', str(path), '--mock-replicates', '2', '--output-prefix', str(prefix)]
    if command != 'reference-structure':
        args += ['--geometry', 'levenshtein', '--nboot', '3']
    assert main(args) == 0
    if command == 'reference-coverage':
        identity = json.loads((tmp_path/'out.manifest.json').read_text())['control']['pgen_features']
    else:
        identity = json.loads((tmp_path/'out.manifest.json').read_text())['motif_information']['control_resampling']['pgen_features']
    assert identity['sha256'] == sha256(path) and identity['options']['mismatches'] == 1


def test_feature_bank_can_move_with_its_manifest(tmp_path, monkeypatch):
    path, manifest_path, frame, _ = _bank(tmp_path, monkeypatch)
    destination = tmp_path/'relocated'
    destination.mkdir()
    path.rename(destination/path.name)
    manifest_path.rename(destination/manifest_path.name)
    result, _ = load_pgen_features(destination/path.name, frame, frame)
    assert result.equals(frame.select('junction_aa', 'pgen'))
