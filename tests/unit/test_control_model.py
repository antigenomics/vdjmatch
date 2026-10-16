"""Synthetic round-trip and integrity checks for stored observed-control models."""
import json

import numpy as np
import polars as pl
import pytest

from vdjmatch.cli.control_model import load_control_model
from vdjmatch.db.cache import sha256


def _prepare(tmp_path):
    from vdjmatch.cli.__main__ import main
    raw, directory = tmp_path/'raw.tsv', tmp_path/'model'
    pl.DataFrame({'junction_aa': ['CAAAAAAF', 'CAAAAAAF', 'CAGGGGGF', 'CAAAAAAF'],
        'v_call': ['TRBV1*01', 'TRBV1*01', 'TRBV2*01', 'TRBV1*02'],
        'j_call': ['TRBJ1-1']*4, 'occurrences': [3, 2, 5, 7]}).write_csv(raw, separator='\t')
    assert main(['prepare-control', '--control', str(raw), '--species', 'human',
        '--locus', 'TRB', '--threads', '1', '--output-dir', str(directory)]) == 0
    return raw, directory


def test_prepared_control_round_trip_matches_existing_reducers(tmp_path, monkeypatch):
    from vdjmatch.cli.combined import _raw_control
    from vdjmatch.evalue.mock_reference import control_events
    from vdjmatch.aggregate.structure import control_gene_usage
    from vdjmatch.aggregate.motif_information import motif_population
    raw, directory = _prepare(tmp_path)
    frame, _ = _raw_control(raw)
    expected = control_events(frame)
    expected_usage = control_gene_usage(expected, locus='TRB')
    modes = ('exact3', 'imgt3', 'wild4', 'imgt_wild4')
    expected_populations = {mode: motif_population(expected.rename({'junction_aa': 'cdr3'}),
        locus='TRB', mode=mode) for mode in modes}
    def forbidden(*args, **kwargs):
        pytest.fail('loading a stored control model must not rebuild its reducers')
    monkeypatch.setattr('vdjmatch.cli.combined._raw_control', forbidden)
    monkeypatch.setattr('vdjmatch.evalue.mock_reference.control_events', forbidden)
    monkeypatch.setattr('vdjmatch.aggregate.structure.control_gene_usage', forbidden)
    monkeypatch.setattr('vdjmatch.aggregate.motif_information.motif_population', forbidden)
    events, usage, populations, info = load_control_model(directory, raw,
        species='HomoSapiens', locus='TRB', modes=modes)
    assert events.equals(expected) and usage.equals(expected_usage)
    assert info['event_keys'] == 3 and info['total_events'] == 17
    assert info['control_model']['manifest_sha256'] == sha256(directory/'manifest.json')
    manifest = json.loads((directory/'manifest.json').read_text())
    assert manifest['raw']['sha256_before'] == manifest['raw']['sha256_after'] == sha256(raw)
    assert manifest['producer']['sources_before'] == manifest['producer']['sources_after']
    for mode, population in populations.items():
        assert population['vocabulary'].equals(expected_populations[mode]['vocabulary'])
        for scope in ('exact_length', 'all_lengths'):
            assert population[scope][0].equals(expected_populations[mode][scope][0])
            np.testing.assert_array_equal(population[scope][1], expected_populations[mode][scope][1])


def test_events_only_does_not_read_unused_motif_or_gene_files(tmp_path):
    raw, directory = _prepare(tmp_path)
    for path in directory.iterdir():
        if path.name not in ('events.parquet', 'manifest.json'):
            path.unlink()
    events, usage, populations, _ = load_control_model(directory, raw, species='human', locus='TRB', modes=())
    assert events.height == 3 and usage is None and populations == {}


@pytest.mark.parametrize('mismatch', ['raw', 'species', 'locus', 'mode', 'tamper', 'path'])
def test_prepared_control_rejects_mismatched_inputs(tmp_path, mismatch):
    raw, directory = _prepare(tmp_path)
    species, locus, modes = 'human', 'TRB', ('exact3',)
    if mismatch == 'raw':
        raw.write_text(raw.read_text()+'\n')
    elif mismatch == 'species':
        species = 'mouse'
    elif mismatch == 'locus':
        locus = 'TRA'
    elif mismatch == 'mode':
        modes = ('bad',)
    elif mismatch == 'tamper':
        with (directory/'events.parquet').open('ab') as stream:
            stream.write(b'changed')
    else:
        original = directory/'events.parquet'
        outside = tmp_path/'elsewhere.parquet'
        original.rename(outside)
        original.symlink_to(outside)
    with pytest.raises(ValueError):
        load_control_model(directory, raw, species=species, locus=locus, modes=modes)


@pytest.mark.parametrize('corruption', ['shape', 'schema'])
def test_prepared_control_rejects_rehashed_structural_corruption(tmp_path, corruption):
    raw, directory = _prepare(tmp_path)
    manifest_path = directory/'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if corruption == 'shape':
        path = directory/'exact3.counts.npz'
        np.savez_compressed(path, exact_length=np.zeros((1, 1), dtype=np.int64),
            all_lengths=np.zeros((1, 1), dtype=np.int64))
    else:
        path = directory/'exact3.exact_length.parquet'
        pl.read_parquet(path).drop('keys').write_parquet(path)
    manifest['files'][path.name]['sha256'] = sha256(path)
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        load_control_model(directory, raw, species='human', locus='TRB', modes=('exact3',))


def test_prepare_control_refuses_existing_directory_before_reading_raw(tmp_path):
    from vdjmatch.cli.__main__ import main
    existing = tmp_path/'existing'
    existing.mkdir()
    with pytest.raises(SystemExit) as exc:
        main(['prepare-control', '--control', str(tmp_path/'absent.tsv'), '--locus', 'TRB',
            '--output-dir', str(existing)])
    assert exc.value.code == 2 and not list(existing.iterdir())


def test_prepare_control_detects_raw_change_and_leaves_no_complete_manifest(tmp_path, monkeypatch):
    from vdjmatch.cli.__main__ import main
    import vdjmatch.aggregate.motif_information as motifs
    raw, directory = _prepare(tmp_path)
    destination = tmp_path/'drifted-model'
    original = motifs.motif_population
    def mutate(*args, **kwargs):
        raw.write_text(raw.read_text()+'\n')
        return original(*args, **kwargs)
    monkeypatch.setattr(motifs, 'motif_population', mutate)
    with pytest.raises(SystemExit) as exc:
        main(['prepare-control', '--control', str(raw), '--locus', 'TRB',
            '--output-dir', str(destination)])
    assert exc.value.code == 2 and not (destination/'manifest.json').exists()
    assert (directory/'manifest.json').exists()


def test_control_model_parser_import_keeps_numerical_kernels_lazy():
    import subprocess
    import sys
    subprocess.run([sys.executable, '-c',
        "import sys; from argparse import ArgumentParser; from vdjmatch.cli.control_model import register; "
        "register(ArgumentParser().add_subparsers()); assert 'numpy' not in sys.modules; "
        "assert 'polars' not in sys.modules"], check=True, timeout=10)
