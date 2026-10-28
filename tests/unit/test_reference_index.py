"""Prepared-reference assets preserve observations, native ordering and source identity."""
import json

import polars as pl
import pytest
from seqtree import SearchParams


def reference_file(path):
    pl.DataFrame({'cdr3': ['CASSLGQAYEQYF', 'CASSLGQAYEQYF', 'CASSLGCAYEQYF', 'CAVRDSNYQLIW'],
        'v': ['TRBV19', 'TRBV19', 'TRBV19', 'TRAV1-1'],
        'j': ['TRBJ2-7', 'TRBJ2-7', 'TRBJ2-7', 'TRAJ33'],
        'gene': ['TRB', 'TRB', 'TRB', 'TRA'], 'species': ['HomoSapiens']*4,
        'epitope': ['DEMO', 'OTHER', 'DEMO', 'DEMO'], 'record_id': ['r1', 'r2', 'r3', 'r4'],
        'complex_id': ['pair1', '0', '0', 'pair1'], 'mhc_a': ['HLA-A*02:01']*4,
        'mhc_b': ['B2M']*4, 'mhc_class': ['MHCI']*4}).write_csv(path, separator='\t')
    return path


def test_missing_source_does_not_create_sidecars(tmp_path):
    from vdjmatch.db.index import load_or_build_index
    source = tmp_path/'missing.tsv'
    with pytest.raises(FileNotFoundError):
        load_or_build_index(source)
    assert not source.with_name(source.name+'.vdjmatch').exists()


@pytest.mark.parametrize('allocation', [1, 2])
@pytest.mark.parametrize('entrypoint', ['index', 'api'])
def test_index_defaults_resolve_allocation(tmp_path, monkeypatch, allocation, entrypoint):
    from vdjmatch.api import Annotator
    from vdjmatch.db.index import load_or_build_index
    import vdjtools.cores as cores

    monkeypatch.setattr(cores, 'available_cores', lambda: allocation)
    resolved, resolver = [], cores.kernel_threads
    def observed(threads=0):
        result = resolver(threads)
        resolved.append((threads, result))
        return result
    monkeypatch.setattr(cores, 'kernel_threads', observed)
    load = Annotator.from_path if entrypoint == 'api' else load_or_build_index
    source = reference_file(tmp_path/'reference.tsv')
    load(source)
    assert resolved[0] == (0, allocation)
    with pytest.raises(ValueError, match='threads'):
        load(source, threads=allocation + 1)


def test_saved_native_roundtrip_preserves_keys_observations_and_skips_producers(tmp_path, monkeypatch):
    from vdjmatch.db.index import load_or_build_index
    from vdjmatch import db
    from vdjmatch.match.engine import VdjdbIndex
    import vdjmatch.aggregate.fast_profiles as profiles_module
    source = reference_file(tmp_path/'reference.tsv')
    first, profiles, manifest = load_or_build_index(source, threads=1)
    from pathlib import Path
    assert pl.read_csv(Path(manifest['artifact_path'])/'profiles.tsv', separator='\t').height == profiles.height
    assert first.reference.height == 4
    assert first.reference['record_id'].to_list() == ['r1', 'r2', 'r3', 'r4']
    def forbidden(*args, **kwargs):
        pytest.fail('unchanged prepared reference must not rerun a producer')
    monkeypatch.setattr(db, 'load', forbidden)
    monkeypatch.setattr(VdjdbIndex, 'build', forbidden)
    monkeypatch.setattr(profiles_module, 'reference_profiles', forbidden)
    second, restored, manifest2 = load_or_build_index(source, threads=1)
    assert manifest2['artifact_path'] == manifest['artifact_path']
    assert restored.equals(profiles)
    assert second.reference.equals(first.reference)
    parameters = SearchParams(max_subs=1, max_ins=0, max_dels=0)
    for gene in first.genes:
        before, keys, records = first._by_gene[gene]
        after, restored_keys, restored_records = second._by_gene[gene]
        assert restored_keys.equals(keys) and restored_records.equals(records)
        assert [after.ref_seq(i) for i in range(len(after))] == keys['cdr3'].to_list()
        assert [(h.ref_id, h.score) for h in before.search('CASSLGQAYEQYF', parameters)] == [
            (h.ref_id, h.score) for h in after.search('CASSLGQAYEQYF', parameters)]


@pytest.mark.parametrize('damage', ['missing', 'corrupt', 'manifest', 'native-format', 'profiles-tsv'])
def test_incomplete_or_corrupt_assets_rebuild(tmp_path, monkeypatch, damage):
    from pathlib import Path
    from vdjmatch.db.index import load_or_build_index
    from vdjmatch import db
    from vdjmatch.db.cache import sha256
    source = reference_file(tmp_path/'reference.tsv')
    expected, _, manifest = load_or_build_index(source, threads=1)
    bundle = Path(manifest['artifact_path'])
    member = bundle/'normalized-reference.parquet'
    if damage == 'missing':
        member.unlink()
    elif damage == 'corrupt':
        member.write_bytes(b'corrupt')
    elif damage == 'manifest':
        (bundle/'manifest.json').unlink()
    elif damage == 'profiles-tsv':
        (bundle/'profiles.tsv').unlink()
    else:
        member = bundle/'TRB.seqtree'
        member.write_bytes(b'bad native format')
        entry = json.loads((bundle/'manifest.json').read_text())
        entry['files'][member.name]['sha256'] = sha256(member)
        (bundle/'manifest.json').write_text(json.dumps(entry))
        (bundle/'manifest.json.sha256').write_text(sha256(bundle/'manifest.json'))
    calls, loader = [], db.load
    def tracked(*args, **kwargs):
        calls.append(1)
        return loader(*args, **kwargs)
    monkeypatch.setattr(db, 'load', tracked)
    rebuilt, _, _ = load_or_build_index(source, threads=1)
    assert calls == [1] and rebuilt.reference.equals(expected.reference)


def test_source_bytes_filters_and_directory_members_invalidate(tmp_path):
    from vdjmatch.db.index import load_or_build_index
    source = reference_file(tmp_path/'reference.tsv')
    _, _, original = load_or_build_index(source, threads=1)
    source.write_text(source.read_text().replace('OTHER', 'CHANGED'))
    _, _, changed = load_or_build_index(source, threads=1)
    selected, _, filtered = load_or_build_index(source, epitope='DEMO', threads=1)
    assert len({m['artifact_path'] for m in (original, changed, filtered)}) == 3
    assert selected.reference.height == 3
    directory = tmp_path/'tables'
    directory.mkdir()
    pl.read_csv(source, separator='\t').write_parquet(directory/'vdjdb.parquet')
    _, _, first = load_or_build_index(directory, threads=1)
    (directory/'source-notes.txt').write_text('source sidecar A')
    _, _, second = load_or_build_index(directory, threads=1)
    assert first['artifact_path'] != second['artifact_path']
    _, _, unchanged = load_or_build_index(directory, threads=1)
    assert second['artifact_path'] == unchanged['artifact_path']


def test_source_checksum_changes_do_not_bypass_existing_validation(tmp_path):
    from vdjmatch.db.index import load_or_build_index
    source = reference_file(tmp_path/'reference.tsv')
    load_or_build_index(source, threads=1)
    source.with_name(source.name+'.sha256').write_text('0'*64)
    with pytest.raises(ValueError, match='checksum mismatch'):
        load_or_build_index(source, threads=1)
