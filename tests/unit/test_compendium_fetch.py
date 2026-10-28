"""Reference-table bootstrap integrity; all downloads are mocked."""
import sys
from types import SimpleNamespace

import pytest
import polars as pl

from vdjmatch.db.cache import valid
from vdjmatch.precursor.compendium import _fetch, _NUMERIC, _HF_REVISION


@pytest.fixture
def downloads(tmp_path, monkeypatch):
    calls = []

    def download(**kwargs):
        calls.append(kwargs)
        source = tmp_path / 'download.tsv'
        pl.DataFrame({**{name: [1.] for name in _NUMERIC if name != 'n_ball'},
            'freq_per_1e6_naive': [1.], 'freq_per_1e6_published': [1.], 'naive': [1],
            'design': ['measured'], 'species': ['human'], 'epitope_seq': ['PEP'],
            'lab': ['synthetic'], 'source_type': ['primary'], 'chain': ['TRB'],
            'source': [kwargs['repo_id']]}).write_csv(source, separator='\t')
        return str(source)

    monkeypatch.setitem(sys.modules, 'huggingface_hub', SimpleNamespace(hf_hub_download=download))
    return calls


def test_fetch_separates_repositories_and_reuses_verified_file(tmp_path, downloads):
    cache = tmp_path / 'cache'
    first = _fetch('compendium', cache, repo='owner/first')
    second = _fetch('compendium', cache, repo='owner/second')
    assert first != second
    assert pl.read_csv(first, separator='\t')['source'].to_list() == ['owner/first']
    assert pl.read_csv(second, separator='\t')['source'].to_list() == ['owner/second']
    assert valid(first) and valid(second)
    assert _fetch('compendium', cache, repo='owner/first') == first
    assert len(downloads) == 2
    assert downloads[0]['revision'] == 'main'
    assert downloads[0]['cache_dir'] == str(cache / 'hf_hub')


def test_fetch_pins_default_revision_and_separates_overrides(tmp_path, downloads):
    first = _fetch('compendium', tmp_path)
    second = _fetch('compendium', tmp_path, revision='other-revision')
    assert first != second
    assert [call['revision'] for call in downloads] == [_HF_REVISION, 'other-revision']


def test_malformed_refresh_preserves_complete_file(tmp_path, downloads, monkeypatch):
    path = _fetch('compendium', tmp_path)
    previous = path.read_bytes()
    malformed = tmp_path / 'malformed.tsv'
    malformed.write_text('wrong\nvalue\n')
    monkeypatch.setitem(sys.modules, 'huggingface_hub', SimpleNamespace(
        hf_hub_download=lambda **kwargs: malformed))
    with pytest.raises(ValueError, match='missing columns'):
        _fetch('compendium', tmp_path, force=True)
    assert path.read_bytes() == previous
    assert valid(path)


@pytest.mark.parametrize('loader', ['load_compendium', 'load_estimates'])
def test_public_loaders_forward_revision(tmp_path, downloads, loader):
    from vdjmatch.precursor import compendium
    getattr(compendium, loader)(cache=tmp_path, revision='explicit-revision')
    assert downloads[0]['revision'] == 'explicit-revision'


def test_estimates_accepts_missing_optional_n_ball(tmp_path, downloads):
    from vdjmatch.precursor.compendium import load_estimates
    frame = load_estimates(cache=tmp_path)
    assert frame.height == 1
    assert 'n_ball' not in frame.columns
    assert valid(_fetch('estimates', tmp_path))


@pytest.mark.parametrize('damage', ['content', 'checksum'])
def test_fetch_retries_corrupted_file(tmp_path, downloads, damage):
    path = _fetch('estimates', tmp_path / 'cache')
    if damage == 'content':
        path.write_text('partial')
    else:
        path.with_name(path.name + '.sha256').unlink()
    assert not valid(path)
    assert _fetch('estimates', tmp_path / 'cache') == path
    assert valid(path)
    assert len(downloads) == 2
