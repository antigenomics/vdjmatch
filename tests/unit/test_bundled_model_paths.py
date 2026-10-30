"""Released upstream aliases must resolve to the files used by the loader."""
import json

import pytest

from vdjmatch.precursor.model import bundled_model_path


@pytest.mark.parametrize('source,family', [('olga', 'olga'), ('olga-legacy', 'olga'),
    ('arda', 'isalgo'), ('learned', 'isalgo'), ('isalgo', 'isalgo')])
def test_bundled_source_files(source, family):
    path = bundled_model_path('TRA', source)
    assert path.parent.name == family
    assert json.loads((path / 'manifest.json').read_text())['locus'] == 'TRA'


def test_mouse_source_files_and_missing_contexts():
    assert bundled_model_path('TRB', 'arda', 'mouse').name == 'mouse_TRB'
    with pytest.raises((ValueError, FileNotFoundError)):
        bundled_model_path('TRG', 'isalgo', 'mouse')
    with pytest.raises(ValueError):
        bundled_model_path('TRA', 'unknown')
