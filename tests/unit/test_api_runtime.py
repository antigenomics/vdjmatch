"""Public annotation startup and outer-worker native budget regressions."""
import json
import os
import subprocess
import sys

import pytest


def test_api_initializes_within_scheduler_budget_without_loading_models():
    env = dict(os.environ, SLURM_CPUS_PER_TASK='2', SLURM_CPUS_ON_NODE='2')
    env.pop('POLARS_MAX_THREADS', None)
    code = '''
import json, sys
from vdjmatch import Annotator
import polars as pl
print(json.dumps({'threads': pl.thread_pool_size(),
                  'model': 'vdjtools.model' in sys.modules,
                  'numpy': 'numpy' in sys.modules}))
'''
    result = subprocess.run([sys.executable, '-c', code], env=env, text=True,
                            capture_output=True, check=True, timeout=20)
    observed = json.loads(result.stdout)
    assert observed['threads'] <= 2
    assert not observed['model']
    assert not observed['numpy']


def test_unique_search_resolves_outer_worker_budget_and_keeps_one_batch():
    from vdjtools.cores import serial_kernels
    from vdjmatch.match.engine import _search_unique_junctions
    calls = []

    class Index:
        def search_batch(self, keys, params, threads):
            calls.append((keys, threads))
            return [[key] for key in keys]

    with serial_kernels():
        assert _search_unique_junctions(Index(), ['CASSF', 'CASSF', 'CAVW'], None, 0) == [
            ['CASSF'], ['CASSF'], ['CAVW']]
        with pytest.raises(ValueError, match='threads'):
            _search_unique_junctions(Index(), ['CASSF'], None, 2)
    assert calls == [(['CASSF', 'CAVW'], 1)]


def test_public_default_uses_allocation_and_preserves_serial_calls(monkeypatch):
    import polars as pl
    from polars.testing import assert_frame_equal
    from vdjtools.cores import kernel_threads
    from vdjmatch import Annotator
    from vdjmatch.match.engine import VdjdbIndex

    reference = pl.DataFrame({'gene': ['TRB'], 'cdr3': ['CASSF'], 'v': ['TRBV1'],
                              'j': ['TRBJ1'], 'epitope': ['EXAMPLE']})
    ann = Annotator.from_frame(reference)
    observed = []
    original = VdjdbIndex.annotate

    def record(self, *args, **kwargs):
        observed.append(kwargs['threads'])
        return original(self, *args, **kwargs)

    monkeypatch.setattr(VdjdbIndex, 'annotate', record)
    queries = pl.DataFrame({'cdr3': ['CASSF', 'CASSF', 'CWWWF']})
    serial = ann.annotate(queries, scope='0', threads=1)
    automatic = ann.annotate(queries, scope='0')
    assert observed == [1, kernel_threads(0)]
    assert_frame_equal(serial, automatic)


def test_preinitialized_oversized_table_pool_fails_before_reference_build():
    code = '''
import polars as pl
assert pl.thread_pool_size() == 4
from vdjmatch import Annotator
try:
    Annotator.from_frame(pl.DataFrame({'gene': ['TRB'], 'cdr3': ['CASSF'], 'epitope': ['EXAMPLE']}))
except RuntimeError as exc:
    assert 'Polars' in str(exc) and 'fresh process' in str(exc)
else:
    raise AssertionError('oversized immutable Polars pool accepted')
'''
    env = dict(os.environ, POLARS_MAX_THREADS='4', SLURM_CPUS_PER_TASK='2', SLURM_CPUS_ON_NODE='2')
    subprocess.run([sys.executable, '-c', code], env=env, text=True,
                   capture_output=True, check=True, timeout=20)
