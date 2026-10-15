"""The example profiler must measure a process and surface command failure."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


def test_example_profile(tmp_path):
    pytest.importorskip('psutil')
    spec = importlib.util.spec_from_file_location('example_profile',
        Path(__file__).parents[2]/'examples/profile_cli.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.profile([sys.executable, '-c', 'import time; time.sleep(.2)'],
                            tmp_path/'success.log', timeout=5)
    assert result['seconds'] >= .2
    assert result['peak_tree_rss_GB'] > 0
    assert result['peak_live_threads'] >= 1
    assert result['completed'] is True and result['exit_code']==0
    assert result['cpu_peak_window_seconds']==1.0 and result['peak_sampled_cpu_cores'] is None
    assert (tmp_path/'success.resources.json').is_file()
    with pytest.raises(subprocess.CalledProcessError):
        module.profile([sys.executable, '-c', 'raise SystemExit(3)'],
                       tmp_path/'failure.log', timeout=5)
    with pytest.raises(subprocess.TimeoutExpired):
        module.profile([sys.executable, '-c', 'import time; time.sleep(10)'],
                       tmp_path/'timeout.log', timeout=.1)
    with pytest.raises(MemoryError,match='exceeds budget'):
        module.profile([sys.executable,'-c','import time; time.sleep(10)'],
                       tmp_path/'memory.log',timeout=5,max_rss_bytes=1)
    for name,kind in [('failure','CalledProcessError'),('timeout','TimeoutExpired'),('memory','MemoryError')]:
        record=json.loads((tmp_path/f'{name}.resources.json').read_text())
        assert record['completed'] is False and record['failure_kind']==kind
        assert record['exit_code']!=0
    observed=module.profile([sys.executable,'-c','import time; time.sleep(1.2)'],
                            tmp_path/'cpu-window.log',timeout=5)
    assert observed['peak_sampled_cpu_cores'] is not None
