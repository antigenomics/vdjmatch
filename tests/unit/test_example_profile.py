"""The example profiler must measure a process and surface command failure."""
import importlib.util
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
    assert (tmp_path/'success.resources.json').is_file()
    with pytest.raises(subprocess.CalledProcessError):
        module.profile([sys.executable, '-c', 'raise SystemExit(3)'],
                       tmp_path/'failure.log', timeout=5)
    with pytest.raises(subprocess.TimeoutExpired):
        module.profile([sys.executable, '-c', 'import time; time.sleep(10)'],
                       tmp_path/'timeout.log', timeout=.1)
