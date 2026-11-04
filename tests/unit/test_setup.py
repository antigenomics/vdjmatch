"""Bootstrap flags select declared extras and the package-owned runtime entrypoint."""
import os
from pathlib import Path
import subprocess

import pytest


@pytest.mark.parametrize('flag', [None, '--test', '--tests', '--unknown'])
def test_setup_flags_and_runtime_initialization(tmp_path, flag):
    script = Path(__file__).resolve().parents[2] / 'setup.sh'
    (tmp_path / 'setup.sh').write_text(script.read_text())
    tools = tmp_path / 'tools'
    tools.mkdir()
    python = tmp_path / '.venv/bin/python'
    python.parent.mkdir(parents=True)
    stub = '#!/bin/sh\nprintf "%s\\n" "$*" >> "$SETUP_LOG"\n'
    for path in [tools / 'uv', python]:
        path.write_text(stub)
        path.chmod(0o755)
    log = tmp_path / 'commands.txt'
    env = dict(os.environ, PATH=str(tools) + ':' + os.environ['PATH'], SETUP_LOG=str(log))
    result = subprocess.run(['sh', str(tmp_path / 'setup.sh'), *([flag] if flag else [])],
                            env=env, capture_output=True, text=True, timeout=10)
    if flag == '--unknown':
        assert result.returncode == 2 and not log.exists()
        return
    assert result.returncode == 0, result.stderr
    commands = log.read_text()
    tests = flag is not None
    assert ('test,control,bench,diagnostics' in commands) == tests
    assert ('import vdjtools' in commands) == tests
