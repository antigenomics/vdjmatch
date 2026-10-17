"""Notebook-backed runtimes retain imports, workers and installed tool files."""
import importlib.util
import json
from pathlib import Path
import stat
import subprocess
import sys

import pytest


def test_wrapper_imports_and_portable_export(tmp_path, monkeypatch):
    source_file = Path(__file__).parents[2]/'bench/comparison_loader.py'
    if not source_file.is_file():
        pytest.skip('comparison notebook sources are checkout-only')
    spec = importlib.util.spec_from_file_location('comparison_loader_test', source_file)
    monkeypatch.setattr(sys, 'meta_path', sys.meta_path.copy())
    loader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loader)
    root = tmp_path/'source'
    (root/'comparisons').mkdir(parents=True)
    tool = root/'external'/'tool'
    tool.parent.mkdir()
    tool.write_text('tool payload')
    tool.chmod(0o755)
    for name, source in [('entry', 'import helper\ndef run(): return helper.read()\n'),
            ('helper', 'from pathlib import Path\ndef read(): return (Path(__file__).parent/"external"/"tool").read_text()\n'),
            ('worker', 'print(42)\n')]:
        (root/'comparisons'/(name+'.ipynb')).write_text(json.dumps({
            'metadata': {'legacy_filename':name+'.py', 'program_cell_type':'code'},
            'cells':[{'cell_type':'code','metadata':{'program':True},'source':[source]}]}))
    monkeypatch.setattr(loader, 'ROOT', root)
    monkeypatch.setattr(sys, 'path', sys.path.copy())
    namespace = {'__name__':'entry', '__file__':str(root/'entry.py')}
    try:
        loader.load(namespace, 'entry')
        assert namespace['run']() == 'tool payload'
        worker = loader.program_path('worker.py')
        result = subprocess.run([sys.executable, str(worker)], capture_output=True,
            text=True, check=True, timeout=5)
        assert result.stdout.strip() == '42'
        with pytest.raises(ValueError):
            loader.program_path('../worker.py')
        destination = loader.export(tmp_path/'portable')
        assert not (destination/'external').is_symlink()
        assert (destination/'external'/'tool').stat().st_mode & stat.S_IXUSR
        tool.unlink()
        result = subprocess.run([sys.executable, '-c', 'import entry; print(entry.run())'],
            cwd=destination, capture_output=True, text=True, check=True, timeout=5)
        assert result.stdout.strip() == 'tool payload'
    finally:
        sys.modules.pop('helper', None)
        if loader._worker_directory is not None:
            loader._worker_directory.cleanup()
