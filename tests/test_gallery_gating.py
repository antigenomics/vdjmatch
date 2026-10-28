"""Check gallery selection without executing CLI workflows or notebook cells."""
import contextlib
import io
import os
from pathlib import Path
import runpy
import sys
import tempfile
from types import SimpleNamespace
from unittest.mock import patch


def test_gallery_scientific_opt_in():
    script = Path(__file__).resolve().parents[1] / 'examples/check_gallery.py'
    for scientific in (False, True):
        seen = []
        nbformat = SimpleNamespace(
            read=lambda path, **kwargs: seen.append(str(path)) or SimpleNamespace(metadata={}),
            write=lambda *args: None)
        nbclient = SimpleNamespace(NotebookClient=lambda *args, **kwargs: SimpleNamespace(execute=lambda: None))
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ), patch.dict(
                sys.modules, nbformat=nbformat, nbclient=nbclient), patch(
                'subprocess.run') as command, patch.object(
                sys, 'argv', [str(script), str(Path(folder)/'gallery')] + (['--scientific'] if scientific else [])), contextlib.redirect_stdout(io.StringIO()):
            try:
                runpy.run_path(str(script), run_name='__main__')
            except SystemExit as error:
                assert error.code == 0
            assert len(seen) == (6 if scientific else 2)
            assert command.call_count == (4 if scientific else 3)
            assert any('precursor.ipynb' in path for path in seen) == scientific
            assert any('vaccination.ipynb' in path for path in seen) == scientific
