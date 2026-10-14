"""Competitor branding is confined to canonical comparison notebook sources."""
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def test_repository_naming_boundary():
    policy = json.loads((ROOT/'examples/comparison-provenance.ipynb').read_text())
    pattern = re.compile('|'.join(policy['metadata']['forbidden_package_patterns']), re.I)
    names = subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard'],
        cwd=ROOT,text=True).splitlines()
    failures = []
    text_suffixes = {'.py','.md','.rst','.sh','.json','.js','.yml','.yaml','.toml','.txt','.ipynb'}
    for name in names:
        path = ROOT/name
        if not path.is_file():
            continue
        canonical = path.resolve().relative_to(ROOT)
        comparison = (canonical.parts[:2] == ('bench','comparisons') or
            str(canonical) in {'examples/comparison.ipynb','examples/comparison-provenance.ipynb'})
        if comparison and path.suffix == '.ipynb':
            notebook = json.loads(path.read_text())
            assert notebook['metadata'].get('comparison_notebook') is True
            continue
        if pattern.search(str(canonical)):
            failures.append(name+': filename')
        if path.suffix not in text_suffixes:
            continue
        if path.suffix == '.ipynb':
            notebook = json.loads(path.read_text())
            text = json.dumps({k:v for k,v in notebook.items() if k != 'cells'})
            for cell in notebook['cells']:
                text += ''.join(cell['source'])
                for output in cell.get('outputs',[]):
                    text += ''.join(output.get('text',[]))
                    data = output.get('data',{})
                    text += ''.join(data.get('text/plain',[]))
                    text += ''.join(data.get('text/html',[]))
        else:
            text = path.read_text()
        match = pattern.search(text)
        if match:
            failures.append(name+': '+match.group())
    assert not failures, '\n'.join(failures)


def test_cli_help_uses_only_package_names(capsys):
    import pytest
    from vdjmatch.cli.__main__ import main
    policy = json.loads((ROOT/'examples/comparison-provenance.ipynb').read_text())
    pattern = re.compile('|'.join(policy['metadata']['forbidden_package_patterns']), re.I)
    with pytest.raises(SystemExit) as stopped:
        main(['--help'])
    assert stopped.value.code == 0
    root_help=capsys.readouterr().out
    assert pattern.search(root_help) is None
    commands=re.search(r'\{([^{}]+)\}',root_help).group(1).split(',')
    for command in commands:
        with pytest.raises(SystemExit) as stopped:
            main([command, '--help'])
        assert stopped.value.code == 0
        assert pattern.search(capsys.readouterr().out) is None
