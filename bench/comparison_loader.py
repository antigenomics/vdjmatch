"""Read comparison programs from their canonical notebook sources."""
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parent


def programs():
    for path in sorted((ROOT/'comparisons').glob('*.ipynb')):
        notebook = json.loads(path.read_text())
        yield path, notebook


def load(namespace, name):
    """Keep unbranded benchmark helper imports working from notebook-owned code."""
    for path, notebook in programs():
        if notebook['metadata'].get('legacy_filename') == name+'.py':
            source = ''.join(''.join(c['source'])+'\n' for c in notebook['cells']
                             if c['cell_type'] == 'code' and c['metadata'].get('program'))
            exec(compile(source, str(path), 'exec'), namespace)
            return
    raise ImportError('comparison notebook not found: '+name)


def export(destination):
    """Export a fresh executable comparison bundle; no prediction is computed."""
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    for path in ROOT.iterdir():
        if path.is_file() and path.suffix in ('.py','.js','.md'):
            shutil.copyfile(path,destination/path.name)
    for _, notebook in programs():
        metadata = notebook['metadata']
        original = metadata['legacy_filename']
        if Path(original).name != original:
            raise ValueError('comparison program filename must have no directory')
        kind = metadata['program_cell_type']
        source = ''.join(''.join(c['source'])+'\n' for c in notebook['cells']
                         if c['cell_type'] == kind and c['metadata'].get('program'))
        (destination/original).write_text(source)
    return destination


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, help='new comparison runtime directory')
    args = parser.parse_args()
    print(export(args.output))
