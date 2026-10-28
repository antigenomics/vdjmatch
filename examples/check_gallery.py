"""Run shell examples and gallery notebooks with this checkout's installed package.

From the repository root: python examples/check_gallery.py NEW_OUTPUT_DIRECTORY
Install examples/requirements-notebooks.txt in an isolated environment first.
"""
import json
import argparse
import os
from pathlib import Path
import subprocess
import sys
import time

# Budget numerical pools before imports and before starting any notebook kernel.
for name in ("POLARS_MAX_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
             "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

import nbformat
from nbclient import NotebookClient

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('output', type=Path)
parser.add_argument('--scientific', action='store_true', help='also run biological notebooks from local raw mirrors')
args = parser.parse_args()
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=False)
os.environ["VDJMATCH_CACHE"] = str(output / "reference")
os.environ["HF_HOME"] = str(output / "huggingface")
os.environ["MPLCONFIGDIR"] = str(output / "matplotlib")
os.environ["VDJMATCH"] = str(Path(sys.executable).with_name("vdjmatch"))
os.environ["PATH"] = str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]
# A local kernel specification ensures that Jupyter uses the same installed package.
kernel = output / "jupyter/kernels/gallery"
kernel.mkdir(parents=True)
(kernel / "kernel.json").write_text(json.dumps({
    "argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
    "display_name": "Gallery verification", "language": "python",
}))
os.environ["JUPYTER_PATH"] = str(output / "jupyter")
results = []


def run(name, action):
    start = time.perf_counter()
    record = {"workflow": name, "status": "passed"}
    try:
        action()
    except Exception as error:
        record.update(status="failed", error=str(error))
    record["seconds"] = round(time.perf_counter() - start, 3)
    results.append(record)
    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(name, record["status"], record["seconds"], flush=True)


for name in ("annotate", "compare", "database"):
    run(name, lambda name=name: subprocess.run(
        ["sh", str(root / f"examples/{name}.sh"), str(output / name)],
        cwd=root, check=True, timeout=120))


def notebook(path):
    nb = nbformat.read(root / path, as_version=4)
    try:
        timeout=nb.metadata.get('execution_timeout_seconds',300)
        if not isinstance(timeout,int) or not 0<timeout<=1800:
            raise ValueError('notebook execution budget must be between 1 and 1800 seconds')
        NotebookClient(nb, timeout=timeout, kernel_name="gallery",
                       resources={"metadata": {"path": str(root)}}, record_timing=True,
                       on_cell_start=lambda cell, cell_index: print(
                           path, "cell", cell_index + 1, flush=True)).execute()
    finally:
        nbformat.write(nb, output / (Path(path).stem + ".ipynb"))


for path in ("examples/annotation.ipynb", "examples/comparison.ipynb"):
    run(path, lambda path=path: notebook(path))
if args.scientific:
    for name in ('accuracy', 'vaccination', 'precursor_estimates'):
        run(f'examples/{name}.ipynb', lambda name=name: notebook(f'examples/{name}.ipynb'))
    run("docs/notebooks/precursor.ipynb", lambda: notebook("docs/notebooks/precursor.ipynb"))
    # Export a source-identical copy so marimo session files stay with this fresh run.
    compendium = output / "precursor_compendium.py"
    compendium.write_bytes((root / "docs/notebooks/precursor_compendium.py").read_bytes())
    run("precursor_compendium", lambda: subprocess.run(
        [sys.executable, "-m", "marimo", "export", "html",
         str(compendium),
         "-o", str(output / "precursor_compendium.html")],
        cwd=root, check=True, timeout=300))
sys.exit(any(record["status"] == "failed" for record in results))
