# CLI-first worked examples

Run from a fresh installation with `vdjmatch` on PATH:

Browse the [examples gallery](https://docs.isalgo.dev/vdjmatch/examples.html) and
[notebooks gallery](https://docs.isalgo.dev/vdjmatch/notebooks.html).

| Workflow | Read online | Notebook source | Shell script |
| --- | --- | --- | --- |
| Cohort and paired annotation | [Preview](https://docs.isalgo.dev/vdjmatch/notebooks/annotation.html) | [annotation.ipynb](annotation.ipynb) | [annotate.sh](annotate.sh) |
| Geometry comparison | [Preview](https://docs.isalgo.dev/vdjmatch/notebooks/comparison.html) | [comparison.ipynb](comparison.ipynb) | [compare.sh](compare.sh) |

```sh
sh examples/annotate.sh NEW_ANNOTATION_DIRECTORY
sh examples/compare.sh NEW_COMPARISON_DIRECTORY
```

Keep this complete directory together. The scripts find their raw inputs beside
themselves, use the offline reference, default to four native threads and fail if
the output directory exists. Every invocation calculates fresh CLI tables;
there are no feature caches, predictions as inputs, skip or resume paths.

The annotation notebook covers a four-row metadata cohort and two wide pairs,
including an incomplete pair. The comparison notebook calls the CLI with exact
exclusion and detailed default stem4 / advanced vcdr12 settings, then joins its
continuous evidence on source identity and complete pMHC keys. Open either
notebook from the repository root or this directory. No outputs are committed.

All inputs were written specifically as deterministic software fixtures, not
sampled from a biological dataset. `DEMO_A`/`DEMO_B` are demonstration labels,
not experimental binding truth. Literal junctions include conserved anchors.
There is no RNG or fit and no accuracy or calibration claim.

Offline reference SHA256: `3ed6c1baafe2ec51bac658a02fa806d4a866ffcffc1ff7e1a391188050bcb76d`.
The Git version of this directory pins the raw inputs; CLI manifests record the
reference and scoring configuration of each fresh invocation.

## Verify the complete gallery against PyPI

From the repository root, create an isolated environment and run every shell
example, all three Jupyter notebooks and the marimo compendium:

```sh
uv venv .venv-gallery
uv pip install --python .venv-gallery/bin/python --index-url https://pypi.org/simple --no-cache -r examples/requirements-notebooks.txt
.venv-gallery/bin/python examples/check_gallery.py NEW_GALLERY_OUTPUT_DIRECTORY
```

The checker creates a kernel for that interpreter, keeps executed notebooks and
the HTML export in the new output directory, and exits nonzero on any failure.
Shell examples have a 120-second timeout; notebook cells and the marimo export
have a 300-second timeout. The precursor tutorial downloads the pinned public
Hugging Face VDJdb snapshot and uses four native threads. The compendium downloads
its documented public tables and reference. Reference downloads stay in the
output directory; computed results are rebuilt on every invocation.

## Biological examples

| Notebook | Raw inputs | Figures |
| --- | --- | --- |
| [accuracy.ipynb](accuracy.ipynb) | Complete dataset 1/2 TRB panels, pinned VDJdb | ROC/PR, score ECDFs, exact TP/FP trajectories |
| [vaccination.ipynb](vaccination.ipynb) | One baseline `airr_yfv19` repertoire | Abundance/coverage and read-weighted LLW evidence |
| [precursor_estimates.ipynb](precursor_estimates.ipynb) | Complete TRB sets for GIL, NLV and LLW | Model frequency, expected cells, occupancy |

Set `AIRR_BENCHMARK` to the complete benchmark mirror (default
`~/hf/airr_benchmark`) and `AIRR_YFV19` to the vaccination data directory
(default `~/hf/airr_yfv19`). `YFV_INPUT` chooses a specific raw vaccination
file; without it the notebook chooses the first table in sorted filename order.
The selected table must be baseline for the baseline interpretation. Original
sample identifiers and source mappings stay private. No biological inputs are
copied into this repository. The vaccination lesson takes the 10,000 most
abundant unique V/junction/J clonotypes and reports its complete-read denominator;
that subset is not an unbiased repertoire sample.

The notebooks use the installed CLI for every core calculation, one native
batch per population, four native threads and one-thread numerical helper pools.
Their final tables are regenerated in a new temporary directory on every run.
`profile_cli.py` measures whole-CLI elapsed time, sampled process-tree RSS, CPU
seconds and live threads; install the `bench` extra for its existing `psutil`
dependency. On macOS, process inspection requires an ordinary terminal/kernel
with permission to inspect its own processes. A restricted sandbox may reject it.
See the [resource measurements](https://docs.isalgo.dev/vdjmatch/performance.html).

The committed biological notebooks retain reviewed figure outputs and aggregate counts. To
execute these alongside the original gallery in a fresh isolated environment:

```sh
.venv-gallery/bin/python examples/check_gallery.py NEW_OUTPUT_DIRECTORY --scientific
```

Use a new directory for each verification. The checker records cell timing and
keeps complete execution outputs private; docs builds never execute notebooks.
