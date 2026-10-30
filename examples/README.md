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
sh examples/database.sh NEW_DATABASE_DIRECTORY
# Use an existing VDJdb primary build without downloading or rebuilding it:
sh examples/annotate.sh NEW_BUILD_DIRECTORY /path/to/build/tables/vdjdb.parquet
```

Keep this complete directory together. The scripts find their raw inputs beside
themselves, use the offline reference unless a second annotation argument supplies a
local primary or legacy table (or an AIRR bundle with receptor species metadata), default to four native threads and fail if
the output directory exists. Every invocation calculates fresh CLI tables;
there are no feature caches, predictions as inputs, skip or resume paths.

Annotation defaults to the curated manuscript reference population in 0.6.0.
Pass `--reference-population all-valid` for the accepted all-valid baseline;
keep that population label with every output, curve and interval. The stem4
formula is unchanged, and these tiny fixtures do not establish manuscript accuracy.

The annotation notebook covers a four-row metadata cohort and two wide pairs,
including an incomplete pair. The comparison notebook calls the CLI with exact
exclusion and detailed default stem4 / advanced vcdr12 settings, then joins its
continuous evidence on source identity and complete pMHC keys. Open either
notebook from the repository root or this directory. No outputs are committed.

`database.sh` checks three synthetic queries against the same offline reference
and one supplied synthetic raw control. Exact/one-edit/no-hit queries return
`p_both` values 1.0/1.0/1.0; exact presence is a separate boolean. These
small fixtures explain the calculation, not a biological null law.

All inputs were written specifically as deterministic software fixtures, not
sampled from a biological dataset. `DEMO_A`/`DEMO_B` are demonstration labels,
not experimental binding truth. Literal junctions include conserved anchors.
There is no RNG or fit and no accuracy or calibration claim.

Offline reference SHA256: `3ed6c1baafe2ec51bac658a02fa806d4a866ffcffc1ff7e1a391188050bcb76d`.
The Git version of this directory pins the raw inputs; CLI manifests record the
reference and scoring configuration of each fresh invocation.

## Verify the gallery against this checkout

From the repository root, install the matching source in an isolated environment
and run the offline shell workflows and annotation/comparison notebooks:

```sh
uv venv .venv-gallery
uv pip install --python .venv-gallery/bin/python -r examples/requirements-notebooks.txt
.venv-gallery/bin/python examples/check_gallery.py NEW_GALLERY_OUTPUT_DIRECTORY
```

The checker creates a kernel for that interpreter, keeps executed notebooks and
the HTML export in the new output directory, and exits nonzero on any failure.
Use the same Git checkout for source, scripts and documentation. Released
packages need the examples from their matching release tag.
Shell examples have a 120-second timeout; ordinary notebook cells and the marimo
export have a 300-second timeout. The full-cohort vaccination calculation has
an explicit 1,800-second cell budget. With `--scientific`, the precursor tutorial downloads the pinned public
Hugging Face VDJdb snapshot and uses four native threads. With `--scientific`, the compendium downloads
its documented public tables and reference. The default check requires no downloads. Reference downloads stay in the
output directory; computed results are rebuilt on every invocation.

## Biological examples

Run the three-score accuracy comparison entirely from the terminal:

```sh
python examples/accuracy.py --benchmark /path/to/airr_benchmark \
  --reference reference.zip --germline-background raw-prior.airr.tsv \
  --sorted-input original-sorted.tsv --threads 4 --output accuracy-run
```

The optional sorted input preserves original LLW/LLL labels and V/J calls.
Without it, the contrast uses explicitly labelled mirror gene calls. The raw
germline background supplies original scoring's junction/V/J prior population.
Outputs include all three source-row scores, exact counts, ROC-AUC/AP, individual
opposite-class rank losses, TP/FP trajectories, figures and per-stage resources.
The complete mirror control tasks remain separate from the sorted contrast.

| Notebook | Raw inputs | Figures |
| --- | --- | --- |
| [accuracy.ipynb](accuracy.ipynb) | Complete dataset 1/2 TRB panels, pinned VDJdb | ROC/PR, score ECDFs, exact TP/FP trajectories |
| [vaccination.ipynb](vaccination.ipynb) | All local `airr_yfv19` repertoires and metadata | LLW evidence ranks, complete read coverage and neighbourhood support |
| [precursor_estimates.ipynb](precursor_estimates.ipynb) | Complete TRB sets for GIL, NLV and LLW | Model frequency, expected cells, occupancy |

Set `AIRR_BENCHMARK` to the complete benchmark mirror (default
`~/hf/airr_benchmark`) and `AIRR_YFV19` to the vaccination data directory
(default `~/hf/airr_yfv19`). The vaccination runner reads every local AIRR
repertoire named in `metadata.txt`, preserving donors, days and technical
replicates. All donors have the author-declared A*02 restriction. It scores every
eligible TRB V/junction/J key in one batch and restores complete replicate read
counts. Original sample identifiers and mappings stay private. No biological
inputs are copied into this repository.
The producer's full ``junction_aa`` is retained whenever its amino acids are
valid; noncanonical anchor residues are flagged rather than discarded.
The inputs do not report receptor functionality.

```sh
python examples/vaccination.py --dataset /path/to/airr_yfv19 \
  --reference reference.zip --threads 4 --memory-gb 24 --output vaccination-run
```

`coverage.tsv` records retained/dropped rows and read denominators;
`observations.tsv` preserves replicate counts, LLW scores and query dispositions.
`support.tsv` counts radius and exact-junction support. The evidence-ranked
coverage curves use complete score-tied groups; no abundance cutoff or binding
threshold is fitted. The CLI has a 30-minute timeout and a sampled RAM ceiling
(24 decimal GB by default, at most 32 GB).

The notebooks use the installed CLI for every core calculation, one native
batch per population, four native threads and one-thread numerical helper pools.
Their final tables are regenerated in a new temporary directory on every run.
`profile_cli.py` measures whole-CLI elapsed time, sampled process-tree RSS, CPU
seconds and live threads; install the `bench` extra for its existing `psutil`
dependency. On macOS, process inspection requires an ordinary terminal/kernel
with permission to inspect its own processes. A restricted sandbox may reject it.
See the [resource measurements](https://docs.isalgo.dev/vdjmatch/performance.html).

The accuracy and precursor notebooks retain reviewed aggregate outputs. The
full-cohort vaccination notebook is a workflow draft: its four-worker run has
not met the 30-minute budget and has no completed cohort outputs. To execute
the biological workflows alongside the original gallery in a fresh environment:

```sh
.venv-gallery/bin/python examples/check_gallery.py NEW_OUTPUT_DIRECTORY --scientific
```

Use a new directory for each verification. The checker records cell timing and
keeps complete execution outputs private; docs builds never execute notebooks.
