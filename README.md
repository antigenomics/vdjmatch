<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/vdjmatch_dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="assets/vdjmatch_light.svg">
    <!-- Absolute PNG fallback: PyPI strips <picture>/<source> and cannot render a relative or
         raw-served SVG, so the logo must be an absolute-URL raster here. GitHub uses the SVG sources.
         Use the light (dark-on-transparent) raster so it reads on PyPI's white page. -->
    <img alt="vdjmatch" src="https://raw.githubusercontent.com/antigenomics/vdjmatch/master/assets/vdjmatch_light.png" width="340">
  </picture>
</p>

<h1 align="center">vdjmatch — TCR peptide–MHC annotation</h1>

<p align="center">
  <a href="https://pypi.org/project/vdjmatch/"><img alt="PyPI" src="https://img.shields.io/pypi/v/vdjmatch"></a>
  <a href="https://github.com/antigenomics/vdjmatch/actions/workflows/tests.yml"><img alt="tests" src="https://github.com/antigenomics/vdjmatch/actions/workflows/tests.yml/badge.svg"></a>
  <a href="https://antigenomics.github.io/vdjmatch/"><img alt="docs" src="https://github.com/antigenomics/vdjmatch/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-GPLv3-green"></a>
</p>

`vdjmatch` annotates selected T-cell receptors and repertoires against
[VDJdb](https://github.com/antigenomics/vdjdb-db), preserving query identities and
reporting peptide–MHC candidates, competing evidence and optional background calibration.
Native sequence search uses [`seqtree`](https://github.com/antigenomics/seqtree);
repertoire ingestion reuses published `vdjtools` readers.

The Python implementation follows the legacy Java/Groovy package, preserved on the
`legacy-java` branch. See the [tutorial](https://antigenomics.github.io/vdjmatch/tutorial.html),
[workflow guide](https://antigenomics.github.io/vdjmatch/how-to.html) and
[score explanation](https://antigenomics.github.io/vdjmatch/explanation.html).

## Features

- Latest/pinned GitHub releases and offline legacy, rich or AIRR reference inputs.
- Source-aware junction ingestion, stable row identities and explicit ingestion counts.
- Native batched searches by locus, with retained reference observations.
- All peptide–MHC candidates, NED-v1 ranking and competing-label support.
- Optional finite-background calibration with a stated counting unit and search predicate.
- Paired-reference evidence requiring both chains in the same reference complex.
- CLI hit/candidate/call/summary/ingestion tables and a reproducibility manifest.
- Separate precursor-frequency APIs and their established worked examples.

### Precursor frequency

How much repertoire mass can see an epitope, and how many cognate clonotypes exist — seen and unseen.
Use the documented precursor setup (`pip install 'vdjmatch[precursor]'`) for these examples:

```console
$ vdjmatch precursor --vdjdb --mhc-class MHCI --min-junctions 10 \
      --n-eff 1e8 --selection auto -o precursor.txt
```

```python
from vdjmatch import precursor as P

model = P.load_model("TRB")
junctions = ["CASSIRSSYEQYF", "CASSLGQAYEQYF", ...]     # junctions, not IMGT CDR3s

P.union_mass(model, junctions)      # exact union of the 1-mm balls, and the overlap a sum invents
P.occupancy(model, junctions, n_eff=1e8, selection=P.SELECTION_BY_CHAIN["TRB"])
# -> {"S": ..., "F": ..., "n_seen": ..., "n_unseen": ..., "seen_fraction": ...}
```

`S` is the effective cognate-set size, `F` the fraction of the naive repertoire, and `n_seen`
saturates in depth — a cognate set concentrated on a few high-`Pgen` junctions is exhausted at
shallow sequencing, a broad one keeps accumulating clonotypes. `SELECTION_BY_CHAIN` carries the
measured per-chain depth factors (TRB 4.62, TRA 1.07); they differ and should not be pooled.

Two things worth knowing before trusting a number:

- **`--source` picks the recombination model, and it matters for a mass over a set.** The default
  `olga` is a bit-faithful import of OLGA's published models — faithful to OLGA's deletion-bin grid
  too, under which **5.7% of human TRA junctions in VDJdb score `Pgen` exactly zero** (up to 15% for
  some epitopes) and vanish from the total without an error. `learned` and `arda` are refits that do
  not inherit it; `arda` is the only set carrying mouse. The `n_zero_pgen` output column reports the
  loss per group and the CLI warns above 1%.
- **`--species` selects records, `--organism` selects the model, and they must agree.** A mismatch
  used to run and return a plausible number that was wrong by an epitope-dependent factor; since
  0.3.0 it is refused.

#### Checking it against what laboratories measured

`F(e)` estimates a **naive repertoire mass**, and that quantity has been measured directly — by
tetramer and dextramer enrichment of unexposed donors, cord blood and naive mice. Both sides of that
comparison ship, and are deliberately kept apart:

```python
from vdjmatch.precursor import load_compendium, load_estimates

comp = load_compendium()                                  # experimental: 147 published measurements,
                                                          # 14 studies, 31 originating laboratories
est = load_estimates(species="human", chain="TRB", model="arda")   # derived: what the model predicts
est.select("epitope_seq", (est["union"] * 1e6).alias("predicted"), "measured_per_1e6")
```

Both fetch from `isalgo/airr_benchmark` (`vdjmatch/precursor_freq/`) and cache next to the VDJdb
releases. Merging them into one column would conflate a measurement with a prediction of it.

Beyond a point estimate, the same two inputs — the `Pgen` spectrum and the catalogued junction count
— answer more than one question, and the answers are not interchangeable:

```python
u = P.union_mass(model, junctions, r=1, count_members=True)
u["n_union"] - u["n_seqs"]                     # uncatalogued neighbours: a census, finite and exact
P.observed_mass(model, junctions) / u["union"] # share of the neighbourhood's mass already catalogued
u["overlap"]                                   # what naive per-sequence summing would double-count
```

**Unseen mass converges; unseen count does not.** The Horvitz–Thompson weight `p/π → 1/N` as
`p → 0`, which is what makes `coverage_corrected_mass` work without guessing the shape of the tail.
But `1/π` diverges, so a richness extrapolation on real TCR data returns 10¹¹–10¹⁴ unseen junctions
— not a number. `unseen_junctions` marks `richness_reliable` false when the rarest observed junction
was captured too seldom, and the ball census above is the finite alternative.

An interactive [marimo](https://marimo.io) notebook walks the whole comparison — the scatter against
measurement, the `Pgen` spectrum inside one cognate set, the ball census, and the step from a mass to
a number of cells:

```console
$ pip install 'vdjmatch[precursor,notebook]'
$ marimo edit docs/notebooks/precursor_compendium.py
```

The method, its benchmarks and the paper live outside this repository:
[repseq/2026-precursor-freq](https://github.com/repseq/2026-precursor-freq) (benchmarks, result
tables, the literature compendium) and
[repseq/2026-precursor-freq-ms](https://github.com/repseq/2026-precursor-freq-ms) (the manuscript).
This repository holds only the software, its CLI, tests, docs and the worked examples
(`docs/notebooks/`).

## Install

```console
python -m pip install vdjmatch
```

For development, use an isolated environment:

```console
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[test,bench]"
```

## Python API

```python
import polars as pl
from vdjmatch import Annotator

ann = Annotator.latest()  # GitHub latest; .version(tag) pins a release
# ann = Annotator.from_path("reference.zip")  # existing local input, offline
queries = pl.DataFrame({
    "sequence_id": ["query-a", "query-b"],
    "junction_aa": ["CASSIRSSYEQYF", "CWWWWF"],
    "locus": ["TRB", "TRB"],
})
calls = ann.annotate(queries, threads=1)  # every input row remains
candidates = ann.candidates(queries, threads=1)  # all peptide–MHC alternatives
calibrated = ann.annotate(queries, calibrate=True, species="human", threads=1)
```

AIRR `junction_aa` contains the conserved anchors; AIRR `cdr3_aa` excludes them.
True AIRR input without a junction is rejected. Legacy VDJdb `cdr3` can already
contain the junction. Custom ingestion requires an explicit junction convention;
the software never fabricates anchors.

The API uses uncalibrated ranking unless `calibrate=True` or an explicit control
index is supplied. NED ranks candidates; `E` is an expected background match count,
and `p_enrichment` tests the stated reference set. None is a posterior probability
of the winning specificity. Exact ranking ties abstain as `ambiguous`.

`ann.annotate_paired(cell_df, cdr3a="cdr3a", cdr3b="cdr3b")` requires both chains
to match the same reference complex and peptide–MHC label. Paired background
calibration states an independent-chain null explicitly. Agreement between
separate alpha and beta top labels alone does not establish paired-reference evidence.

## Command line

```console
vdjmatch update --asset default --cache reference-inputs
vdjmatch match --input-format airr --threads 1 -o results/run sample.tsv
vdjmatch match --vdjdb reference.zip --no-evalue -o results/offline sample.tsv
vdjmatch match --vdjdb reference.zip --paired --link cell_id -o results/cells cells.tsv
vdjmatch match --vdjdb reference.zip --search-mode ball -o results/balls sample.tsv
vdjmatch match --vdjdb reference.zip --search-mode ranked --top-k 10 -o results/ranked sample.tsv
```

CLI matching requests calibration by default; `--no-evalue` selects uncalibrated
ranking. The command emits hits, candidates, row-preserving calls, descriptive
summaries and ingestion diagnostics as TSV, plus a JSON manifest.

Use `--input-format custom --sequence-convention junction` for an ambiguous custom
table. `--scope s,i,d,t` sets substitution/insertion/deletion/total budgets. Hard
`--match-v` or `--match-j` requests require `--no-evalue` with sequence-only backgrounds.
Native search defaults to one thread. See `vdjmatch match -h` and the
[CLI reference](https://antigenomics.github.io/vdjmatch/cli.html).

Ranked mode scans the complete eligible reference with a native restricted gap-block scorer,
retaining at most K unique junctions or genuine paired keys before expanding observations.
Its global fixed-K enrichment test is separate from descriptive epitope ranks. Ball mode
reports radii1–5 separately. Neither selects the best K/radius after seeing the results.
These new modes require seqtree1.0.2; see the CLI reference for scoring and null assumptions.

## Scoring: what works (and what doesn't)

An empirical study on VDJdb (see `appendix/vdjmatch_scoring.tex`; regenerated on the **2026-06-11-ZENODO**
release with composition controls and balanced metrics) settles the scoring question honestly:

- **Hamming distance 1 is the signal:noise optimum** — macro purity (per-epitope mean) falls
  0.49 → 0.07 across edit distance 1–5, a 56× → 2.5× enrichment over chance (reproducing the original
  VDJdb observation, Shugay et al. NAR 2018). The search-ball radius, not the substitution matrix, is
  the dominant lever. (On the dense 2026 release this only shows up once a few 10× mega-studies are
  capped and the random tail is treated as an admixed control — naive pooled purity reads a flat ~0.9.)
- **Central (NDN) substitutions carry the specificity signal** — a mismatch in the CDR3 core most
  often changes specificity (P(same epitope) ≈ 0.31) while near-anchor mismatches are germline noise
  (≈ 0.75); the NDN core is also ~31–34% glycine (insertion/D-gene signature).
- **No amino-acid matrix clearly beats BLOSUM62 — and a genetic-code null ties it.** BLOSUM62 ≈ PAM250 ≈
  structural > Hamming > data-derived VDJAM. Strikingly, **VDJAMr** — a matrix built from the *genetic
  code alone* (how mutationally accessible one AA is from another; `loo_vdjam.codon_dissim`) — matches
  BLOSUM62 (0.581 vs 0.564 @≤2), so TCR CDR3 substitution structure is **generative, not chemical**.
  A published TCR-specific matrix, [tcrBLOSUM](https://doi.org/10.1093/bib/bbae602), does *not* beat
  BLOSUM62 either (0.555 vs 0.567; `bench/tcrblosum_refute.py`) — its same-epitope counts don't transfer.
- **Position-weighting BLOSUM62 does beat it.** Encoding the central-substitution finding as a
  seqtree positional matrix (`PositionalMatrix.from_weights(BLOSUM62, …)`, centre ~2× the borders)
  raises leave-one-out retrieval (balanced PR-AUC) above flat BLOSUM62 in **7/8 held-out epitopes**
  (0.564 → 0.598 @≤2; 6/8 @≤4). The weight is **end-anchored** (offset from each germline anchor,
  Beta-Binomial-smoothed) and BLOSUM severity discriminates exactly where the weight is high (the core).
  For CDR3, *where* a mismatch falls matters more than *which* residue it is; the first-order statistic
  is still the control-calibrated E-value.
- **The V gene is a strong, near-binary prior — partly recovered at near-exact germline identity.**
  Same-V neighbours share the epitope **~45–64%** of the time vs **~6–17%** cross-V (ratio up to ~7×).
  *Loose* CDR1/CDR2 similarity barely predicts it (point-biserial r ≈ 0), but cross-V co-specificity
  **rises monotonically as germline CDR1+CDR2 approach identity** — ~11% at ≥6 mismatches to ~32% at
  edit-0: a real ~3× lift, but only **~half** the same-V level (32% vs 64%), recovered *whole-loop*
  rather than via a sparse pseudosequence (per-position lift flat except CDR2 pos 5; edit-0 bin n≈25,
  noisy). So the germline loops carry part of the V prior at near-exact tolerance; the rest is
  gene-identity-specific (`bench/vregion_decompose.py`, `bench/vpseudo.py`; `vdjmatch.match.vgene`).

## Benchmark

The standing benchmark is the **2026-06-11-ZENODO VDJdb release** (mirrored on the
[`isalgo/airr_benchmark`](https://huggingface.co/datasets/isalgo/airr_benchmark) HF dataset, fetched via
`db.fetch_hf`). For the scoring studies it is **composition-controlled** (`_bench.long_list`: keep
epitopes ≥30 clonotypes, cap mega-epitopes to a random 3000, drop spectratype-anomalous spike studies),
with imbalance-robust metrics (`bench/metrics.py`: ROC-AUC + balanced PR/F1). Its highest-confidence
subset — the **shortlist** of clonotype–epitope pairs in **≥2 independent references** (`db.replicated`,
~3000 TRB + ~2000 TRA) — is the gold standard (as in `mhcmatch`), kept separate from the long-list.
Leave-one-out NN annotation (`bench/shortlist_accuracy.py`, subs 1, exact self excluded) reaches
**~44–47% top-1**, modestly above single-reference controls (the dramatic 3× gap on a sparse older
export was an artefact — the dense release makes controls findable too).

## License

GPL-3.0-or-later (it builds on `seqtree`, which is GPL-3.0-or-later).
