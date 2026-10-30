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
  <a href="https://docs.isalgo.dev/vdjmatch/"><img alt="docs" src="https://github.com/antigenomics/vdjmatch/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-Academic%20%26%20Non--Commercial-green"></a>
</p>

`vdjmatch` annotates selected T-cell receptors and repertoires against
[VDJdb](https://github.com/antigenomics/vdjdb-db), preserving query identities and
reporting peptide–MHC candidates, competing evidence and optional background calibration.
Native sequence search uses [`seqtree`](https://github.com/antigenomics/seqtree);
repertoire ingestion reuses published `vdjtools` readers.
Requires released `seqtree>=1.3.2` and `vdjtools>=5.4.0`, including allocation-aware
native thread validation.

The Python implementation follows the legacy Java/Groovy package, preserved on the
`legacy-java` branch and [Java tag `1.3.1`](https://github.com/antigenomics/vdjmatch/tree/1.3.1).
Existing Java consumers use the separate, pinned JitPack artifact
`com.github.antigenomics:vdjmatch:1.3.1`; Python releases do not replace that artifact.
See the [tutorial](https://docs.isalgo.dev/vdjmatch/tutorial.html),
[workflow guide](https://docs.isalgo.dev/vdjmatch/how-to.html) and
[score explanation](https://docs.isalgo.dev/vdjmatch/explanation.html).

## Pick an entry point

| Task | Entry point | Guide |
| --- | --- | --- |
| Annotate a repertoire or cohort | `vdjmatch annotate` | [Tutorial](https://docs.isalgo.dev/vdjmatch/tutorial.html) |
| Compare retained geometries from raw inputs | `sh examples/compare.sh NEW_DIRECTORY` | [CLI examples](https://docs.isalgo.dev/vdjmatch/examples.html) |
| Calibrate the combined score from one control | `vdjmatch annotate --control CONTROL` | [Calibration guide](https://docs.isalgo.dev/vdjmatch/control-calibration.html) |
| Check literal database presence and unusual resemblance | `vdjmatch database-evidence` | [Score explanation](https://docs.isalgo.dev/vdjmatch/explanation.html) |
| Choose a pinned or offline reference | `--vdjdb LOCAL_REFERENCE` or `vdjmatch update` | [Workflow guide](https://docs.isalgo.dev/vdjmatch/how-to.html) |

## Runnable offline examples

Browse the [examples gallery](https://docs.isalgo.dev/vdjmatch/examples.html) or
the [notebooks gallery](https://docs.isalgo.dev/vdjmatch/notebooks.html) for
rendered previews, source downloads and raw inputs.
The [versioned Hugging Face asset mirror](https://huggingface.co/datasets/isalgo/airr_benchmark/tree/314b72d2955e9578157c6767845741444d89f540/vdjmatch/releases/0.6.3)
also provides the synthetic offline examples and bundled static models with a checksum manifest.

The [examples directory](examples/) contains tiny raw synthetic inputs, a pinned
offline reference, [annotation](examples/annotate.sh) and
[comparison](examples/compare.sh) shell scripts, and output-free
[annotation](examples/annotation.ipynb) and
[comparison](examples/comparison.ipynb) notebooks. Activate your installed
environment and run from the repository root:

```sh
sh examples/annotate.sh annotation-demo
sh examples/compare.sh comparison-demo
```

Both use four native threads, reject existing output directories and fail on CLI
errors. The comparison explicitly runs exact-excluded detailed stem4 and advanced
vcdr12; notebooks analyze only fresh CLI tables. `DEMO_A`/`DEMO_B` labels are
synthetic software examples, not binding truth. See the
[example guide](https://docs.isalgo.dev/vdjmatch/examples.html).

## Scientific worked examples

Three longer notebooks connect raw biological inputs to fresh CLI tables and
publication-size figures, with editable PDF/SVG and 300-dpi PNG exports:

| Question | Notebook |
| --- | --- |
| How do common epitopes separate supplied positives and controls? | [ROC/PR, score overlap and TP/FP counts](examples/accuracy.ipynb) |
| Draft: LLW neighbourhood evidence across donors (validation deferred) | [Vaccination dataset walkthrough](examples/vaccination.ipynb) |
| How do model mass, precursor frequency and expected cells differ? | [Precursor estimates](examples/precursor_estimates.ipynb) |

Install the matching checkout from its repository root with
`python -m pip install -r examples/requirements-notebooks.txt`; keep the complete
`examples` directory together. Inputs remain in the local `airr_benchmark` and
`airr_yfv19` mirrors, configured through `AIRR_BENCHMARK` and `AIRR_YFV19`.
The notebooks preserve supplied labels and never reuse predictions as inputs.
Browse the [rendered gallery](https://docs.isalgo.dev/vdjmatch/notebooks.html)
and [measured laptop resource table](https://docs.isalgo.dev/vdjmatch/performance.html).
Four threads sped up the 2,670-query example from 13.71 to 4.69 seconds
(2.92×), with identical evidence and about 0.85 GB peak CLI process-tree RSS.

## Features

- Latest/pinned GitHub releases and offline legacy, rich or AIRR reference inputs.
- Source-aware junction ingestion, stable row identities and explicit ingestion counts.
- Native batched searches by locus, with retained reference observations.
- Targeted stem4 kernel annotation; frozen alignment scores, PSSM-density, fixed-scope unit-edit distance
  and neighbor counts, and other retained scorers remain comparison options.
- All peptide–MHC candidates, NED-v1 ranking and competing-label support.
- Optional finite-background calibration with a stated counting unit and search predicate.
- Paired evidence from linked reference complexes plus unlinked same-target chain marginals.
- CLI hit/candidate/call/summary/ingestion tables and a reproducibility manifest.
- Separate precursor-frequency APIs and their established worked examples.

### Targeted annotation

```console
$ vdjmatch annotate sample.tsv --threads 4 -o annotation
$ vdjmatch annotate airr_directory/ --threads 4 -o cohort
$ vdjmatch annotate --metadata samples.tsv --vdjdb reference.zip --threads 4 -o cohort
```

The reference is fetched automatically; `--vdjdb` selects a pinned offline table
or archive. A metadata TSV has `path` and optional `sample_id` columns; relative
paths resolve beside the metadata file. Directory discovery accepts TSV, CSV,
text, gzip and Parquet tables. Sample names must be unique.

AIRR input uses anchor-inclusive `junction_aa` and V calls. Version 0.6.0 defaults
to `--reference-population curated-benchmark`, matching the manuscript primary
reference view. Use `--reference-population all-valid` to reproduce the accepted
all-valid baseline. Curation excludes KLGGALQAK and records sourced to PMID
40498839; natural observations from other sources remain eligible. Exact hits
remain included unless `--exclude-exact` is supplied. The default scores all
peptide–MHC targets with the stem4 BLOSUM62 kernel ladder. `--locus paired` also
retains independent TRA/TRB evidence from unlinked references. Scores are relative
evidence; whole-selection null calibration is required for false-call recommendations.
The original PSSM-density scorer remains available through `historical-density`;
`trained-score` uses the bundled frozen alignment coefficients and hit informativeness.
Advanced alternatives remain available through `historical-density` and `neighbours`.

| Query format | Required sequence/gene fields | Identity/linkage |
| --- | --- | --- |
| Single AIRR | `junction_aa`, `v_call`, `locus` (or locus inferred from V) | Optional `sequence_id`, `query_id`, `clonotype_id`, `clone_id` |
| Simple repertoire | `junction_aa`, `v`, `locus`; optional `j` | Optional `clonotype_id`, `clone_id` or `cloneid` |
| Linked AIRR | One TRA/TRB row per linked receptor, `junction_aa`, `v_call`, `locus` | Shared `cell_id`, `pair_id`, `clone_id` or `cloneid`; cell/pair linkage takes priority |
| Wide paired repertoire | `junction_aa_alpha`, `junction_aa_beta`, `v_alpha`, `v_beta`; optional `j_alpha`, `j_beta` | Optional `clone_id` or `cloneid`; otherwise row identity links the pair |

`--locus auto` infers single or paired populations; `--locus TRA`, `TRB` or
`paired` selects explicitly. Missing chains and invalid rows retain dispositions;
multiple chains in one clone/locus require resolution before annotation.
`clonotype_id` and `clone_id` are retained as separate identities. Bare AIRR
`cdr3_aa` excludes conserved anchors and cannot replace `junction_aa`.
For amino-acid clonotype inputs, collapse by the supplied V, junction and J
calls, summing counts before annotation. Distinct J calls remain distinct
clonotypes even when their search geometry is shared.

The default writes `.queries.tsv`, `.selection.tsv` and `.manifest.json`, retaining
sample and source identities, competing scores, ties and no-support rows.
`--details` additionally writes candidates and per-target evidence; its memory
grows with queries × targets. `--epitope PEPTIDE` restricts the requested targets.
`--exclude-exact` is for honest benchmarks; ordinary annotation retains exact hits.
V-loop lookup accepts model-backed dual-use gene spellings such as
`TRAV14_DV4*02` and `TRAV14/DV4*02` as the same allele, while retaining
the submitted call in outputs. Unknown alleles are not silently renamed.

Advanced scorer and geometry recipes remain in the [CLI reference](https://docs.isalgo.dev/vdjmatch/cli.html).
The manuscript's curated primary comparison has higher equal-task mean average
precision than each of seven available rivals on that rival's fixed matched task
set. This is a matched-mean result, not a claim of winning every individual task.
The accepted all-valid baseline and its curves and intervals remain separate;
changing the reference population does not transfer calibration between views.
Calibration is specific to the complete score and selection procedure.
`--threads N` budgets native compute threads; searches batch queries by population
and do not create a worker process per receptor or input file.
Mixed single/paired cohorts share marginal chain searches. V-loop tables are
prepared once per invocation, including the selected matrix and V regions.
The library also enforces this ceiling for the Polars table pool. If Polars was
already imported with a larger pool, start a fresh process with
`POLARS_MAX_THREADS=N` set before importing Polars; the API fails explicitly
rather than exceeding the requested budget. Input metadata whose names conflict
with generated output fields is retained under `input_` prefixes.

```python
from vdjmatch import annotate_repertoires

result = annotate_repertoires("airr_directory", reference="reference.zip",
                             output_prefix="cohort", threads=4)
queries, selections = result["queries"], result["selection"]
```

### Complete-score calibration from one control

```sh
vdjmatch annotate sample.tsv --vdjdb reference.zip --control control.tsv --threads 4 -o calibrated
```

With the default `--background vdjdb-other`, supplying one raw control selects the
opt-in combined score and computes both conditional randomization tails on demand.
It writes `.queries.tsv`, `.database.tsv`, `.epitopes.tsv` and `.calibration.tsv`;
no fitted calibration file or JSON input is needed. `--mock-replicates` defaults
to 200 and `--seed` to 0. Reference selectors (`--epitope`, MHC,
`--evidence-type`, `--min-score`) apply before calibration.

Benchmark `--exclude-exact` filters exact hits only; reference models, counts,
priors and graph weights remain fixed across queries.

Your own reference can be a `junction_aa,locus,epitope` table, with `pair_id`
linking long paired rows, or a wide `junction_aa_alpha/beta,v_alpha/beta,j_alpha/beta,epitope`
table. Supplied pairing and metadata are retained; calibration uses that selected table.

The control requires full junctions and unambiguous V/J calls. Productive flags
and valid amino-acid junctions define eligibility; excluded counts are reported.
Rearrangement-event weights, rather than read abundance, define sampling.
For paired queries, the same control table supplies TRA and TRB, drawn independently.
Missing partners retain their raw rows and use eligible single-chain evidence.

Database tails include label and paired-axis selection; epitope tails test a
fixed full peptide–MHC group. `p_both = max(p_reference, p_query)` requires evidence
in both directions. Null match rates are separate `k/B` estimates; none of these
quantities is a binding probability. This candidate recipe derives the original
component prior from the matching chain control (TRB for paired queries), which
differs from historical benchmark model assets. It remains opt-in pending full
nonregression checks. The no-control stem4 default and older real/generative
control routes remain available. See the [calibration guide](https://docs.isalgo.dev/vdjmatch/control-calibration.html).

The original `Annotator` and `annotate` APIs remain available. Agent usage is
documented in the repository [vdjmatch skill](skills/vdjmatch/SKILL.md).
The synthetic library/CLI batch resource gate runs with test/bench extras:
`python -m pytest tests/test_annotation_resources.py -q`.

Offline call-policy experiments remain available through `reference-calls`,
`calibrate-calls` and `apply-call-calibration`; see the
[CLI reference](https://docs.isalgo.dev/vdjmatch/cli.html). Their pinned operating
points belong to the exact scorer, reference, null population and paired selection
event used to fit them. They are not inputs to the one-control workflow above.

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
  some epitopes) and vanish from the total without an error. `isalgo` selects current gene-level fits; `learned` and `arda` alias those same
  weights. Mouse TRA/TRB/IGH/IGK/IGL are available; mouse TRG/TRD fail explicitly.
  Historical learned/ARDA results require their original pinned weights. The `n_zero_pgen` output column reports the
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
calls = ann.annotate(queries, threads=4)  # every input row remains
candidates = ann.candidates(queries, threads=4)  # all peptide–MHC alternatives
calibrated = ann.annotate(queries, calibrate=True, species="human", threads=4)
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
vdjmatch match --input-format airr --threads 4 -o results/run sample.tsv
vdjmatch match --vdjdb reference.zip --no-evalue -o results/offline sample.tsv
vdjmatch match --vdjdb reference.zip --paired --link cell_id -o results/cells cells.tsv
vdjmatch match --vdjdb reference.zip --search-mode ball -o results/balls sample.tsv
```

CLI matching requests calibration by default; `--no-evalue` selects uncalibrated
ranking. The command emits hits, candidates, row-preserving calls, descriptive
summaries and ingestion diagnostics as TSV, plus a JSON manifest.

Use `--input-format custom --sequence-convention junction` for an ambiguous custom
table. `--scope s,i,d,t` sets substitution/insertion/deletion/total budgets. Hard
`--match-v` or `--match-j` requests require `--no-evalue` with sequence-only backgrounds.
Ordinary CLI search defaults to allocated CPUs (`--threads 0`); use `--threads 1` for serial execution. See `vdjmatch match -h` and the
[CLI reference](https://docs.isalgo.dev/vdjmatch/cli.html).

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

ISALGO LAB Academic and Non-Commercial Research License v1.0. See [`LICENSE`](LICENSE).

### Repertoire-level enrichment (M3)

`vdjmatch repertoire-enrichment samples.tsv --null reference --vdjdb reference.zip --control control.tsv.gz --threads 4 --output-prefix enrichment`
holds each repertoire fixed and compares its unique receptor support against fresh
Pgen-matched mock epitope corpora. `--mock-policy pgen-vj` preserves V/J usage;
`pgen-different-v` changes V and junction identities. Paired inputs retain orphan
chains and count a receptor once per target. `--significance nominal|bh|both`
selects explicit inference options; no binary significance call is the default.
All no-hit targets remain in the outputs. See the
[CLI contract](https://docs.isalgo.dev/vdjmatch/cli.html) for matching, Pgen,
replicate resolution, test families and input schemas. The historical
`--null repertoire` mode remains separate.

An independent counts-only comparison filters both the sample and each pMHC
reference with the existing TCRnet neighbourhood test, then runs flexible
VDJmatch matching:

```sh
vdjmatch repertoire-tcrnet samples.tsv --vdjdb reference.zip \
  --control control.tsv.gz --tcrnet-fdr 0.1 --scope 3,1,1,3 \
  --threads 4 --output-prefix filtered
```

This mode exports filter evidence, matching hits and raw-denominator counts;
it does not produce repertoire P-values or assume random backgrounds filter to
an empty set. Paired inputs retain orphan chains and count each receptor once.


For inputs too large to retain as one sample, `vdjmatch.io.iter_rearrangement`
yields normalized batches with file-wide positional IDs. The underlying
`vdjtools.io.iter_table_batches` supports projected raw CSV/TSV/gzip/Parquet
batches. Prefer the simpler one-sample-at-a-time path when profiling shows it
fits; within-sample batching does not by itself bound reference/control memory.
