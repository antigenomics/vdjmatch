---
name: vdjmatch
description: Annotate single or paired TCR repertoires against VDJdb with the vdjmatch library or CLI, including AIRR directories and sample-path metadata.
---

# VDJmatch repertoire annotation

Use the installed `vdjmatch annotate` command or Python `annotate_repertoires`.
Read [README.md](../../README.md) for input schemas and [CLI documentation](../../docs/cli.rst)
for advanced scorers. Run `vdjmatch annotate --help` to verify installed options.
Use `vdjtools>=5.0.1` for ingestion and allocation-aware native thread validation.

```sh
vdjmatch annotate airr_directory/ --threads 4 -o annotation
vdjmatch annotate --metadata samples.tsv --vdjdb reference.zip --threads 4 -o annotation
```

Metadata has `path` and optional `sample_id`; relative paths resolve beside the
metadata TSV. Directory entries must be repertoire tables and have unique names.
Omit `--vdjdb` to fetch the current reference; use an explicit offline reference
for reproducible comparisons.

```python
from vdjmatch import annotate_repertoires
result = annotate_repertoires("airr_directory", reference="reference.zip",
                             threads=4, output_prefix="annotation")
queries, selections = result["queries"], result["selection"]
```

- AIRR sequences are anchor-inclusive `junction_aa`, with `v_call` and `locus`.
  Never fabricate anchors from bare AIRR `cdr3_aa`.
- Simple tables use `locus`, `v`, optional `j`, and `junction_aa`.
  Unique amino-acid clonotype panels collapse the supplied V/junction/J key
  and sum counts; retain different J calls and do not deduplicate by nucleotides.
- Long pairs link TRA/TRB rows by the supplied cell/pair/clone identity. Wide pairs
  use `junction_aa_alpha/beta`, `v_alpha/beta` and optional `j_alpha/beta`.
  Keep clone and clonotype IDs separate; never select one ambiguous chain silently.
- Auto locus routing is the default. Use `--locus TRA`, `TRB` or `paired` explicitly
  when the intended population cannot be inferred, including incomplete pair cohorts.
- Keep `sample_id` and source identities when joining outputs. Retain no-support,
  invalid and missing-chain dispositions; these are not experimentally negative labels.
- Default stem4 kernel scores and support shares are relative evidence, not binding
  probabilities, E-values or calibrated false-call rates. Linked and unlinked marginal
  paired evidence remain separate; absence of linked references cannot veto marginals.
- Ordinary annotation retains exact hits and conflicting-label ambiguity.
  `--exclude-exact` is the explicit benchmark option.
- Compact output is the default. `--details` exports continuous evidence per target
  and candidates; its memory grows with queries × targets. `--epitope` restricts
  requested targets. `--kernel-scorer` selects retained advanced kernel variants;
  `--v-regions cdr12 --stem-depth 0 --distance-matrix clipped
  --junction-ends trim3` selects the retained `vcdr12` geometry.
  Matrix extras use `--distance-matrix pam250|mj|vdjam` with stem depth 0 and trim3.
  `sig` uses compact `--position-weighting significance-trimmed
  --distance-matrix clipped --stem-depth 0 --junction-ends trim3`; `sigb62` uses
  blosum62. Preserve their pre-V junction rounding; calibration is geometry specific.
  `--geometry-extra flank|framev|framemid|framej|freegap|jw|stem4j|jlen|vjlen|norm|jnorm|sw|nw`
  selects a fixed retained recipe and sets its matrix, trim3 ends and stem term.
  Require compact uniform full-alphabet geometry and all three V regions.
  Frame recipes preserve the recorded callback ordering. `norm`/`jnorm` retain
  pre-exclusion same-V target counts in their denominator, so exact-excluded
  runs are historical reproduction diagnostics, not fully punctured models.
  Do not transfer targeted stem4 calibration to any alternative.
  `historical-density` retains the original scorer and its advanced options.
  For original held-out components, use `historical-density --legacy-components
  --reference-population all-valid --germline-background RAW --control RAW`;
  retain the owning full reference union and supplied prior/control populations.
- Supply an explicit `--threads` budget. Keep queries in population batches; do not
  wrap annotation/search in a per-receptor or per-file multiprocessing pool.
  For Python, set `POLARS_MAX_THREADS` before importing Polars if its existing pool
  would exceed the requested budget; the API rejects that mismatch before loading.

For library/CLI engineering, the synthetic serial/parallel and process-tree memory
gate is `python -m pytest tests/test_annotation_resources.py -q` with test/bench
extras installed. Scientific comparisons require matched populations and exact
exclusion policies; preserve established results and report every task delta.

Single-chain complete-score calibration uses `reference-calls --scorer stem4
--reference-population curated-benchmark --leave-epitope-out`, followed by
`calibrate-calls` on its `.calibration-cases.tsv`. Pin the reference, model,
scorer, exclusion policy and generated null source before recommending a cutoff.
Use `--input-manifest` for fitting and application to enforce the producer binding.
For genuine paired withheld nulls, add `--locus paired`; fit with
`--score-column kernel_any_axis_score`. This is the maximum of final selected
linked, alpha and beta scores; labels and ties stay separate on each axis.
Include matched generated paired nulls before recommending repertoire calls.
Do not transfer a single-chain or linked-only cutoff to this any-axis event,
or a withheld benchmark recipe to an ordinary all-valid exact-retained model.
The original `reference-calls` mode remains available and is its default.

## CLI-first examples

Use [examples](../../docs/examples.rst) for raw synthetic offline inputs,
annotation/comparison shell commands and fresh-kernel notebooks. Run
`sh examples/annotate.sh NEW_DIRECTORY` or
`sh examples/compare.sh NEW_DIRECTORY` with the installed CLI on PATH.
All computed results come from CLI tables; preserve source identity and full
pMHC keys in downstream joins. Scripts default to threads4, reject an existing
output directory and propagate errors. DEMO labels are not experimental truth.
Do not add computed caches or use old notebook outputs as prediction inputs.

## CLI-first scientific examples and measured resources

Use [the notebook gallery](../../examples/README.md) for labelled TRB ROC/PR,
TP/FP trajectories, baseline vaccination repertoire abundance, and precursor
frequency. Keep the complete examples bundle, including `profile_cli.py`.
Configure local raw mirrors with `AIRR_BENCHMARK` and `AIRR_YFV19`;
`YFV_INPUT` chooses a specific baseline table. Do not publish original sample
names, donor identifiers or source-identifiable raw/execution tables.

Run core calculations through the installed CLI; plot only its final tables.
Assayed-target ROC/AP uses full pMHC joins, supplied labels, exact exclusion,
and all source observations including unsupported scores. It is distinct from
all-label selected-call calibration. Never infer experimental negatives or
vaccine expansion from unlabelled baseline annotation. Keep read-count abundance,
model precursor frequency and expected cells as separate quantities.

Four native threads are the worked-example setting; the runtime default remains
one. Limit numerical helper pools before imports, run one native batch per
population and no outer search pool. Profile whole-invocation wall time and
summed concurrent CLI process-tree RSS, with CPU seconds and live-thread counts;
live idle threads alone do not demonstrate oversubscription. Compare serial and
parallel final evidence for equality. Report host, versions, raw/reference hashes,
query/target counts and timing scope; do not extrapolate laptop measurements to
other hardware or a full manuscript gate. See [performance](../../docs/performance.rst).

For notebook figures, state the estimand, class balance/counts and caption before
polishing. Keep colour/line meanings consistent, use log axes for broad positive
score ranges, export editable PDF/SVG and 300-dpi PNG, and inspect the final-size
render. Execute from a clean kernel, then commit only reviewed aggregate figures and counts;
complete logs and identifier mappings stay private. Documentation builds never
execute the notebooks.

## Reference morphology for structure-guided score investigation

`reference-structure --vdjdb REFERENCE --species human --locus TRA|TRB --threads 4
--output-prefix OUT` exports full-pMHC Hamming-one graph, junction-length and V-usage
features. `--targets-from-sample RAW_AIRR` uses declared compatible pMHC unions.
Each vertex is a unique full junction. The puncture table removes that complete
junction, including all associated V keys/records, and computes articulation
fragments exactly. Keep the graph-node, clonotype-key and source-record units
separate; entropies are in nats and lengths include anchors. Empty structure is
unavailable, not negative recognition evidence. One native population batch feeds
all target graphs; no per-query/target native pool or persisted feature cache.
These features do not fit fusion weights or calibrate recognition. Any harmonized
score sum still needs matched component identities and final-selection calibration.
