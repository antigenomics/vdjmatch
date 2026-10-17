---
name: vdjmatch
description: Annotate single or paired TCR repertoires against VDJdb with the vdjmatch library or CLI, including AIRR directories and sample-path metadata.
---

# VDJmatch repertoire annotation

`reference-structure --control` exports length KL(reference||control) and
equal-mixture Jensen–Shannon divergence in bits per junction/V/J key in
`.spectratype.tsv`. Unsupported control lengths give infinite KL; JS remains
bounded by one bit. `.degrees.tsv` includes self-inclusive degree, and
`.graph-acquisition.tsv` interpolates expected edges and connected sampled
vertices of the fixed Hamming<=1 graph without replacement. These are
descriptive effects, not class support tests or unseen-family coverage.
Add `--publicity` to search all reference junctions against unique control
junctions in one native Hamming<=1 batch. Summary fields retain valid source
record counts, matched records and record/control pair denominators.
Structure summaries include uniform-vertex binary-similarity Hill0/1/2,
Shannon-equivalent entropy and N/Hill enrichment from existing graph degrees.
Degree+1 Chao1/coverage is explicitly a synthetic degree-mass proxy; actual
occurrence-based coverage remains in `reference-coverage`.
`reference-structure --motif-mode full20 --epitope-idf` emits full20 exact3
and single-X wild4 views, with distinct V/length-conditional epitope-document
IDF and control rarity. Index profiles preserve core3 support and add the
separate full20 union's selection-aware class and marginal any-chain support.
`score-confusion candidates.tsv --truth truth.tsv --profiles profiles.tsv`
exports source-eligible per-epitope TP/FN/FP/TN, MCC and balanced accuracy from
complete explicit scorer dispositions. Ties, multiple truth labels and no-hits
are retained. Exact-hit filtering belongs to the producing scorer.
Its `.class-groups.tsv` uses equal-epitope means/medians/ranges and explicit
defined denominators for alpha/TRA, beta/TRB and joint/paired panels. Paired
eligibility requires one chain>=30 source records, never a sum of small chains.
New degree/Hill/Chao/acquisition indices use allele-pooled junction/V-gene/J-gene
clonotypes and keep unresolved categories; `.clonotype-degrees.tsv` reports
their self-inclusive degree histogram. Legacy `hamming_*`/`.degrees.tsv`
statistics remain junction-only.

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
  see the CLI documentation for recorded geometry recipes. Calibration does not
  transfer between scorers, geometries or selection procedures. Preserve the
  original `historical-density` route and its explicitly supplied model assets.
  PSSM-density and naive fixed-scope scores are current alternatives, not obsolete
  methods. `fixed-scope` reports unit-edit distances and neighbor counts separately;
  use a query table with declared full pMHC targets. State the exact operation
  budgets and counting unit, and preserve every no-hit and supplied negative.
  `trained-score` uses the bundled frozen alignment model and hit informativeness;
  its candidate scope and aggregate indel cap precede NW alignment. Preserve its
  trained marginal outputs without inventing a paired probability or transferring
  another scorer's calibration.
- Supply an explicit `--threads` budget. Keep queries in population batches; do not
  wrap annotation/search in a per-receptor or per-file multiprocessing pool.
  For Python, set `POLARS_MAX_THREADS` before importing Polars if its existing pool
  would exceed the requested budget; the API rejects that mismatch before loading.
- Use matrix-relative distance settings for comparisons. The shortlist API's
  omitted radius and stem4 reference-call geometry derive their scales from the
  selected matrix. Historical absolute overrides remain explicit reproduction
  options; a censored sparse distance is not an uncensored nearest distance.

## One raw control, complete-score calibration

```sh
vdjmatch annotate sample.tsv --vdjdb reference.zip --control control.tsv --threads 4 -o calibrated
```

Pass the reporting threshold with `--fpr` (default `0.01`); the calibration
library accepts `report_p=0.01`. Compact reports retain all supported labels
passing the inclusive any-call threshold, plus declared targets regardless of
support. Calibration metadata records the threshold as `report_p`.
Query-null storage keeps finite supported scores and per-draw all-label maxima;
no-hits remain in every denominator. The 2 GiB auxiliary guard includes finite
record growth and sorting scratch; full raw-component limits remain separate.

With default `--background vdjdb-other`, this explicitly opts into the combined
score. No fitted JSON input/output is required. Outputs are `.queries.tsv`
(every raw chain row with calibration identity/status), `.database.tsv`
(selected-label/axis tails), `.epitopes.tsv` (fixed full-pMHC tails and score
components; continuous `score` versus support-gated `call_score`), and `.calibration.tsv` (population, seed/model/native provenance).
Use `--mock-replicates` (default199) and `--seed` (default0) explicitly when pinning
reproduction. Reference species/locus/evidence/score/pMHC filtering precedes mocks.
Custom labelled references use long `junction_aa,locus,epitope` rows with `pair_id`,
or wide `junction_aa_alpha/beta,v_alpha/beta,j_alpha/beta,epitope` rows; retain their
links and metadata and rebuild calibration against the selected table.

Controls need valid anchor-inclusive junctions and unambiguous V/J calls.
Productive flags and invalid-AA exclusions are counted. Sampling weights describe
rearrangement events, not reads. Paired controls require both loci with independent
randomization; an incomplete genuine pair keeps its missing row and uses eligible
marginal evidence. Retain generated calibration IDs when joining outputs across
populations; original query IDs need not be globally unique.

`p_both=max(p_reference,p_query)` tests both conditional directions; match rates
are separate `k/B` estimates, never binding probabilities. The candidate prior
comes from the matching-chain control (TRB for paired scoring), distinct from
historical benchmark assets. Keep this route opt-in until nonregression gates
pass; ordinary no-control annotation keeps the accepted stem4 scorer.
Nondefault geometries/target-selection options fail before computation. Legacy
`--background real/generative --control` routes retain their existing meaning.
See [the calibration guide](../../docs/control-calibration.rst).

For library/CLI engineering, the synthetic serial/parallel and process-tree memory
gate is `python -m pytest tests/test_annotation_resources.py -q` with test/bench
extras installed. Scientific comparisons require matched populations and exact
exclusion policies; preserve established results and report every task delta.

For separate offline call-policy experiments, use the documented
`reference-calls` / `calibrate-calls` / `apply-call-calibration` commands. Preserve
producer manifests, null identity and exact-hit policy. A single-chain or
linked-only operating point cannot calibrate the paired any-axis event. These
advanced recipes do not supply files to the one-control workflow above.

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
TP/FP trajectories, complete vaccination-repertoire evidence, and precursor
frequency. Keep the complete examples bundle, including `profile_cli.py`.
Configure local raw mirrors with `AIRR_BENCHMARK` and `AIRR_YFV19`.
The vaccination driver retains every locally available repertoire named in
`metadata.txt`, with donor/day/technical-replicate identities and complete counts. Do not publish original sample
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

## Hard rule: benchmark exact exclusion

Exclude exact junctions only from search hits. Never remove a query from model
populations, V/J/length priors, motifs, graph structure, weights or fitted parameters.
Build each real/mock reference model once. Historical punctured diagnostic exports
are not benchmark scoring inputs and must not define an “honest” scorer pool.

## Reference morphology for structure-guided score investigation

`reference-structure --vdjdb REFERENCE --species human --locus TRA|TRB --threads 4
--output-prefix OUT` exports full-pMHC Hamming-one graph, junction-length and V-usage
features. `--targets-from-sample RAW_AIRR` uses declared compatible pMHC unions.
Each vertex is a unique full junction. The command exports the complete fixed
reference graph, degree histogram and `.lengths.tsv` exact-length coverage table;
the latter separates valid source records, modeled records and unique graph
junctions, retaining short junctions and unresolved V genes. No query-specific graph is rebuilt or
exported. Keep graph-node, clonotype-key and source-record units separate;
entropies are in nats and lengths include anchors. Empty reference structure is
unavailable. One native population batch feeds all target graphs, without a
per-query/target native pool or persisted feature cache.
The `.genes.tsv` export gives V/J gene usage overall and at each exact junction
length, with explicit record and junction/gene-key denominators. Alleles are
pooled; unresolved genes remain a counted null category. Usage is composition,
not background enrichment or proof of germline contacts.
These features do not fit fusion weights or calibrate recognition. Any harmonized
score sum still needs matched component identities and final-selection calibration.

With `--control`, choose `--motif-mode exact3` for full-AA 3-mers or `wild4`
for full-AA 4-mers with all four single-wildcard masks. `imgt3` and the default
`imgt_wild4` use eleven IMGT classes. `.motifs.tsv` labels the mode/alphabet;
counts retain once-per-junction/V-gene presence and fixed V/length denominators.
Keep these binary selection bits separate from TF-IDF, bootstrap uncertainty
and final specificity probabilities. Each invocation uses the same one native
graph batch; the existing 256 MiB count-buffer guard applies to every mode.

Motif rows now include global-control document frequency, `idf` and `tfidf`:
documents are distinct junction/V-gene keys, IDF is `1+ln((M+1)/(df+1))`, and TF
is row-stratum motif presence divided by all target keys (no vector norm).
`--mock-replicates 50 --seed 0` adds matched-control TF-IDF mean/95% envelopes.
Fresh vdjtools Pgen feeds the existing joint-event sampler; record counts are
matched before key deduplication, with no V/J/length matching. Keep the fixed
IDF, absent mock features as zero, and these descriptive envelopes distinct
from uncertainty on selection bits or final specificity probabilities.

`reference-coverage --vdjdb REFERENCE --locus TRB --geometry levenshtein
--scope 1,1,1,1 --nboot 50 --threads 4 --output-prefix OUT` uses optional
`vdjmatch[diagnostics]` (seeded Leiden). `--geometry sw --sw-radius R` requests
full SW/BLOSUM62 self-score distances over the complete junction and CDR1/2/2.5;
`--cdr-loops` explicitly chooses a smaller loop set. No trimmed/stem surrogate.
One population search batch supplies every real/mock graph; pair/membership/
expansion/bootstrap guards fail explicitly. Supply `--control RAW_CONTROL
--mock-replicates 50` for matched control populations with their own partitions.
VDJdb occurrence counts are retained at exact junction/V/J-call and Leiden
community resolution. `.coverage.tsv` separates Chao1 richness, coverage/unseen
mass and unseen-type fraction; `.rarefaction.tsv` has effort/coverage curves
and uncertainty, reusing vdjtools statistics/augmented bootstrap. Partitions
are fixed within each population's rarefaction; no claim of unseen-family
exhaustion, donor/study incidence inference or automatic adequacy cutoff.

`combine-evidence EVIDENCE --null-evidence NULLS --three-channel --output SCORES`
puts `original_score`, `gapblock_score`, `stem4_score` on conservative matched
null-tail scales and sums their explicit `original_weight`, `gapblock_weight`,
`stem4_weight` contributions. Weights are nonnegative and sum to one; retain
them in output and derive them from the declared reference morphology, never
assay labels or epitope lookup. Missing components contribute zero without weight
redistribution. `combined_score` retains continuous ranking evidence;
`combined_call_score` and `combined_status` separately describe weighted hit
support. Finite prior evidence does not by itself authorize a call. The
low-level command consumes explicit weights; `annotate --control` derives them
from the reference and recomputes them for mock references. Input
scorer/reference/exclusion/null identities must match. Component tails are not
the final fused P-value.

## Database-level diagnostic

For a separate scoped nearest/count diagnostic, use
`database-evidence SAMPLE --vdjdb REFERENCE --control RAW_CONTROL --locus TRB
--species human --scope 1,0,0,1 --statistic nearest --threads 4 --output-prefix OUT`.
It computes fresh native marginal Pgen, samples reference-size/per-pMHC Pgen
strata on demand and tests both directions. `p_reference` holds query fixed and
randomizes reference; `p_query` holds reference fixed and randomizes query within
its one-bit Pgen stratum. Ties count; all draws/no-hits remain denominators.
`p_both=max(p_reference,p_query)` tests evidence in both directions without
independence, conditional on each null being valid. No calibration files or JSON
are needed; output is three flat `.database/.epitopes/.calibration.tsv` files.
Scope and reference evidence filters are applied before recalibration. `count`
uses distinct matching junctions. Null hit rates `k/B` are estimates of any match
within scope, distinct from plus-one tails at the observed score. Event weights
ignore read abundance and use occurrences/NT events/pooled NT incidence or the
available convergence unit. Incidental exact control/reference overlaps remain.
The offline `examples/database.sh` uses an explicitly synthetic single control
and native Pgen; it is not a validated biological null. This adapter is junction
scope/count only: never transfer these P-values to V loops, combined scores,
selected labels or paired receptors. Reference-derived score weights must be
recomputed on mock references when implementing those adapters.

The following older route remains available as a legacy diagnostic:

`database-evidence SAMPLE --vdjdb REFERENCE --locus TRA|TRB --species human
--null-sample RAW_NULL_DRAWS --threads 4 --output-prefix OUT` separates a literal
`database_exact_junction` lookup from `database_p`: the matching-null probability
of being at least as close to any database junction. It ignores antigen labels,
V/J and source-study support cutoffs, retains ordinary exact records, and uses
one native population batch for query/null nearest full-junction edit distances.
Null repetitions and no-hits remain in the denominator; invalid/mismatched null
rows are rejected. For multiple supplied sources, retain the largest tail.
An invalid query has unavailable diagnostics; a valid no-hit has P=1. State the
fixed reference, search limits and null law/exchangeability assumption. This is
a junction-level diagnostic, not epitope confidence, V/J-clonotype membership or
a paired P-value. Its artificial null draws demonstrate software behavior, not
a biological null law. Do not substitute the legacy `p_any` or a component tail.

`database_exact_gene_clonotype` separately checks a complete junction/V-gene/
J-gene key with alleles pooled; missing/ambiguous query calls remain unavailable.
The P-value still ignores V/J. `--max-edits` bounds full-junction Levenshtein
distance (default5); native workers follow `--threads`, Polars/BLAS helpers use1.

- `reference-structure --control` exports gene-consensus CDR1/CDR2/CDR2.5 residue factors and KL plus `residue_bits_lower`, allowing all compatible background alleles. Keep full-junction/V-key mass after projection, allele ambiguity and one-based model alignment columns explicit; do not sum correlated residue views or interpret the bound as a sampling interval/P-value. See `docs/selection-information.rst`.

### Experimental frozen reference policy

`vdjmatch score-policy components.tsv --policy reference-policy.tsv --output-prefix selected`
selects existing continuous V/conditional-rank scores using a frozen reference-only
model. Input columns: species/locus/epitope, original_score (or --original-column),
rule_v_score/rule_rank_score. Requires the policy companion manifest and table hash.
It preserves cases, falls back explicitly, does not fit or change defaults, and
provides no transferred component P/E/posterior. See docs/scoring.rst.
