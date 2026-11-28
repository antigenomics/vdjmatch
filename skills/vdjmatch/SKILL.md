---
name: vdjmatch
description: Annotate single or paired TCR repertoires against VDJdb with the vdjmatch library or CLI, including AIRR directories and sample-path metadata.
---

# VDJmatch repertoire annotation

## Choose the reference

Use `Annotator.latest()` for official vdjdb-db, `from_path(path)` for custom
VDJdb-format data, and `from_labels(rows_or_table, locus="TRA"/"TRB")` for simple
`junction_aa,label` references. Three-value rows or wide alpha/beta junction
columns define real pairs. `match --reference PATH --reference-format labels`
uses the same direct index; add `--locus TRA/TRB` or `--paired` as appropriate.
Labels need no V/J, MHC or species boilerplate and skip reference diagnostics.
Exact hits are retained in ordinary queries. Exclusion is explicit classifier
assessment/record reliability, never an ordinary-query default. See
[reference interfaces](https://docs.isalgo.dev/vdjmatch/api.html#user-owned-labeled-references).

## Read before selecting a scorer

Read [the scoring evidence and definitions](https://docs.isalgo.dev/vdjmatch/scoring.html) before
choosing a fuzzy mode.
Preserve one-mismatch matching; prioritize BLOSUM62 and the established
end-anchored positional-BLOSUM62 recipe. Data-derived VDJAM underperformed
BLOSUM62. VDJAMr means genetic-code substitution accessibility, a different
comparator; neither is a reason to introduce a chemistry-driven replacement.
Reuse `match.regions.significance_weights/significance_pssm` and the existing
historical matcher. Mode availability is not evidence of improved accuracy.

For the frozen **vX research score**, see the
[score and decision equations](https://docs.isalgo.dev/vdjmatch/vx-formalism.html).
Keep TCRvdb, measured NLV, and ELA–LLW–LLL as separate benchmark cohorts.
TCRvdb was compiled without those two experimental series. Their records
may be present in VDJdb, so verify reference and model exposure separately.
Order target hits by the weighted vX distance, with the fixed TRB
reference-ID discount applied before ordering. With first and second hits
`d1,d2`, native competitor distance `c`, `L` observed chains and
maximum weighted junction edit cost `C=12`, use
`X=max(d1+L*m-c,0)`, `G=min(2*L*C,max((d1+d2)/2-d1,0))`, and
`S=-d1-gamma*X+b*G`. Evaluate the midpoint before subtraction. The fixed cap bounds the second-hit effect. For paired-chain receptors, sum the two chain distances within each reference pair before ranking hits. Repeat
a singleton hit, use zero penalty for an infinite background, and give
unsupported targets negative infinity. Apply benchmark exact exclusion to
hit lists with the model held fixed.

For vX classification, require a unique finite winner. The score-only assay
diagnostic adds `log1p(N_e)/400`, where `N_e` counts fixed target keys.
The coefficient lies inside the verified three-class assignment plateau.
The separate within-assay diagnostic uses
`H=-log(pK)+(S-median_assay(S))/40+0.95*log(assay_class_fraction)`.
The fixed size offset cancels when centering by target. Here `pK` is the
inclusive add-one kernel tail. The generated
null evaluates raw `S` after selection. Count all null queries in the FPR
denominator. At a common cutoff, `Q*FPR` is the empirical expected number
of false calls for `Q` query opportunities under that null. Sample
enrichment uses the same endpoint in observed and Pgen-matched mock panels.
The accepted yellow-fever curve-mean decision uses its positive-hit mock
tail and donor/day correction.

For epitope structure, define the vertices, radius and gene gates before
computing a degree histogram. Use the existing graph acquisition and
vdjtools subset-rarefaction APIs. The vX clustering benchmark uses reciprocal
k-NN on weighted pairwise distances with all boundary ties retained. The vX
target score instead uses the fixed target and competitor panels. A common
vX cutoff reproduces exact or Hamming calls only when the score ordering
strictly separates the two membership sets. See the
[structure and rarefaction equations](https://docs.isalgo.dev/vdjmatch/vx-structure.html).

For contact-derived geometry research, align native structural positions in
IMGT order and smooth frequencies using resolved-stratum counts. Use the
valleys to define windows and retain measured positive sites. For weights
`f` and clipped site costs `delta`, compare the fixed-denominator mean
`sum(f*delta)/sum(f)` with the effective-width sum using
`sum(f)^2/sum(f^2)` as the width. Check native wins and every frozen endpoint
on the same cases before carrying a geometry through the remaining benchmarks.
When testing a coupled TRA–TRB term, estimate edge weights on structures
that resolve both sites. Normalize the weighted product of the two clipped
site costs by the edge-weight sum and the maximum single-site cost.
For the benchmarked human research extension, add one thirty-second of
that result to each paired-chain target distance before ordering hits. Keep
the native competitor distances and frozen fallback for missing templates.
Use the species-specific observed pair weights. The class-I atlas supports
54 human and 31 mouse V-domain pairs, including 29 shared pairs. The
paired performance checks use the human weights and human queries.

Matrix derivation is research-only. Use the companion analysis repository’s
`scripts/vdjam/derive_matrices.py` and `matrix_models.py`, pinned to vdjtools>=5.9.0
for mixed-origin work. Read the
[origin-aware formalism](https://docs.isalgo.dev/vdjmatch/rearrangement-origin.html) before interpreting.
Keep penalties, full kernel scores and sample-level P-values separate.
Positive kernel support alone is not a calibrated precision threshold, and
nominal-radius pruning drops its outer-ladder contributions.

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
`--neighbor-draws 200` exports uniform same-N clonotype background draws from
VDJdb and optional control. One native graph batch per sampled population;
shared nested prefixes across N. These are background distributions, not
bootstrap CIs or knee-selected statistics.
For a fresh classification rebuild use `reference-profiles --vdjdb RAW
--output-prefix NEW --threads 4`. It rebuilds raw graphs/shuffles on every call;
`index` remains the reusable annotation-model preparation route. The fresh
manifest checks package source, actual native bytes and raw identities.
New profile nulls default to999label shuffles (`--shuffles`, minimum200).
Diagnostic matched-control and Chao-bootstrap defaults are200replicates.
Permutation/matched-control envelopes are not Efron bootstrap CIs.

Use the installed `vdjmatch annotate` command or Python `annotate_repertoires`.
Read [README.md](../../README.md) for input schemas and [CLI documentation](https://docs.isalgo.dev/vdjmatch/cli.html)
for advanced scorers. Run `vdjmatch annotate --help` to verify installed options.
Use released `seqtree>=1.3.2` and `vdjtools>=5.7.0` for native search, ingestion and abundance statistics.
Version 0.6.0 defaults to `--reference-population curated-benchmark`, the manuscript
primary population. Preserve `--reference-population all-valid` for the accepted
all-valid baseline. Curation excludes KLGGALQAK and PMID 40498839 source records;
exact hits remain included unless benchmark exclusion is explicitly requested.
The stem4 scoring formula is unchanged; neither the resolved-first repair nor
paired-score diagnostic is promoted by this reference-population change.

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
- Ordinary human and mouse annotation defaults to vX. Use one matching
  `--vx-null-dir` per primary locus for calibrated calls. Without a bank,
  `score` ranks evidence and `selection_status` marks ambiguous or unsupported
  queries. Paired inputs report marginal TRA/TRB evidence in separate ranking
  tables. Apply null cutoffs only to the primary selection. Use
  `--score-model kernel` to reproduce stem4 scores.
- Ordinary annotation retains exact hits and conflicting-label ambiguity.
  `--exclude-exact` is the explicit benchmark option.
- Compact output is the default. `--details` selects the kernel workflow and
  exports per-target evidence and candidates. Its memory grows with queries ×
  targets. `--epitope` restricts requested targets in the kernel workflow.
  `--kernel-scorer` selects retained advanced kernel variants;
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
- Ordinary annotation and public `Annotator` searches default to `threads=0`
  (allocated CPUs); a positive count selects a smaller budget. Import the API
  before initializing Polars. VDJtools owns shared runtime setup and resolves
  native batches to one thread inside its sample workers; excessive explicit
  counts and oversized preinitialized table pools raise. Keep queries in population
  batches; do not wrap search in a per-receptor multiprocessing pool.
- Use matrix-relative distance settings for comparisons. The shortlist API's
  omitted radius and stem4 reference-call geometry derive their scales from the
  selected matrix. Historical absolute overrides remain explicit reproduction
  options; a censored sparse distance is not an uncensored nearest distance.

## Further reading

Read [advanced workflows and contracts](references/advanced-workflows.md) for calibration, reference diagnostics, scope evidence and repertoire inference.

For paired-score producer digests and equal-tail recipe compatibility, see
[advanced workflows](references/advanced-workflows.md).
