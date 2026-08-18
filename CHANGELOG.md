# Changelog

## [Unreleased]

### Added
- **`load_compendium()` and `load_estimates()`** — the naive precursor frequencies other
  laboratories measured, and what this library predicts for the same epitopes. 147 published
  records from 14 studies tracing to **31 originating laboratories**, harmonized to cells per
  million of the *naive* compartment, human and mouse, CD4 and CD8. Both bootstrap from
  `isalgo/airr_benchmark` (`vdjmatch/precursor_freq/`) and cache beside the VDJdb releases, so an
  installed user can reproduce the validation without checking out the analysis repository. They
  stay two tables on purpose: one is what a laboratory counted, the other a mass under a generative
  model, and merging them would conflate a measurement with a prediction of it.
- **`docs/notebooks/precursor_compendium.py`** — a marimo notebook walking the comparison:
  prediction against measurement per epitope, the `Pgen` spectrum inside one cognate set, the ball
  census of uncatalogued neighbours, and the step from a mass to a number of cells. Install with
  `pip install 'vdjmatch[precursor,notebook]'`.
- **`docs/compendium.rst`** — why unseen *mass* converges where unseen *count* does not, what the
  ball overlap measures, and why a set total is the wrong functional.

### Fixed
- **`--selection auto` help quoted the superseded per-chain factors.** It named TRB 4.42 / TRA 1.05
  while `SELECTION_BY_CHAIN` has shipped **4.62 / 1.07** since 0.3.0, so the one place a user reads
  the constants off disagreed with the one place the code uses them.
- **`ALICE_Q`'s docstring quoted a retired comparison.** It put the model-free event ratio at a
  median 14.8x; against `union_mass` at radius 1 — the comparable quantity — that is **8.8x and
  9.9x** on two beta cohorts and **4.5x** on alpha, so beta is the same order as ALICE's 9.41 and
  alpha is not.
- **`test_cli_precursor_switches_source_for_non_human_organism` needed the network.** It exercised
  the source switch through `--vdjdb`, which downloads the database, so on a machine with no
  network it failed with `URLError` before reaching the switch it was testing. It now reads a
  sample path instead and is hermetic.

## [0.3.1] - 2026-08-18

### Fixed
- **`vdjmatch match` crashed on its own default flags.** `--asset full` fetches `vdjdb_full.txt`,
  which is one row per **complex** — the chains sit side by side in `cdr3.alpha`/`cdr3.beta` and
  there is no `gene` column — while `slim`/`default` are one row per chain. `schema.normalize()`
  assumed a single layout and filled the four missing columns with nulls "so downstream code can
  rely on the schema", so the whole table arrived with **every `cdr3`, `v`, `j` and `gene` null**:
  176,130 human records and 14,852 mouse ones, all unusable. The matcher then died inside seqtree
  with `TypeError: build(): incompatible function arguments ... Invoked with: [None]`, which names
  neither the asset nor the schema. The paired layout is now unpivoted back to one row per chain on
  load, so `--asset full` yields 286,013 usable records (29,670 alpha-only + 69,823 beta-only +
  93,260 paired, counted independently of the loader) and `complex_id` is synthesized VDJdb's own
  way — shared and non-zero across the two chains of a complex, `0` for a single-chain record — so
  `paired_only=True` works against `full` for the first time.

## [0.3.0] - 2026-08-18

### Added
- **`n_zero_pgen`** in the `precursor` output — junctions that pass the anchor check and still score
  `Pgen` **exactly zero**. Such a junction contributes nothing to any mass and raises nothing, so
  until now a cognate set could be quietly short. It is not rare: under the default `olga` model
  **5.7% of human TRA junctions in VDJdb** score zero, reaching **13.4%** for `YVLDHLIVV` and
  **15.5%** for `GLCTLVAML`. The CLI warns when the overall rate exceeds 1%.

### Fixed
- **`vdjmatch precursor` refuses a species/organism mismatch.** `--species MusMusculus` against the
  default `--organism human` used to run to completion and emit a plausible number for every group:
  a median **0.157x** of the right answer with an **11x spread across epitopes**, so the ranking was
  corrupted too and nothing in the output said so. Now an error naming the correct invocation.
- **`--organism mouse` selects `arda` automatically.** `olga` ships no mouse model, so the default
  source could only ever fail; there is nothing else the user could have meant.
- **`load_model`'s warning about `source="learned"` was stale and wrong.** It claimed 68 of 89 TRB V
  alleles have `P(V) = 0`; the current bundled set has **zero** such alleles, and `learned` loses
  fewer junctions than the default `olga` on both human loci. The docstring now describes what the
  three sets actually are — `olga` a bit-faithful import of OLGA's published models (defects
  included, hence the deletion-grid zeros), `learned` and `arda` refits from real 5'RACE reads, and
  `arda` the only set carrying mouse — and says which to prefer for a mass over a set.

### Changed
- `SELECTION_BY_CHAIN` refined to **TRB 4.62 / TRA 1.07** (was 4.42 / 1.05) after adding
  `airr_covid19_vacc` as a third cohort — TRB now rests on three independent cohorts (4.75, 4.92,
  3.79) and TRA on two (1.05, 1.09), with residual RMS 0.054 and 0.067 decades over n=777 and
  n=336. At matched depth the replication is 1.04x on both chains. The shift is inside the residual
  spread, so 0.2.0's constants were not wrong, just fitted on less data. Stale copies of the old
  values in `README.md` and `docs/cli.rst` corrected.

## [0.2.0] - 2026-08-17

### Added
- **`vdjmatch.precursor`** — T-cell precursor frequency for an epitope, moved here from
  `mhcmatch.precursor` (which now re-exports it) and given a CLI. Estimators: `event_ratio`
  (model-free, events over events), `observed_mass`, `coverage_corrected_mass`, `union_mass`,
  `shell_profile`, `motif_mass`, `unseen_junctions`, `cross_check`, `precursor_frequency`.
  Needs the new optional extra: `pip install 'vdjmatch[precursor]'`.
- **`union_mass` — the exact union of Hamming-`r` balls, without enumerating it.** The naive sum
  counts every ball member `cov(x)` times, so
  `m(union) = sum_a m(B_r(a)) - sum_{x: cov(x)>=2} (cov(x)-1) Pgen(x)`, exactly. No
  inclusion–exclusion, hence no truncation error — truncating I–E at pairs and triples is unsafe,
  a four-junction component has a non-empty four-way term. Only centres inside one connected
  component of the `2r` graph can contribute, and only their multiply-covered members need a `Pgen`
  call, so singleton components (41.8% of VDJdb human TRB junctions) cost nothing. Measured against
  the enumeration oracle at rel ~1e-15; 13–186x faster at `r=2`, and 20 s for NLVPMVATV's 13,338
  TRB junctions where enumeration needs 3.3M `Pgen` calls.
- **`closed_ball_mass` — closed-form ball mass at any radius**, by an alternating sum over
  wildcarded motifs against vdjtools' masked transfer-matrix DP:
  `m(B_r(a)) = sum_k (-1)^(r-k) C(L-k-1, r-k) sum_{|S|=k} m(W_S)`. 106 DP passes against 33,117
  enumerated sequences at `L=14, r=2`.
- **`unseen_junctions`** — the richness counterpart of the coverage correction: how many cognate
  junctions were never catalogued and how rare they are. Reported two ways, because they answer
  differently: the **ball census** (`n_unseen_ball`, finite and exact given the ball) and the
  **Horvitz–Thompson** extrapolation (`n_unseen_ht`, which diverges in the tail — `richness_reliable`
  flags when the rarest observed junction was captured too seldom for the count to mean anything).
  The *mass* converges either way; the *count* does not.
- `precursor_frequency` / `expected_cells` / `p_at_least` / `paired_frequency` — from a set of
  `Pgen` values to a frequency, an expected precursor count at `n_cells`, and `P(>=k precursors)`.
  `F(e)` answers "is there a precursor"; `p_ge_k` is what a detectable response actually needs.
- **`occupancy` — seen and unseen clonotype counts.** `S = sum_k alpha^k n_k`,
  `n_seen(N) = sum_k alpha^k n_k E_k[1 - exp(-N Pgen)]`, `n_unseen = S - n_seen`. Replaces the summed
  mass for anything count-shaped: under the size bias a set total has expectation
  `n * (sum pi^2) / F`, i.e. proportional to how many TCRs were catalogued, proportional to
  concentration, and *inversely* proportional to what it estimates. The occupancy form has no such
  term and **saturates** in depth. Validated out-of-sample with zero free parameters against how many
  distinct junctions real cohorts show — `airr_covid19` TRA observed/predicted median 1.015
  (IQR 0.93–1.10, n=168), TRB 2.043 (n=259), `airr_hip` TRB 1.510 (n=259).
- **`SELECTION_BY_CHAIN`** — measured per-chain depth factors, TRB 4.42 and TRA 1.05, fitted through
  the saturating curve. The two TRB cohorts agree at 4.75 and 3.79 despite differing in protocol and
  8.5x in depth. The chains do not share a value and a joint fit suits neither. Distinct from
  `ALICE_Q`, which scales a probability to a frequency rather than an observer's depth.
- `vdjmatch precursor` CLI, `--vdjdb` to run natively on the whole database, `--selection auto` for
  the measured constants.
- Registered the `slow` pytest marker that CI's `-m "not slow"` had been filtering on unregistered.

### Changed
- `shell_profile` obtains shells by **differencing unions** rather than enumerating them, so it has
  no memory ceiling: the `r=2` profile that cost ~9.9M materialised strings for 300 junctions is now
  `r+1` closed-form calls. It no longer takes `max_members`; that guard moved to `union_mass`, where
  it applies per connected component. Shell `n` is `None` when the union is too large to census —
  the masses never are.
- `seqtree>=0.7.0` (was `>=0.6.1`) for `neighbourhood` / `neighbourhood_union` / `union_size`.
- `__version__` now reads the installed distribution metadata instead of a hand-maintained literal,
  which had drifted (`0.1.0` in code against `0.1.2` in `pyproject.toml`).

## [0.1.2] - 2026-07-30

### Fixed
- **`cluster.overlap()` raised `SchemaError` whenever a query had zero fuzzy matches** —
  regression from 0.1.1's index-only hit rows: with no hits, `pl.DataFrame([], schema=[names])`
  has no data to infer a dtype from and defaulted `a_idx`/`b_idx` to `Null`, which then failed to
  join against the lookup frames' `Int64` idx columns. `hits` now declares an explicit
  `{name: dtype}` schema so the empty and non-empty paths agree. Added a regression test
  (`test_overlap_no_matches_returns_empty_not_schema_error`) — the prior release shipped without
  one, since the existing `test_overlap_within_drops_self` fixture always has at least one hit.

## [0.1.1] - 2026-07-30

### Fixed
- `cluster.overlap()` no longer copies full CDR3 strings into every hit row; hit rows are now
  built from indices only and expanded to `a_cdr3`/`b_cdr3` via a vectorized join afterward
  (same pattern as `match.engine.VdjdbIndex.annotate()`). Cuts peak memory on dense fuzzy-match
  results (e.g. BCR light-chain overlap queries against large reference sets).

### Changed
- Requires `seqtree>=0.6.1` (was `>=0.4.0`), which corrects the MJ `structural` matrix's A-N
  contact energy (`0.00` → `0.15`); scores from `structural()` change for sequences containing
  A, W, or B.

## [0.1.0] - prior release

- Initial PyPI release; requires `seqtree>=0.4.0`.
