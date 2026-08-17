# Changelog

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
- `vdjmatch precursor` CLI, `--vdjdb` to run natively on the whole database.
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
