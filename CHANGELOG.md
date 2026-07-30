# Changelog

## [0.1.1] - 2026-07-30

### Fixed
- `cluster.overlap()` no longer copies full CDR3 strings into every hit row; hit rows are now
  built from indices only and expanded to `a_cdr3`/`b_cdr3` via a vectorized join afterward
  (same pattern as `match.engine.VdjdbIndex.annotate()`). Cuts peak memory on dense fuzzy-match
  results (e.g. BCR light-chain overlap queries against large reference sets).

### Changed
- Requires `seqtree>=0.6.0` (was `>=0.4.0`), which corrects the MJ `structural` matrix's A-N
  contact energy (`0.00` → `0.15`); scores from `structural()` change for sequences containing
  A, W, or B.

## [0.1.0] - prior release

- Initial PyPI release; requires `seqtree>=0.4.0`.
