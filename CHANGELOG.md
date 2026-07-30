# Changelog

## [0.1.1] - 2026-07-30

### Fixed
- `cluster.overlap()` no longer copies full CDR3 strings into every hit row; hit rows are now
  built from indices only and expanded to `a_cdr3`/`b_cdr3` via a vectorized join afterward
  (same pattern as `match.engine.VdjdbIndex.annotate()`). Cuts peak memory on dense fuzzy-match
  results (e.g. BCR light-chain overlap queries against large reference sets).

## [0.1.0] - prior release

- Initial PyPI release; requires `seqtree>=0.4.0`.
