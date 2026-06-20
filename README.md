# vdjmatch

Fast, control-calibrated annotation of **T-cell receptor antigen specificity**.

`vdjmatch` annotates clonotypes in large AIRR repertoires against [VDJdb](https://github.com/antigenomics/vdjdb-db)
by fuzzy CDR3 search, reporting a **control-calibrated E-value** (BLAST-style significance against a
background repertoire) and enriched antigen-specificity labels. It is a Python rewrite of the legacy
Java/Groovy vdjmatch, built on the [`seqtree`](https://github.com/antigenomics/seqtree) search core.

> **Status:** 2.0 alpha, under active development on `dev`. The legacy Java tool is preserved on the
> `legacy-java` branch (tags `1.1.4`–`1.3.1`).

## Features (target)

- Fetch the latest VDJdb release and annotate AIRR Rearrangement / Cell (paired α/β) samples.
- Extremely fast, multithreaded search of million-scale repertoires (via `seqtree`).
- Control-calibrated **E-values** (single-chain now; paired α/β and single-chain-paired estimates).
- Custom substitution matrices, including segment-specific (V / NDN / J) scoring; the TCR-specific
  **VDJAM** matrix is bundled.
- Rich per-hit output: ranked hits, CIGAR + alignment match/gap, alignment scores, E-values.
- Epitope-level enrichment summaries; pairwise sample overlap.
- `polars` throughout for I/O.

## Install (development)

```fish
python -m venv .venv
source .venv/bin/activate.fish
pip install -e .[test,bench]
```

`seqtree` (the search engine) is installed from PyPI as a dependency.

## License

GPL-3.0-or-later (it builds on `seqtree`, which is GPL-3.0-or-later).
