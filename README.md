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

## Scoring: what works (and what doesn't)

An empirical study on VDJdb (see `appendix/scoring.tex`) settles the scoring question honestly:

- **Hamming distance 1 is the signal:noise optimum** — neighbour purity falls 0.53 → 0.13 across edit
  distance 1–5 (reproducing the original VDJdb observation, Shugay et al. NAR 2018). The search-ball
  radius, not the substitution matrix, is the dominant lever.
- **Central (NDN) substitutions carry the specificity signal** — a mismatch in the CDR3 core most
  often changes specificity (P(same epitope) ≈ 0.41) while V/J-border mismatches are germline noise
  (≈ 0.70); the NDN core is also ~30–36% glycine (insertion/D-gene signature).
- **No substitution matrix beats BLOSUM62** for epitope retrieval (BLOSUM62 ≈ PAM250 > structural >
  Hamming > the data-derived VDJAM; region weighting helps VDJAM a little but not past BLOSUM). The
  matrix is a second-order lever; the first-order statistic is the control-calibrated E-value.

## License

GPL-3.0-or-later (it builds on `seqtree`, which is GPL-3.0-or-later).
