# vdjmatch

<p>
  <a href="https://pypi.org/project/vdjmatch/"><img alt="PyPI" src="https://img.shields.io/pypi/v/vdjmatch"></a>
  <a href="https://github.com/antigenomics/vdjmatch/actions/workflows/tests.yml"><img alt="tests" src="https://github.com/antigenomics/vdjmatch/actions/workflows/tests.yml/badge.svg"></a>
  <a href="https://docs.isalgo.dev/vdjmatch/"><img alt="docs" src="https://github.com/antigenomics/vdjmatch/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-GPLv3-green"></a>
</p>

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

An empirical study on VDJdb (see `appendix/vdjmatch_scoring.tex`; regenerated on the **2026-06-11-ZENODO**
release with composition controls and balanced metrics) settles the scoring question honestly:

- **Hamming distance 1 is the signal:noise optimum** — macro purity (per-epitope mean) falls
  0.44 → 0.07 across edit distance 1–5, a 44× → 2.3× enrichment over chance (reproducing the original
  VDJdb observation, Shugay et al. NAR 2018). The search-ball radius, not the substitution matrix, is
  the dominant lever. (On the dense 2026 release this only shows up once a few 10× mega-studies are
  capped and the random tail is treated as an admixed control — naive pooled purity reads a flat ~0.9.)
- **Central (NDN) substitutions carry the specificity signal** — a mismatch in the CDR3 core most
  often changes specificity (P(same epitope) ≈ 0.31) while V/J-border mismatches are germline noise
  (≈ 0.90); the NDN core is also ~31–34% glycine (insertion/D-gene signature).
- **No amino-acid matrix clearly beats BLOSUM62 — and a genetic-code null ties it.** BLOSUM62 ≈ PAM250 ≈
  structural > Hamming > data-derived VDJAM. Strikingly, **VDJAMr** — a matrix built from the *genetic
  code alone* (how mutationally accessible one AA is from another; `loo_vdjam.codon_dissim`) — matches
  BLOSUM62 (0.585 vs 0.572 @≤2), so TCR CDR3 substitution structure is **generative, not chemical**.
  A published TCR-specific matrix, [tcrBLOSUM](https://doi.org/10.1093/bib/bbae602), does *not* beat
  BLOSUM62 either (0.557 vs 0.570; `bench/tcrblosum_refute.py`) — its same-epitope counts don't transfer.
- **Position-weighting BLOSUM62 does beat it.** Encoding the central-substitution finding as a
  seqtree positional matrix (`PositionalMatrix.from_weights(BLOSUM62, …)`, centre ~2× the V/J borders)
  raises leave-one-out retrieval (balanced PR-AUC) above flat BLOSUM62 in **7/8 held-out epitopes**
  (0.572 → 0.623 @≤2; 0.589 → 0.617 @≤4). For CDR3, *where* a mismatch falls matters more than *which*
  residue it is. The first-order statistic is still the control-calibrated E-value.
- **The V gene is a strong, near-binary prior — recovered only at near-exact germline identity.**
  Same-V neighbours share the epitope **~44–64%** of the time vs **~6–17%** cross-V (ratio up to 8×).
  *Loose* CDR1/CDR2 similarity barely predicts it (point-biserial r ≈ 0), but cross-V co-specificity
  **rises monotonically as germline CDR1+CDR2 approach identity** — 11% at ≥6 mismatches up to ~60% at
  edit-0, approaching the same-V level — at a tight, CDR3-like tolerance, recovered *whole-loop* rather
  than via a sparse pseudosequence (per-position lift flat except CDR2 pos 5). So the V contribution
  lives in the germline contacting loops but needs near-exact match; a useful soft-V match must demand
  it (`bench/vregion_decompose.py`, `bench/vpseudo.py`; helpers in `vdjmatch.match.vgene`).

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

GPL-3.0-or-later (it builds on `seqtree`, which is GPL-3.0-or-later).
