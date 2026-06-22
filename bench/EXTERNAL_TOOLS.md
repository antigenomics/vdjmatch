# External comparators — sandboxing & open questions

Each non-`vdjmatch` method runs in its **own conda env** (isolated deps, often old TF/torch), writes a
standard predictions file, and `bench/compare.py` reads it. Standard output contract per method/dataset:

```
predictions/<method>/<dataset>.tsv     # columns: query_id, locus, epitope, score, significant(0/1)
```
`score` higher = more likely binder; `significant` is the method's own call (omit if it has none → we
threshold the score). One row per (query, candidate epitope) the method scores; queries with no call
may be omitted (counted as negatives).

Convention: `conda create -n cmp-<tool> ...`; a `bench/external/<tool>/run.sh` wrapper produces the
predictions file. **Where install is non-obvious, the questions below are for @mikessh.**

---

## TCRMatch  (PMID 33777034)
- **Install:** bioconda — `conda create -n cmp-tcrmatch -c bioconda -c conda-forge tcrmatch`. Self-
  contained C++; reference = IEDB CDR3β set it ships, or feed our VDJdb β as the reference.
- **I/O:** input = CDR3β list; output = nearest-reference matches + score (k-mer similarity 0–1).
- **Q1:** Use TCRMatch's bundled IEDB reference, or rebuild its reference from our leakage-removed VDJdb
  β so it's apples-to-apples? (Recommend the latter — fairness.) β-only, so TRA is N/A for it.

## tcrdist3  (Mayer-Blackwell 2021)
- **Install:** `conda create -n cmp-tcrdist python=3.8 && pip install tcrdist3` (pulls parasail, numba).
- **I/O:** V + CDR3 (α and/or β) → pairwise TCRdist; we do nearest-neighbour label transfer over the
  leakage-removed reference. Handles TRA, TRB, and paired.
- **Q2:** OK to define its "prediction" as the epitope of the nearest reference TCR (1-NN) with score =
  −TCRdist, or majority over k-NN? (Recommend 1-NN + a k-NN variant.)

## GLIPH2  (Huang 2020)
- **Install:** standalone binary (irepertoire) or the `gliph2` source; needs a reference CDR3 set + HLA.
  Not on conda — likely a downloaded binary.
- **Q3:** Do you have the GLIPH2 binary / reference set, or should we use the web tool? GLIPH clusters
  query+reference and we label-transfer within clusters — confirm that scoring scheme is acceptable.

## ERGO-II  (PMID 33981311)
- **Install:** `conda create -n cmp-ergo python=3.7` + their repo (PyTorch); ships pretrained VDJdb/McPAS
  models. Supervised TCR–peptide binding: input (CDR3, peptide) → P(bind).
- **Q4:** Use their pretrained VDJdb model (risk: trained on overlapping VDJdb → leakage we can't fully
  remove on their side), or retrain on our leakage-removed split? Retraining is more work but the only
  clean comparison — your call.

## ImmuneWatch DETECT  (commercial)
- **Install:** none — we only have the provided outputs in `test_data/immunedetect_results/*.ods`
  (read via the `ods` extra). Accuracy-only (no throughput).
- **Q5:** Which samples do the provided DETECT `.ods` cover, and what is their score/decision column?

## Throughput / RAM (for the speed half of Deliverable 3)
- **Q6:** Run all methods on the same machine (this M3, 32 GB) for wall-time + peak-RSS, or is there a
  shared cluster spec you want reported? sample5 (300k) is the scaling point — some tools may OOM/timeout;
  we record where.

---

## Methods you may want to add (placeholders)
NetTCR-2.2, pMTnet, MixTCRpred, STAPLER, TCRBert, DeepTCR, clusTCR, ATM-TCR. Tell me which to include and
I'll add a sandbox stanza + a `run.sh` for each.

## What I need from you
- Answers to Q1–Q6 above (or "use defaults" and I'll pick the recommended option).
- GLIPH2 binary + HLA reference, and the DETECT `.ods` mapping (Q3, Q5), if you have them.
- The final method list (I'll wire each into `compare.py` via its predictions file).
