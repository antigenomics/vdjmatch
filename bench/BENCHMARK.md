# vdjmatch comparison benchmark

Plan + protocol for the cross-method epitope-annotation benchmark (manuscript Deliverable 3). The
companion `EXTERNAL_TOOLS.md` covers sandboxing the non-`vdjmatch` methods; the CLI is `bench/compare.py`.

## Goal

For each method, report **per-locus (TRA / TRB), per-epitope** annotation quality as a *distribution*
over epitopes — **F1, PR-AUC, retention, purity** — drawn as **boxplots + beeswarm** across methods, on
the latest VDJdb release and its high-confidence shortlist, plus four held-out samples. We are
deliberately strict with ourselves: `vdjmatch` may **never** use an exact CDR3 match to a benchmark
query (those reference rows are removed; see *Leakage*).

## Datasets

Reference DB = **VDJdb 2026-06-11-ZENODO** (`db.fetch_hf`), with its **≥2-reference shortlist**
(`db.replicated`) as the gold standard. Query sets (in `2026-vdjmatch/test_data/`):

| set | file | chain(s) | role | truth |
|---|---|---|---|---|
| full VDJdb | (HF release) | TRA, TRB | leave-one-out-by-epitope retrieval | epitope label |
| shortlist | (HF, ≥2 ref) | TRA, TRB | gold-standard LOO | epitope label |
| sample1 | `sample1_airr.txt` (489) | TRB | NLVPMVATV (CMV) sorted+ | **positives** (all NLV) |
| sample2 | `sample2_airr.txt` (319) | TRB | LLWNGPMAV (YF) + LLL… (BST-2 neoag) sorted+ | **positives** |
| sample3 | `sample3_airr.txt` / `sample3_vdjdb.txt` | TRA,TRB | VDJdb 2025 (older release) | epitope label |
| sample5 | `sample5_olga_airr.txt` (300k) | TRB | random OLGA repertoire | **negatives** (no hits expected) |
| sample6 | `sample6_TCRvdb.csv` (3693, paired) | TRA+TRB | TCRvdb, 2 epitopes | **pos/neg by `padj < 1e-5`** |

`sample6_vdb_airr.txt` is the single-chain (β) view of the same TCRvdb data; `sample6_TCRvdb.csv` is the
paired (α+β) table with `cdr3_alpha_aa`, `cdr3_beta_aa`, `epitope_aa`, `log2FoldChange`, `pvalue`,
`padj` — used for the **paired-chain E-value** test (positives = enriched binders at `padj < 1e-5`).

## Leakage removal (be fair to ourselves)

Before annotating any benchmark query with `vdjmatch`, drop from the reference DB **every row whose
CDR3 (per locus) exactly matches a query CDR3**, and for full-VDJdb / shortlist LOO also exclude the
held-out epitope's own clonotypes. So a correct call must come from *other* TCRs, never a database
copy of the query. (sample5 OLGA negatives need no removal — they shouldn't be in VDJdb at all; if any
are, that is itself a measured false-positive.)

## Methods (algorithm list — expand later)

Pluggable `Method` in `compare.py`; external ones read a standard predictions file produced in their
own conda sandbox (`EXTERNAL_TOOLS.md`).

- **vdjmatch** (ours): control-calibrated E-value over the leakage-removed reference; paired E-value for
  sample6.
- **TCRMatch** — CDR3β k-mer similarity to IEDB/VDJdb (PMID 33777034).
- **tcrdist3** — TCRdist nearest-neighbour over V+CDR3 (Mayer-Blackwell 2021).
- **GLIPH2** — motif/global clustering, then label transfer (Huang 2020).
- **ERGO-II** — supervised TCR–peptide binding (LSTM/AE; PMID 33981311).
- **ImmuneWatch DETECT** — commercial; accuracy-only from provided outputs (`immunedetect_results/`).
- *(placeholders for the user to add: NetTCR-2.2, pMTnet, STAPLER, TCRBert, MixTCRpred, …)*

A "method" maps (query CDR3[s], V/J, locus) → per-epitope score (higher = more likely binder), plus a
significance/decision flag where it has one.

## Metrics (per locus × per epitope)

For epitope *E* and locus *L*, restrict to queries of locus *L*; the binary task is "does this query
bind *E*". Let the method give each query a score for *E* and a binary call (its significant set).

- **PR-AUC** — area under precision–recall, ranking queries by the *E*-score (prevalence-aware;
  `metrics.pr_auc`; balanced variant for cross-method comparability where prevalence differs).
- **F1** — at the method's own decision threshold (E-value significant / score cutoff): `2PR/(P+R)`.
- **purity** — precision of the method's *E*-set: fraction of queries it calls *E* that truly bind *E*
  (set homogeneity).
- **retention** — recall after the significance filter: fraction of true-*E* queries that keep ≥1
  significant same-*E* hit. (Purity = "are the calls clean", retention = "did we keep the true ones".)

Only epitopes with ≥ N (default 30) labelled queries in that locus enter the distribution. sample5 has
no true epitope → scored as **specificity / FPR** (fraction of OLGA queries wrongly called significant).

## Plots

Per metric: one panel per locus, x = method, y = per-epitope value, **boxplot (quartiles) + beeswarm
(one dot per epitope)** overlaid (seaborn `boxplot` + `swarmplot`). Paired Wilcoxon vdjmatch-vs-each
across the shared epitopes for the stats line.

## Paired-chain (sample6 TCRvdb)

Run `vdjmatch` paired α+β E-values; positives = `padj < 1e-5`, negatives = the rest. Report ROC-AUC /
PR-AUC of the paired E-value as a binder classifier, and check paired E < min(single-chain E) for true
pairs. Single-chain β baseline from `sample6_vdb_airr.txt`.

## CLI

```fish
# vdjmatch arm on the annotated samples (leakage-removed reference)
python bench/compare.py --methods vdjmatch --datasets sample1 sample2 sample5 --out bench/out
# add external methods once their predictions/<method>/<dataset>.tsv files exist
python bench/compare.py --methods vdjmatch tcrmatch --datasets shortlist --out bench/out
```
Emits `bench/out/results.tsv` (tidy: method, dataset, locus, epitope, metric, value) + the boxplot/
beeswarm PNGs.

## Figures/tables selected for the appendix (now)

Keep the small set that carries the argument; cut the rest (text update later, from `docs/`):

- **Fig — purity vs distance** (Hamming-1 optimum) — keep.
- **Fig — matrices bar** (no AA matrix beats BLOSUM62; VDJAMr ties; possig wins) — keep.
- **Fig — position significance** (end-anchored profile) — keep.
- **Fig — vgene_scan** (V prior universal) — keep.
- **Table — matrices** (balanced PR-AUC) and **Table — vregion** (near-exact recovery) — keep.
- **This benchmark's boxplot+beeswarm** (cross-method, per-epitope) — the headline of Deliverable 3.
- Cut/condense: tra_trb_gly, region_corr, retention, loo_prauc, vgene/vgene-sim/vgene-scan tables,
  shortlist table → fold into one or two summary panels.

## Open questions

See `EXTERNAL_TOOLS.md` for per-tool install/run questions.
