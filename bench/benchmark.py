#!/usr/bin/env python3
"""Comprehensive cross-method TCR-specificity benchmark (see bench/BENCHMARK.md).

Conditions (classification -> ROC / PR / FP / F1, at each method's own significance call):
  C1 sample-LLW/LLL : sample2 LLW(pos) vs LLL(neg) and reverse;  ref = vdjdb2026 (exact-filtered); TRB
  C2 sample-NLV     : sample1 cmv(NLV+) vs control(NLV-);        ref = vdjdb2026;  TRA, TRB
  C3 tcrvdb         : padj<1e-5 (pos) vs >=1e-5 (neg), GLC+YLQ;  ref = vdjdb2026;  TRA, TRB, paired
  C4 vdjdb2025-LOO  : epitope X(pos) vs other(neg), within 2025; ref = 2025 \\ self; TRA, TRB, paired
  C5 new2026-LOO    : new-in-2026 X(pos) vs other-new(neg);      ref = vdjdb2025 (temporal); TRA, TRB, paired
  C6 olga-FP        : OLGA TRA(sample5)+TRB(sample4), 10k;       ref = vdjdb2025;  any hit = FP
Clustering (retention / purity): vdjdb2026 self-cluster (separate, methods that cluster).

Per-method per-condition: a by-epitope table and a summary row (mean +- sd across epitopes).
Phase A = framework + vdjmatch; external methods (tcrmatch / tcrdist / ERGO-II / GLIPH2 / DETECT) plug
in via predictions/<method>/<condition>_<locus>.tsv (query_id, epitope, score[, significant]).

    python bench/benchmark.py --conditions C1 C2 C6 --methods vdjmatch
"""
from __future__ import annotations

import argparse
import math
import statistics as st
import sys
from collections import Counter, defaultdict
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench
from compare import TESTDATA, _big_epis, _top_held, load_sample
from metrics import pr_auc_balanced, roc_auc
from vdjmatch import db
from vdjmatch.evalue import background, first_hit
from seqtree import Index

EPI = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV",
       "GLC": "GLCTLVAML", "YLQ": "YLQPRTFLL"}
_REL = {}  # release cache


def release(which: str, species="HomoSapiens") -> pl.DataFrame:
    if which not in _REL:
        if which == "vdjdb2026":
            _REL[which] = db.load(_bench.source(), species=species)
        else:                                                     # vdjdb2025: prefer cached file (no API call)
            from vdjmatch.db.vdjdb import cache_dir
            p = cache_dir(None) / "vdjdb-2025-12-29.default.txt"
            _REL[which] = (db.load(str(p), species=species) if p.exists()
                           else db.load(pin="2025-12-29", asset="default", species=species))
    return _REL[which]


def ref_index(df: pl.DataFrame, locus: str):
    """VDJdb frame -> (Index, ref_epitopes, N_per_epitope) for a single locus."""
    r = (_bench.valid_cdr3(df.filter(pl.col("gene") == locus)).group_by("cdr3")
         .agg(pl.col("epitope").first()))
    epi = r["epitope"].to_list()
    return Index.build(r["cdr3"].to_list(), "aa"), epi, Counter(epi)


# ---- vdjmatch classification (first-hit E-value) --------------------------------------------------
def vdjmatch_classify(tgt, ref_epi, n_epi, ctrl, queries, epitopes, alpha, exclude_exact):
    """queries -> per query {epitope: (-log10 p_enrichment, significant)} for the candidate epitopes,
    plus an overall-significant flag (any-epitope first-hit p < alpha) for FP estimation."""
    N, M = len(tgt), len(ctrl)
    th, cc = first_hit.scan(tgt, ref_epi, ctrl, queries, exclude_exact=exclude_exact)
    scores, overall = {}, {}
    for q, t, c in zip(queries, th, cc):
        overall[q] = first_hit.pvalue(t, c, N, M)["p_enrichment"] < alpha
        scores[q] = {}
        for e in epitopes:
            p = first_hit.pvalue(t, c, n_epi.get(e, 1), M, epitope=e)["p_enrichment"]
            scores[q][e] = (-math.log10(max(p, 1e-300)), p < alpha)
    return scores, overall


def classify_metrics(scores, pos, neg, epi):
    """ROC / PR (threshold-free) + FP / F1 (at the method's significance call) for one epitope."""
    pairs = [(1, scores[q][epi][0]) for q in pos] + [(0, scores[q][epi][0]) for q in neg]
    tp = sum(scores[q][epi][1] for q in pos)
    fp = sum(scores[q][epi][1] for q in neg)
    fn, tn = len(pos) - tp, len(neg) - fp
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    return {"roc_auc": roc_auc(pairs), "pr_auc": pr_auc_balanced(pairs),
            "fp": fp / (fp + tn) if fp + tn else 0.0,
            "f1": 2 * prec * rec / (prec + rec) if prec + rec else 0.0,
            "n_pos": len(pos), "n_neg": len(neg)}


# ---- conditions: each yields (locus, ref_df, queries, [(epitope, pos_ids, neg_ids)], exclude_exact)-
def cond_sample_pair(name, sample, epi_keys, loci):
    """C1/C2: discriminate two query groups by matching against vdjdb2026."""
    for locus in loci:
        if name == "C1":                                          # sample2: LLW vs LLL by antigen.epitope
            d = (pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t")
                 .rename({"antigen.epitope": "label"}).select("cdr3", "label"))
            d = _bench.valid_cdr3(d).unique("cdr3")
            groups = {EPI["LLW"]: d.filter(pl.col("label") == EPI["LLW"])["cdr3"].to_list(),
                      EPI["LLL"]: d.filter(pl.col("label") == EPI["LLL"])["cdr3"].to_list()}
            tasks = [(EPI["LLW"], groups[EPI["LLW"]], groups[EPI["LLL"]]),
                     (EPI["LLL"], groups[EPI["LLL"]], groups[EPI["LLW"]])]
        else:                                                     # C2 sample1: cmv(NLV+) vs control(NLV-)
            d = (pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
                 .filter(pl.col("gene") == locus).select("cdr3", label="type"))
            d = _bench.valid_cdr3(d).unique("cdr3")
            pos = d.filter(pl.col("label") == "cmv")["cdr3"].to_list()
            neg = d.filter(pl.col("label") == "control")["cdr3"].to_list()
            tasks = [(EPI["NLV"], pos, neg)]
        yield locus, release("vdjdb2026"), tasks, True


def cond_olga_fp(loci, n=10000):
    """C6: OLGA negatives (sample5=TRA, sample4=TRB) vs vdjdb2025 -> any significant hit is a FP."""
    files = {"TRA": "sample5", "TRB": "sample4"}
    for locus in loci:
        q = load_sample(files[locus])
        if q.height > n:
            q = q.sample(n, seed=0)
        yield locus, release("vdjdb2025"), q["cdr3"].to_list()


def cond_tcrvdb(loci):
    """C3: TCRvdb padj<1e-5 (pos) vs >=1e-5 (neg) per epitope, annotated vs vdjdb2026."""
    t = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv").with_columns(pos=pl.col("padj") < 1e-5)
    chain_col = {"TRA": "cdr3_alpha_aa", "TRB": "cdr3_beta_aa"}
    for locus in loci:
        d = (t.select(cdr3=chain_col[locus], epitope="epitope_aa", pos="pos")
             .pipe(_bench.valid_cdr3).unique("cdr3"))
        tasks = []
        for e in (EPI["GLC"], EPI["YLQ"]):
            de = d.filter(pl.col("epitope") == e)
            pos = de.filter("pos")["cdr3"].to_list()
            neg = de.filter(~pl.col("pos"))["cdr3"].to_list()
            if pos and neg:
                tasks.append((e, pos, neg))
        yield locus, release("vdjdb2026"), tasks, True


def cond_loo(which, loci, top, min_epi, max_q):
    """C4 (within vdjdb2025) / C5 (new-in-2026 vs 2025 reference): per-epitope X(pos) vs other(neg)."""
    v25 = release("vdjdb2025")
    for locus in loci:
        ref_cell = _bench.valid_cdr3(v25.filter(pl.col("gene") == locus))
        if which == "C4":
            test_cell, excl = _bench.long_list(v25.filter(pl.col("gene") == locus), cap=max_q,
                                               min_n=min_epi), True
        else:                                                     # C5: 2026 clonotypes new vs 2025
            v26 = _bench.long_list(release("vdjdb2026").filter(pl.col("gene") == locus), cap=max_q,
                                   min_n=min_epi)
            test_cell = v26.join(ref_cell.select("cdr3", "epitope").unique(),
                                 on=["cdr3", "epitope"], how="anti")
            excl = False
        held = _top_held(test_cell, _big_epis(ref_cell, min_epi), top, min_epi)
        per = {e: test_cell.filter(pl.col("epitope") == e).unique("cdr3")["cdr3"].to_list()[:max_q]
               for e in held}
        tasks = [(e, per[e], [q for e2 in held if e2 != e for q in per[e2]]) for e in held]
        yield locus, v25, tasks, excl


def summarize(rows):
    """rows: (condition, locus, method, epitope, metric, value) -> by-epitope df + mean+-sd summary."""
    by = pl.DataFrame(rows, schema=["condition", "locus", "method", "epitope", "metric", "value"],
                      orient="row")
    summ = (by.group_by(["condition", "locus", "method", "metric"])
            .agg(pl.col("value").mean().round(3).alias("mean"),
                 pl.col("value").std().round(3).alias("sd"),
                 pl.col("epitope").n_unique().alias("n_epi"))
            .sort(["condition", "locus", "method", "metric"]))
    return by, summ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conditions", nargs="+", default=["C1", "C2", "C3", "C4", "C5", "C6"])
    ap.add_argument("--methods", nargs="+", default=["vdjmatch"])
    ap.add_argument("--loci", nargs="+", default=["TRA", "TRB"])
    ap.add_argument("--alpha", type=float, default=1e-3)
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--min-epi", type=int, default=30)
    ap.add_argument("--max-queries", type=int, default=300)
    ap.add_argument("--out", default="bench/out/benchmark")
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rows, fp_rows = [], []

    def classify_cells(cond):
        if cond == "C1":
            return cond_sample_pair("C1", None, None, ["TRB"])
        if cond == "C2":
            return cond_sample_pair("C2", None, None, args.loci)
        if cond == "C3":
            return cond_tcrvdb(args.loci)
        if cond in ("C4", "C5"):
            return cond_loo(cond, args.loci, args.top, args.min_epi, args.max_queries)
        return None

    for cond in args.conditions:
        if cond == "C6":
            for locus, ref_df, queries in cond_olga_fp(args.loci):
                tgt, ref_epi, n_epi = ref_index(ref_df, locus)
                if "vdjmatch" in args.methods:
                    _, overall = vdjmatch_classify(tgt, ref_epi, n_epi, background(locus), queries,
                                                   [], args.alpha, True)
                    fp = sum(overall.values()) / len(queries)
                    fp_rows.append(("C6", locus, "vdjmatch", fp, len(queries)))
                    print(f"C6/{locus}: FP={fp*100:.2f}% over {len(queries)} OLGA negatives")
            continue
        cells = classify_cells(cond)
        if cells is None:
            print(f"  ({cond} not yet implemented)")
            continue
        for locus, ref_df, tasks, excl in cells:
            if not tasks:
                print(f"{cond}/{locus}: no tasks (skipped)")
                continue
            tgt, ref_epi, n_epi = ref_index(ref_df, locus)
            epitopes = sorted({e for e, _, _ in tasks})
            allq = sorted({q for _, p, n in tasks for q in (*p, *n)})
            if "vdjmatch" in args.methods:
                scores, _ = vdjmatch_classify(tgt, ref_epi, n_epi, background(locus), allq,
                                              epitopes, args.alpha, excl)
                for epi, pos, neg in tasks:
                    m = classify_metrics(scores, pos, neg, epi)
                    for k in ("roc_auc", "pr_auc", "fp", "f1"):
                        rows.append((cond, locus, "vdjmatch", epi, k, m[k]))
                print(f"{cond}/{locus}: {len(tasks)} epitope-task(s), {len(allq)} queries scored")

    if rows:
        by, summ = summarize(rows)
        by.write_csv(out / "by_epitope.tsv", separator="\t")
        summ.write_csv(out / "summary.tsv", separator="\t")
        with pl.Config(tbl_rows=100):
            print("\n=== summary (mean +- sd across epitopes) ===\n", summ)
    if fp_rows:
        fpdf = pl.DataFrame(fp_rows, schema=["condition", "locus", "method", "fp_rate", "n"], orient="row")
        fpdf.write_csv(out / "fp_rates.tsv", separator="\t")
        print("\n=== C6 OLGA false-positive rates ===\n", fpdf)


if __name__ == "__main__":
    main()
