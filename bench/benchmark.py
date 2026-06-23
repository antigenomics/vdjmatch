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


def vgene(v) -> str:
    """Normalise a V gene to its allele-stripped gene name for V-matching (e.g. TRBV20-1*01 -> TRBV20-1)."""
    return (v or "").split("*")[0]


def shortlist(df: pl.DataFrame, min_refs: int = 2) -> pl.DataFrame:
    """Clonotype-pMHC pairs in >= min_refs distinct references (keeps ALL columns incl complex_id)."""
    key = ["gene", "cdr3", "v", "j", "epitope"]
    keep = (df.group_by(key).agg(pl.col("reference_id").n_unique().alias("nr"))
            .filter(pl.col("nr") >= min_refs).select(key))
    return df.join(keep, on=key, how="semi")


def ref_index(df: pl.DataFrame, locus: str):
    """VDJdb frame -> (Index, ref_epi, ref_v, n_epi, n_epi_v, n_v) for a single locus. One entry per
    unique CDR3 (representative V/epitope); the V vector enables the V+CDR3 joint E-value."""
    r = (_bench.valid_cdr3(df.filter(pl.col("gene") == locus)).group_by("cdr3")
         .agg(pl.col("epitope").first(), pl.col("v").first()))
    epi = r["epitope"].to_list()
    v = [vgene(x) for x in r["v"].to_list()]
    return (Index.build(r["cdr3"].to_list(), "aa"), epi, v,
            Counter(epi), Counter(zip(v, epi)), Counter(v))


# ---- vdjmatch classification (first-hit E-value; optional V+CDR3 joint null) ----------------------
def vdjmatch_classify(tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v, ctrl, queries, query_v, epitopes,
                      alpha, exclude_exact, v_mode="none", params=None):
    """queries -> per query {epitope: (-log10 p_enrichment, significant)} for the candidate epitopes,
    plus an overall-significant flag (any-epitope first-hit p < alpha) for FP estimation. ``v_mode``:
    ``none`` = CDR3-only; ``match_v`` = V+CDR3 joint E-value (same-V-restricted Poisson tail)."""
    M, Ntot = len(ctrl), len(ref_epi)
    match_v = v_mode == "match_v"
    t2 = lambda hits: [(c, e) for c, e, *_ in hits]            # drop the V tag for the V-agnostic pvalue
    th, cc = first_hit.scan(tgt, ref_epi, ctrl, queries, target_v=ref_v, params=params,
                            exclude_exact=exclude_exact)
    scores, overall = {}, {}
    for q, qv, t, c in zip(queries, query_v, th, cc):
        # FP/significance always V-AGNOSTIC (the same-V null is tight for rare V -> over-calls); the
        # V+CDR3 prior sharpens the per-epitope SCORE (ranking), not the significance threshold.
        overall[q] = first_hit.pvalue(t2(t), c, Ntot, M)["p_enrichment"] < alpha
        scores[q] = {}
        for e in epitopes:
            N_e = n_epi_v.get((qv, e), 1) if match_v else n_epi.get(e, 1)
            p_score = first_hit.pvalue_v(t, c, qv, N_e, M, epitope=e, match_v=match_v)["p_enrichment"]
            p_sig = first_hit.pvalue(t2(t), c, n_epi.get(e, 1), M, epitope=e)["p_enrichment"]
            scores[q][e] = (-math.log10(max(p_score, 1e-300)), p_sig < alpha)
    return scores, overall


def read_predictions(path: Path, allq, epitopes):
    """predictions/<method>/<cond>_<locus>.tsv (query_id, epitope, score[, significant]) -> scores dict
    {query: {epitope: (score, significant)}} with 0/False defaults for unscored (query, epitope)."""
    scores = {q: {e: (0.0, False) for e in epitopes} for q in allq}
    if not path.exists():
        return None
    df = pl.read_csv(path, separator="\t")
    has_sig = "significant" in df.columns
    for r in df.iter_rows(named=True):
        q, e = r["query_id"], r["epitope"]
        if q in scores and e in epitopes:
            scores[q][e] = (float(r["score"]), bool(r["significant"]) if has_sig else r["score"] > 0)
    return scores


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
def cond_sample_pair(name, loci):
    """C1/C2: discriminate two query groups by matching against vdjdb2026."""
    for locus in loci:
        if name == "C1":                                          # sample2: LLW vs LLL by antigen.epitope
            d = (pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t")
                 .rename({"antigen.epitope": "label"}).select("cdr3", "label", v="v.segm"))
            d = _bench.valid_cdr3(d).unique("cdr3")
            groups = {EPI["LLW"]: d.filter(pl.col("label") == EPI["LLW"])["cdr3"].to_list(),
                      EPI["LLL"]: d.filter(pl.col("label") == EPI["LLL"])["cdr3"].to_list()}
            tasks = [(EPI["LLW"], groups[EPI["LLW"]], groups[EPI["LLL"]]),
                     (EPI["LLL"], groups[EPI["LLL"]], groups[EPI["LLW"]])]
        else:                                                     # C2 sample1: cmv(NLV+) vs control(NLV-)
            d = (pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
                 .filter(pl.col("gene") == locus).select("cdr3", label="type", v="v.segm"))
            d = _bench.valid_cdr3(d).unique("cdr3")
            pos = d.filter(pl.col("label") == "cmv")["cdr3"].to_list()
            neg = d.filter(pl.col("label") == "control")["cdr3"].to_list()
            tasks = [(EPI["NLV"], pos, neg)]
        qv = {c: vgene(v) for c, v in zip(d["cdr3"], d["v"])}
        yield locus, release("vdjdb2026"), tasks, True, qv


def cond_olga_fp(loci, n=1000):
    """C6: OLGA negatives (sample5=TRA, sample4=TRB) vs vdjdb2025 -> any significant hit is a FP."""
    files = {"TRA": "sample5", "TRB": "sample4"}
    for locus in loci:
        d = (pl.read_csv(TESTDATA / f"{files[locus]}_olga_airr.txt", separator="\t")
             .select(cdr3="junction_aa", v="v_gene").pipe(_bench.valid_cdr3).unique("cdr3"))
        if d.height > n:
            d = d.sample(n, seed=0)
        qv = {c: vgene(v) for c, v in zip(d["cdr3"], d["v"])}
        yield locus, release("vdjdb2025"), d["cdr3"].to_list(), qv


def cond_tcrvdb(loci):
    """C3: TCRvdb padj<1e-5 (pos) vs >=1e-5 (neg) per epitope, annotated vs vdjdb2026."""
    t = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv").with_columns(pos=pl.col("padj") < 1e-5)
    chain_col = {"TRA": "cdr3_alpha_aa", "TRB": "cdr3_beta_aa"}
    v_col = {"TRA": "TRAV", "TRB": "TRBV"}
    for locus in loci:
        d = (t.select(cdr3=chain_col[locus], epitope="epitope_aa", pos="pos", v=v_col[locus])
             .pipe(_bench.valid_cdr3).unique("cdr3"))
        qv = {c: vgene(v) for c, v in zip(d["cdr3"], d["v"])}
        tasks = []
        for e in (EPI["GLC"], EPI["YLQ"]):
            de = d.filter(pl.col("epitope") == e)
            pos = de.filter("pos")["cdr3"].to_list()
            neg = de.filter(~pl.col("pos"))["cdr3"].to_list()
            if pos and neg:
                tasks.append((e, pos, neg))
        yield locus, release("vdjdb2026"), tasks, True, qv


def cond_loo(which, loci, top, min_epi, max_q):
    """C4 (within vdjdb2025) / C5 (new-in-2026 vs 2025 reference): per-epitope X(pos) vs other(neg).

    Both releases are restricted to the SHORTLIST: clonotype-pMHC pairs seen in >=2 references
    (``db.replicated``), per request — the public/learnable subset.
    """
    v25 = shortlist(release("vdjdb2025"))
    for locus in loci:
        ref_cell = _bench.valid_cdr3(v25.filter(pl.col("gene") == locus))
        if which == "C4":
            test_cell, excl = _bench.long_list(v25.filter(pl.col("gene") == locus), cap=max_q,
                                               min_n=min_epi), True
        else:                                                     # C5: 2026 shortlist clonotypes new vs 2025
            v26 = _bench.long_list(shortlist(release("vdjdb2026"))
                                   .filter(pl.col("gene") == locus), cap=max_q, min_n=min_epi)
            test_cell = v26.join(ref_cell.select("cdr3", "epitope").unique(),
                                 on=["cdr3", "epitope"], how="anti")
            excl = False
        held = _top_held(test_cell, _big_epis(ref_cell, min_epi), top, min_epi)
        per = {e: test_cell.filter(pl.col("epitope") == e).unique("cdr3")["cdr3"].to_list()[:max_q]
               for e in held}
        tasks = [(e, per[e], [q for e2 in held if e2 != e for q in per[e2]]) for e in held]
        qv = {c: vgene(v) for c, v in zip(test_cell["cdr3"], test_cell["v"])}
        yield locus, v25, tasks, excl, qv


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
    ap.add_argument("--conditions", nargs="+", default=["C1", "C2", "C3", "C6"])
    ap.add_argument("--methods", nargs="+", default=["vdjmatch"])
    ap.add_argument("--loci", nargs="+", default=["TRA", "TRB"])
    ap.add_argument("--v-mode", default="none", choices=["none", "match_v"],
                    help="vdjmatch V+CDR3 joint E-value (same-V-restricted null)")
    ap.add_argument("--scope", default="5,2,2", help="first-hit scope max_edits,max_ins,max_dels")
    ap.add_argument("--alpha", type=float, default=1e-3)
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--min-epi", type=int, default=30)
    ap.add_argument("--max-queries", type=int, default=300)
    ap.add_argument("--olga-n", type=int, default=2000, help="C6: OLGA negatives per locus")
    ap.add_argument("--pred-dir", default="bench/predictions")
    ap.add_argument("--out", default="bench/out/benchmark")
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rows, fp_rows, n_rows = [], [], []
    params = first_hit.scope(*[int(x) for x in args.scope.split(",")])

    def classify_cells(cond):
        if cond in ("C1", "C2"):
            return cond_sample_pair(cond, ["TRB"] if cond == "C1" else args.loci)
        if cond == "C3":
            return cond_tcrvdb(args.loci)
        if cond in ("C4", "C5"):
            return cond_loo(cond, args.loci, args.top, args.min_epi, args.max_queries)
        return None

    for cond in args.conditions:
        if cond == "C6":
            for locus, ref_df, queries, qv in cond_olga_fp(args.loci, args.olga_n):
                if "vdjmatch" in args.methods:
                    tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v = ref_index(ref_df, locus)
                    _, overall = vdjmatch_classify(tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v,
                                                   background(locus), queries,
                                                   [qv[q] for q in queries], [], args.alpha, True,
                                                   v_mode=args.v_mode, params=params)
                    fp = sum(overall.values()) / len(queries)
                    fp_rows.append(("C6", locus, "vdjmatch", fp, len(queries)))
                    print(f"C6/{locus}: FP={fp*100:.2f}% over {len(queries)} OLGA negatives")
            continue
        cells = classify_cells(cond)
        if cells is None:
            print(f"  ({cond} not yet implemented)")
            continue
        for locus, ref_df, tasks, excl, qv in cells:
            if not tasks:
                print(f"{cond}/{locus}: no tasks (skipped)")
                continue
            epitopes = sorted({e for e, _, _ in tasks})
            allq = sorted({q for _, p, n in tasks for q in (*p, *n)})
            n_rows.append((cond, locus, len(allq), len(tasks)))
            tgt = None                                            # lazy: only build the index for vdjmatch
            for method in args.methods:
                if method == "vdjmatch":
                    if tgt is None:
                        tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v = ref_index(ref_df, locus)
                    scores, _ = vdjmatch_classify(tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v,
                                                  background(locus), allq, [qv[q] for q in allq],
                                                  epitopes, args.alpha, excl, v_mode=args.v_mode, params=params)
                else:
                    scores = read_predictions(Path(args.pred_dir) / method / f"{cond}_{locus}.tsv",
                                              allq, epitopes)
                    if scores is None:
                        print(f"  skip {method}/{cond}/{locus}: no predictions")
                        continue
                for epi, pos, neg in tasks:
                    m = classify_metrics(scores, pos, neg, epi)
                    for k in ("roc_auc", "pr_auc", "fp", "f1"):
                        rows.append((cond, locus, method, epi, k, m[k]))
                    if method == "vdjmatch":
                        rows.append((cond, locus, method, epi, "n_pos", float(m["n_pos"])))
                        rows.append((cond, locus, method, epi, "n_neg", float(m["n_neg"])))
                print(f"{cond}/{locus}/{method}: {len(tasks)} epitope-task(s), {len(allq)} queries")

    if n_rows:
        pl.DataFrame(n_rows, schema=["condition", "locus", "n_queries", "n_epitopes"],
                     orient="row").write_csv(out / "dataset_n.tsv", separator="\t")
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
