#!/usr/bin/env python3
"""Cross-method epitope-annotation benchmark CLI (see bench/BENCHMARK.md).

Headline mode (`--datasets shortlist|vdjdb`): leave-one-out-by-epitope over the leakage-removed
2026-06-11-ZENODO reference, computing per-locus (TRA/TRB) per-epitope F1 / PR-AUC / retention / purity,
drawn as boxplot + beeswarm across methods. vdjmatch's arm scores each candidate epitope by its count
of leakage-removed same-epitope neighbours (a query NEVER scores against an exact copy of its own CDR3).
External methods plug in via predictions/<method>/<dataset>.tsv (see EXTERNAL_TOOLS.md).

    python bench/compare.py --methods vdjmatch --datasets shortlist --locus TRB --top 20 --out bench/out
"""
from __future__ import annotations

import argparse
import os
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

import polars as pl
from seqtree import Index, SearchParams

import _bench
from metrics import pr_auc_balanced, roc_auc
from vdjmatch import db


# ---- per-epitope metrics over (query -> {epitope: score}, threshold) -----------------------------
def f1_purity_retention(scores, truth, epi, thresh=1.0):
    """For one epitope: F1, purity (precision), retention (recall) of the score>=thresh call set."""
    tp = fp = fn = 0
    for q, sc in scores.items():
        called = sc.get(epi, 0) >= thresh
        is_pos = truth[q] == epi
        tp += called and is_pos
        fp += called and not is_pos
        fn += (not called) and is_pos
    purity = tp / (tp + fp) if tp + fp else float("nan")     # precision
    retention = tp / (tp + fn) if tp + fn else float("nan")  # recall
    f1 = 2 * purity * retention / (purity + retention) if tp else 0.0
    return f1, purity, retention


def pr_auc_epi(scores, truth, epi):
    pairs = [(1 if truth[q] == epi else 0, sc.get(epi, 0.0)) for q, sc in scores.items()]
    return pr_auc_balanced(pairs)


# ---- vdjmatch arm: leakage-removed neighbour-vote --------------------------------------------------
def run_vdjmatch(uc: pl.DataFrame, held: list[str], subs: int, max_q: int):
    """Returns scores[query_id] = {epitope: same-epitope-neighbour count} and truth[query_id]=epitope,
    over the held epitopes' clonotypes, scored against the whole cell with exact-CDR3 self-hits removed."""
    refs = uc.group_by("cdr3").agg(pl.col("epitope").first())
    ref_cdr3, ref_epi = refs["cdr3"].to_list(), refs["epitope"].to_list()
    index = Index.build(ref_cdr3, "aa")
    params = SearchParams(max_subs=subs, max_total_edits=subs, engine="seqtm")
    scores: dict[str, dict[str, float]] = {}
    truth: dict[str, str] = {}
    for epi in held:
        q = uc.filter(pl.col("epitope") == epi).unique("cdr3").head(max_q)
        qs = q["cdr3"].to_list()
        for qi, hl in zip(qs, index.search_batch(qs, params, 0)):
            qid = f"{epi}|{qi}"
            truth[qid] = epi
            votes = Counter()
            for h in hl:
                r = ref_cdr3[h.ref_id]
                if r == qi:
                    continue                                   # leakage: never score an exact self-copy
                votes[ref_epi[h.ref_id]] += 1
            scores[qid] = dict(votes)
    return scores, truth


def load_predictions(path: Path):
    """External method predictions -> scores/truth (truth must come from the dataset, not the file)."""
    scores = defaultdict(dict)
    df = pl.read_csv(path, separator="\t")
    for qid, epi, sc in zip(df["query_id"], df["epitope"], df["score"]):
        scores[qid][epi] = float(sc)
    return scores


# ---- plotting --------------------------------------------------------------------------------------
def plot(results: pl.DataFrame, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    metrics = ["f1", "pr_auc", "retention", "purity"]
    loci = sorted(results["locus"].unique().to_list())
    fig, axes = plt.subplots(len(metrics), len(loci), figsize=(4 * len(loci), 3 * len(metrics)),
                             squeeze=False)
    for i, m in enumerate(metrics):
        for j, loc in enumerate(loci):
            ax = axes[i][j]
            d = results.filter((pl.col("metric") == m) & (pl.col("locus") == loc))
            if d.height:                                               # seaborn takes x/y vectors (no pandas)
                mx, vy = d["method"].to_list(), d["value"].to_list()
                sns.boxplot(x=mx, y=vy, ax=ax, color="#cbd5e1", fliersize=0, width=0.5)
                sns.swarmplot(x=mx, y=vy, ax=ax, color="#2563eb", size=3, alpha=0.7)
            ax.set_title(f"{m} — {loc}"); ax.set_xlabel(""); ax.set_ylabel(m if j == 0 else "")
            ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(out / "compare_boxplots.png", dpi=150)
    print("wrote", out / "compare_boxplots.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--methods", nargs="+", default=["vdjmatch"])
    ap.add_argument("--datasets", nargs="+", default=["shortlist"])
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--locus", nargs="+", default=["TRA", "TRB"])
    ap.add_argument("--subs", type=int, default=1)
    ap.add_argument("--top", type=int, default=20, help="epitopes per locus entering the distribution")
    ap.add_argument("--min-epi", type=int, default=30)
    ap.add_argument("--max-queries", type=int, default=300)
    ap.add_argument("--pred-dir", default="bench/predictions")
    ap.add_argument("--out", default="bench/out")
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    rows = []
    vdj_all = db.load(_bench.source(), species=args.species)
    for dataset in args.datasets:
        for locus in args.locus:
            cell = vdj_all.filter(pl.col("gene") == locus)
            if dataset == "shortlist":
                short = db.replicated(vdj_all, min_refs=2)
                keep = short.filter(pl.col("gene") == locus)["epitope"].unique()
                cell = cell.filter(pl.col("epitope").is_in(keep))
            uc = _bench.long_list(cell, cap=3000, min_n=args.min_epi)
            sizes = (uc.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
                       .filter(pl.col("n") >= args.min_epi).sort(["n", "epitope"], descending=[True, False]))
            held = sizes["epitope"].to_list()[:args.top]
            if not held:
                continue
            for method in args.methods:
                if method == "vdjmatch":
                    scores, truth = run_vdjmatch(uc, held, args.subs, args.max_queries)
                else:
                    p = Path(args.pred_dir) / method / f"{dataset}_{locus}.tsv"
                    if not p.exists():
                        print(f"  skip {method}/{dataset}/{locus}: no predictions at {p}")
                        continue
                    scores = load_predictions(p)
                    truth = {q: e for q, e in ((q, q.split("|", 1)[0]) for q in scores)}
                for epi in held:
                    f1, pur, ret = f1_purity_retention(scores, truth, epi)
                    pra = pr_auc_epi(scores, truth, epi)
                    for metric, v in (("f1", f1), ("pr_auc", pra), ("retention", ret), ("purity", pur)):
                        rows.append((method, dataset, locus, epi, metric, v))
            print(f"{dataset}/{locus}: {len(held)} epitopes scored")

    res = pl.DataFrame(rows, schema=["method", "dataset", "locus", "epitope", "metric", "value"],
                       orient="row")
    res.write_csv(out / "results.tsv", separator="\t")
    print("wrote", out / "results.tsv", f"({res.height} rows)")
    summ = (res.group_by(["method", "locus", "metric"]).agg(pl.col("value").median().alias("median"))
               .sort(["metric", "locus", "method"]))
    print(summ)
    if res.height:
        plot(res, out)


if __name__ == "__main__":
    main()
