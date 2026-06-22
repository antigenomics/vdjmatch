#!/usr/bin/env python3
"""Cross-method epitope-annotation benchmark CLI (see bench/BENCHMARK.md).

Headline mode (`--datasets shortlist|vdjdb`): leave-one-out-by-epitope over the leakage-removed
2026-06-11-ZENODO reference, computing per-locus (TRA/TRB) per-epitope F1 / PR-AUC / retention / purity,
drawn as boxplot + beeswarm across methods. vdjmatch's arm scores each candidate epitope by its count
of leakage-removed same-epitope neighbours (a query NEVER scores against an exact copy of its own CDR3).
External methods plug in via predictions/<method>/<dataset>_<locus>.tsv (see EXTERNAL_TOOLS.md).

Other datasets (`--datasets ...`):
  temporal : temporal holdout — 2025 release (`--vdjdb2025-pin`) as reference, test = 2026 clonotypes
             NEW vs 2025 (held out by time); exact matches KEPT (cross-release match is a legitimate
             annotation, not leakage). The real-world 'annotate tomorrow's data' benchmark.
  tcrvdb   : validated TCRvdb paired clonotypes (`--tcrvdb-padj`) as test, annotated per chain against
             the 2026 reference (exact-match removed — TCRvdb overlaps vdjdb).
  samples  : NLV/LLW/LLL ROC+PR vs OLGA negatives (see run_samples).

    python bench/compare.py --methods vdjmatch --datasets shortlist --locus TRB --top 20 --out bench/out
"""
from __future__ import annotations

import argparse
import math
import os
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

import polars as pl
from seqtree import Index, SearchParams

import _bench
from metrics import pr_auc_balanced, roc_auc
from vdjmatch import db
from vdjmatch.evalue import background, first_hit

TESTDATA = Path(os.environ.get("VDJMATCH_TESTDATA",
                               "/Users/mikesh/vcs/manuscripts/2026-vdjmatch/test_data"))
SAMPLE_EPI = {"NLVPMVATV": "NLV (CMV)", "LLWNGPMAV": "LLW (YF)", "LLLGIGILV": "LLL (BST-2)"}


def load_sample(name: str) -> pl.DataFrame:
    """External-validation query set -> df[cdr3, true_epitope] (TRB; OLGA negatives have null epitope)."""
    if name == "sample1":                                          # CMV NLV (labeled-positive)
        d = (pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
               .filter(pl.col("gene") == "TRB").select("cdr3")
               .with_columns(true_epitope=pl.lit("NLVPMVATV")))
    elif name == "sample2":                                        # YF LLW + BST-2 LLL (per-row labeled)
        d = (pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t")
               .rename({"antigen.epitope": "true_epitope"}).select("cdr3", "true_epitope"))
    elif name == "sample5":                                        # random OLGA: negatives
        d = (pl.read_csv(TESTDATA / "sample5_olga_airr.txt", separator="\t")
               .select(cdr3="junction_aa").with_columns(true_epitope=pl.lit(None, pl.Utf8)))
    else:
        raise ValueError(name)
    return _bench.valid_cdr3(d).unique("cdr3")


def _roc_pr(pairs):
    """(label,score) -> ((fpr,tpr,auc),(recall,prec,ap)) point arrays for plotting."""
    s = sorted(pairs, key=lambda x: -x[1])
    P = sum(l for l, _ in s) or 1
    N = len(s) - sum(l for l, _ in s) or 1
    tp = fp = 0
    fpr, tpr, rec, prec = [0.0], [0.0], [], []
    for lab, _ in s:
        tp += lab; fp += 1 - lab
        fpr.append(fp / N); tpr.append(tp / P)
        rec.append(tp / P); prec.append(tp / (tp + fp))
    return (fpr, tpr, roc_auc(pairs)), (rec, prec, pr_auc_balanced(pairs))


def run_samples(args, out: Path):
    """NLV / LLW / LLL ROC+PR (positives vs OLGA negatives) and OLGA spurious-hit filtering, per method."""
    v = not args.quiet
    log = (lambda *a: print(*a)) if v else (lambda *a: None)
    log("loading VDJdb-beta reference + control...")
    vdj = db.load(_bench.source(), species=args.species).filter(pl.col("gene") == "TRB")
    ref = _bench.valid_cdr3(vdj).group_by("cdr3").agg(pl.col("epitope").first())
    ref_cdr3, ref_epi = ref["cdr3"].to_list(), ref["epitope"].to_list()
    tgt = Index.build(ref_cdr3, "aa")
    ctrl = background("TRB")
    N, M, N_epi = len(tgt), len(ctrl), Counter(ref_epi)
    log(f"  target N={N} unique CDR3; control M={M}")
    q1, q2, q5 = load_sample("sample1"), load_sample("sample2"), load_sample("sample5")
    if args.olga_n and q5.height > args.olga_n:
        q5 = q5.sample(args.olga_n, seed=0)
    queries = pl.concat([q1, q2, q5]).unique("cdr3")
    qlist = queries["cdr3"].to_list()
    truth = dict(zip(queries["cdr3"], queries["true_epitope"]))
    epis = list(SAMPLE_EPI)
    log(f"  queries: {q1.height} NLV + {q2.height} LLW/LLL + {q5.height} OLGA = {len(qlist)} unique")

    # vdjmatch: one wide first-hit scan; score(query, epitope) = -log10 p_enrichment at the first E-hit;
    # OLGA "significant" = overall first-hit p_enrichment < alpha (any epitope).
    log(f"first-hit scan (scope 5 edits, <=2 ins, <=2 del) over {len(qlist)} queries...")
    th, cc = first_hit.scan(tgt, ref_epi, ctrl, qlist, exclude_exact=True, progress=v)
    vdj_score, olga_sig = {}, {"vdjmatch": {}}
    for q, t, c in zip(qlist, th, cc):
        olga_sig["vdjmatch"][q] = first_hit.pvalue(t, c, N, M)["p_enrichment"] < args.alpha
        vdj_score[q] = {e: -math.log10(max(first_hit.pvalue(t, c, N_epi.get(e, 1), M, epitope=e)
                                           ["p_enrichment"], 1e-300)) for e in epis}
    methods = {"vdjmatch": vdj_score}
    # external methods: predictions/<m>/samples.tsv (query_id, epitope, score [, significant])
    for m in args.methods:
        if m == "vdjmatch":
            continue
        p = Path(args.pred_dir) / m / "samples.tsv"
        if p.exists():
            sc, sig = defaultdict(dict), {}
            for row in pl.read_csv(p, separator="\t").iter_rows(named=True):
                sc[row["query_id"]][row["epitope"]] = float(row["score"])
                sig[row["query_id"]] = sig.get(row["query_id"], False) or bool(row.get("significant"))
            methods[m] = dict(sc)
            olga_sig[m] = {q: sig.get(q, max(sc.get(q, {}).values(), default=0) > 0) for q in qlist}

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    epis = list(SAMPLE_EPI)
    fig, axes = plt.subplots(2, len(epis), figsize=(4 * len(epis), 7), squeeze=False)
    rows = []
    for j, epi in enumerate(epis):
        for m, sc in methods.items():
            pairs = [(1 if truth[c] == epi else 0, sc.get(c, {}).get(epi, 0.0))
                     for c in qlist if truth[c] == epi or truth[c] is None]   # positives vs OLGA
            (fpr, tpr, auc), (rec, prec, ap) = _roc_pr(pairs)
            axes[0][j].plot(fpr, tpr, label=f"{m} ({auc:.2f})")
            axes[1][j].plot(rec, prec, label=f"{m} ({ap:.2f})")
            rows.append((m, epi, "roc_auc", auc)); rows.append((m, epi, "pr_auc", ap))
        axes[0][j].plot([0, 1], [0, 1], "k--", lw=0.5); axes[0][j].set_title(SAMPLE_EPI[epi])
        axes[0][j].set_xlabel("FPR"); axes[0][j].set_ylabel("TPR" if j == 0 else ""); axes[0][j].legend(fontsize=7)
        axes[1][j].set_xlabel("recall"); axes[1][j].set_ylabel("precision" if j == 0 else ""); axes[1][j].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out / "samples_roc_pr.png", dpi=150)
    print("wrote", out / "samples_roc_pr.png")

    print(f"\nOLGA spurious-hit filtering (fraction called significant at p<{args.alpha}; lower=better):")
    olga = [c for c in qlist if truth[c] is None]
    for m in methods:
        frac = sum(1 for c in olga if olga_sig[m].get(c, False)) / len(olga)
        print(f"  {m:12} {frac*100:6.3f}%  (n={len(olga)})")
        rows.append((m, "OLGA", "spurious_rate", frac))
    pl.DataFrame(rows, schema=["method", "epitope", "metric", "value"], orient="row").write_csv(
        out / "samples_metrics.tsv", separator="\t")


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


# ---- vdjmatch arm: neighbour-vote (ref index, test queries) -----------------------------------------
def run_vdjmatch(ref_df: pl.DataFrame, test_df: pl.DataFrame, held: list[str], subs: int, max_q: int,
                 exclude_exact: bool = True):
    """Score each held-epitope test clonotype against the ``ref_df`` index by same-epitope neighbour
    vote. Returns scores[query_id]={epitope: vote count} and truth[query_id]=epitope. ``exclude_exact``
    drops exact-CDR3 self-hits (within-release leakage control); set False for a temporal holdout where
    the reference predates the query, so an exact match is a legitimate cross-time annotation."""
    refs = ref_df.group_by("cdr3").agg(pl.col("epitope").first())
    ref_cdr3, ref_epi = refs["cdr3"].to_list(), refs["epitope"].to_list()
    index = Index.build(ref_cdr3, "aa")
    params = SearchParams(max_subs=subs, max_total_edits=subs, engine="seqtm")
    scores: dict[str, dict[str, float]] = {}
    truth: dict[str, str] = {}
    for epi in held:
        q = test_df.filter(pl.col("epitope") == epi).unique("cdr3").head(max_q)
        qs = q["cdr3"].to_list()
        for qi, hl in zip(qs, index.search_batch(qs, params, 0)):
            qid = f"{epi}|{qi}"
            truth[qid] = epi
            votes = Counter()
            for h in hl:
                r = ref_cdr3[h.ref_id]
                if exclude_exact and r == qi:
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


# ---- dataset cell builders (ref / test / held epitopes per locus) ----------------------------------
_REL_CACHE: dict[str, pl.DataFrame] = {}


def _release(key: str, species: str) -> pl.DataFrame:
    """Cached VDJdb release loader. ``key`` = '2026' (HF-pinned benchmark) or a GitHub release tag."""
    if key not in _REL_CACHE:
        _REL_CACHE[key] = (db.load(_bench.source(), species=species) if key == "2026"
                           else db.load(pin=key, asset="default", species=species))
    return _REL_CACHE[key]


def _top_held(df: pl.DataFrame, ref_epis, top: int, min_n: int) -> list[str]:
    """Held epitopes: present in the reference, >= min_n unique test CDR3, largest first."""
    sizes = (df.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
               .filter((pl.col("n") >= min_n) & pl.col("epitope").is_in(list(ref_epis)))
               .sort(["n", "epitope"], descending=[True, False]))
    return sizes["epitope"].to_list()[:top]


def load_tcrvdb(padj: float) -> pl.DataFrame:
    """TCRvdb validated paired clonotypes -> long df[gene, cdr3, epitope] (one row per chain)."""
    src = TESTDATA / "sample6_TCRvdb.csv"
    t = pl.read_csv(src).filter(pl.col("padj") < padj)
    a = t.select(cdr3="cdr3_alpha_aa", epitope="epitope_aa").with_columns(gene=pl.lit("TRA"))
    b = t.select(cdr3="cdr3_beta_aa", epitope="epitope_aa").with_columns(gene=pl.lit("TRB"))
    return _bench.valid_cdr3(pl.concat([a, b])).unique(["gene", "cdr3", "epitope"])


def dataset_cells(dataset: str, args):
    """Yield (locus, ref_df, test_df, held, exclude_exact) for a dataset.

    shortlist/vdjdb: within-release LOO on the 2026 benchmark (exact-self removed).
    temporal: 2025 release reference, test = 2026 clonotypes new vs 2025 (held out by time); exact
              matches KEPT (cross-release match is a legitimate annotation, not leakage).
    tcrvdb: validated TCRvdb (padj<thresh) test annotated against the 2026 reference, per chain.
    """
    if dataset in ("shortlist", "vdjdb"):
        vdj = _release("2026", args.species)
        for locus in args.locus:
            cell = vdj.filter(pl.col("gene") == locus)
            if dataset == "shortlist":
                keep = db.replicated(vdj, min_refs=2).filter(pl.col("gene") == locus)["epitope"].unique()
                cell = cell.filter(pl.col("epitope").is_in(keep))
            uc = _bench.long_list(cell, cap=3000, min_n=args.min_epi)
            held = _top_held(uc, set(uc["epitope"].unique()), args.top, args.min_epi)
            yield locus, uc, uc, held, True
    elif dataset == "temporal":
        v25, v26 = _release(args.vdjdb2025_pin, args.species), _release("2026", args.species)
        for locus in args.locus:
            ref = _bench.valid_cdr3(v25.filter(pl.col("gene") == locus))
            test26 = _bench.long_list(v26.filter(pl.col("gene") == locus), cap=3000, min_n=args.min_epi)
            new = test26.join(ref.select("cdr3", "epitope").unique(),       # held out by time
                              on=["cdr3", "epitope"], how="anti")
            held = _top_held(new, set(ref["epitope"].unique()), args.top, args.min_epi)
            yield locus, ref, new, held, False                 # keep exact matches (temporal)
    elif dataset == "tcrvdb":
        vdj, tcr = _release("2026", args.species), load_tcrvdb(args.tcrvdb_padj)
        for locus in args.locus:
            ref = _bench.valid_cdr3(vdj.filter(pl.col("gene") == locus))
            test = tcr.filter(pl.col("gene") == locus)
            held = _top_held(test, set(ref["epitope"].unique()), args.top, args.min_epi)
            yield locus, ref, test, held, True                 # TCRvdb overlaps vdjdb -> leakage control
    else:
        raise ValueError(f"unknown dataset {dataset!r}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--methods", nargs="+", default=["vdjmatch"])
    ap.add_argument("--datasets", nargs="+", default=["shortlist"])
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--locus", nargs="+", default=["TRA", "TRB"])
    ap.add_argument("--subs", type=int, default=1)
    ap.add_argument("--alpha", type=float, default=1e-3, help="first-hit E-value significance cutoff")
    ap.add_argument("--olga-n", type=int, default=0, help="subsample OLGA negatives (0 = all ~240k)")
    ap.add_argument("--quiet", action="store_true", help="suppress progress bars / stage logging")
    ap.add_argument("--top", type=int, default=20, help="epitopes per locus entering the distribution")
    ap.add_argument("--min-epi", type=int, default=30)
    ap.add_argument("--max-queries", type=int, default=300)
    ap.add_argument("--vdjdb2025-pin", default="2025-12-29", help="VDJdb release tag for the temporal-holdout reference")
    ap.add_argument("--tcrvdb-padj", type=float, default=1e-5, help="TCRvdb validated-clonotype padj cutoff")
    ap.add_argument("--pred-dir", default="bench/predictions")
    ap.add_argument("--out", default="bench/out")
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    if "samples" in args.datasets:
        run_samples(args, out)
        args.datasets = [d for d in args.datasets if d != "samples"]
        if not args.datasets:
            return

    rows = []
    for dataset in args.datasets:
        for locus, ref_df, test_df, held, exclude_exact in dataset_cells(dataset, args):
            if not held:
                print(f"  {dataset}/{locus}: no held epitopes")
                continue
            for method in args.methods:
                if method == "vdjmatch":
                    scores, truth = run_vdjmatch(ref_df, test_df, held, args.subs, args.max_queries,
                                                 exclude_exact)
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
