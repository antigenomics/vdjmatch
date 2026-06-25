"""Which cases does imw-DETECT rescue on GLC-TRB, and is it memorisation or signal?

imw-DETECT is trained on public TCRs (incl. VDJdb GLC binders); vdjmatch EXCLUDES exact matches (no
leakage). Hypothesis: imw's GLC-TRB edge is the exact-match positives it memorised. Test: split positives
into exact (== a VDJdb GLC reference CDR3) vs fuzzy; compare NED vs imw ROC overall and on fuzzy-only.

    .venv/bin/python bench/hardcase_imw.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _feat_probe as FP                                             # noqa: E402
from metrics import roc_auc, pr_auc_balanced                        # noqa: E402

PRED = Path(__file__).resolve().parent / "predictions"


def analyze(task, locus, methods=("imw-detect", "ergo", "nettcr")):
    d = FP.task_table(task, locus)
    ned = FP.baseline_scores(task, locus)
    ref = set(FP.ref_table(task, locus)["cdr3"].to_list())
    rows = []
    for c, lab in zip(d["cdr3"].to_list(), d["label"].to_list()):
        rows.append((c, int(lab), float(ned.get(c, 0.0)), c in ref))
    df = pl.DataFrame(rows, schema=["cdr3", "label", "ned", "exact"], orient="row")
    npos = int(df["label"].sum())
    nexact = df.filter((pl.col("label") == 1) & pl.col("exact")).height
    print(f"\n=== {task}-{locus} ===  n={df.height} pos={npos} (exact-match pos={nexact}, "
          f"{100*nexact/max(npos,1):.0f}%) neg={df.height-npos}")
    y = df["label"].to_numpy()
    fz = df.filter(~((pl.col("label") == 1) & pl.col("exact")))      # drop memorisable exact-match positives
    yf = fz["label"].to_numpy()
    print(f"  {'method':12s} {'ROC(all)':>9} {'ROC(fuzzy)':>11}   score: exact-pos / fuzzy-pos / neg")
    print(f"  {'NED':12s} {roc_auc(list(zip(y, df['ned'].to_list()))):>9.3f} "
          f"{roc_auc(list(zip(yf, fz['ned'].to_list()))):>11.3f}")
    for m in methods:
        f = PRED / m / f"{task}_{locus}.tsv"
        if not f.exists():
            continue
        sc = {r["query_id"]: float(r["score"]) for r in pl.read_csv(f, separator="\t").iter_rows(named=True)}
        dd = df.with_columns(s=pl.col("cdr3").replace_strict(sc, default=None)).drop_nulls("s")
        ff = dd.filter(~((pl.col("label") == 1) & pl.col("exact")))
        ra = roc_auc(list(zip(dd["label"].to_numpy(), dd["s"].to_list())))
        rf = roc_auc(list(zip(ff["label"].to_numpy(), ff["s"].to_list())))
        med = lambda s: f"{s.median():.2f}" if s.len() else "--"
        rff = f"{rf:.3f}" if not (rf != rf) else "  (no fuzzy pos)"
        ep = dd.filter((pl.col("label") == 1) & pl.col("exact"))["s"]
        fp = dd.filter((pl.col("label") == 1) & ~pl.col("exact"))["s"]
        ng = dd.filter(pl.col("label") == 0)["s"]
        drop = f"{ra-rf:+.3f}" if not (rf != rf) else "n/a"
        print(f"  {m:12s} {ra:>9.3f} {rff:>11}   {med(ep)} / {med(fp)} / {med(ng)}"
              f"   (drop {drop} when exact-match positives removed)")


if __name__ == "__main__":
    for t, l in (("GLC", "TRB"), ("GLC", "TRA"), ("NLV", "TRB"), ("LLW", "TRB")):
        analyze(t, l)
