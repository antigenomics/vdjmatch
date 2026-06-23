#!/usr/bin/env python3
"""MixTCRpred on condition C3 (paired TCRvdb): padj<1e-5 (pos) vs >=1e-5 (neg), epitopes GLC + YLQ.

MixTCRpred is a paired-alpha/beta, per-epitope binding classifier (one pretrained model per pMHC).
We score each TCRvdb paired clonotype with the GLC and YLQ models and ask whether the validated
(padj<1e-5) clonotypes score higher than the non-validated ones. Data is read from the manuscript
test_data at runtime (never copied); MixTCRpred runs in the cmp-mixtcrpred conda env on a derived,
gitignored input CSV. Emits predictions/mixtcrpred/C3_paired.tsv + prints ROC/PR/FP/F1.

    python bench/mixtcrpred_c3.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import TESTDATA
from metrics import pr_auc_balanced, roc_auc

ENV = "cmp-mixtcrpred"
MIX = Path(__file__).resolve().parent / "external" / "MixTCRpred"
EPI = {"GLCTLVAML": "A0201_GLCTLVAML", "YLQPRTFLL": "A0201_YLQPRTFLL"}
WORK = Path("bench/out/_mix_work")


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    t = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv").with_columns(pos=pl.col("padj") < 1e-5)
    pred_rows, metric_rows = [], []
    for epi, model in EPI.items():
        d = t.filter(pl.col("epitope_aa") == epi).select(
            cdr3_TRA="cdr3_alpha_aa", cdr3_TRB="cdr3_beta_aa",
            TRAV="TRAV", TRAJ="TRAJ", TRBV="TRBV", TRBJ="TRBJ", pos="pos").drop_nulls(
            ["cdr3_TRA", "cdr3_TRB", "pos"])           # null padj -> no pos/neg label, exclude
        inp, outp = WORK / f"{epi}_in.csv", WORK / f"{epi}_out.csv"
        d.select("cdr3_TRA", "cdr3_TRB", "TRAV", "TRAJ", "TRBV", "TRBJ").write_csv(inp)
        subprocess.run(["conda", "run", "-n", ENV, "python", str(MIX / "MixTCRpred.py"),
                        "--model", model, "--input", str(inp.resolve()),
                        "--output", str(outp.resolve())], check=True, cwd=MIX,
                       capture_output=True, text=True)
        out = pl.read_csv(outp, comment_prefix="#")
        score_col = "score" if "score" in out.columns else out.columns[-1]
        out = out.with_columns(pos=d["pos"])
        labels = out["pos"].to_list()
        scores = out[score_col].cast(pl.Float64).to_list()
        pairs = list(zip([int(x) for x in labels], scores))
        # MixTCRpred call: %rank <= 5 (strong binder) if a rank column exists, else score > 0
        rank_col = next((c for c in out.columns if "rank" in c.lower()), None)
        if rank_col:
            sig = (out[rank_col].cast(pl.Float64) <= 5.0).to_list()
        else:
            sig = [s > 0 for s in scores]
        tp = sum(s for s, p in zip(sig, labels) if p)
        fp = sum(s for s, p in zip(sig, labels) if not p)
        npos, nneg = sum(labels), len(labels) - sum(labels)
        fn, tn = npos - tp, nneg - fp
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        m = {"roc_auc": roc_auc(pairs), "pr_auc": pr_auc_balanced(pairs),
             "fp": fp / (fp + tn) if fp + tn else 0.0,
             "f1": 2 * prec * rec / (prec + rec) if prec + rec else 0.0}
        for q, s, sg in zip(d["cdr3_TRB"], scores, sig):
            pred_rows.append((q, epi, s, int(bool(sg))))
        for k in ("roc_auc", "pr_auc", "fp", "f1"):
            metric_rows.append(("C3", "paired", "mixtcrpred", epi, k, round(m[k], 3)))
        print(f"C3/paired mixtcrpred {epi}: ROC {m['roc_auc']:.3f} PR {m['pr_auc']:.3f} "
              f"FP {m['fp']:.3f} F1 {m['f1']:.3f} (n_pos={npos} n_neg={nneg})")

    pdir = Path("bench/predictions/mixtcrpred"); pdir.mkdir(parents=True, exist_ok=True)
    preds = pl.DataFrame(pred_rows, schema=["query_id", "epitope", "score", "significant"], orient="row")
    preds.write_csv(pdir / "C3_paired.tsv", separator="\t")
    # also emit per-task competitor files (query_id = CDR3beta) so the manuscript scores MixTCRpred on
    # the SAME query set as the other methods (read_predictions defaults unscored queries to 0).
    task_of = {"YLQPRTFLL": "YLQ", "GLCTLVAML": "GLC"}
    for epi, task in task_of.items():
        preds.filter(pl.col("epitope") == epi).write_csv(pdir / f"{task}_TRB.tsv", separator="\t")
    pl.DataFrame(metric_rows, schema=["condition", "locus", "method", "epitope", "metric", "value"],
                 orient="row").write_csv(Path("bench/out/benchmark") / "mixtcrpred_C3.tsv", separator="\t")


if __name__ == "__main__":
    main()
