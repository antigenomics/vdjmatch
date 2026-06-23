#!/usr/bin/env python3
"""ERGO-II (PMID 33981311) predictions -> detection-benchmark prediction files.

ERGO-II is a supervised TCR-peptide binding classifier: input (CDR3β [+optional α/V/J/MHC], peptide)
-> P(bind). We use their **pretrained VDJdb model** (leakage on their side is accepted per the plan).
For each detection task we score every query CDR3β against the task's target peptide and emit
predictions/ergo/<TASK>_TRB.tsv (query_id = CDR3β, score = P(bind), significant = 1 iff P>=0.5).

ERGO runs in the cmp-ergo conda env (py3.7, torch 1.4 / pytorch-lightning 0.7.1) on a derived,
gitignored input CSV; data is read from the manuscript test_data at runtime via _feat_probe.task_table.

    python bench/ergo_detect.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from _feat_probe import task_table
from metrics import roc_auc

ENV = "cmp-ergo"
ERGO = Path(__file__).resolve().parent / "external" / "ERGO-II"
WORK = Path("bench/out/_ergo_work")
THRESH = 0.5
MHC = "HLA-A*02"                                                   # all 5 targets are HLA-A*02 restricted
EPI = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV",
       "YLQ": "YLQPRTFLL", "GLC": "GLCTLVAML"}


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    pdir = Path("bench/predictions/ergo"); pdir.mkdir(parents=True, exist_ok=True)
    for task, pep in EPI.items():
        d = task_table(task)                                       # cdr3,v,j,label[,a_cdr3,a_v,a_j],query_id
        # ERGO input CSV: TRA,TRB,TRAV,TRAJ,TRBV,TRBJ,T-Cell-Type,Peptide,MHC. We supply TRB + V/J + MHC;
        # alpha is left empty (-> UNK inside ERGO). Unknown V/J/MHC degrade to UNK (index 0) gracefully.
        inp = WORK / f"{task}_in.csv"
        outp = WORK / f"{task}_out.csv"
        ergo_in = d.select(
            TRA=pl.lit(None, dtype=pl.Utf8), TRB="cdr3",
            TRAV=pl.lit(None, dtype=pl.Utf8), TRAJ=pl.lit(None, dtype=pl.Utf8),
            TRBV="v", TRBJ="j").with_columns(
            **{"T-Cell-Type": pl.lit("CD8")}, Peptide=pl.lit(pep), MHC=pl.lit(MHC))
        ergo_in.write_csv(inp)
        subprocess.run(["conda", "run", "-n", ENV, "python", "_run_predict.py",
                        "vdjdb", str(inp.resolve()), str(outp.resolve())],
                       check=True, cwd=ERGO, capture_output=True, text=True)
        out = pl.read_csv(outp)
        # ERGO drops rows with invalid TRB/peptide; map score back by TRB string.
        score_by_cdr3 = dict(zip(out["TRB"].to_list(), out["Score"].cast(pl.Float64).to_list()))
        rows = []
        for cdr3 in d["cdr3"].to_list():
            if cdr3 in score_by_cdr3:
                s = score_by_cdr3[cdr3]
                rows.append((cdr3, pep, float(s), int(s >= THRESH)))
        pl.DataFrame(rows, schema=["query_id", "epitope", "score", "significant"],
                     orient="row").write_csv(pdir / f"{task}_TRB.tsv", separator="\t")
        # ROC sanity (positives should outscore negatives)
        sc = {c: s for c, _, s, _ in rows}
        pairs = [(int(lab), sc[c]) for c, lab in zip(d["cdr3"], d["label"]) if c in sc]
        auc = roc_auc(pairs) if pairs else float("nan")
        npos = sum(p for p, _ in pairs)
        print(f"ergo {task}/TRB: {len(rows)} scored | ROC-AUC {auc:.3f} (n_pos={npos} n={len(pairs)})")


if __name__ == "__main__":
    main()
