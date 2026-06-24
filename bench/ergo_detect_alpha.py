#!/usr/bin/env python3
"""ERGO-II (PMID 33981311) ALPHA-chain detection predictions -> predictions/ergo/{YLQ,GLC}_TRA.tsv.

bench/ergo_detect.py runs ERGO on the CDR3beta (alpha left empty -> UNK). ERGO's pretrained VDJdb model
is beta-conditioned: read_input_file() SKIPS any row whose TRB is invalid and defaults TRA to 'UNK'
when invalid. So there is no "alpha-only" code path - to obtain an alpha-based score we feed the real
alpha CDR3 (+ TRAV/TRAJ) through ERGO's TCR encoder via the TRB/TRBV/TRBJ slots, leaving the true TRA
empty. This is the alpha sequence run through the same model, not a fabricated score. Same peptide +
MHC=HLA-A*02 + T-Cell-Type=CD8 as the beta run. Only YLQ/GLC carry paired alpha (sample6).

Runs in the cmp-ergo conda env; data read from manuscript test_data at runtime via task_table(.,'TRA').
query_id = the ALPHA CDR3 (so paired_product.py can join on it). significant = P >= 0.5.

    python bench/ergo_detect_alpha.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _feat_probe import task_table
from metrics import roc_auc

ENV = "cmp-ergo"
ERGO = Path(__file__).resolve().parent / "external" / "ERGO-II"
WORK = Path("bench/out/_ergo_alpha_work")
THRESH = 0.5
MHC = "HLA-A*02"
EPI = {"YLQ": "YLQPRTFLL", "GLC": "GLCTLVAML"}                      # only these have paired alpha


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    pdir = Path("bench/predictions/ergo"); pdir.mkdir(parents=True, exist_ok=True)
    for task, pep in EPI.items():
        d = task_table(task, "TRA")        # query_id=cdr3=alpha CDR3; v/j=TRAV/TRAJ; label from padj
        inp = WORK / f"{task}_in.csv"
        outp = WORK / f"{task}_out.csv"
        # Feed alpha CDR3 + TRAV/TRAJ through ERGO's TCR encoder via the TRB slots; true TRA left empty.
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
        score_by_cdr3 = dict(zip(out["TRB"].to_list(), out["Score"].cast(pl.Float64).to_list()))
        rows = []
        for cdr3 in d["cdr3"].to_list():                           # cdr3 here is the alpha CDR3
            if cdr3 in score_by_cdr3:
                s = score_by_cdr3[cdr3]
                rows.append((cdr3, pep, float(s), int(s >= THRESH)))
        pl.DataFrame(rows, schema=["query_id", "epitope", "score", "significant"],
                     orient="row").write_csv(pdir / f"{task}_TRA.tsv", separator="\t")
        sc = {c: s for c, _, s, _ in rows}
        pairs = [(int(lab), sc[c]) for c, lab in zip(d["cdr3"], d["label"]) if c in sc]
        auc = roc_auc(pairs) if pairs else float("nan")
        npos = sum(p for p, _ in pairs)
        print(f"ergo {task}/TRA: {len(rows)} scored | ROC-AUC {auc:.3f} (n_pos={npos} n={len(pairs)})")


if __name__ == "__main__":
    main()
