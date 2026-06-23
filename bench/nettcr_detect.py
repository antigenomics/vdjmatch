#!/usr/bin/env python3
"""NetTCR-2.0 (Montemurro 2021, PMID 34508155) predictions -> detection-benchmark prediction files.

NetTCR-2.0 is a BLOSUM-CNN TCR-peptide binding classifier (CDR3β, or paired CDR3α+β, + peptide ->
binding score). The repo ships no saved weights, so we reproduce "their pretrained model" by training
once on their bundled published data (train_beta_90 for β, train_ab_90_alphabeta for paired) with
their architecture, then predicting our queries. NetTCR encodes C/F-trimmed CDR3s, so we strip the
flanking C and F/W before encoding; predictions/nettcr/<TASK>_TRB.tsv keys on the ORIGINAL full CDR3β.

NOTE: NetTCR-2.0 (2021) is per-peptide-conditioned and only trained on a fixed peptide panel. Of our 5
targets only NLVPMVATV and GLCTLVAML are in its training set; LLW/LLL/YLQ are absent, so the model has
learned nothing for them and scores there are expected to be ~uninformative (ROC~0.5). We still emit
all 5 files so the benchmark can show this coverage limitation. Runs in the cmp-nettcr conda env
(TF 2.13 ARM); data read from manuscript test_data at runtime via _feat_probe.task_table.

    python bench/nettcr_detect.py
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

ENV = "cmp-nettcr"
NET = Path(__file__).resolve().parent / "external" / "NetTCR-2.0"
WORK = Path("bench/out/_nettcr_work")
THRESH = 0.5
TRAIN_B = "data/train_beta_90.csv"                                 # their primary published β set (90% partition)
TRAIN_AB = "data/train_ab_90_alphabeta.csv"                        # their primary published paired set
EPI = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV",
       "YLQ": "YLQPRTFLL", "GLC": "GLCTLVAML"}
PAIRED = ("YLQ", "GLC")


def trim(cdr3: str) -> str:
    """NetTCR encodes IMGT CDR3 with the flanking C and F/W removed (e.g. CASS..EQYF -> ASS..EQY)."""
    s = cdr3
    if s.startswith("C"):
        s = s[1:]
    if s.endswith(("F", "W")):
        s = s[:-1]
    return s


def _roc(d: pl.DataFrame, score_by_cdr3: dict) -> tuple[float, int, int]:
    pairs = [(int(lab), score_by_cdr3[c]) for c, lab in zip(d["cdr3"], d["label"]) if c in score_by_cdr3]
    auc = roc_auc(pairs) if pairs else float("nan")
    return auc, sum(p for p, _ in pairs), len(pairs)


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    pdir = Path("bench/predictions/nettcr"); pdir.mkdir(parents=True, exist_ok=True)
    tt = {t: task_table(t) for t in EPI}

    # --- beta-only: one model, all 5 tasks ---
    beta_specs = []
    for task, pep in EPI.items():
        d = tt[task]
        inp = WORK / f"{task}_b_in.csv"
        d.select(CDR3b=pl.col("cdr3").map_elements(trim, return_dtype=pl.Utf8)).with_columns(
            peptide=pl.lit(pep)).write_csv(inp)
        beta_specs.append(f"{task}:{inp.resolve()}:{(WORK / f'{task}_b_out.csv').resolve()}")
    cmd = ["conda", "run", "-n", ENV, "python", "_run_nettcr.py", "--chain", "b", "--train", TRAIN_B]
    for s in beta_specs:
        cmd += ["--pred", s]
    subprocess.run(cmd, check=True, cwd=NET, capture_output=True, text=True)

    for task, pep in EPI.items():
        d = tt[task]
        out = pl.read_csv(WORK / f"{task}_b_out.csv")
        # map trimmed -> score, then back to original full cdr3
        score_trim = dict(zip(out["CDR3b"].to_list(), out["prediction"].cast(pl.Float64).to_list()))
        rows, score_by_cdr3 = [], {}
        for cdr3 in d["cdr3"].to_list():
            s = score_trim.get(trim(cdr3))
            if s is not None:
                rows.append((cdr3, pep, float(s), int(s >= THRESH)))
                score_by_cdr3[cdr3] = float(s)
        pl.DataFrame(rows, schema=["query_id", "epitope", "score", "significant"],
                     orient="row").write_csv(pdir / f"{task}_TRB.tsv", separator="\t")
        auc, npos, n = _roc(d, score_by_cdr3)
        print(f"nettcr {task}/TRB: {len(rows)} scored | ROC-AUC {auc:.3f} (n_pos={npos} n={n})")

    # --- paired alpha+beta: one model, YLQ + GLC ---
    paired_specs = []
    for task in PAIRED:
        d = tt[task]
        inp = WORK / f"{task}_ab_in.csv"
        d.select(CDR3a=pl.col("a_cdr3").map_elements(trim, return_dtype=pl.Utf8),
                 CDR3b=pl.col("cdr3").map_elements(trim, return_dtype=pl.Utf8)).with_columns(
            peptide=pl.lit(EPI[task])).write_csv(inp)
        paired_specs.append(f"{task}:{inp.resolve()}:{(WORK / f'{task}_ab_out.csv').resolve()}")
    cmd = ["conda", "run", "-n", ENV, "python", "_run_nettcr.py", "--chain", "ab", "--train", TRAIN_AB]
    for s in paired_specs:
        cmd += ["--pred", s]
    subprocess.run(cmd, check=True, cwd=NET, capture_output=True, text=True)

    for task in PAIRED:
        d, pep = tt[task], EPI[task]
        out = pl.read_csv(WORK / f"{task}_ab_out.csv")
        # key on the (trimmed) beta to recover original beta cdr3
        score_trim = dict(zip(out["CDR3b"].to_list(), out["prediction"].cast(pl.Float64).to_list()))
        rows, score_by_cdr3 = [], {}
        for cdr3 in d["cdr3"].to_list():
            s = score_trim.get(trim(cdr3))
            if s is not None:
                rows.append((cdr3, pep, float(s), int(s >= THRESH)))
                score_by_cdr3[cdr3] = float(s)
        pl.DataFrame(rows, schema=["query_id", "epitope", "score", "significant"],
                     orient="row").write_csv(pdir / f"{task}_paired.tsv", separator="\t")
        auc, npos, n = _roc(d, score_by_cdr3)
        print(f"nettcr {task}/paired: {len(rows)} scored | ROC-AUC {auc:.3f} (n_pos={npos} n={n})")


if __name__ == "__main__":
    main()
