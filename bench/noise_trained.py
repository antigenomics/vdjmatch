#!/usr/bin/env python3
"""OLGA-noise false-positive rate for the TRAINED classifiers (extends Fig 3a of the benchmark).

A "noise FP" = the method confidently calls a RANDOM OLGA receptor a binder of SOME benchmark epitope.
Random receptors are the same 1000-per-locus OLGA set used by benchmark.cond_olga_fp (sample4=TRB,
sample5=TRA): valid_cdr3 -> unique(cdr3) -> sample(1000, seed=0), read from the manuscript test_data
at runtime (never copied). The 5-epitope detection panel is all HLA-A*02:01 / CD8:

    NLVPMVATV  LLWNGPMAV  LLLGIGILV  YLQPRTFLL  GLCTLVAML

Per (method, locus): FP rate = fraction of the 1000 randoms whose score crosses the method's standard
call threshold (0.5) for ANY panel epitope.

  imw-DETECT  (TRB+TRA): top-1 epitope + Score per TCR from DETECT's OWN vocabulary; FP if that top
              Score >= 0.5 (a confident call on a random). DETECT was already run on sample4/sample5
              by ImmuneWatch -> read predictions_sample{4,5}.tsv.ods directly (no env, no model call).
  ERGO-II     (TRB; beta-trained): score each random x each of the 5 peptides (MHC=HLA-A*02, CD8); FP
              if max-over-peptides P(bind) >= 0.5. Reuses ERGO-II/_run_predict.py in the cmp-ergo env.
              TRA: ERGO's VDJdb model is beta-conditioned (no alpha-only path), so per ergo_detect_alpha.py
              we feed the alpha CDR3 (+TRAV/TRAJ) through ERGO's TCR encoder via the TRB/TRBV/TRBJ slots.
  NetTCR-2.0  (TRB; beta model trained once on data/train_beta_90.csv): score each random x each of the 5
              peptides; FP if max >= 0.5. Reuses NetTCR-2.0/_run_nettcr.py in the cmp-nettcr env. TRA is
              SKIPPED: the beta model cannot accept an alpha chain and the paired model needs both chains
              (no single-chain alpha random set exists).
  MixTCRpred  SKIPPED: paired-only model; no single-chain random set exists.

Writes per-method detail files + the summary table under bench/predictions/_noise/. Run with the
vdjmatch venv python (it shells into the cmp-* conda envs like the other producers):

    .venv/bin/python bench/noise_trained.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench
from benchmark import vgene
from compare import TESTDATA

THRESH = 0.5
MHC = "HLA-A*02"
N = 1000
SEED = 0
PANEL = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV",
         "YLQ": "YLQPRTFLL", "GLC": "GLCTLVAML"}
PEPTIDES = list(PANEL.values())
FILES = {"TRB": "sample4", "TRA": "sample5"}                      # cond_olga_fp: 4=TRB, 5=TRA

OUT = Path("bench/predictions/_noise"); OUT.mkdir(parents=True, exist_ok=True)
WORK = OUT / "_work"; WORK.mkdir(parents=True, exist_ok=True)
ERGO = Path(__file__).resolve().parent / "external" / "ERGO-II"
NET = Path(__file__).resolve().parent / "external" / "NetTCR-2.0"


def noise_set(locus: str) -> pl.DataFrame:
    """The canonical 1000-receptor OLGA noise set for a locus (cdr3, v, j) — identical sampling to
    benchmark.cond_olga_fp so the trained-method FP rate is comparable to the matching methods."""
    d = (pl.read_csv(TESTDATA / f"{FILES[locus]}_olga_airr.txt", separator="\t")
         .select(cdr3="junction_aa", v="v_gene", j="j_gene")
         .pipe(_bench.valid_cdr3).unique("cdr3"))
    if d.height > N:
        d = d.sample(N, seed=SEED)
    return d.with_columns(v=pl.col("v").map_elements(vgene, return_dtype=pl.Utf8),
                          j=pl.col("j").map_elements(vgene, return_dtype=pl.Utf8))


# --------------------------------------------------------------------------- imw-DETECT (no model call)
def run_imw(rows: list):
    ODS = TESTDATA / "immunedetect_results"
    for locus in ("TRB", "TRA"):
        canon = set(noise_set(locus)["cdr3"].to_list())
        od = pd.read_excel(ODS / f"predictions_{FILES[locus]}.tsv.ods", engine="odf")
        od = od.rename(columns={"junction_aa": "cdr3", "Epitope": "epitope", "Score": "score"})
        sub = od[od["cdr3"].astype(str).isin(canon)].copy()
        # one top-1 prediction per receptor; FP iff that confident call crosses threshold
        sub = sub.sort_values("score", ascending=False).drop_duplicates("cdr3", keep="first")
        sub["significant"] = (sub["score"].astype(float) >= THRESH).astype(int)
        det = sub[["cdr3", "epitope", "score", "significant"]].rename(columns={"cdr3": "query_id"})
        det.to_csv(OUT / f"imw-detect_{locus}.tsv", sep="\t", index=False)
        n, n_fp = len(canon), int(det["significant"].sum())
        rows.append(("imw-DETECT", locus, n, n_fp))
        print(f"imw-DETECT {locus}: {n_fp}/{n} FP (top-1 Score>={THRESH}) | scored {len(det)}")


# --------------------------------------------------------------------------- ERGO-II (cmp-ergo env)
def run_ergo(rows: list):
    for locus in ("TRB", "TRA"):
        d = noise_set(locus)
        # one ERGO input file: each random x each panel peptide. For TRA, the alpha CDR3 (+TRAV/TRAJ)
        # goes through the TRB/TRBV/TRBJ slots (ERGO has no alpha-only path); true TRA left empty.
        recs = []
        for cdr3, v, j in zip(d["cdr3"], d["v"], d["j"]):
            for pep in PEPTIDES:
                recs.append((cdr3, v, j, pep))
        ergo_in = pl.DataFrame(recs, schema=["TRB", "TRBV", "TRBJ", "Peptide"], orient="row").with_columns(
            TRA=pl.lit(None, dtype=pl.Utf8), TRAV=pl.lit(None, dtype=pl.Utf8),
            TRAJ=pl.lit(None, dtype=pl.Utf8)).with_columns(
            **{"T-Cell-Type": pl.lit("CD8")}, MHC=pl.lit(MHC)).select(
            "TRA", "TRB", "TRAV", "TRAJ", "TRBV", "TRBJ", "T-Cell-Type", "Peptide", "MHC")
        inp = WORK / f"ergo_{locus}_in.csv"
        outp = WORK / f"ergo_{locus}_out.csv"
        ergo_in.write_csv(inp)
        subprocess.run(["conda", "run", "-n", "cmp-ergo", "python", "_run_predict.py",
                        "vdjdb", str(inp.resolve()), str(outp.resolve())],
                       check=True, cwd=ERGO, capture_output=True, text=True)
        out = pl.read_csv(outp)
        # max P(bind) over the 5 peptides, per (TRB,Peptide) -> per receptor
        best = (out.group_by("TRB").agg(pl.col("Score").cast(pl.Float64).max().alias("score")))
        score_by_cdr3 = dict(zip(best["TRB"].to_list(), best["score"].to_list()))
        det = []
        for cdr3 in d["cdr3"].to_list():
            s = score_by_cdr3.get(cdr3)
            if s is not None:
                det.append((cdr3, float(s), int(s >= THRESH)))
        pl.DataFrame(det, schema=["query_id", "max_score", "significant"], orient="row").write_csv(
            OUT / f"ergo_{locus}.tsv", separator="\t")
        n, n_fp = len(det), sum(sig for _, _, sig in det)
        rows.append(("ERGO-II", locus, n, n_fp))
        note = " (alpha via beta-encoder)" if locus == "TRA" else ""
        print(f"ERGO-II {locus}{note}: {n_fp}/{n} FP (max P(bind)>={THRESH}) | scored {n} of {d.height}")


# --------------------------------------------------------------------------- NetTCR-2.0 (cmp-nettcr env)
def trim(cdr3: str) -> str:
    s = cdr3
    if s.startswith("C"):
        s = s[1:]
    if s.endswith(("F", "W")):
        s = s[:-1]
    return s


def run_nettcr(rows: list):
    locus = "TRB"                                                  # beta model only; alpha has no path
    d = noise_set(locus)
    recs = []
    for cdr3 in d["cdr3"].to_list():
        for pep in PEPTIDES:
            recs.append((trim(cdr3), pep))
    inp = WORK / "nettcr_TRB_in.csv"
    outp = WORK / "nettcr_TRB_out.csv"
    pl.DataFrame(recs, schema=["CDR3b", "peptide"], orient="row").write_csv(inp)
    subprocess.run(["conda", "run", "-n", "cmp-nettcr", "python", "_run_nettcr.py",
                    "--chain", "b", "--train", "data/train_beta_90.csv",
                    "--pred", f"TRB:{inp.resolve()}:{outp.resolve()}"],
                   check=True, cwd=NET, capture_output=True, text=True)
    out = pl.read_csv(outp)
    best = (out.group_by("CDR3b").agg(pl.col("prediction").cast(pl.Float64).max().alias("score")))
    score_trim = dict(zip(best["CDR3b"].to_list(), best["score"].to_list()))
    det = []
    for cdr3 in d["cdr3"].to_list():
        s = score_trim.get(trim(cdr3))
        if s is not None:
            det.append((cdr3, float(s), int(s >= THRESH)))
    pl.DataFrame(det, schema=["query_id", "max_score", "significant"], orient="row").write_csv(
        OUT / "nettcr_TRB.tsv", separator="\t")
    n, n_fp = len(det), sum(sig for _, _, sig in det)
    rows.append(("NetTCR-2.0", locus, n, n_fp))
    print(f"NetTCR-2.0 {locus}: {n_fp}/{n} FP (max score>={THRESH}) | scored {n} of {d.height}")
    print("NetTCR-2.0 TRA: SKIPPED (beta model takes no alpha; paired model needs both chains, "
          "no single-chain alpha random set)")


# method -> (label, [(locus, detail_file)]); each detail file has a 'significant' 0/1 column
DETAIL = {
    "imw-DETECT": [("TRB", "imw-detect_TRB.tsv"), ("TRA", "imw-detect_TRA.tsv")],
    "ERGO-II": [("TRB", "ergo_TRB.tsv"), ("TRA", "ergo_TRA.tsv")],
    "NetTCR-2.0": [("TRB", "nettcr_TRB.tsv")],
}


def assemble():
    """Build noise_fp_trained.tsv by tallying the per-method detail files (single source of truth)."""
    rows = []
    for method, specs in DETAIL.items():
        for locus, fname in specs:
            f = OUT / fname
            if not f.exists():
                print(f"{method} {locus}: detail file missing ({fname}); skipping in summary")
                continue
            d = pl.read_csv(f, separator="\t")
            rows.append((method, locus, d.height, int(d["significant"].sum())))
    summary = pl.DataFrame(rows, schema=["method", "locus", "n", "n_fp"], orient="row").with_columns(
        fp_pct=(pl.col("n_fp") / pl.col("n") * 100).round(3))
    summary.write_csv(OUT / "noise_fp_trained.tsv", separator="\t")
    print("\n=== noise_fp_trained.tsv ===")
    with pl.Config(tbl_rows=-1):
        print(summary)


def main():
    rows: list = []
    run_imw(rows)
    run_ergo(rows)
    run_nettcr(rows)
    print("MixTCRpred: SKIPPED (paired-only model; no single-chain random set exists)")
    assemble()


if __name__ == "__main__":
    if "--assemble" in sys.argv:                                  # rebuild summary from existing detail
        assemble()
    else:
        main()
