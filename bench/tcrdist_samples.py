#!/usr/bin/env python3
"""Produce predictions/{tcrdist,tcrdist-knn}/samples.tsv via tcrdist3 (samples mode, EXTERNAL_TOOLS Q2).

Data prep runs here (.venv): build the same query universe as compare.run_samples (sample1/2/5 + V/J),
plus the VDJdb-beta reference (cdr3, v, j, epitope); export to TSVs; then run the tcrdist3 compute in
the `cmp-tcrdist` conda env (bench/_tcrdist_compute.py), which writes the prediction files. Two arms per
"use both": 1-NN (tcrdist) and k-NN (tcrdist-knn). Run from repo root with the project venv.

    python bench/tcrdist_samples.py --olga-n 5000
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench
from compare import SAMPLE_EPI, TESTDATA
from vdjmatch import db

CONDA_ENV = "cmp-tcrdist"


def norm_gene(col: pl.Expr) -> pl.Expr:
    """IMGT-ify: strip, default allele *01 when missing (tcrdist3 wants e.g. TRBV19*01)."""
    g = col.cast(pl.Utf8).str.strip_chars()
    return pl.when(g.str.contains(r"\*")).then(g).otherwise(g + "*01")


def load_sample_vj(name: str) -> pl.DataFrame:
    """Query set -> df[cdr3, v, j] (TRB), matching compare.load_sample's cdr3 universe."""
    if name == "sample1":
        d = (pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
               .filter(pl.col("gene") == "TRB").select("cdr3", v="v.segm", j="j.segm"))
    elif name == "sample2":
        d = (pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t")
               .select("cdr3", v="v.segm", j="j.segm"))
    elif name in ("sample4", "sample5"):                           # OLGA AIRR (4=TRB, 5=TRA)
        d = (pl.read_csv(TESTDATA / f"{name}_olga_airr.txt", separator="\t")
               .select(cdr3="junction_aa", v="v_gene", j="j_gene"))
    else:
        raise ValueError(name)
    d = d.with_columns(norm_gene(pl.col("v")).alias("v"), norm_gene(pl.col("j")).alias("j"))
    return _bench.valid_cdr3(d).unique("cdr3")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--olga-n", type=int, default=0)
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--radius", type=int, default=90)
    ap.add_argument("--knn", type=int, default=5)
    ap.add_argument("--pred-dir", default="bench/predictions")
    ap.add_argument("--work", default="bench/predictions/_tcrdist_work")
    args = ap.parse_args()

    q1, q2, q5 = load_sample_vj("sample1"), load_sample_vj("sample2"), load_sample_vj("sample4")  # TRB OLGA
    if args.olga_n and q5.height > args.olga_n:
        q5 = q5.sample(args.olga_n, seed=0)
    queries = pl.concat([q1, q2, q5]).unique("cdr3")
    vdj = db.load(_bench.source(), species=args.species).filter(pl.col("gene") == "TRB")
    ref = (_bench.valid_cdr3(vdj).select("cdr3", v=norm_gene(pl.col("v")), j=norm_gene(pl.col("j")),
                                         epitope="epitope").unique("cdr3"))
    print(f"queries: {queries.height} unique CDR3; reference: {ref.height} VDJdb-beta clonotypes")

    work = Path(args.work); work.mkdir(parents=True, exist_ok=True)
    queries.select("cdr3", "v", "j").write_csv(work / "queries.tsv", separator="\t")
    ref.write_csv(work / "ref.tsv", separator="\t")

    pred = Path(args.pred_dir)
    (pred / "tcrdist").mkdir(parents=True, exist_ok=True)
    (pred / "tcrdist-knn").mkdir(parents=True, exist_ok=True)
    cmd = ["conda", "run", "-n", CONDA_ENV, "python",
           str(Path(__file__).resolve().parent / "_tcrdist_compute.py"),
           "--ref", str(work / "ref.tsv"), "--queries", str(work / "queries.tsv"),
           "--targets", ",".join(SAMPLE_EPI), "--radius", str(args.radius), "--knn", str(args.knn),
           "--out-1nn", str(pred / "tcrdist" / "samples.tsv"),
           "--out-knn", str(pred / "tcrdist-knn" / "samples.tsv")]
    print("running tcrdist3 compute in conda env", CONDA_ENV, "...")
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
