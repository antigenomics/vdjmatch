#!/usr/bin/env python3
"""tcrdist3 -> benchmark prediction files for the sample conditions C1/C2/C3 (TRB).

Data prep here (.venv): per condition, build sample CDR3-beta + V/J queries and the VDJdb2026-beta
reference (cdr3, v, j, epitope); run the parallel tcrdist3 compute (cmp-tcrdist env, _tcrdist_compute.py)
which writes 1-NN (tcrdist) and k-NN (tcrdist-knn) predictions per candidate epitope. Read from
manuscript test_data at runtime (never copied).

    python bench/tcrdist_bench.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench
from compare import TESTDATA
from tcrdist_samples import norm_gene
from vdjmatch import db

EPI = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV", "GLC": "GLCTLVAML", "YLQ": "YLQPRTFLL"}
PLAN = {"C1": [EPI["LLW"], EPI["LLL"]], "C2": [EPI["NLV"]], "C3": [EPI["GLC"], EPI["YLQ"]]}
ENV, COMPUTE = "cmp-tcrdist", Path(__file__).resolve().parent / "_tcrdist_compute.py"
WORK = Path("bench/out/_tcrdist_bench")


def queries(cond: str) -> pl.DataFrame:
    if cond == "C1":
        d = pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t").select(
            "cdr3", v="v.segm", j="j.segm")
    elif cond == "C2":
        d = (pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
             .filter(pl.col("gene") == "TRB").select("cdr3", v="v.segm", j="j.segm"))
    else:
        d = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv").select(cdr3="cdr3_beta_aa", v="TRBV", j="TRBJ")
    return _bench.valid_cdr3(d.with_columns(v=norm_gene(pl.col("v")), j=norm_gene(pl.col("j")))).unique("cdr3")


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    vdj = db.load(_bench.source(), species="HomoSapiens").filter(pl.col("gene") == "TRB")
    ref = (_bench.valid_cdr3(vdj).select("cdr3", v=norm_gene(pl.col("v")), j=norm_gene(pl.col("j")),
                                         epitope="epitope").unique("cdr3"))
    reff = WORK / "ref.tsv"; ref.write_csv(reff, separator="\t")
    for d1, d2 in (("tcrdist", "tcrdist-knn"),):
        (Path("bench/predictions") / d1).mkdir(parents=True, exist_ok=True)
        (Path("bench/predictions") / d2).mkdir(parents=True, exist_ok=True)
    for cond, cand in PLAN.items():
        q = queries(cond); qf = WORK / f"{cond}_q.tsv"; q.select("cdr3", "v", "j").write_csv(qf, separator="\t")
        subprocess.run(["conda", "run", "-n", ENV, "python", str(COMPUTE),
                        "--ref", str(reff.resolve()), "--queries", str(qf.resolve()),
                        "--targets", ",".join(cand), "--radius", "90",
                        "--out-1nn", str(Path("bench/predictions/tcrdist") / f"{cond}_TRB.tsv"),
                        "--out-knn", str(Path("bench/predictions/tcrdist-knn") / f"{cond}_TRB.tsv")],
                       check=True)
        print(f"tcrdist {cond}/TRB: {q.height} queries done")


if __name__ == "__main__":
    main()
