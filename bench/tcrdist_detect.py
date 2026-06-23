#!/usr/bin/env python3
"""tcrdist3 detection producer (beta tasks): match each query against ONE epitope's A*02 VDJdb2026
reference (exact removed) -> 1-NN / k-NN scores. Emits predictions/{tcrdist,tcrdist-knn}/<task>_TRB.tsv
and prints detection ROC/PR/recall so we can compare to vdjmatch (esp. NLV recall).

    python bench/tcrdist_detect.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench
import benchmark as B
from compare import TESTDATA
from metrics import pr_auc_balanced, roc_auc
from tcrdist_samples import norm_gene
from vdjmatch import db

ENV, COMPUTE = "cmp-tcrdist", Path(__file__).resolve().parent / "_tcrdist_compute.py"
WORK = Path("bench/out/_tcrdist_detect"); WORK.mkdir(parents=True, exist_ok=True)
RADIUS, SIGR = 90, 24


def sample_vj(task):
    """query cdr3 -> (v,j) for a TRB task, from the source sample."""
    if task == "NLV":
        d = pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t").filter(pl.col("gene") == "TRB")
        d = d.select("cdr3", v="v.segm", j="j.segm")
    elif task in ("LLW", "LLL"):
        d = pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t").select("cdr3", v="v.segm", j="j.segm")
    else:
        d = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv").select(cdr3="cdr3_beta_aa", v="TRBV", j="TRBJ")
    d = _bench.valid_cdr3(d.with_columns(v=norm_gene(pl.col("v")), j=norm_gene(pl.col("j")))).unique("cdr3")
    return {c: (v, j) for c, v, j in zip(d["cdr3"], d["v"], d["j"])}


def main():
    v26 = db.load(_bench.source(), species="HomoSapiens")
    for task, estr, locus, pos_qv, neg_qv in B.detection_tasks(["TRB"]):
        vj = sample_vj(task)
        ref = (_bench.valid_cdr3(v26.filter((pl.col("epitope") == estr) & pl.col("mhc_a").str.contains(B.A02)
                                            & (pl.col("gene") == "TRB")))
               .select("cdr3", v=norm_gene(pl.col("v")), j=norm_gene(pl.col("j")),
                       epitope="epitope").unique("cdr3"))
        reff = WORK / f"{task}_ref.tsv"; ref.write_csv(reff, separator="\t")
        allq = [q for q in (list(pos_qv) + list(neg_qv)) if q in vj]
        qf = WORK / f"{task}_q.tsv"
        pl.DataFrame([(q, vj[q][0], vj[q][1]) for q in allq], schema=["cdr3", "v", "j"],
                     orient="row").write_csv(qf, separator="\t")
        for d in ("tcrdist", "tcrdist-knn"):
            (Path("bench/predictions") / d).mkdir(parents=True, exist_ok=True)
        subprocess.run(["conda", "run", "-n", ENV, "python", str(COMPUTE), "--ref", str(reff.resolve()),
                        "--queries", str(qf.resolve()), "--targets", estr, "--radius", str(RADIUS),
                        "--sig-radius", str(SIGR),
                        "--out-1nn", str(Path("bench/predictions/tcrdist") / f"{task}_TRB.tsv"),
                        "--out-knn", str(Path("bench/predictions/tcrdist-knn") / f"{task}_TRB.tsv")],
                       check=True)
        # read back, compute detection metrics
        for arm in ("tcrdist", "tcrdist-knn"):
            sc = B.read_predictions(Path("bench/predictions") / arm / f"{task}_TRB.tsv",
                                    list(pos_qv) + list(neg_qv), [estr])
            m = B.classify_metrics(sc, list(pos_qv), list(neg_qv), estr)
            print(f"{task}/TRB/{arm}: ROC={m['roc_auc']:.3f} PR={m['pr_auc']:.3f} "
                  f"rec={m['recall']:.2f} prec={m['precision']:.2f} TP={m['tp']} FP={m['fp']}")


if __name__ == "__main__":
    main()
