"""Full held-out + clustering run: emit result TSVs for the manuscript with all adopted metrics.

  holdout_annotation.tsv : leakage-robust scorer (NED+V+J+len), TRB/TRA x full/shortlist, per-epitope
                           ROC + AUC0.1 + confusion.
  holdout_clustering.tsv : full/shortlist x {TRA,TRB,paired} x {none,apex} with purity, retention, and
                           internal metrics HS/CTS/silhouette/CH (single-chain).

    .venv/bin/python bench/run_full.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path.home() / "vcs/manuscripts/2026-vdjmatch/benchmarks/scripts"))
import benchmark as B                                                # noqa: E402
import cluster_metrics as CM                                         # noqa: E402
import cluster_results as CR                                         # noqa: E402
import cluster_trim as CT                                            # noqa: E402
import holdout_features as HF                                        # noqa: E402

RES = Path.home() / "vcs/manuscripts/2026-vdjmatch/benchmarks/results"
CK = ["ref", "chain", "trim", "n", "purity", "retention", "HS", "CTS", "silhouette", "CH"]


def annotation():
    rows = []
    for locus in ("TRB", "TRA"):
        for ref in ("full", "shortlist"):
            rows += HF.robust_eval(locus, ref)
    pl.DataFrame(rows).write_csv(RES / "holdout_annotation.tsv", separator="\t")
    print(f"[wrote] holdout_annotation.tsv  ({len(rows)} rows)", file=sys.stderr)


def clustering():
    d = B.release("vdjdb2026"); sl = B.shortlist(d)
    rows = []
    for refn, df in (("shortlist", sl), ("full", d)):
        for chain in ("TRA", "TRB"):
            cdr3, epi, v = CR.single_clonotypes(df, chain)
            for trim in (False, True):
                m = CM.metrics(cdr3, v, epi, trim)
                rows.append({"ref": refn, "chain": chain, "trim": "apex" if trim else "none",
                             **{k: round(m[k], 4) if isinstance(m[k], float) else m[k] for k in
                                ("n", "purity", "retention", "HS", "CTS", "silhouette", "CH")}})
        ca, va, cb, vb, epi = CR.paired_clonotypes(d, df)
        for trim in (False, True):
            pur, ret, nc, n = CT.cluster_pr_paired(ca, va, cb, vb, epi, trim)
            rows.append({"ref": refn, "chain": "paired", "trim": "apex" if trim else "none", "n": n,
                         "purity": round(pur, 4), "retention": round(ret, 4),
                         "HS": None, "CTS": None, "silhouette": None, "CH": None})
    pl.DataFrame(rows, schema={k: (pl.Float64 if k in ("purity", "retention", "HS", "CTS", "silhouette", "CH")
                                   else pl.Int64 if k == "n" else pl.Utf8) for k in CK}) \
        .write_csv(RES / "holdout_clustering.tsv", separator="\t")
    print(f"[wrote] holdout_clustering.tsv  ({len(rows)} rows)", file=sys.stderr)


if __name__ == "__main__":
    annotation()
    clustering()
