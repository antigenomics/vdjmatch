"""Merge vdjmatch (leakage-robust) + competitor held-out results into one comparison table.

vdjmatch from holdout_annotation.tsv; competitors from competitors_holdout.tsv. Emits
holdout_comparison.tsv (locus, epitope, then ROC and AUC0.1 per method) + macro rows. Competitor cells
are blank where the method does not cover that (epitope, locus). All trained competitors carry a VDJdb
training-leakage caveat (their held-out positives are in their training data); vdjmatch is leakage-robust.

    .venv/bin/python bench/compete_compare.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

RES = Path.home() / "vcs/manuscripts/2026-vdjmatch/benchmarks/results"
METHODS = ["vdjmatch", "tcrbert", "pmtnet", "tcrgp"]
EPI_ORDER = ["NLV", "LLW", "LLL", "ELA", "YLQ", "GLC"]


def main():
    vdj = (pl.read_csv(RES / "holdout_annotation.tsv", separator="\t")
           .filter(pl.col("ref") == "full")
           .select("locus", "epitope", roc="roc", auc01="auc01")
           .with_columns(tool=pl.lit("vdjmatch")))
    comp = pl.read_csv(RES / "competitors_holdout.tsv", separator="\t") if (RES / "competitors_holdout.tsv").exists() \
        else pl.DataFrame(schema={"tool": pl.Utf8, "locus": pl.Utf8, "epitope": pl.Utf8, "n_pos": pl.Int64,
                                  "roc": pl.Float64, "auc01": pl.Float64})
    allm = pl.concat([vdj.select("tool", "locus", "epitope", "roc", "auc01"),
                      comp.select("tool", "locus", "epitope", "roc", "auc01")], how="vertical")
    rows = []
    for locus in ("TRB", "TRA"):
        for sh in EPI_ORDER:
            r = {"locus": locus, "epitope": sh}
            for m in METHODS:
                cell = allm.filter((pl.col("tool") == m) & (pl.col("locus") == locus) & (pl.col("epitope") == sh))
                r[f"{m}_roc"] = round(cell["roc"][0], 3) if cell.height and cell["roc"][0] is not None else None
                r[f"{m}_auc01"] = round(cell["auc01"][0], 3) if cell.height and cell["auc01"][0] is not None else None
            rows.append(r)
    df = pl.DataFrame(rows)
    df.write_csv(RES / "holdout_comparison.tsv", separator="\t")
    # macro per method/locus (mean over epitopes where the method has a value)
    print("=== held-out comparison (ROC; full A*02 reference) ===")
    with pl.Config(tbl_rows=20, tbl_cols=20):
        print(df.select("locus", "epitope", *[f"{m}_roc" for m in METHODS]))
    print("\n=== macro ROC / AUC0.1 per method (mean over covered epitopes) ===")
    for locus in ("TRB", "TRA"):
        sub = allm.filter(pl.col("locus") == locus)
        for m in METHODS:
            mm = sub.filter((pl.col("tool") == m) & pl.col("roc").is_not_nan() & pl.col("roc").is_not_null())
            if mm.height:
                print(f"  {locus} {m:9} macroROC {mm['roc'].mean():.3f}  macroAUC0.1 {mm['auc01'].mean():.3f}  "
                      f"(n_epi={mm.height})")
    return df


if __name__ == "__main__":
    main()
