"""Read query repertoires: single-chain (AIRR Rearrangement / VDJtools) and paired (AIRR Cell
long-form linked by cell/clone id, or TCRvdb wide-form with α/β columns in one row).
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl

from . import columns

_SEP = {".csv": ",", ".tsv": "\t", ".txt": "\t"}


def _read_table(path: str | os.PathLike) -> pl.DataFrame:
    p = Path(path)
    sep = _SEP.get(p.suffix.lower(), "\t")
    return pl.read_csv(p, separator=sep, quote_char='"', infer_schema_length=0)


def read_rearrangement(path: str | os.PathLike, dedup: bool = True,
                       valid_aa: bool = True) -> pl.DataFrame:
    """Single-chain query → schema ``cdr3, v, j, locus, count, pair_id``. If ``dedup`` (default),
    collapse to unique (cdr3, v, j, locus) summing ``count`` (required for E-values)."""
    df = columns.normalize_query(_read_table(path), valid_aa=valid_aa)
    if dedup:
        df = (df.group_by("cdr3", "v", "j", "locus")
                .agg(pl.col("count").sum(), pl.col("pair_id").first()))
    return df


def read_cell(path: str | os.PathLike, link: str | None = None,
              valid_aa: bool = True) -> pl.DataFrame:
    """Paired long-form: one row per chain, α/β linked by ``link`` (default: auto = cell_id /
    clone_id, resolved to ``pair_id``). Pivots to one row per cell with ``cdr3a/va/ja`` and
    ``cdr3b/vb/jb`` (locus prefix TRA→a, TRB→b)."""
    df = columns.normalize_query(_read_table(path), valid_aa=valid_aa)
    if df["pair_id"].null_count() == df.height:
        raise ValueError("no cell/clone link column found for paired (cell) input")
    df = df.with_columns(
        pl.when(pl.col("locus") == "TRA").then(pl.lit("a"))
          .when(pl.col("locus") == "TRB").then(pl.lit("b"))
          .otherwise(pl.col("locus").str.slice(2, 1).str.to_lowercase()).alias("ch"))
    wide = df.pivot(values=["cdr3", "v", "j"], index="pair_id", on="ch",
                    aggregate_function="first")
    rename = {f"{c}_{ch}": f"{c}{ch}" for c in ("cdr3", "v", "j") for ch in ("a", "b")}
    return wide.rename({k: v for k, v in rename.items() if k in wide.columns})


def read_tcrvdb(path: str | os.PathLike, valid_aa: bool = True) -> pl.DataFrame:
    """TCRvdb wide-form CSV: α/β CDR3s in one row. → ``cdr3a, va, ja, cdr3b, vb, jb,
    epitope, mhc`` (epitope/mhc kept as ground-truth labels for benchmarking)."""
    df = _read_table(path)
    cols = {c.lower(): c for c in df.columns}
    g = lambda *names: next((cols[n] for n in names if n in cols), None)  # noqa: E731
    sel = {
        "cdr3a": g("cdr3_alpha_aa"), "va": g("trav"), "ja": g("traj"),
        "cdr3b": g("cdr3_beta_aa"), "vb": g("trbv"), "jb": g("trbj"),
        "epitope": g("epitope_aa", "antigen.epitope", "epitope"),
        "mhc": g("hla_short", "hla_long", "mhc.a"),
    }
    out = df.select([pl.col(src).alias(dst) for dst, src in sel.items() if src is not None])
    out = out.with_columns(
        pl.col("cdr3a").str.strip_chars().str.to_uppercase(),
        pl.col("cdr3b").str.strip_chars().str.to_uppercase())
    if valid_aa:
        rx = r"^[ACDEFGHIKLMNPQRSTVWY]+$"
        out = out.filter(pl.col("cdr3a").str.contains(rx) & pl.col("cdr3b").str.contains(rx))
    return out
