"""Column alias resolution for query repertoire files.

Accepts AIRR Rearrangement canonical headers and the VDJtools / short forms seen in the wild
(``v_gene``, ``j_gene``, ``junction_aa``, ``clone_id`` …). Normalizes to the vdjmatch query
schema: ``cdr3, v, j, locus, count, pair_id``.
"""
from __future__ import annotations

import polars as pl

AA = set("ACDEFGHIKLMNPQRSTVWY")

# canonical query column -> accepted source headers (case-insensitive), most specific first
ALIASES: dict[str, tuple[str, ...]] = {
    "cdr3": ("cdr3_aa", "junction_aa", "cdr3"),     # note: AIRR `junction`/`cdr3` (nt) excluded by suffix
    "v": ("v_call", "v_gene", "v.segm", "v"),
    "j": ("j_call", "j_gene", "j.segm", "j"),
    "locus": ("locus", "chain"),
    "count": ("duplicate_count", "reads", "count"),
    "pair_id": ("cell_id", "clone_id", "complex.id"),
}


def _resolve(df: pl.DataFrame) -> dict[str, str]:
    lower = {c.lower(): c for c in df.columns}
    out = {}
    for canon, srcs in ALIASES.items():
        for s in srcs:
            if s in lower:
                out[canon] = lower[s]
                break
    return out


def normalize_query(df: pl.DataFrame, valid_aa: bool = True) -> pl.DataFrame:
    """Rename to the query schema, derive ``locus`` from the V gene when absent, default
    ``count`` to 1, and (if ``valid_aa``) drop CDR3s containing non-standard amino acids."""
    found = _resolve(df)
    if "cdr3" not in found:
        raise ValueError(f"no CDR3 column found; have {df.columns}")
    df = df.rename({src: canon for canon, src in found.items()})
    df = df.with_columns(pl.col("cdr3").str.strip_chars().str.to_uppercase())

    for col in ("v", "j"):
        if col not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias(col))
    if "locus" not in df.columns:
        df = df.with_columns(pl.col("v").str.slice(0, 3).alias("locus"))
    if "count" not in df.columns:
        df = df.with_columns(pl.lit(1, dtype=pl.Int64).alias("count"))
    else:
        df = df.with_columns(pl.col("count").cast(pl.Int64, strict=False).fill_null(1))
    if "pair_id" not in df.columns:
        df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("pair_id"))
    else:
        df = df.with_columns(pl.col("pair_id").cast(pl.Utf8))

    df = df.filter(pl.col("cdr3").is_not_null() & (pl.col("cdr3").str.len_chars() > 0))
    if valid_aa:
        df = df.filter(pl.col("cdr3").str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$"))
    return df.select("cdr3", "v", "j", "locus", "count", "pair_id")
