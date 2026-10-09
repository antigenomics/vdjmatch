"""Source-aware normalization of query junctions and their original identities."""

from __future__ import annotations

import polars as pl

AA = set("ACDEFGHIKLMNPQRSTVWY")
ALIASES: dict[str, tuple[str, ...]] = {
    "cdr3": ("junction_aa", "cdr3_aa", "cdr3"),
    "v": ("v_call", "v_gene", "v.segm", "v"),
    "j": ("j_call", "j_gene", "j.segm", "j"),
    "locus": ("locus", "chain"),
    "count": (
        "duplicate_count",
        "reads",
        "count",
        "clonecount",
        "clone count",
        "readcount",
        "uniquetagcountmolecule",
    ),
    "pair_id": ("pair_id", "cell_id", "clone_id", "complex.id"),
}


def _resolve(df: pl.DataFrame) -> dict[str, str]:
    lower = {c.lower(): c for c in df.columns}
    return {
        canon: next(lower[s] for s in srcs if s in lower)
        for canon, srcs in ALIASES.items()
        if any(s in lower for s in srcs)
    }


def _integer_count(name: str) -> pl.Expr:
    value = pl.col(name).cast(pl.String).str.strip_chars()
    return (
        pl.when(value.str.contains(r"^[+]?[0-9]+(?:\.0+)?$"))
        .then(value.str.replace(r"\.0+$", "").cast(pl.Int64, strict=False))
        .otherwise(pl.lit(None, dtype=pl.Int64))
    )


def validate_counts(df: pl.DataFrame) -> None:
    """Reject missing, malformed, fractional, negative or overflowing supplied counts."""
    lower = {c.lower().lstrip("#"): c for c in df.columns}
    name = _resolve(df).get("count") or next(
        (
            lower[c]
            for c in (
                "count",
                "clone count",
                "clonecount",
                "readcount",
                "uniquetagcountmolecule",
                "read.count",
                "number of reads",
                "count (templates/reads)",
                "templates",
            )
            if c in lower
        ),
        None,
    )
    if name is None:
        return
    if df.select(_integer_count(name).is_null().any()).item():
        raise ValueError(
            f"{name} must contain non-negative integer counts without missing values"
        )


def normalize_query(
    df: pl.DataFrame,
    valid_aa: bool = True,
    *,
    source: str = "auto",
    sequence_convention: str | None = None,
    return_report: bool = False,
):
    """Normalize junctions, preserving supplied query/sequence identities and row order.

    ``source='airr'`` requires ``junction_aa``. ``source='custom'`` requires an
    explicit ``sequence_convention='junction'`` for ambiguous sequence columns.
    Auto accepts legacy short-form producers, but recognizes AIRR ``v_call`` or
    ``sequence_id`` headers. Bare IMGT CDR3s cannot be converted without anchors.
    ``return_report`` returns ``(frame, counts)`` with retained/dropped row counts.
    """
    if source not in {"auto", "airr", "legacy", "custom"}:
        raise ValueError(f"unknown query source {source!r}")
    if sequence_convention not in {None, "junction", "cdr3"}:
        raise ValueError("sequence_convention must be 'junction' or 'cdr3'")
    found = _resolve(df)
    if "cdr3" not in found:
        raise ValueError(f"no junction amino-acid column found; have {df.columns}")
    lower = {c.lower() for c in df.columns}
    is_airr = source == "airr" or (
        source == "auto" and bool(lower & {"v_call", "sequence_id"})
    )
    junction = found["cdr3"].lower() == "junction_aa"
    if not junction:
        if is_airr:
            raise ValueError(
                "AIRR input requires junction_aa; cdr3_aa excludes conserved anchors"
            )
        if sequence_convention == "cdr3":
            raise ValueError(
                "bare CDR3 input cannot be matched without a supplied junction_aa"
            )
        if source == "custom" and sequence_convention != "junction":
            raise ValueError("custom input requires sequence_convention='junction'")
    validate_counts(df)
    n_input = df.height
    # Select aliases instead of renaming: both junction_aa and legacy cdr3 may coexist.
    identity = [c for c in ("query_id", "sequence_id") if c in df.columns]
    df = df.select(
        [pl.col(src).alias(canon) for canon, src in found.items()] + identity
    )
    df = df.with_columns(
        pl.col("cdr3").cast(pl.String).str.strip_chars().str.to_uppercase()
    )
    for col in ("v", "j"):
        if col not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=pl.String).alias(col))
        else:
            df = df.with_columns(pl.col(col).cast(pl.String))
    if "locus" not in df.columns:
        df = df.with_columns(pl.col("v").str.slice(0, 3).alias("locus"))
    else:
        df = df.with_columns(pl.col("locus").cast(pl.String).str.to_uppercase())
    if "count" not in df.columns:
        df = df.with_columns(pl.lit(1, dtype=pl.Int64).alias("count"))
    else:
        df = df.with_columns(_integer_count("count").alias("count"))
    if "pair_id" not in df.columns:
        df = df.with_columns(pl.lit(None, dtype=pl.String).alias("pair_id"))
    else:
        df = df.with_columns(pl.col("pair_id").cast(pl.String))
    missing = pl.col("cdr3").is_null() | (pl.col("cdr3").str.len_chars() == 0)
    n_missing = df.select(missing.sum()).item()
    df = df.filter(~missing)
    n_nonmissing = df.height
    if valid_aa:
        df = df.filter(pl.col("cdr3").str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$"))
    report = {
        "input_rows": n_input,
        "missing_sequences": n_missing,
        "invalid_sequences": n_nonmissing - df.height,
        "retained_rows": df.height,
        "dropped_rows": n_input - df.height,
    }
    df = df.select("cdr3", "v", "j", "locus", "count", "pair_id", *identity)
    return (df, report) if return_report else df
