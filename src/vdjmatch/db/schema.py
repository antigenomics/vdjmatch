"""Canonical VDJdb schema and column normalization.

VDJdb ships dot-named columns (``v.segm``, ``antigen.epitope``, ``mhc.a`` …), which we normalize
to snake_case canonical names so the rest of the package never special-cases the source format.

The exports do not share a *layout*, though. ``vdjdb.slim.txt`` and ``vdjdb.txt`` are one row per
chain, with ``cdr3``/``v.segm``/``j.segm``/``gene``; ``vdjdb_full.txt`` is one row per **complex**,
with the chains side by side in ``cdr3.alpha``/``cdr3.beta`` and no ``gene`` column at all. The
paired layout is unpivoted back to one row per chain on load (:func:`_unpivot_paired`), so both
reach the rest of the package in the same shape.
"""
from __future__ import annotations

import polars as pl

# canonical name -> set of accepted source headers (VDJdb full + slim variants)
_VDJDB_MAP: dict[str, tuple[str, ...]] = {
    "complex_id": ("complex.id",),
    "gene": ("gene",),
    "cdr3": ("cdr3",),
    "v": ("v.segm", "v"),
    "j": ("j.segm", "j"),
    "species": ("species",),
    "mhc_a": ("mhc.a",),
    "mhc_b": ("mhc.b",),
    "mhc_class": ("mhc.class",),
    "epitope": ("antigen.epitope",),
    "antigen_gene": ("antigen.gene",),
    "antigen_species": ("antigen.species",),
    "reference_id": ("reference.id",),
    "vdjdb_score": ("vdjdb.score",),
}

# the columns vdjmatch actually uses downstream (others are dropped on load)
CANONICAL: tuple[str, ...] = tuple(_VDJDB_MAP.keys())


# vdjdb_full.txt keeps the two chains side by side; canonical name -> its alpha/beta source column
_PAIRED: dict[str, dict[str, str]] = {
    "TRA": {"cdr3": "cdr3.alpha", "v": "v.alpha", "j": "j.alpha"},
    "TRB": {"cdr3": "cdr3.beta", "v": "v.beta", "j": "j.beta"},
}


def _unpivot_paired(df: pl.DataFrame) -> pl.DataFrame:
    """Split VDJdb's complex-per-row ``full`` export into one row per chain.

    The chain columns are dropped in favour of ``cdr3``/``v``/``j`` plus a synthesized ``gene``,
    and rows whose chain is absent (the export carries alpha-only and beta-only records as well as
    paired ones) are dropped rather than kept as empty strings. ``complex.id`` is synthesized the
    way VDJdb itself uses it: a shared non-zero id for the two chains of one complex, ``0`` for a
    record that has only one chain.
    """
    chain_cols = {c for cols in _PAIRED.values() for c in cols.values()}
    shared = [c for c in df.columns if c not in chain_cols]
    both = (pl.col("cdr3.alpha").fill_null("") != "") & (pl.col("cdr3.beta").fill_null("") != "")
    df = (df.with_row_index("_i")
            .with_columns(pl.when(both).then(pl.col("_i") + 1).otherwise(0).alias("complex.id")))
    parts = [df.select(*shared, "complex.id",
                       *[pl.col(src).alias(dst) for dst, src in cols.items()],
                       pl.lit(gene).alias("gene"))
               .filter(pl.col("cdr3").fill_null("") != "")
             for gene, cols in _PAIRED.items()]
    return pl.concat(parts, how="vertical")


def normalize(df: pl.DataFrame) -> pl.DataFrame:
    """Rename VDJdb source columns to canonical names; keep only recognized columns.

    Missing canonical columns are filled with nulls so downstream code can rely on the
    schema regardless of whether a ``slim`` or ``full`` export was loaded.
    """
    if "cdr3" not in df.columns and "cdr3.alpha" in df.columns:
        df = _unpivot_paired(df)
    present = set(df.columns)
    renames = {src: canon for canon, srcs in _VDJDB_MAP.items()
               for src in srcs if src in present}
    df = df.rename(renames)
    have = set(df.columns)
    out_cols = []
    for canon in CANONICAL:
        if canon in have:
            out_cols.append(pl.col(canon))
        else:
            out_cols.append(pl.lit(None).alias(canon))
    return df.select(out_cols)
