"""Canonical VDJdb schema and column normalization.

VDJdb ships dot-named columns (``v.segm``, ``antigen.epitope``, ``mhc.a`` …) in both the
``full`` and ``slim`` exports. We normalize to snake_case canonical names so the rest of the
package never special-cases the source format.
"""
from __future__ import annotations

import polars as pl

# canonical name -> set of accepted source headers (VDJdb full + slim variants)
_VDJDB_MAP: dict[str, tuple[str, ...]] = {
    "complex_id": ("complex.id",),
    "gene": ("gene",),
    "cdr3": ("cdr3",),
    "v": ("v.segm", "v.end", "v"),
    "j": ("j.segm", "j.start", "j"),
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


def normalize(df: pl.DataFrame) -> pl.DataFrame:
    """Rename VDJdb source columns to canonical names; keep only recognized columns.

    Missing canonical columns are filled with nulls so downstream code can rely on the
    schema regardless of whether a ``slim`` or ``full`` export was loaded.
    """
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
