"""Normalize reference names while retaining every source metadata column.

Matching uses anchored junctions: legacy VDJdb ``cdr3`` and AIRR ``junction_aa``.
AIRR ``cdr3_aa`` alone is never interpreted as a junction.
"""

from __future__ import annotations

import polars as pl

_VDJDB_MAP = {
    "complex_id": ("complex.id",),
    "gene": ("gene", "locus"),
    "cdr3": ("cdr3", "junction_aa"),
    "v": ("v.segm", "v_segm", "v_call", "v"),
    "j": ("j.segm", "j_segm", "j_call", "j"),
    "species": ("species",),
    "mhc_a": ("mhc.a", "mhc_allele_1"),
    "mhc_b": ("mhc.b", "mhc_allele_2"),
    "mhc_class": ("mhc.class",),
    "epitope": ("antigen.epitope", "antigen_epitope", "peptide_sequence_aa"),
    "antigen_gene": ("antigen.gene", "antigen"),
    "antigen_species": ("antigen.species", "antigen_source_species"),
    "reference_id": ("reference.id", "reactivity_refs"),
    "vdjdb_score": ("vdjdb.score",),
}
CANONICAL = tuple(_VDJDB_MAP)
_PAIRED = {
    "TRA": {"cdr3": "cdr3.alpha", "v": "v.alpha", "j": "j.alpha"},
    "TRB": {"cdr3": "cdr3.beta", "v": "v.beta", "j": "j.beta"},
}


def _unpivot_paired(df: pl.DataFrame) -> pl.DataFrame:
    both = (pl.col("cdr3.alpha").fill_null("") != "") & (
        pl.col("cdr3.beta").fill_null("") != ""
    )
    # Retain the original wide fields, including repairs and submitted calls.
    df = df.with_row_index("__full_row")
    source_link = next((c for c in ("complex_id", "complex.id") if c in df.columns), None)
    if source_link:
        key = pl.col(source_link).cast(pl.String)
        declared = key.is_not_null() & ~key.is_in(["", "0"])
        # Explicit links win; physical wide rows still link otherwise missing
        # partners. A separate row namespace avoids merging legacy numeric links.
        fallback = pl.concat_str(pl.lit("__wide_row:"), pl.col("__full_row"))
        generated = df.filter(~declared.fill_null(False) & both).select(fallback.alias("key"))
        existing = df.filter(declared).select(key.alias("key"))
        if generated.join(existing, on="key", how="inner").height:
            raise ValueError("generated wide complex identity collides with a declared link")
        link = pl.when(declared).then(key).otherwise(
            pl.when(both).then(fallback).otherwise(pl.lit("0"))
        )
    else:
        link = pl.when(both).then(pl.col("__full_row") + 1).otherwise(0)
    df = df.with_columns(link.alias("complex_id"))
    parts = [
        df.with_columns(
            *[
                pl.col(src).alias(dst) if src in df.columns else pl.lit(None).alias(dst)
                for dst, src in cols.items()
            ],
            pl.lit(gene).alias("gene"),
        )
        .filter(pl.col("cdr3").fill_null("") != "")
        .drop("__full_row")
        for gene, cols in _PAIRED.items()
    ]
    return pl.concat(parts, how="vertical")


def normalize(df: pl.DataFrame) -> pl.DataFrame:
    """Add canonical matching fields without discarding source fields or identifiers."""
    if "cdr3" not in df.columns and "cdr3.alpha" in df.columns:
        if "cdr3.beta" not in df.columns:
            df = df.with_columns(pl.lit("").alias("cdr3.beta"))
        df = _unpivot_paired(df)
    if not any(c in df.columns for c in ("cdr3", "junction_aa")):
        raise ValueError(
            "reference requires VDJdb cdr3 or AIRR junction_aa; bare cdr3_aa has no anchors"
        )
    if not any(c in df.columns for c in ("gene", "locus")):
        raise ValueError("reference requires gene/locus")
    if not any(
        c in df.columns
        for c in (
            "epitope",
            "antigen.epitope",
            "antigen_epitope",
            "peptide_sequence_aa",
        )
    ):
        raise ValueError("reference requires peptide annotation (epitope)")
    exprs = []
    for canon, aliases in _VDJDB_MAP.items():
        if canon not in df.columns:
            src = next((s for s in aliases if s in df.columns), None)
            exprs.append(
                pl.col(src).alias(canon)
                if src
                else pl.lit(None, dtype=pl.String).alias(canon)
            )
    df = df.with_columns(exprs)
    if df.filter(
        pl.any_horizontal(
            *(
                pl.col(c).is_null() | (pl.col(c).cast(pl.String) == "")
                for c in ("gene", "epitope")
            )
        )
    ).height:
        raise ValueError("reference has missing junction, locus or peptide annotation")
    if df.filter(~pl.col("gene").is_in(["TRA", "TRB"])).height:
        raise ValueError("VDJdb reference supports only TRA/TRB loci")
    if (
        "reactivity_value" in df.columns
        and "vdjdb.score" not in df.columns
        and df["vdjdb_score"].null_count() == df.height
    ):
        df = df.with_columns(
            pl.when(
                (pl.col("reactivity_readout") == "confidence")
                & (pl.col("reactivity_unit") == "vdjdb.score")
            )
            .then(pl.col("reactivity_value"))
            .otherwise(None)
            .alias("vdjdb_score")
        )
    df = df.with_columns(
        (
            pl.col("cdr3").is_not_null()
            & pl.col("cdr3").cast(pl.String).str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
        )
        .fill_null(False)
        .alias("reference_valid")
    )
    df = df.with_columns(
        pl.col("mhc_class").replace({"MHC-I": "MHCI", "MHC-II": "MHCII"})
    )
    # Declared pair links outrank inferred observation links. Record IDs may be
    # independent chain identities, particularly after VdjdbIndex expansion.
    identity = (
        "record_id" if "record_id" in df.columns
        else "cell_id" if "cell_id" in df.columns else None
    )
    complex_key = pl.col("complex_id").cast(pl.String)
    declared = complex_key.is_not_null() & ~complex_key.is_in(["", "0"])
    link = pl.when(declared).then(complex_key).otherwise(pl.lit("0"))
    if identity:
        if (
            df[identity].null_count()
            or df.filter(pl.col(identity).cast(pl.String) == "").height
        ):
            raise ValueError(f"reference has missing {identity}")
        if df.select(identity, "gene").n_unique() != df.height:
            raise ValueError(f"reference has duplicate ({identity}, gene) identities")
        known = pl.when(declared).then(complex_key).otherwise(None).drop_nulls()
        declared_count = known.n_unique().over(identity)
        if df.filter(declared_count > 1).height:
            raise ValueError(f"reference has conflicting complex links within {identity}")
        paired = pl.col("gene").eq("TRA").any().over(identity) & pl.col("gene").eq(
            "TRB"
        ).any().over(identity)
        missing_pair = ~declared.fill_null(False) & paired
        # Zero is a valid source observation ID, but the complex field reserves
        # it for unpaired rows. Keep the source ID and namespace its pair link.
        observation_key = pl.col(identity).cast(pl.String)
        inferred_key = pl.when(observation_key == "0").then(
            pl.lit(f"__{identity}:0")
        ).otherwise(observation_key)
        # Inferred record/cell keys must not silently merge another declared pair.
        unresolved = df.filter(missing_pair & (declared_count == 0)).select(
            inferred_key.alias("key"), observation_key.alias("observation")
        ).unique()
        if unresolved.group_by("key").agg(pl.col("observation").n_unique()).filter(
            pl.col("observation") > 1
        ).height:
            raise ValueError("inferred complex identities collide between observations")
        unresolved = unresolved.select("key").unique()
        existing = df.filter(declared).select(complex_key.alias("key")).unique()
        if unresolved.join(existing, on="key", how="inner").height:
            raise ValueError("inferred complex identity collides with a declared link")
        link = pl.when(declared).then(complex_key).otherwise(
            pl.when(paired)
            .then(pl.coalesce(known.first().over(identity), inferred_key))
            .otherwise(pl.lit("0"))
        )
        if "record_id" not in df.columns:
            df = df.with_columns(pl.col(identity).alias("record_id"))
    links = df.select(link.alias("complex_id"))["complex_id"]
    # Legacy integer IDs stay integers when conversion is lossless. Meaningful
    # string IDs (including e.g. "01") are preserved rather than coerced to zero.
    integers = links.cast(pl.Int64, strict=False)
    if integers.null_count() == 0 and integers.cast(pl.String).equals(links):
        links = integers
    df = df.with_columns(links)
    if df.filter(
        pl.col("vdjdb_score").is_not_null()
        & (pl.col("vdjdb_score").cast(pl.String) != "")
        & pl.col("vdjdb_score").cast(pl.Int64, strict=False).is_null()
    ).height:
        raise ValueError("reference has invalid VDJdb confidence score")
    return df.with_columns(
        pl.col("vdjdb_score").cast(pl.Int64, strict=False).fill_null(0)
    )
