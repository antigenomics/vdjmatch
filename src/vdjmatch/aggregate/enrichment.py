"""Descriptive pMHC match counts and shared candidate selection."""

from __future__ import annotations

import polars as pl
from .candidates import PMHC, candidates


def epitope_summary(hits):
    """Unique junction/V/J/locus clonotypes, query rows and their abundances, without database-record duplication."""
    groups = [*PMHC, "antigen_species"]
    if not hits.height:
        return pl.DataFrame(
            schema={
                **{c: pl.String for c in groups},
                "unique": pl.UInt32,
                "n_query_rows": pl.UInt32,
                "reads": pl.Int64,
                "best_score": pl.Int32,
            }
        )
    per = hits.group_by(*groups, "query_id").agg(
        pl.col("count").first(),
        pl.col("score").min().alias("best_score"),
        pl.struct("query_cdr3", "query_v", "query_j", "query_locus")
        .first()
        .alias("clonotype"),
    )
    return (
        per.group_by(groups)
        .agg(
            pl.col("clonotype").n_unique().alias("unique"),
            pl.len().alias("n_query_rows"),
            pl.col("count").sum().cast(pl.Int64).alias("reads"),
            pl.col("best_score").min(),
        )
        .sort("unique", descending=True)
    )


def best_call(hits, evals=None):
    """Compatibility reducer using NED-v1; optional global statistics stay separately named."""
    if not hits.height:
        return candidates(hits, pl.DataFrame())
    ref = hits.rename({"db_cdr3": "cdr3", "db_v": "v", "db_j": "j"})
    c = candidates(hits, ref).filter(pl.col("rank") == 1)
    if evals is not None:
        keys = ["query_id"] if "query_id" in evals.columns else ["query_cdr3"]
        if evals.select(keys).n_unique() != evals.height:
            raise ValueError("global E-value join keys are not unique")
        stat = evals.select(
            *keys,
            *[s for s in ["E", "p_enrichment", "rule_of_three"] if s in evals.columns],
        )
        stat = stat.rename({s: "global_" + s for s in stat.columns if s not in keys})
        c = c.join(stat, on=keys, how="left", validate="m:1")
    return c
