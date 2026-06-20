"""Aggregate per-hit annotations into epitope-level enrichment and per-query best calls."""
from __future__ import annotations

import polars as pl


def epitope_summary(hits: pl.DataFrame) -> pl.DataFrame:
    """Epitope-level summary: unique matched query clonotypes, reads, and best alignment score,
    grouped by (epitope, mhc_class, antigen_species). Sorted by unique clonotypes desc."""
    if hits.height == 0:
        return pl.DataFrame(schema={"epitope": pl.Utf8, "mhc_class": pl.Utf8,
                                    "antigen_species": pl.Utf8, "unique": pl.UInt32,
                                    "reads": pl.Int64, "best_score": pl.Int32})
    per_q = (hits.group_by("epitope", "mhc_class", "antigen_species", "query_cdr3")
                 .agg(pl.col("count").first(), pl.col("score").min().alias("best_score")))
    return (per_q.group_by("epitope", "mhc_class", "antigen_species")
                 .agg(pl.col("query_cdr3").n_unique().alias("unique"),
                      pl.col("count").sum().alias("reads"),
                      pl.col("best_score").min().alias("best_score"))
                 .sort("unique", descending=True))


def best_call(hits: pl.DataFrame, evals: pl.DataFrame | None = None) -> pl.DataFrame:
    """Per query clonotype, the predicted epitope = most-supported (most DB records, then highest
    VDJdb confidence, then best alignment). If ``evals`` (query_evalues output) is supplied, attach
    the query's E / p_enrichment as the annotation confidence."""
    if hits.height == 0:
        return hits
    support = (hits.group_by("query_cdr3", "epitope")
                   .agg(pl.len().alias("n_records"),
                        pl.col("vdjdb_score").max().alias("db_score"),
                        pl.col("score").min().alias("best_score"),
                        pl.col("mhc_class").first(), pl.col("antigen_species").first()))
    best = (support.sort(["query_cdr3", "n_records", "db_score", "best_score"],
                         descending=[False, True, True, False])
                   .group_by("query_cdr3", maintain_order=True).first())
    if evals is not None:
        best = best.join(evals.select("query_cdr3", "E", "p_enrichment", "rule_of_three"),
                         on="query_cdr3", how="left")
    return best
