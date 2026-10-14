"""Population-relative ordering of fixed-radius linked-pair evidence, without calibration."""
from __future__ import annotations

import polars as pl

from .candidates import PMHC

GROUP = ["cohort", "species", *PMHC]
SCORES = ["paired_score", "fallback_score", "alpha_score", "beta_score"]
COUNTS = ["n_reference_pair", "n_reference_alpha", "n_reference_beta"]
OUTPUTS = ["evidence_block", "evidence_status", "evidence_score", "within_block_rank"]


def paired_evidence_order(frame: pl.DataFrame) -> pl.DataFrame:
    """Rank a complete cohort separately for each species/full pMHC.

    Higher scores are better. Block 1 uses fixed-radius paired rank evidence;
    block 2 sums average ranks of original alpha/beta scores *within block 2*;
    block 3 uses the original paired score. Null counts mean unavailable support.
    Missing pairs/chains enter block 3 via paired_available=False. Null scores
    tie below finite scores, including observed chains with no reference hits.
    Scores are population-relative, not comparable across batches, and carry no
    P-value, E-value, posterior or biological confidence interpretation.
    """
    required = ["query_id", *GROUP, *SCORES, *COUNTS, "paired_available"]
    if missing := set(required) - set(frame.columns):
        raise ValueError("missing columns: " + ", ".join(sorted(missing)))
    if set(OUTPUTS) & set(frame.columns):
        raise ValueError("input already contains paired-evidence output columns")
    for name in ["query_id", "cohort", "species", "epitope"]:
        if frame[name].null_count() or frame[name].cast(pl.String).str.strip_chars().eq("").any():
            raise ValueError(f"{name} must be non-null and nonempty")
    if frame.select(pl.struct([*GROUP, "query_id"]).is_duplicated().any()).item():
        raise ValueError("query_id must be unique within cohort/species/full-pMHC")
    if frame.schema["paired_available"] != pl.Boolean or frame["paired_available"].null_count():
        raise ValueError("paired_available must be a non-null Boolean")
    for name in SCORES + COUNTS:
        dtype = frame.schema[name]
        if dtype != pl.Null and not dtype.is_numeric():
            raise ValueError(f"{name} must be numeric")
        values = frame[name].drop_nulls()
        if not values.is_finite().all():
            raise ValueError(f"{name} must be finite; use null for unavailable evidence")
        if name in COUNTS and len(values) and (values.min() < 0 or values.max() > 2**32 - 1 or (values != values.floor()).any()):
            raise ValueError(f"{name} must contain integer counts in UInt32 range")
    # Work separately so every original column (including its dtype) is preserved.
    work = frame.select(required).with_columns(
        pl.col(SCORES).cast(pl.Float64), pl.col(COUNTS).fill_null(0),
    )
    available = pl.col("paired_available")
    common = available & (pl.col("n_reference_pair") > 0)
    if work.filter(common & pl.col("paired_score").is_null()).height:
        raise ValueError("available common-pair support requires paired_score")
    work = work.with_columns(
        pl.when(common).then(1)
        .when(available & ((pl.col("n_reference_alpha") > 0) | (pl.col("n_reference_beta") > 0))).then(2)
        .otherwise(3).cast(pl.UInt8).alias("evidence_block")
    )
    groups = [*GROUP, "evidence_block"]
    work = work.with_columns(
        pl.when(pl.col("evidence_block") == 1).then(pl.col("paired_score"))
        .when(pl.col("evidence_block") == 2).then(
            pl.col("alpha_score").fill_null(float("-inf")).rank("average").over(groups)
            + pl.col("beta_score").fill_null(float("-inf")).rank("average").over(groups)
        ).otherwise(pl.col("fallback_score").fill_null(float("-inf"))).alias("_statistic")
    ).with_columns(pl.col("_statistic").rank("average").over(groups).alias("within_block_rank"))
    work = work.with_columns(
        (3 - pl.col("evidence_block") + pl.col("within_block_rank") / (pl.len().over(groups) + 1)).alias("evidence_score"),
        pl.when(pl.col("evidence_block") == 1).then(pl.lit("common_reference_pair"))
        .when(pl.col("evidence_block") == 2).then(pl.lit("single_chain_only"))
        .when(~available).then(pl.lit("paired_evidence_unavailable"))
        .otherwise(pl.lit("no_reference_neighbour")).alias("evidence_status"),
    )
    return frame.with_columns(pl.col(COUNTS).cast(pl.UInt32)).hstack(work.select(OUTPUTS))
