"""Per-pMHC evidence; NED ranks candidates, enrichment tests the stated reference set.

NED-v1 is the size-invariant +1-control neighbour-density formula in the existing
holdout scorer. It is a ranking statistic, not a posterior probability of specificity.
"""

from __future__ import annotations

from difflib import SequenceMatcher
import math

import polars as pl
from seqtree.evalue import evalue_result

from ..match.engine import _genefam
from ..match.vgene import load_v_cdr12

PMHC = ["epitope", "mhc_a", "mhc_b", "mhc_class"]
CLONE = ["db_cdr3", "db_v", "db_j"]
ESTIMATOR = "ned-v1"


def candidate_schema(query_id_type=pl.UInt32):
    return {
        "query_id": query_id_type,
        **{
            c: pl.String
            for c in PMHC
            + [
                "query_cdr3",
                "query_v",
                "query_j",
                "query_locus",
                "estimator",
                "calibration",
            ]
        },
        **{
            c: pl.UInt32
            for c in [
                "rank",
                "n_target",
                "n_records",
                "n_clonotypes",
                "n_studies",
                "n_competing_clonotypes",
                "n_reference",
                "n_control",
                "nearest_edits",
                "n_exact",
                "n_hits",
            ]
        },
        "distance_sum": pl.UInt64,
        "nearest_score": pl.Int32,
        "db_score": pl.Int64,
        **{
            c: pl.Float64
            for c in [
                "ned_score",
                "competing_score",
                "score_margin",
                "E",
                "p_enrichment",
                "p_any",
            ]
        },
        "rule_of_three": pl.Boolean,
    }


def candidates(
    hits,
    reference,
    controls=None,
    *,
    control_size=None,
    score_scale=400.0,
    soft_v=True,
    match_v=False,
    match_j=False,
    control_counts=None,
    distance_column="edits",
):
    """Reduce independent observations to candidate evidence for each stable query_id.

    ``controls``: hit table query_id/edits, one row per unique control sequence match.
    Counts use unique (junction,V,J) keys; study support excludes empty identifiers.
    ``control_size=None`` requests explicitly uncalibrated ranking. Soft V uses existing
    germline-loop similarity (weight 0.25 across families); unknown/missing calls are neutral.
    """
    if distance_column not in {"edits", "score"}:
        raise ValueError("distance_column must be edits or score")
    if controls is not None and control_counts is not None:
        raise ValueError("supply control hits or cumulative control counts, not both")
    if not math.isfinite(score_scale) or score_scale <= 0:
        raise ValueError("score_scale must be positive")
    if control_size is not None and control_size <= 0:
        raise ValueError("control index is empty")
    dtype = hits.schema.get("query_id", pl.UInt32)
    if hits.height == 0:
        return pl.DataFrame(schema=candidate_schema(dtype))
    h = hits.with_columns(
        (pl.col("n_subs").cast(pl.UInt32) + pl.col("n_ins") + pl.col("n_dels")).alias(
            "edits"
        ),
        _genefam(pl.col("query_v")).alias("_qv"),
        _genefam(pl.col("db_v")).alias("_dv"),
    )
    groups = ["query_id", *PMHC]
    h = h.with_columns(
        pl.col("reference_id").cast(pl.String),
        pl.col("vdjdb_score").cast(pl.Int64, strict=False).fill_null(0),
    )
    obs = h.group_by(groups).agg(
        pl.len().alias("n_records"),
        pl.col("reference_id")
        .filter(pl.col("reference_id").is_not_null() & (pl.col("reference_id") != ""))
        .n_unique()
        .alias("n_studies"),
        pl.col("vdjdb_score").max().alias("db_score"),
    )
    # Search-key evidence is independent of observation multiplicity.
    u = h.sort("score").unique(subset=[*groups, *CLONE], maintain_order=True)
    if control_counts is not None:
        u = u.join(
            control_counts.select("query_id", distance_column, "_nc"),
            on=["query_id", distance_column], how="left", validate="m:1",
        )
        if u["_nc"].null_count() or u.filter(pl.col("_nc") < 0).height:
            raise ValueError("cumulative control counts must cover every target score and be nonnegative")
    elif controls is not None and controls.height:
        hist = (
            controls.group_by("query_id", distance_column)
            .len()
            .sort(["query_id", distance_column])
            .with_columns(pl.col("len").cum_sum().over("query_id").alias("_nc"))
            .drop("len")
        )
        u = u.sort(["query_id", distance_column]).join_asof(
            hist, on=distance_column, by="query_id", strategy="backward", check_sortedness=False
        )
    else:
        u = u.with_columns(pl.lit(0, dtype=pl.UInt32).alias("_nc"))
    u = u.with_columns(pl.col("_nc").fill_null(0))
    if soft_v:
        pairs = u.select("_qv", "_dv").drop_nulls().unique()
        loops = load_v_cdr12()
        weights = []
        for a, b in pairs.iter_rows():
            weights.append(
                1.0
                if a == b
                else 0.25 * SequenceMatcher(None, loops[a], loops[b]).ratio()
                if a in loops and b in loops
                else 1.0
            )
        pairs = pairs.with_columns(pl.Series("_vweight", weights, dtype=pl.Float64))
        u = u.join(pairs, on=["_qv", "_dv"], how="left", validate="m:1").with_columns(
            pl.col("_vweight").fill_null(1.0)
        )
    else:
        u = u.with_columns(pl.lit(1.0).alias("_vweight"))
    u = u.with_columns(
        (
            pl.col("_vweight")
            * (-pl.col("score") / score_scale).exp()
            / (pl.col("_nc") + 1)
        ).alias("_w")
    )
    out = u.group_by(groups).agg(
        pl.len().alias("n_clonotypes"),
        pl.col("db_cdr3").n_unique().alias("n_target"),
        pl.col("edits").min().alias("nearest_edits"),
        pl.col("edits").cast(pl.UInt64).sum().alias("distance_sum"),
        pl.col("score").min().alias("nearest_score"),
        (pl.col("query_cdr3") == pl.col("db_cdr3")).sum().cast(pl.UInt32).alias("n_exact"),
        pl.col("_w").sum().alias("ned_score"),
        *[
            pl.col(c).first()
            for c in ["query_cdr3", "query_v", "query_j", "query_locus"]
        ],
    )
    if distance_column == "score":
        out = out.with_columns(
            pl.lit(None, dtype=pl.UInt32).alias("nearest_edits"),
            pl.lit(None, dtype=pl.UInt64).alias("distance_sum"),
        )
    # Competitors may share a clone with this candidate. Subtract only exclusive clones.
    assignments = u.group_by(["query_id", *CLONE]).len().rename({"len": "_labels"})
    total = assignments.group_by("query_id").len().rename({"len": "_total"})
    exclusive = (
        u.join(assignments, on=["query_id", *CLONE], nulls_equal=True, validate="m:1")
        .group_by(groups)
        .agg((pl.col("_labels") == 1).sum().alias("_exclusive"))
    )
    out = (
        out.join(obs, on=groups, nulls_equal=True, validate="1:1")
        .join(exclusive, on=groups, nulls_equal=True)
        .join(total, on="query_id")
    )
    out = out.with_columns(
        (pl.col("_total") - pl.col("_exclusive"))
        .cast(pl.UInt32)
        .alias("n_competing_clonotypes")
    )
    # Reference sizes obey the same V/J restrictions as target hits.
    r = reference.rename({"cdr3": "db_cdr3", "v": "db_v", "j": "db_j"})
    covkeys = PMHC.copy()
    joinkeys = PMHC.copy()
    if match_v:
        r = r.with_columns(_genefam(pl.col("db_v")).alias("_qv"))
        out = out.with_columns(_genefam(pl.col("query_v")).alias("_qv"))
        covkeys.append("_qv")
        joinkeys.append("_qv")
    if match_j:
        r = r.with_columns(_genefam(pl.col("db_j")).alias("_qj"))
        out = out.with_columns(_genefam(pl.col("query_j")).alias("_qj"))
        covkeys.append("_qj")
        joinkeys.append("_qj")
    coverage = (
        r.unique(subset=[*covkeys, "db_cdr3"])
        .group_by(covkeys)
        .len()
        .rename({"len": "n_reference"})
    )
    out = out.join(coverage, on=joinkeys, how="left", nulls_equal=True, validate="m:1")
    if controls is not None:
        counts = controls.group_by("query_id").len().rename({"len": "n_control"})
        out = out.join(counts, on="query_id", how="left").with_columns(
            pl.col("n_control").fill_null(0)
        )
    else:
        out = out.with_columns(pl.lit(0, dtype=pl.UInt32).alias("n_control"))
    # Public upstream statistics; no target re-search. One reduction per candidate.
    if control_size is not None:
        stats = [
            evalue_result(nt, nc, nr, control_size)
            for nt, nc, nr in out.select(
                "n_target", "n_control", "n_reference"
            ).iter_rows()
        ]
        out = out.with_columns(
            *[
                pl.Series(c, [s[c] for s in stats])
                for c in ["E", "p_enrichment", "p_any", "rule_of_three"]
            ]
        )
    else:
        out = out.with_columns(
            *[
                pl.lit(None, dtype=pl.Float64).alias(c)
                for c in ["E", "p_enrichment", "p_any"]
            ],
            pl.lit(False).alias("rule_of_three"),
        )
    out = out.sort(
        ["query_id", "ned_score", "nearest_score", *PMHC],
        descending=[False, True, False, *([False] * len(PMHC))],
        nulls_last=True,
    )
    out = out.with_columns(
        pl.int_range(1, pl.len() + 1).over("query_id").cast(pl.UInt32).alias("rank"),
        pl.col("n_records").sum().over("query_id").cast(pl.UInt32).alias("n_hits"),
        pl.lit(ESTIMATOR).alias("estimator"),
        pl.lit("fixed_ball" if control_size is not None else "uncalibrated").alias(
            "calibration"
        ),
    )
    # Best alternative to each candidate; margin can be negative for non-top candidates.
    best = out.group_by("query_id").agg(
        pl.col("ned_score").first().alias("_best"),
        pl.col("ned_score").get(1, null_on_oob=True).fill_null(0.0).alias("_second"),
    )
    out = out.join(best, on="query_id").with_columns(
        pl.when(pl.col("rank") == 1)
        .then(pl.col("_second"))
        .otherwise(pl.col("_best"))
        .alias("competing_score")
    )
    out = out.with_columns(
        (pl.col("ned_score") - pl.col("competing_score")).alias("score_margin")
    )
    return out.select([pl.col(c).cast(t) for c, t in candidate_schema(dtype).items()])
