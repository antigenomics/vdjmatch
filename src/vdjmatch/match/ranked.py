"""Full-reference gap-block top-K evidence; global rank tests are separate from labels.

Search keys are unique junctions (or genuine junction pairs). Observation metadata
is expanded after the native batch. Candidate weights are descriptive and do not
inherit fixed-ball epitope calibration.
"""
from __future__ import annotations

import math

import polars as pl
from seqtree import gapblock

from ..aggregate.candidates import PMHC, candidates as reduce_candidates
from ..evalue.paired import build_paired_ref
from ..evalue.ranked import order_statistic, paired_order_statistic

ESTIMATOR = "gapblock-topk-v1"


def _refs(control):
    if control is None:
        return None
    refs = control.ref_seqs()
    if not refs:
        raise ValueError("control index is empty")
    if len(set(refs)) != len(refs):
        raise ValueError("ranked controls must contain unique junctions")
    return refs


def _global(qid, k, n, scores, calibrated, result=None):
    complete = len(scores) == k
    return {
        "query_id": qid, "requested_k": k, "n_retained": len(scores),
        "n_reference": n, "threshold": max(scores) if complete else None,
        "estimator": ESTIMATOR,
        "calibration": "order-statistic" if calibrated else "uncalibrated",
        "p_global": 1.0 if not complete else None,
        "p_upper": None, "E_raw": None, "E": None,
        "finite_control_delta": None,
        **(result or {}),
    }


def ranked_evidence(annotator, q, *, k=10, threads=1, matrix=None,
                    gap_positions=(3, 4, -4, -3), control=None,
                    exclude_exact=False, paired=False):
    """Return ``(observation_hits, pMHC_candidates, global_statistics)``.

    K is fixed before searching. All control ties at the Kth score contribute to
    its global test. Fewer than K eligible target keys have no finite threshold
    and yield p_global=1. Pair radius is max(chain scores), with independent-chain
    controls; candidate weights use exp(-(scoreA+scoreB)/temperature)/(CDF_A*CDF_B+1).
    Temperature is 1 for unit costs and 400 for matrix penalties. Matrix=None
    selects unit costs; this does not load a scoring matrix implicitly.
    Missing/invalid query rows remain the responsibility of the caller's calls table.
    """
    from ..api import _prepare

    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise ValueError("k must be a positive integer")
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 0:
        raise ValueError("threads must be a nonnegative integer")
    kwargs = dict(matrix=matrix, gap_open=2 * matrix.scale() if matrix is not None else 2, gap_extend=1,
                  gap_prior=gapblock.positions_prior(gap_positions), threads=threads)
    if paired:
        return _paired(annotator, q, k, control, exclude_exact, kwargs)
    _, prepared = _prepare(q)
    hits, candidates, statistics = [], [], []
    for locus in sorted(prepared["locus"].drop_nulls().unique()):
        if locus not in annotator._index.genes:
            continue
        if isinstance(control, dict) and control.get(locus) is None:
            raise ValueError("ranked calibration requires a control for every active locus")
        queries = prepared.filter(pl.col("locus") == locus)
        records = annotator._index.records_for(locus).filter(pl.col("reference_valid"))
        refs = records["cdr3"].unique().sort().to_list()
        ctrl = _refs(control.get(locus) if isinstance(control, dict) else control)
        seqs = queries["cdr3"].to_list()
        result = gapblock.topk_batch(seqs, refs, k=k, exclude_exact=exclude_exact, **kwargs)
        result = [[h for h in row if h.score < gapblock.UNREACHABLE] for row in result]
        thresholds = [[h.score for h in row] for row in result]
        counts = gapblock.count_batch(seqs, ctrl, thresholds, exclude_exact=exclude_exact, **kwargs) if ctrl is not None else [[0] * len(row) for row in result]
        rows, cdf = [], []
        for query, row, nc in zip(queries.iter_rows(named=True), result, counts):
            scores = [h.score for h in row]
            stat = order_statistic(k, len(refs), nc[-1], len(ctrl)) if ctrl is not None and len(row) == k else None
            statistics.append(_global(query["query_id"], k, len(refs), scores, ctrl is not None, stat))
            for hit, count in zip(row, nc):
                rows.append({"query_id": query["query_id"], "query_cdr3": query["cdr3"],
                             "query_v": query["v"], "query_j": query["j"], "query_locus": locus,
                             "db_cdr3": refs[hit.ref_id], "score": hit.score})
                cdf.append((query["query_id"], hit.score, count))
        if not rows:
            continue
        detail = pl.DataFrame(rows, schema_overrides={"query_id": prepared.schema["query_id"], "score": pl.Int32, "query_v": pl.String, "query_j": pl.String}).join(
            records.rename({"cdr3": "db_cdr3", "v": "db_v", "j": "db_j"}), on="db_cdr3", how="inner")
        detail = detail.with_columns(*(pl.lit(None, dtype=pl.UInt16).alias(c) for c in ("n_subs", "n_ins", "n_dels")))
        cframe = pl.DataFrame(cdf, orient="row", schema={"query_id": prepared.schema["query_id"], "score": pl.Int32, "_nc": pl.UInt64}).unique()
        reduced = reduce_candidates(detail, records, control_counts=cframe, distance_column="score", soft_v=False,
                                    score_scale=400.0 if matrix is not None else 1.0)
        reduced = reduced.with_columns(pl.lit(ESTIMATOR).alias("estimator"), pl.lit("uncalibrated").alias("calibration"), pl.lit("gapblock_score").alias("distance_name"))
        hits.append(detail)
        candidates.append(reduced)
    return (pl.concat(hits, how="diagonal_relaxed") if hits else annotator._index.empty_hits(prepared.schema["query_id"]),
            pl.concat(candidates, how="diagonal_relaxed") if candidates else _empty_candidates(prepared.schema["query_id"]),
            _statistics(statistics, prepared.schema["query_id"]))


def _empty_candidates(dtype):
    from ..aggregate.candidates import candidate_schema
    return pl.DataFrame(schema={**candidate_schema(dtype), "distance_name": pl.String})


def _statistics(rows, dtype):
    schema = {"query_id": dtype, "requested_k": pl.UInt32, "n_retained": pl.UInt32,
              "n_reference": pl.UInt64, "threshold": pl.Int32, "estimator": pl.String,
              "calibration": pl.String, **{c: pl.Float64 for c in ("p_global", "p_upper", "E_raw", "E", "finite_control_delta")}}
    return pl.DataFrame(rows, schema=schema, strict=False)


def _paired(annotator, q, k, control, exclude_exact, kwargs):
    from ..api import _prepare
    if not isinstance(q, pl.DataFrame) or not {"cdr3a", "cdr3b"} <= set(q.columns):
        raise ValueError("paired ranked queries require cdr3a and cdr3b")
    _, qa = _prepare(q, "cdr3a", locus="TRA")
    _, qb = _prepare(q, "cdr3b", locus="TRB")
    shared = qa.select("query_id").join(qb.select("query_id"), on="query_id")
    qa = shared.join(qa, on="query_id", maintain_order="left")
    qb = shared.join(qb, on="query_id", maintain_order="left")
    observations = build_paired_ref(annotator._index.reference)
    keys = observations.select("alpha", "beta").unique().sort("alpha", "beta")
    ra, rb = keys["alpha"].to_list(), keys["beta"].to_list()
    if control is not None and (not isinstance(control, dict) or any(control.get(c) is None for c in ("TRA", "TRB"))):
        raise ValueError("paired controls require both TRA and TRB")
    ca, cb = (_refs(control[c]) if control is not None else None for c in ("TRA", "TRB"))
    sa, sb = qa["cdr3"].to_list(), qb["cdr3"].to_list()
    result = gapblock.paired_topk_batch(sa, sb, ra, rb, k=k, exclude_exact=exclude_exact, **kwargs)
    result = [[r for r in row if max(r[1], r[2]) < gapblock.UNREACHABLE] for row in result]
    ta, tb = [[r[1] for r in row] for row in result], [[r[2] for r in row] for row in result]
    # Include the max-chain threshold for the global rectangular CDF, separately
    # from per-hit chain thresholds used by candidate ranking.
    gta = [row + ([max(max(r[1], r[2]) for r in hits)] if hits else []) for row, hits in zip(ta, result)]
    gtb = [row + ([max(max(r[1], r[2]) for r in hits)] if hits else []) for row, hits in zip(tb, result)]
    na = gapblock.count_batch(sa, ca, gta, **kwargs) if ca is not None else [[0] * len(row) for row in gta]
    nb = gapblock.count_batch(sb, cb, gtb, **kwargs) if cb is not None else [[0] * len(row) for row in gtb]
    aset, bset = set(ca or []), set(cb or [])
    rows, stats = [], []
    for i, (query, row) in enumerate(zip(qa.iter_rows(named=True), result)):
        za, zb = int(sa[i] in aset), int(sb[i] in bset)
        stat = paired_order_statistic(k, len(ra), na[i][-1], nb[i][-1], za, zb, len(ca), len(cb), exclude_exact) if ca is not None and len(row) == k else None
        stats.append(_global(query["query_id"], k, len(ra), [max(a, b) for _, a, b in row], ca is not None, stat))
        for j, (rid, a, b) in enumerate(row):
            nc = na[i][j] * nb[i][j] - (za * zb if exclude_exact else 0)
            rows.append({"query_id": query["query_id"], "query_cdr3_alpha": sa[i], "query_cdr3_beta": sb[i],
                         "alpha": ra[rid], "beta": rb[rid], "score_alpha": a, "score_beta": b,
                         "score": max(a, b), "_weight": math.exp(-(a+b)/(400.0 if kwargs["matrix"] is not None else 1.0)) / (nc+1),
                         "_exact": sa[i] == ra[rid] and sb[i] == rb[rid]})
    unique = pl.DataFrame(rows, schema={"query_id": qa.schema["query_id"],
        **{c: pl.String for c in ("query_cdr3_alpha", "query_cdr3_beta", "alpha", "beta")},
        **{c: pl.Int32 for c in ("score_alpha", "score_beta", "score")},
        "_weight": pl.Float64, "_exact": pl.Boolean})
    detail = unique.join(observations, on=["alpha", "beta"], how="inner")
    raw = annotator._index.reference
    if "record_id" not in raw.columns:
        raw = raw.with_row_index("record_id").with_columns(pl.col("record_id").cast(pl.String))
    raw = raw.join(observations.select("complex_id"), on="complex_id", how="semi")
    for chain, prefix in (("TRA", "alpha_"), ("TRB", "beta_")):
        metadata = raw.filter(pl.col("gene") == chain).unique(subset=["complex_id", "cdr3", "v", "j"]).select(
            "complex_id", *(pl.col(c).alias(prefix+c) for c in raw.columns if c != "complex_id"))
        detail = detail.join(metadata, on="complex_id", how="left", validate="m:1")
    if not rows:
        return detail.drop("_weight", "_exact"), _empty_candidates(qa.schema["query_id"]), _statistics(stats, qa.schema["query_id"])
    groups = ["query_id", *[c for c in PMHC if c in detail.columns]]
    votes = detail.unique(subset=[*groups, "alpha", "beta"])
    out = votes.group_by(groups).agg(pl.len().alias("n_clonotypes"), pl.col("_weight").sum().alias("ned_score"), pl.col("score").min().alias("nearest_score"), pl.col("_exact").sum().alias("n_exact"))
    support = detail.group_by(groups).agg(
        pl.col("complex_id").n_unique().alias("n_records"),
        pl.max_horizontal("alpha_vdjdb_score", "beta_vdjdb_score").max().alias("db_score"))
    ref = annotator._index.reference
    if "reference_id" in ref.columns:
        study = ref.select("complex_id", "reference_id").unique()
        support = support.join(detail.join(study, on="complex_id").group_by(groups).agg(pl.col("reference_id").filter(pl.col("reference_id").is_not_null() & (pl.col("reference_id") != "")).n_unique().alias("n_studies")), on=groups, nulls_equal=True)
    else:
        support = support.with_columns(pl.lit(0).alias("n_studies"))
    out = out.join(support, on=groups, nulls_equal=True).sort(["query_id", "ned_score", *groups[1:]], descending=[False, True, *[False]*(len(groups)-1)], nulls_last=True)
    coverage = observations.unique(subset=[*groups[1:], "alpha", "beta"]).group_by(groups[1:]).len().rename({"len": "n_reference"})
    assignments = votes.group_by("query_id", "alpha", "beta").len().rename({"len": "_labels"})
    total = assignments.group_by("query_id").len().rename({"len": "_total"})
    exclusive = votes.join(assignments, on=["query_id", "alpha", "beta"]).group_by(groups).agg((pl.col("_labels") == 1).sum().alias("_exclusive"))
    out = out.join(coverage, on=groups[1:], nulls_equal=True).join(exclusive, on=groups, nulls_equal=True).join(total, on="query_id").with_columns(
        (pl.col("_total")-pl.col("_exclusive")).alias("n_competing_clonotypes"), pl.col("n_clonotypes").alias("n_target")).drop("_total", "_exclusive")
    out = out.sort(["query_id", "ned_score", *groups[1:]], descending=[False, True, *[False]*(len(groups)-1)], nulls_last=True)
    out = out.with_columns(pl.int_range(1, pl.len()+1).over("query_id").alias("rank"))
    best = out.group_by("query_id").agg(pl.col("ned_score").first().alias("_best"),
        pl.col("ned_score").get(1, null_on_oob=True).fill_null(0.0).alias("_second"))
    out = out.join(best, on="query_id").with_columns(
        pl.when(pl.col("rank") == 1).then(pl.col("_second")).otherwise(pl.col("_best")).alias("competing_score"),
        pl.col("n_records").sum().over("query_id").alias("n_hits"),
        pl.lit(False).alias("rule_of_three")).drop("_best", "_second")
    out = out.with_columns((pl.col("ned_score")-pl.col("competing_score")).alias("score_margin"),
                          pl.lit(ESTIMATOR).alias("estimator"), pl.lit("uncalibrated").alias("calibration"), pl.lit("gapblock_score").alias("distance_name"),
                          *(pl.lit(None, dtype=pl.Float64).alias(c) for c in ("E", "p_enrichment", "p_any")),
                          pl.lit(None, dtype=pl.UInt32).alias("nearest_edits"), pl.lit(None, dtype=pl.UInt64).alias("distance_sum"))
    return detail.drop("_weight", "_exact"), out.sort("query_id", "rank"), _statistics(stats, qa.schema["query_id"])
