"""Fixed geometry and empirical background-mass ranking across length routes."""
from __future__ import annotations

import numpy as np
import polars as pl
from seqtree import SubstitutionMatrix, gapblock

from ..aggregate.candidates import candidate_schema, candidates
from . import search_params


def unified_evidence(ann, q, *, threads=1, control=None, species=None, exclude_exact=False, distance="edit"):
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be a positive integer")
    if distance not in {"edit", "gapblock"}:
        raise ValueError("distance must be edit or gapblock")
    active = sorted(set(q["locus"].drop_nulls()) & set(ann.loci))
    controls = {}
    for locus in active:
        ctrl = control.get(locus) if isinstance(control, dict) else control
        if isinstance(control, dict) and ctrl is None:
            raise ValueError("unified scoring requires a control for every active locus")
        if ctrl is None:
            from ..evalue.control import background
            ctrl = background(locus, ann._control_species([locus], species))
        refs_control = ctrl.ref_seqs()
        if not refs_control or len(set(refs_control)) != len(refs_control):
            raise ValueError("unified controls must contain nonempty unique junctions")
        controls[locus] = ctrl
    control = controls
    short = q if distance == "edit" else q.head(0)
    long = q if distance == "gapblock" else q.head(0)
    hs, cs = [], []
    if short.height:
        h, c = ann._evidence(short, search_params("5,2,2,5"), threads=threads,
                             control=control, calibrate=True, species=species,
                             score_scale=1.0, soft_v=False, exclude_exact=exclude_exact,
                             background_mass=True)
        hs.append(h.with_columns(pl.lit("edit").alias("distance_name")))
        cs.append(c.with_columns(pl.lit("edit").alias("distance_name")))
    matrix = SubstitutionMatrix.blosum62()
    scale = matrix.scale()
    cutoff = 5 * scale
    kwargs = dict(matrix=matrix, gap_open=2 * scale, gap_extend=scale,
                  gap_prior=gapblock.positions_prior((3, 4, -4, -3)), threads=threads)
    for locus in sorted(set(long["locus"].drop_nulls()) & set(ann.loci)):
        queries = long.filter(pl.col("locus") == locus)
        refs_control = controls[locus].ref_seqs()
        records = ann._index.records_for(locus)
        refs = records.select("cdr3").unique().sort("cdr3").with_row_index("_rid")
        refseqs = refs["cdr3"].to_list()
        unique_queries = queries.select("cdr3").unique().sort("cdr3")
        # Bound the native all-versus-all int32 matrix at 64 MiB. Reduce immediately.
        batch_size = max(1, (64 * 1024**2) // (4 * len(refseqs)))
        for batch in unique_queries.iter_slices(batch_size):
            seqs = batch["cdr3"].to_list()
            scores = np.asarray(gapblock.score_matrix(seqs, refseqs, **kwargs))
            mask = scores <= cutoff
            if exclude_exact:
                mask &= np.asarray(seqs)[:, None] != np.asarray(refseqs)[None, :]
            qi, ri = np.nonzero(mask)
            thresholds = [sorted(set(scores[i, mask[i]].tolist()) | {cutoff}) for i in range(len(seqs))]
            counts = gapblock.count_batch(seqs, refs_control, thresholds,
                                         exclude_exact=exclude_exact, **kwargs)
            detail = pl.DataFrame({"query_cdr3": np.asarray(seqs)[qi], "_rid": ri,
                                   "score": scores[qi, ri]}).join(refs, on="_rid").drop("_rid").rename({"cdr3": "db_cdr3"})
            del scores, mask
            detail = detail.join(queries.select("query_id", "count",
                pl.col("cdr3").alias("query_cdr3"), pl.col("v").alias("query_v"),
                pl.col("j").alias("query_j"), pl.col("locus").alias("query_locus"),
                *[pl.col(c).alias("query_"+c) for c in ("sequence_id", "pair_id") if c in queries.columns]), on="query_cdr3")
            detail = detail.with_columns(*(pl.lit(None, dtype=pl.UInt16).alias(c) for c in ("n_subs", "n_ins", "n_dels")))
            metadata = records.rename({"cdr3": "db_cdr3", "v": "db_v", "j": "db_j"})
            shared = (set(detail.columns) & set(metadata.columns)) - {"db_cdr3"}
            metadata = metadata.rename({c: "db_metadata_"+c for c in shared})
            detail = detail.join(metadata, on="db_cdr3")
            detail = detail.with_columns(pl.col("score").cast(pl.Int32),
                *(pl.lit(None, dtype=pl.UInt16).alias(c) for c in ("n_subs", "n_ins", "n_dels")))
            cdf = pl.DataFrame([(s, d, n) for s, ds, ns in zip(seqs, thresholds, counts) for d, n in zip(ds, ns)],
                schema={"cdr3": pl.String, "score": pl.Int32, "_nc": pl.UInt64}, orient="row").join(queries.select("query_id", "cdr3"), on="cdr3")
            totals = cdf.filter(pl.col("score") == cutoff).select("query_id", pl.col("_nc").alias("n_control"))
            c = candidates(detail, records, control_counts=cdf, control_totals=totals,
                           control_size=len(refs_control), distance_column="score", soft_v=False,
                           background_mass=True)
            hs.append(detail.with_columns(pl.lit("blosum62_gapblock").alias("distance_name")))
            cs.append(c.with_columns(pl.lit("blosum62_gapblock").alias("distance_name")))
    return (pl.concat(hs, how="diagonal_relaxed") if hs else ann._index.empty_hits(q.schema["query_id"]),
            pl.concat(cs, how="diagonal_relaxed").sort("query_id", "rank") if cs else pl.DataFrame(schema=candidate_schema(q.schema["query_id"])))
