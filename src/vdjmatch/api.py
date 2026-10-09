"""Standalone selected-receptor annotation with explicit reference and candidate evidence."""

from __future__ import annotations

import polars as pl

from . import db
from .aggregate.candidates import (
    candidates as reduce_candidates,
    candidate_schema,
    PMHC,
)
from .io.columns import normalize_query
from .match.engine import VdjdbIndex, _empty_hits, _search_unique_junctions
from .match.scope import search_params
from .match.scoring import load_vdjam, DEFAULT_SCALE

DEFAULT_SCOPE = "1,0,0,1"
_LOCUS = {
    "TRA": "TRA",
    "TRB": "TRB",
    "A": "TRA",
    "B": "TRB",
    "ALPHA": "TRA",
    "BETA": "TRB",
    "TCRA": "TRA",
    "TCRB": "TRB",
}


def _norm_locus(x):
    return _LOCUS.get(str(x).upper())


def _prepare(data, cdr3="cdr3", v=None, j=None, locus=None, sequence_convention=None):
    if isinstance(data, str):
        data = [data]
    if not isinstance(data, pl.DataFrame):
        data = pl.DataFrame({cdr3: pl.Series(list(data), dtype=pl.String)})
    # Preserve caller data; computation IDs are positional and cannot collide with caller IDs.
    qcol = "junction_aa" if cdr3 == "cdr3" and "junction_aa" in data.columns else cdr3
    if qcol not in data.columns:
        raise ValueError(
            f'no junction column {qcol!r}; pass cdr3="junction_aa" or the producer-specific junction column'
        )
    if qcol == "cdr3_aa" and any(
        c in data.columns for c in ["sequence_id", "productive", "v_call"]
    ):
        raise ValueError("AIRR cdr3_aa excludes anchors; provide junction_aa")
    if locus is None and "locus" in data.columns:
        locus = "locus"
    if v is None and "v" in data.columns:
        v = "v"
    if j is None and "j" in data.columns:
        j = "j"
    loc_expr = pl.col(locus) if locus in data.columns else pl.lit(locus or "TRB")
    q = data.select(
        pl.col(qcol).cast(pl.String).alias("cdr3"),
        pl.col(v or "v_call").cast(pl.String).alias("v")
        if (v or "v_call") in data.columns
        else pl.lit(None, dtype=pl.String).alias("v"),
        pl.col(j or "j_call").cast(pl.String).alias("j")
        if (j or "j_call") in data.columns
        else pl.lit(None, dtype=pl.String).alias("j"),
        loc_expr.cast(pl.String)
        .str.to_uppercase()
        .replace_strict(_LOCUS, default=None)
        .alias("locus"),
        pl.col("count")
        if "count" in data.columns
        else pl.col("duplicate_count").alias("count")
        if "duplicate_count" in data.columns
        else pl.lit(1, dtype=pl.Int64).alias("count"),
    )
    if "query_id" in data.columns:
        q = q.with_columns(data["query_id"])
        if q["query_id"].null_count() or q["query_id"].n_unique() != q.height:
            raise ValueError("query_id must be unique and non-null")
    else:
        q = q.with_row_index("query_id")
    for extra in ("sequence_id", "pair_id"):
        if extra in data.columns:
            q = q.with_columns(data[extra])
    q = normalize_query(q, source="legacy", sequence_convention=sequence_convention)
    return data, q


class Annotator:
    """An explicitly owned reference index, reusable across queries without result memoization."""

    def __init__(self, index: VdjdbIndex, *, provenance=None):
        self._index = index
        self.provenance = provenance or {"source": "custom_frame"}

    @classmethod
    def latest(
        cls,
        *,
        species="HomoSapiens",
        source="github",
        asset="default",
        cache=None,
        pin=None,
    ):
        """Resolve GitHub latest/pinned release; HF benchmark mirror requires source='hf'."""
        if source not in {"github", "hf"}:
            raise ValueError("source must be github or hf")
        path = (
            db.fetch_hf(asset=asset, cache=cache, **({"tag": pin} if pin else {}))
            if source == "hf"
            else db.fetch_latest(asset=asset, cache=cache, pin=pin)
        )
        return cls.from_path(path, species=species, asset=asset)

    @classmethod
    def version(
        cls, tag, *, species="HomoSapiens", source="github", asset="default", cache=None
    ):
        return cls.latest(
            pin=tag, species=species, source=source, asset=asset, cache=cache
        )

    @classmethod
    def from_path(cls, path, *, species=None, asset="default", **filters):
        frame = db.load(path, asset=asset, species=species, **filters)
        return cls(VdjdbIndex.build(frame), provenance=db.provenance(path))

    @classmethod
    def from_frame(cls, vdj, *, species=None):
        return cls(VdjdbIndex.build(vdj, species=species))

    @property
    def reference_report(self):
        """Retained reference observations versus searchable valid junctions."""
        r = self._index.reference
        return {
            "reference_rows": r.height,
            "searchable_rows": int(r["reference_valid"].sum() or 0),
            "unsearchable_rows": r.height - int(r["reference_valid"].sum() or 0),
        }

    @property
    def loci(self):
        return self._index.genes

    def _control_species(self, loci, species):
        from .evalue.control import _organism

        reference = self._index.reference.filter(pl.col("gene").is_in(loci))
        observed = (
            reference["species"].cast(pl.String).str.strip_chars().replace("", None)
        )
        known = {_organism(s) for s in observed.drop_nulls().unique()}
        if len(known) > 1:
            raise ValueError(
                "automatic calibration requires a single reference species"
            )
        if species is not None:
            organism = _organism(species)
            if known and organism not in known:
                raise ValueError("control species does not match reference species")
            return organism
        if not known or observed.null_count():
            raise ValueError(
                "reference species is missing; supply an explicit control species"
            )
        return next(iter(known))

    def _evidence(
        self,
        q,
        sp,
        *,
        threads=1,
        match_v=False,
        match_j=False,
        align=False,
        progress=False,
        control=None,
        calibrate=False,
        species=None,
        score_scale=400.0,
        soft_v=True,
        exclude_exact=False,
        radii=None,
    ):
        if match_v and q.height and q["v"].null_count():
            raise ValueError("match_v requires a V call on every retained query")
        if match_j and q.height and q["j"].null_count():
            raise ValueError("match_j requires a J call on every retained query")
        import math

        if not math.isfinite(score_scale) or score_scale <= 0:
            raise ValueError("score_scale must be finite and positive")
        if (calibrate or control is not None) and (match_v or match_j):
            raise ValueError(
                "V/J-restricted calibration requires V/J-labelled controls; sequence-only controls cannot apply that predicate"
            )
        active = sorted(set(self.loci) & set(q["locus"].drop_nulls()))
        if isinstance(control, dict) and any(
            control.get(gene) is None for gene in active
        ):
            raise ValueError(
                "calibration requires a supplied control for every active locus"
            )
        if calibrate and control is None and active:
            species = self._control_species(active, species)
        hs, cs = [], []
        for gene in self.loci:
            gq = q.filter(pl.col("locus") == gene)
            if not gq.height:
                continue
            h = self._index.annotate(
                gq,
                sp,
                gene=gene,
                threads=threads,
                match_v=match_v,
                match_j=match_j,
                align=align,
                progress=progress,
                exclude_exact=exclude_exact,
            )
            ctrl = control.get(gene) if isinstance(control, dict) else control
            if calibrate and ctrl is None:
                from .evalue.control import background

                ctrl = background(gene, species)
            ctrl_hits = None
            if ctrl is not None:
                if len(ctrl) == 0:
                    raise ValueError("control index is empty")
                results = _search_unique_junctions(ctrl, gq["cdr3"].to_list(), sp, threads)
                rows = [
                    (qid, h.n_subs + h.n_ins + h.n_dels)
                    for qid, junction, hl in zip(gq["query_id"], gq["cdr3"], results)
                    for h in hl
                    if not exclude_exact or ctrl.ref_seq(h.ref_id) != junction
                ]
                ctrl_hits = pl.DataFrame(
                    rows,
                    schema={"query_id": q.schema["query_id"], "edits": pl.UInt32},
                    orient="row",
                )
            for radius in radii if radii is not None else [None]:
                distance = pl.col("n_subs").cast(pl.UInt32) + pl.col("n_ins") + pl.col("n_dels")
                selected = h if radius is None else h.filter(distance <= radius)
                selected_controls = ctrl_hits
                if radius is not None and ctrl_hits is not None:
                    selected_controls = ctrl_hits.filter(pl.col("edits") <= radius)
                c = reduce_candidates(
                    selected, self._index.records_for(gene), selected_controls,
                    control_size=len(ctrl) if ctrl is not None else None,
                    score_scale=score_scale, soft_v=soft_v,
                    match_v=match_v, match_j=match_j,
                )
                if radius is not None:
                    c = c.with_columns(pl.lit(radius, dtype=pl.UInt16).alias("radius"))
                cs.append(c)
            hs.append(h)
        return (
            pl.concat(hs, how="diagonal_relaxed")
            if hs
            else _empty_hits(q.schema["query_id"]),
            pl.concat(cs, how="diagonal_relaxed")
            if cs
            else pl.DataFrame(schema=candidate_schema(q.schema["query_id"])),
        )

    def graded_candidates(self, data, *, cdr3="cdr3", v=None, j=None, locus=None,
                          threads=1, control=None, calibrate=False, species=None,
                          matrix=None, sequence_convention=None, exclude_exact=False,
                          return_hits=False, match_v=False, match_j=False):
        """Five fixed edit balls with radius-constrained weighted rankings.

        Unit-cost balls share one minimum-edit search. Weighted search may choose
        a cheaper alignment with more edits, so each weighted radius is searched
        in one complete query batch. No best-radius P-value is selected.
        Hard V/J predicates are supported for uncalibrated ranking; sequence-only
        controls cannot calibrate these predicates and are rejected.
        """
        _, q = _prepare(data, cdr3, v, j, locus, sequence_convention)
        gap = DEFAULT_SCALE if matrix is not None else 1
        if matrix is None:
            sp = search_params("5,2,2,5", gap_open=gap, gap_extend=gap)
            h, c = self._evidence(q, sp, threads=threads, control=control,
                                  calibrate=calibrate, species=species,
                                  score_scale=1.0, exclude_exact=exclude_exact,
                                  radii=range(1, 6), match_v=match_v, match_j=match_j)
            h = h.with_columns((pl.col("n_subs").cast(pl.UInt32) + pl.col("n_ins") +
                                pl.col("n_dels")).alias("minimum_radius"))
        else:
            hs, cs = [], []
            for radius in range(1, 6):
                sp = search_params(f"{radius},{min(radius,2)},{min(radius,2)},{radius}",
                                   matrix=matrix, gap_open=gap, gap_extend=gap)
                rh, rc = self._evidence(q, sp, threads=threads, control=control,
                                       calibrate=calibrate, species=species,
                                       score_scale=400.0, exclude_exact=exclude_exact,
                                       match_v=match_v, match_j=match_j)
                hs.append(rh.with_columns(pl.lit(radius, dtype=pl.UInt16).alias("radius")))
                cs.append(rc.with_columns(pl.lit(radius, dtype=pl.UInt16).alias("radius")))
            h, c = pl.concat(hs, how="diagonal_relaxed"), pl.concat(cs, how="diagonal_relaxed")
        if "radius" not in c.columns:
            c = c.with_columns(pl.lit(None, dtype=pl.UInt16).alias("radius"))
        return (h, c) if return_hits else c

    def ranked_candidates(self, data, *, cdr3="cdr3", v=None, j=None, locus=None,
                          k=10, threads=1, matrix=None, control=None,
                          sequence_convention=None, exclude_exact=False, paired=False):
        """Global gap-block top-K evidence and a separate global enrichment test.

        K is predeclared; candidate label scores are not calibrated probabilities.
        Paired input uses cdr3a/cdr3b and genuine reference pairs.
        """
        from .match.ranked import ranked_evidence
        if paired:
            q = data
        else:
            _, q = _prepare(data, cdr3, v, j, locus, sequence_convention)
        return ranked_evidence(self, q, k=k, threads=threads, matrix=matrix,
                               control=control, exclude_exact=exclude_exact, paired=paired)

    def hits(
        self,
        cdr3s,
        *,
        locus="TRB",
        scope=DEFAULT_SCOPE,
        match_v=False,
        match_j=False,
        exclude_exact=False,
        threads=1,
        **kw,
    ):
        """Long reference observations; internal cdr3 names denote anchor-inclusive junctions."""
        _, q = _prepare(cdr3s, locus=locus)
        sp = (
            search_params(
                scope,
                matrix=load_vdjam(),
                gap_open=DEFAULT_SCALE,
                gap_extend=DEFAULT_SCALE,
            )
            if isinstance(scope, str)
            else scope
        )
        return self._index.annotate(
            q,
            sp,
            gene=_norm_locus(locus),
            threads=threads,
            match_v=match_v,
            match_j=match_j,
            exclude_exact=exclude_exact,
            **kw,
        )

    def candidates(
        self,
        data,
        *,
        cdr3="cdr3",
        v=None,
        j=None,
        locus=None,
        scope=DEFAULT_SCOPE,
        threads=1,
        match_v=False,
        match_j=False,
        control=None,
        calibrate=False,
        species=None,
        score_scale=400.0,
        soft_v=True,
        sequence_convention=None,
        exclude_exact=False,
    ):
        """All per-pMHC candidates, competing evidence and optional fixed-ball calibration.

        With no control/calibrate request, ranking is explicitly uncalibrated. Supply a seqtree
        Index or locus→Index mapping, or calibrate=True for the versioned standard background.
        NED-v1 uses +1 control counts, temperature 400 and existing soft-V weight 0.25.
        ``exclude_exact`` punctures zero-edit hits in both target and control searches.
        """
        _, q = _prepare(data, cdr3, v, j, locus, sequence_convention)
        sp = (
            search_params(
                scope,
                matrix=load_vdjam(),
                gap_open=DEFAULT_SCALE,
                gap_extend=DEFAULT_SCALE,
            )
            if isinstance(scope, str)
            else scope
        )
        return self._evidence(
            q,
            sp,
            threads=threads,
            match_v=match_v,
            match_j=match_j,
            control=control,
            calibrate=calibrate,
            species=species,
            score_scale=score_scale,
            soft_v=soft_v,
            exclude_exact=exclude_exact,
        )[1]

    def annotate(
        self,
        data,
        *,
        cdr3="cdr3",
        v=None,
        j=None,
        locus=None,
        scope=DEFAULT_SCOPE,
        prefix="vdjmatch_",
        threads=1,
        match_v=False,
        match_j=False,
        control=None,
        calibrate=False,
        species=None,
        score_scale=400.0,
        soft_v=True,
        sequence_convention=None,
        exclude_exact=False,
    ):
        """Append best candidate evidence/status, preserving every input row and its order.

        Exact ranking ties abstain (ambiguous). Confidence is exposed as evidence; NED and
        the fixed-ball enrichment P-value are never a probability of the winning label.
        """
        data, q = _prepare(data, cdr3, v, j, locus, sequence_convention)
        sp = (
            search_params(
                scope,
                matrix=load_vdjam(),
                gap_open=DEFAULT_SCALE,
                gap_extend=DEFAULT_SCALE,
            )
            if isinstance(scope, str)
            else scope
        )
        c = self._evidence(
            q,
            sp,
            threads=threads,
            match_v=match_v,
            match_j=match_j,
            control=control,
            calibrate=calibrate,
            species=species,
            score_scale=score_scale,
            soft_v=soft_v,
            exclude_exact=exclude_exact,
        )[1]
        return _append_calls(data, q, c, self.loci, prefix)

    def paired_candidates(
        self,
        data,
        *,
        cdr3a="cdr3_alpha_aa",
        cdr3b="cdr3_beta_aa",
        scope=DEFAULT_SCOPE,
        threads=1,
        control=None,
        calibrate=False,
        species=None,
        align=False,
        score_scale=400.0,
        return_hits=False,
        progress=False,
        exclude_exact=False,
    ):
        """Genuine same-complex evidence; optional fixed-ball chain-independent calibration.

        ``control`` is a TRA/TRB → seqtree Index mapping. Paired reference similarity ranks
        unique junction pairs. Independent observation counts do not multiply ranking votes.
        Joint enrichment assumes independent background chains; it is not a posterior label
        confidence and does not model repertoire co-occurrence (deferred research).
        ``return_hits=True`` returns ``(joint_hits, candidates)`` from the same searches;
        chain metadata in the detailed hits uses ``alpha_`` and ``beta_`` prefixes.
        ``exclude_exact=True`` removes only the exact/exact pair. Its fixed-ball
        calibration uses simultaneous category bounds and a binomial tail with
        failure budget 1e-6, conditional on IID category/target sampling and
        independent background chains. Controls must index unique junctions;
        duplicate exact memberships fail. Deduplicated controls do not establish IID.
        """
        import math

        from seqtree.evalue import evalue_result

        from .evalue.paired import build_paired_ref, _joint_punctured_results

        if not math.isfinite(score_scale) or score_scale <= 0:
            raise ValueError("score_scale must be finite and positive")
        build_paired_ref(self._index.reference)
        if control is not None and (
            not isinstance(control, dict)
            or any(control.get(chain) is None for chain in ("TRA", "TRB"))
        ):
            raise ValueError("paired calibration requires both TRA and TRB controls")
        if calibrate and control is None:
            species = self._control_species(["TRA", "TRB"], species)
        paired_controls = {}
        for chain in ("TRA", "TRB"):
            ctrl = control.get(chain) if isinstance(control, dict) else None
            if calibrate and ctrl is None:
                from .evalue.control import background

                ctrl = background(chain, species)
            if ctrl is not None and not len(ctrl):
                raise ValueError("control index is empty")
            paired_controls[chain] = ctrl
        _, qa = _prepare(data, cdr3a, locus="TRA")
        _, qb = _prepare(data, cdr3b, locus="TRB")
        sp = (
            search_params(
                scope,
                matrix=load_vdjam(),
                gap_open=DEFAULT_SCALE,
                gap_extend=DEFAULT_SCALE,
            )
            if isinstance(scope, str)
            else scope
        )
        ha = self._index.annotate(
            qa, sp, gene="TRA", threads=threads, align=align, progress=progress
        )
        hb = self._index.annotate(
            qb, sp, gene="TRB", threads=threads, align=align, progress=progress
        )
        keys = ["query_id", "complex_id", *PMHC]

        def projection(h, prefix):
            columns = h.columns if return_hits else ["score", "db_cdr3"] + (["query_cdr3"] if exclude_exact else [])
            return (
                h.with_columns(pl.col("complex_id").cast(pl.String))
                .filter(
                    pl.col("complex_id").is_not_null()
                    & ~pl.col("complex_id").is_in(["0", ""])
                )
                .select(
                    *keys,
                    *[pl.col(c).alias(prefix + c) for c in columns if c not in keys],
                )
            )

        detail = projection(ha, "alpha_").join(
            projection(hb, "beta_"), on=keys, nulls_equal=True
        )
        if exclude_exact:
            detail = detail.filter(
                (pl.col("alpha_db_cdr3") != pl.col("alpha_query_cdr3"))
                | (pl.col("beta_db_cdr3") != pl.col("beta_query_cdr3"))
            )
        if return_hits:
            detail = detail.sort(
                ["query_id", "complex_id", "alpha_record_id", "beta_record_id"]
            )
        joint = detail.select(
            *keys,
            pl.col("alpha_score").alias("score"),
            pl.col("beta_score").alias("score_beta"),
            pl.col("alpha_db_cdr3").alias("db_cdr3"),
            pl.col("beta_db_cdr3").alias("db_cdr3_beta"),
        ).unique()
        dtype = qa.schema["query_id"]
        out_schema = {
            "query_id": dtype,
            **{c: pl.String for c in PMHC},
            "ned_score": pl.Float64,
            "n_hits": pl.UInt32,
            "n_records": pl.UInt32,
            "score_margin": pl.Float64,
            "rank": pl.UInt32,
            "n_reference": pl.UInt32,
            "n_control_alpha": pl.UInt32,
            "n_control_beta": pl.UInt32,
            "E": pl.Float64,
            "p_enrichment": pl.Float64,
            "rule_of_three": pl.Boolean,
            "calibration": pl.String,
            "estimator": pl.String,
        }
        if exclude_exact:
            out_schema.update({
                "n_exact_control_alpha": pl.UInt32, "n_exact_control_beta": pl.UInt32,
                "n_control_joint": pl.UInt64, "control_size_alpha": pl.UInt64,
                "control_size_beta": pl.UInt64, "E_raw": pl.Float64,
                "p_poisson_raw": pl.Float64, "p_upper": pl.Float64,
                "finite_control_delta": pl.Float64,
            })
        if not joint.height:
            empty = pl.DataFrame(schema=out_schema)
            return (detail, empty) if return_hits else empty
        records = joint.group_by(["query_id", *PMHC]).agg(
            pl.col("complex_id").n_unique().alias("n_records")
        )
        unique = joint.unique(subset=["query_id", *PMHC, "db_cdr3", "db_cdr3_beta"])
        out = unique.group_by(["query_id", *PMHC]).agg(
            pl.len().alias("n_hits"),
            ((-(pl.col("score") + pl.col("score_beta")) / score_scale).exp())
            .sum()
            .alias("ned_score"),
        )
        out = out.join(records, on=["query_id", *PMHC], nulls_equal=True)
        # Coverage uses unique paired junctions, rather than publication observations.
        ra = self._index.records_for("TRA")
        rb = self._index.records_for("TRB")
        refkeys = ["complex_id", *PMHC]
        ra = ra.with_columns(pl.col("complex_id").cast(pl.String)).filter(
            ~pl.col("complex_id").is_in(["0", ""])
        )
        rb = rb.with_columns(pl.col("complex_id").cast(pl.String)).filter(
            ~pl.col("complex_id").is_in(["0", ""])
        )
        refs = ra.select(*refkeys, "cdr3").join(
            rb.select(*refkeys, "cdr3"), on=refkeys, nulls_equal=True, suffix="_beta"
        )
        cov = (
            refs.unique(subset=[*PMHC, "cdr3", "cdr3_beta"])
            .group_by(PMHC)
            .len()
            .rename({"len": "n_reference"})
        )
        out = out.join(cov, on=PMHC, nulls_equal=True)
        counts = []
        exact_counts = []
        sizes = []
        for locus, q in [("TRA", qa), ("TRB", qb)]:
            ctrl = paired_controls[locus]
            if ctrl is not None:
                hits = _search_unique_junctions(ctrl, q["cdr3"].to_list(), sp, threads)
                counts.append(dict(zip(q["query_id"], map(len, hits))))
                exact_counts.append(dict(zip(q["query_id"], [
                    sum(ctrl.ref_seq(h.ref_id) == junction for h in group)
                    for junction, group in zip(q["cdr3"], hits)
                ])) if exclude_exact else {})
                sizes.append(len(ctrl))
            else:
                counts.append({})
                exact_counts.append({})
                sizes.append(None)
        nca = [counts[0].get(qid, 0) for qid in out["query_id"]]
        ncb = [counts[1].get(qid, 0) for qid in out["query_id"]]
        out = out.with_columns(
            pl.Series("n_control_alpha", nca), pl.Series("n_control_beta", ncb)
        )
        if exclude_exact:
            out = out.with_columns(
                pl.Series("n_exact_control_alpha", [exact_counts[0].get(qid, 0) for qid in out["query_id"]]),
                pl.Series("n_exact_control_beta", [exact_counts[1].get(qid, 0) for qid in out["query_id"]]),
                pl.lit(sizes[0], dtype=pl.UInt64).alias("control_size_alpha"),
                pl.lit(sizes[1], dtype=pl.UInt64).alias("control_size_beta"),
            )
        if exclude_exact and all(n is not None for n in sizes):
            stats = _joint_punctured_results(
                out["n_hits"].to_list(), nca, ncb, out["n_reference"].to_list(),
                sizes[0], sizes[1], out["n_exact_control_alpha"].to_list(),
                out["n_exact_control_beta"].to_list(),
            )
            out = out.with_columns(*[pl.Series(k, v) for k, v in stats.items()],
                                   pl.lit(False).alias("rule_of_three"))
            calibration = "paired_independent_joint_punctured_binomial_bound"
        elif all(n is not None for n in sizes):
            stats = [
                evalue_result(
                    nt,
                    min(sizes[0], a or 3) * min(sizes[1], b or 3),
                    nr,
                    sizes[0] * sizes[1],
                )
                for nt, a, b, nr in out.select(
                    "n_hits", "n_control_alpha", "n_control_beta", "n_reference"
                ).iter_rows()
            ]
            out = out.with_columns(
                pl.Series("E", [r["E"] for r in stats]),
                pl.Series("p_enrichment", [r["p_enrichment"] for r in stats]),
                (
                    (pl.col("n_control_alpha") == 0) | (pl.col("n_control_beta") == 0)
                ).alias("rule_of_three"),
            )
            calibration = "paired_independent_fixed_ball"
        else:
            out = out.with_columns(
                pl.lit(None, dtype=pl.Float64).alias("E"),
                pl.lit(None, dtype=pl.Float64).alias("p_enrichment"),
                pl.lit(False).alias("rule_of_three"),
            )
            calibration = "uncalibrated"
            if exclude_exact:
                out = out.with_columns(*[
                    pl.lit(None, dtype=out_schema[k]).alias(k)
                    for k in ("n_control_joint", "E_raw", "p_poisson_raw", "p_upper", "finite_control_delta")
                ])
        out = out.sort(
            ["query_id", "ned_score", *PMHC],
            descending=[False, True, False, False, False, False],
            nulls_last=True,
        )
        out = out.with_columns(
            pl.int_range(1, pl.len() + 1).over("query_id").cast(pl.UInt32).alias("rank")
        )
        best = out.group_by("query_id").agg(
            pl.col("ned_score").first().alias("_best"),
            pl.col("ned_score")
            .get(1, null_on_oob=True)
            .fill_null(0.0)
            .alias("_second"),
        )
        out = out.join(best, on="query_id").with_columns(
            (
                pl.col("ned_score")
                - pl.when(pl.col("rank") == 1)
                .then(pl.col("_second"))
                .otherwise(pl.col("_best"))
            ).alias("score_margin"),
            pl.lit(calibration).alias("calibration"),
            pl.lit("paired-reference-v1").alias("estimator"),
        )
        candidates = out.select([pl.col(c).cast(t) for c, t in out_schema.items()])
        return (detail, candidates) if return_hits else candidates

    def annotate_paired(
        self,
        data,
        *,
        cdr3a="cdr3_alpha_aa",
        cdr3b="cdr3_beta_aa",
        scope=DEFAULT_SCOPE,
        prefix="vdjmatch_",
        threads=1,
        control=None,
        calibrate=False,
        species=None,
        align=False,
        score_scale=400.0,
        exclude_exact=False,
    ):
        """Append same-reference-complex evidence with explicit chain-independent calibration."""
        _, qa = _prepare(data, cdr3a, locus="TRA")
        _, qb = _prepare(data, cdr3b, locus="TRB")
        c = self.paired_candidates(
            data,
            cdr3a=cdr3a,
            cdr3b=cdr3b,
            scope=scope,
            threads=threads,
            control=control,
            calibrate=calibrate,
            species=species,
            align=align,
            score_scale=score_scale,
            exclude_exact=exclude_exact,
        )
        q = qa.join(qb.select("query_id"), on="query_id")
        return _append_calls(data, q, c, self.loci, prefix, paired=True)


def _append_calls(data, q, c, loci, prefix, paired=False):
    top = c.filter(pl.col("rank") == 1)
    names = PMHC + ["ned_score", "n_hits", "score_margin"]
    extra = [
        s
        for s in [
            "n_clonotypes",
            "n_studies",
            "n_competing_clonotypes",
            "nearest_edits",
            "E",
            "p_enrichment",
            "rule_of_three",
            "calibration",
            "estimator",
        ]
        if s in c.columns
    ]
    names += extra
    base = (
        data.select("query_id")
        if "query_id" in data.columns
        else pl.DataFrame({"query_id": pl.Series(range(data.height), dtype=pl.UInt32)})
    )
    valid = q.select("query_id", pl.col("locus").alias("_locus"))
    base = (
        base.with_row_index("_order")
        .join(valid, on="query_id", how="left")
        .join(top.select("query_id", *names), on="query_id", how="left")
        .sort("_order")
    )
    base = base.with_columns(
        pl.when(pl.col("_locus").is_null())
        .then(pl.lit("invalid_query"))
        .when(
            ~pl.col("_locus").is_in(loci)
            if not paired
            else pl.lit(not {"TRA", "TRB"} <= set(loci))
        )
        .then(pl.lit("no_reference"))
        .when(pl.col("ned_score").is_null())
        .then(pl.lit("no_joint_hit" if paired else "no_hit"))
        .when(pl.col("score_margin") == 0)
        .then(pl.lit("ambiguous"))
        .otherwise(pl.lit("matched"))
        .alias("status")
    )
    base = base.with_columns(
        *[
            pl.when(pl.col("status") == "matched")
            .then(pl.col(s))
            .otherwise(None)
            .alias(s)
            for s in PMHC
        ]
    )
    rename = {
        c: prefix + ("score" if c == "ned_score" else c) for c in names + ["status"]
    }
    result = base.select(*names, "status").rename(rename)
    overlap = set(data.columns) & set(result.columns)
    if overlap:
        raise ValueError(
            f"annotation output columns already exist: {sorted(overlap)}; choose another prefix"
        )
    return data.hstack(result)


def annotate(data, *, reference=None, **kw):
    """Convenience call; provide an Annotator to reuse reference inputs explicitly."""
    ann = (
        reference
        if isinstance(reference, Annotator)
        else Annotator.from_path(reference)
        if reference is not None
        else Annotator.latest()
    )
    return ann.annotate(data, **kw)
