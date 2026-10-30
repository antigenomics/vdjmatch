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
    for extra in ("sequence_id", "pair_id", "original_v", "v_imputed", "original_j", "j_imputed"):
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
    def from_path(cls, path, *, species=None, asset="default", threads=0, **filters):
        from .db.index import load_or_build_index

        index, profiles, manifest = load_or_build_index(
            path, asset=asset, species=species, threads=threads, **filters)
        result = cls(index, provenance=manifest["identity"]["source"])
        result.reference_profiles = profiles
        result.reference_index_manifest = manifest
        return result

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

    def _matching_queries(self, q, *, match_v=False, match_j=False, species=None):
        """Retain a complete per-query gene-imputation lookup, including no-hits."""
        if not (match_v or match_j) or not q.height:return q
        from .evalue.control import impute_gene_calls
        parts=[]
        for (gene,),part in q.partition_by('locus',as_dict=True,maintain_order=True).items():
            kinds=tuple(k for k,flag in [('v',match_v),('j',match_j)] if flag)
            if gene not in {'TRA','TRB'}:
                parts.append(part.with_columns(*[pl.col(kind).alias('original_'+kind) for kind in kinds if 'original_'+kind not in part.columns],
                    *[pl.lit(False).alias(kind+'_imputed') for kind in kinds if kind+'_imputed' not in part.columns]))
                continue
            ref=self._index.records_for(gene)
            organism=species or next(iter(ref['species'].drop_nulls().unique()),'human')
            parts.append(impute_gene_calls(part,species=organism,locus=gene,kinds=kinds,reference=ref))
        return pl.concat(parts,how='diagonal_relaxed').join(q.select('query_id').with_row_index('_query_order'),on='query_id',validate='1:1').sort('_query_order').drop('_query_order')

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
        background_mass=False,
        exclude_epitopes=None,
        retain_hits=True,
        _prepare_consumer=False,
        _hit_consumer=None,
    ):
        import math

        if not math.isfinite(score_scale) or score_scale <= 0:
            raise ValueError("score_scale must be finite and positive")
        if exclude_epitopes is not None:
            if (exclude_epitopes["query_id"].n_unique() != exclude_epitopes.height
                or exclude_epitopes.select(pl.any_horizontal(pl.all().is_null()).any()).item()):
                raise ValueError("epitope exclusions require one non-null peptide per query")
        if (calibrate or control is not None) and (match_v or match_j) and not (
            isinstance(control, pl.DataFrame) or isinstance(control, dict) and
            all(isinstance(value, pl.DataFrame) for value in control.values())
        ):
            raise ValueError(
                "V/J-restricted calibration requires V/J-labelled controls; sequence-only controls cannot apply that predicate"
            )
        active = sorted(set(self.loci) & set(q["locus"].drop_nulls()))
        if not retain_hits and any(not callable(getattr(self._index.index_for(gene),
                'search_batch_consume', None)) for gene in active):
            raise RuntimeError('bounded candidate reduction requires seqtree Index.search_batch_consume; rebuild/install the current native source')
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
            gene_flags = None
            reference_records = self._index.records_for(gene)
            if match_v or match_j:
                from .evalue.control import impute_gene_calls
                reference_records = impute_gene_calls(reference_records,
                    species=species or next(iter(reference_records['species'].drop_nulls().unique()),'human'),
                    locus=gene, kinds=tuple(k for k,flag in [('v',match_v),('j',match_j)] if flag))
                gq = self._matching_queries(gq, match_v=match_v, match_j=match_j, species=species)
                kinds = tuple(k for k, enabled in [('v', match_v), ('j', match_j)] if enabled)
                gene_flags = gq.select('query_id', *[pl.col(k).alias('query_'+k)
                    for kind in kinds for k in ('original_'+kind, kind+'_imputed')])
            ctrl = control.get(gene) if isinstance(control, dict) else control
            if calibrate and ctrl is None:
                from .evalue.control import background

                ctrl = background(gene, species)
            ctrl_hits = None
            control_cdf = None
            control_totals = control_sizes = None
            if ctrl is not None:
                if len(ctrl) == 0:
                    raise ValueError("control index is empty")
                if isinstance(ctrl, pl.DataFrame):
                    from .evalue.control import conditioned_counts, impute_gene_calls, control_population
                    if not (match_v or match_j):
                        raise ValueError('gene-labelled controls require match_v or match_j')
                    ctrl = control_population(ctrl, species=self._control_species([gene], species), locus=gene)
                    if not ctrl.height:raise ValueError('control has no rows matching requested species/locus')
                    labelled = impute_gene_calls(ctrl, species=self._control_species([gene], species),
                        locus=gene, kinds=kinds)
                    control_cdf, control_totals, control_sizes = conditioned_counts(gq, labelled, sp,
                        match_v=match_v, match_j=match_j, threads=threads, exclude_exact=exclude_exact)
                elif (radii is not None or 0 < sp.max_total_edits <= 5) and sp.engine in {'auto','seqtm'} and sp.mode == 'all' and sp.gap_open > 0 and sp.gap_extend > 0:
                    import numpy as np

                    sequences = gq['cdr3'].unique(maintain_order=True).to_list()
                    histogram = np.asarray(ctrl.edit_histogram_batch(sequences, sp, threads, False), dtype=np.uint64)
                    if exclude_exact:
                        histogram[:, 0] = 0
                    cumulative = histogram.cumsum(axis=1, dtype=np.uint64)
                    control_cdf = pl.DataFrame({
                        'cdr3': np.repeat(sequences, cumulative.shape[1]),
                        'edits': np.tile(np.arange(cumulative.shape[1], dtype=np.uint32), len(sequences)),
                        '_nc': cumulative.ravel(),
                    }).join(gq.select('query_id', 'cdr3'), on='cdr3', validate='m:m').drop('cdr3')
                    control_totals = control_cdf.filter(pl.col('edits') == sp.max_total_edits).select('query_id', pl.col('_nc').alias('n_control'))
                else:
                    results = _search_unique_junctions(ctrl, gq["cdr3"].to_list(), sp, threads)
                    rows = [
                        (qid, h.n_subs + h.n_ins + h.n_dels,
                         exclude_exact and ctrl.ref_seq(h.ref_id) == junction)
                        for qid, junction, hl in zip(gq["query_id"], gq["cdr3"], results)
                        for h in hl
                    ]
                    ctrl_hits = pl.DataFrame(
                        rows,
                        schema={"query_id": q.schema["query_id"], "edits": pl.UInt32, "_self": pl.Boolean},
                        orient="row",
                    )
                    ctrl_hits = ctrl_hits.filter(~pl.col("_self")).drop("_self")
            def consume(h):
                if exclude_epitopes is not None:
                    h = h.join(exclude_epitopes, on=['query_id', 'epitope'], how='anti')
                if _hit_consumer is not None:
                    _hit_consumer(h)
                for radius in radii if radii is not None else [None]:
                    distance = pl.col("n_subs").cast(pl.UInt32) + pl.col("n_ins") + pl.col("n_dels")
                    selected = h if radius is None else h.filter(distance <= radius)
                    selected_controls = ctrl_hits
                    if radius is not None and ctrl_hits is not None:
                        selected_controls = ctrl_hits.filter(pl.col("edits") <= radius)
                    selected_totals = control_totals
                    if control_cdf is not None and radius is not None:
                        selected_totals = control_cdf.filter(pl.col('edits') == radius).select('query_id', pl.col('_nc').alias('n_control'))
                    c = reduce_candidates(
                        selected, reference_records, selected_controls,
                        control_size=len(ctrl) if ctrl is not None else None,
                        control_counts=control_cdf, control_totals=selected_totals, control_sizes=control_sizes,
                        score_scale=score_scale, soft_v=soft_v,
                        match_v=match_v, match_j=match_j, background_mass=background_mass,
                        exclude_exact=exclude_exact,
                    )
                    if gene_flags is not None:
                        c = c.join(gene_flags, on='query_id', how='left', validate='m:1')
                    if radius is not None:
                        c = c.with_columns(pl.lit(radius, dtype=pl.UInt16).alias("radius"))
                    cs.append(c)
                return h
            if _prepare_consumer:
                return consume, cs, reference_records
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
                _reference_records=reference_records if (match_v or match_j) else None,
                **({} if retain_hits else {'consumer':consume, 'chunk':32, '_candidate_only':True}),
            )
            if retain_hits:
                h = consume(h)
                hs.append(h)
        candidates = pl.concat(cs, how='diagonal_relaxed') if cs else pl.DataFrame(schema=candidate_schema(q.schema['query_id']))
        if not retain_hits and candidates.height:
            candidates = candidates.sort('query_locus', *(['radius'] if radii is not None else []), 'query_id', 'rank')
        return (
            pl.concat(hs, how="diagonal_relaxed")
            if hs
            else _empty_hits(q.schema["query_id"]),
            candidates,
        )

    def graded_candidates(self, data, *, cdr3="cdr3", v=None, j=None, locus=None,
                          threads=1, control=None, calibrate=False, species=None,
                          matrix=None, sequence_convention=None, exclude_exact=False,
                          return_hits=False, match_v=False, match_j=False, progress=False,
                          max_radius=5):
        """Fixed edit balls through max_radius, with constrained weighted rankings.

        Unit-cost balls share one minimum-edit search. Weighted search may choose
        a cheaper alignment with more edits, so each weighted radius is searched
        in one complete query batch. No best-radius P-value is selected.
        Hard V/J calibration accepts canonical cdr3/v/j labelled control frames
        (or a locus-to-frame mapping), counting distinct junctions per predicate.
        Sequence-only controls cannot calibrate these predicates and are rejected.
        Missing genes use flagged, label-independent modal imputation.
        """
        if isinstance(max_radius, bool) or not isinstance(max_radius, int) or not 1 <= max_radius <= 5:
            raise ValueError('max_radius must be an integer from 1 through 5')
        _, q = _prepare(data, cdr3, v, j, locus, sequence_convention)
        gap = DEFAULT_SCALE if matrix is not None else 1
        if matrix is None:
            sp = search_params(f"{max_radius},{min(max_radius,2)},{min(max_radius,2)},{max_radius}", gap_open=gap, gap_extend=gap)
            h, c = self._evidence(q, sp, threads=threads, control=control,
                                  calibrate=calibrate, species=species,
                                  score_scale=1.0, exclude_exact=exclude_exact,
                                  radii=range(1, max_radius+1), match_v=match_v, match_j=match_j, retain_hits=return_hits, progress=progress)
            h = h.with_columns((pl.col("n_subs").cast(pl.UInt32) + pl.col("n_ins") +
                                pl.col("n_dels")).alias("minimum_radius"))
        else:
            hs, cs = [], []
            for radius in range(1, max_radius+1):
                sp = search_params(f"{radius},{min(radius,2)},{min(radius,2)},{radius}",
                                   matrix=matrix, gap_open=gap, gap_extend=gap)
                rh, rc = self._evidence(q, sp, threads=threads, control=control,
                                       calibrate=calibrate, species=species,
                                       score_scale=400.0, exclude_exact=exclude_exact,
                                       match_v=match_v, match_j=match_j, retain_hits=return_hits, progress=progress)
                hs.append(rh.with_columns(pl.lit(radius, dtype=pl.UInt16).alias("radius")))
                cs.append(rc.with_columns(pl.lit(radius, dtype=pl.UInt16).alias("radius")))
            h, c = pl.concat(hs, how="diagonal_relaxed"), pl.concat(cs, how="diagonal_relaxed")
        if "radius" not in c.columns:
            c = c.with_columns(pl.lit(None, dtype=pl.UInt16).alias("radius"))
        return (h, c) if return_hits else c

    def unified_candidates(self, data, *, cdr3="cdr3", v=None, j=None, locus=None,
                           threads=1, control=None, species=None, exclude_exact=False,
                           sequence_convention=None, return_hits=False, distance="edit"):
        """Experimental fixed-ball background-mass evidence for a selected distance.

        Both edit and anchored BLOSUM62 gap-block distances cover all query lengths.
        No length switch, label-dependent route or top-K cap.
        Sequence-only controls do not provide V-conditioned or paired calibration.
        """
        from .match.unified import unified_evidence
        _, q = _prepare(data, cdr3, v, j, locus, sequence_convention)
        h, c = unified_evidence(self, q, threads=threads, control=control,
                                species=species, exclude_exact=exclude_exact, distance=distance)
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
        ``exclude_exact`` filters exact hits in target and control searches; model populations stay fixed.
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
            retain_hits=False,
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
        q = self._matching_queries(q, match_v=match_v, match_j=match_j, species=species)
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
            retain_hits=False,
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
        _axis_consumers=None,
        _prepared_queries=None,
        _reference_views=None,
        _linked_radii=None,
        _match_v=False,
        _match_j=False,
    ):
        """Genuine same-complex evidence; optional fixed-ball chain-independent calibration.

        ``control`` is a TRA/TRB → seqtree Index mapping. Paired reference similarity ranks
        unique junction pairs. Independent observation counts do not multiply ranking votes.
        Joint enrichment assumes independent background chains; it is not a posterior label
        confidence and does not model repertoire co-occurrence (deferred research).
        ``return_hits=True`` returns ``(joint_hits, candidates)`` from the same searches;
        chain metadata in the detailed hits uses ``alpha_`` and ``beta_`` prefixes.
        Candidate-only output requires native ``Index.search_paired_batch_consume``:
        one worker pool searches both fixed chain indexes and reduces complete query
        blocks. Hit buffers are bounded; the final candidate table remains in memory.
        ``exclude_exact=True`` removes only the exact/exact hit pair, retaining full
        reference and control population sizes. Its fixed-ball
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
        if _prepared_queries is not None:
            qa, qb = _prepared_queries
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

        groups = ['query_id'] + (['radius'] if _linked_radii is not None else [])
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
        if _linked_radii is not None:
            out_schema['radius'] = pl.UInt16
        def summarize(joint):
            records = joint.group_by([*groups, *PMHC]).agg(
                pl.col("complex_id").n_unique().alias("n_records")
            )
            unique = joint.unique(subset=[*groups, *PMHC, "db_cdr3", "db_cdr3_beta"])
            out = unique.group_by([*groups, *PMHC]).agg(
                pl.len().alias("n_hits"),
                ((-(pl.col("score") + pl.col("score_beta")) / score_scale).exp())
                .sum()
                .alias("ned_score"),
            )
            return out.join(records, on=[*groups, *PMHC], nulls_equal=True)

        def linked(ha, hb):
            detail = projection(ha, "alpha_").join(
                projection(hb, "beta_"), on=keys, nulls_equal=True
            )
            if exclude_exact:
                detail = detail.filter(
                    (pl.col("alpha_db_cdr3") != pl.col("alpha_query_cdr3"))
                    | (pl.col("beta_db_cdr3") != pl.col("beta_query_cdr3"))
                )
            return detail

        if return_hits:
            ha = self._index.annotate(
                qa, sp, gene="TRA", threads=threads, align=align, progress=progress
            )
            hb = self._index.annotate(
                qb, sp, gene="TRB", threads=threads, align=align, progress=progress
            )
            detail = linked(ha, hb).sort(
                ["query_id", "complex_id", "alpha_record_id", "beta_record_id"])
            if not detail.height:
                return detail, pl.DataFrame(schema=out_schema)
            out = summarize(detail.select(
                *keys, pl.col("alpha_score").alias("score"),
                pl.col("beta_score").alias("score_beta"),
                pl.col("alpha_db_cdr3").alias("db_cdr3"),
                pl.col("beta_db_cdr3").alias("db_cdr3_beta")).unique())
        else:
            import numpy as np
            ia, ib = (self._index.index_for(gene) for gene in ("TRA", "TRB"))
            if ia is None or ib is None:
                if _axis_consumers is None:
                    return pl.DataFrame(schema=out_schema)
                from seqtree import Index
                # An absent linked partner cannot veto a supported marginal.
                if ia is None:ia = Index.build([],alphabet='aa')
                if ib is None:ib = Index.build([],alphabet='aa')
            if not callable(getattr(ia, 'search_paired_batch_consume', None)):
                raise RuntimeError('bounded paired reduction requires seqtree Index.search_paired_batch_consume; rebuild/install the current native source')
            valid = lambda c: pl.col(c).is_not_null() & pl.col(c).str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
            population = qa.select('query_id', pl.col('cdr3').alias('_a')).join(
                qb.select('query_id', pl.col('cdr3').alias('_b')), on='query_id', validate='1:1')
            if _axis_consumers is None:
                population = population.filter(valid('_a') & valid('_b'))
            else:
                population = population.with_columns(
                    pl.when(valid('_a')).then(pl.col('_a')).otherwise(pl.lit('')).alias('_a'),
                    pl.when(valid('_b')).then(pl.col('_b')).otherwise(pl.lit('')).alias('_b'))
            population = population.with_row_index('_qi')
            schemas = ['_qi','ref_id','score','n_subs','n_ins','n_dels']
            refs = []
            for axis,(idx,gene) in enumerate(((ia,'TRA'),(ib,'TRB'))):
                records = (_reference_views[axis] if _reference_views is not None else self._index.records_for(gene))
                fields = ['complex_id',*PMHC] + (['v','j','reference_id','vdjdb_score'] if _axis_consumers is not None else [])
                records = records.select(pl.col('cdr3').alias('db_cdr3'),*fields)
                if _axis_consumers is None:
                    records = records.filter(pl.col('complex_id').cast(pl.String).is_not_null() &
                        ~pl.col('complex_id').cast(pl.String).is_in(['0',''])).unique()
                else:
                    records = records.rename({'v':'db_v','j':'db_j'})
                refs.append(pl.DataFrame({'ref_id':pl.Series(range(len(idx)), dtype=pl.UInt32),
                    'db_cdr3':pl.Series(idx.ref_seqs(),dtype=pl.String)}).join(records,on='db_cdr3'))
            outputs = []
            def consume(start, stop, alpha, beta):
                def expand(buffer, axis, records, queries):
                    qfields = queries.select('query_id',pl.col('v').alias('query_v'),
                        pl.col('j').alias('query_j'),pl.col('locus').alias('query_locus'),'count')
                    hits = (pl.DataFrame(np.asarray(buffer), schema=schemas, orient='row')
                        .with_columns(pl.col('score').cast(pl.Int32),pl.col('n_subs','n_ins','n_dels').cast(pl.UInt16))
                        .join(records, on='ref_id')
                        .join(population.slice(start,stop-start).select('_qi','query_id',
                            pl.col(axis).alias('query_cdr3')), on='_qi', validate='m:1')
                        .join(qfields,on='query_id',validate='m:1').drop('_qi','ref_id')
                        .sort('query_id','score','db_cdr3'))
                    if _match_v or _match_j:
                        from .match.engine import _genefam
                        for kind,flag in [('v',_match_v),('j',_match_j)]:
                            if flag:hits = hits.filter(_genefam(pl.col('db_'+kind)) == _genefam(pl.col('query_'+kind)))
                    return hits
                ha,hb = expand(alpha,'_a',refs[0],qa),expand(beta,'_b',refs[1],qb)
                if _axis_consumers is not None:
                    for callback,hits in zip(_axis_consumers,(ha,hb)):
                        marginal = hits.filter(pl.col('query_cdr3')!=pl.col('db_cdr3')) if exclude_exact else hits
                        callback(marginal)
                for radius in _linked_radii if _linked_radii is not None else [None]:
                    def selected(hits):
                        if radius is not None:
                            hits = hits.filter(pl.sum_horizontal('n_subs','n_ins','n_dels')<=radius)
                        # Paired votes are unique complex/pMHC/junction projections;
                        # source multiplicity remains available to marginal reducers above.
                        return hits.select(*keys,'score','db_cdr3','query_cdr3').unique(maintain_order=True)
                    detail = linked(selected(ha), selected(hb))
                    if detail.height:
                        joint = detail.select(*keys,pl.col('alpha_score').alias('score'),
                            pl.col('beta_score').alias('score_beta'),pl.col('alpha_db_cdr3').alias('db_cdr3'),
                            pl.col('beta_db_cdr3').alias('db_cdr3_beta')).unique()
                        if radius is not None:joint = joint.with_columns(pl.lit(radius,dtype=pl.UInt16).alias('radius'))
                        outputs.append(summarize(joint))
            if progress:
                import sys
                print(f'paired: searching {population.height:,} queries in one native population batch', file=sys.stderr)
            ia.search_paired_batch_consume(ib,population['_a'].to_list(),population['_b'].to_list(),
                sp,consume,threads=threads,chunk_size=32,max_hit_bytes=64*1024**2,**({'skip_empty':True} if _axis_consumers is not None else {}))
            if progress:
                print('paired: native search complete',file=sys.stderr)
            if not outputs:
                return pl.DataFrame(schema=out_schema)
            out = pl.concat(outputs)
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
                if (sp.engine in {'auto','seqtm'} and sp.mode == 'all'
                    and sp.gap_open > 0 and sp.gap_extend > 0 and 0 < sp.max_total_edits <= 5
                    and callable(getattr(ctrl, 'edit_histogram_batch', None))):
                    keys = list(dict.fromkeys(q['cdr3']))
                    histogram = ctrl.edit_histogram_batch(keys, sp, threads, False)
                    totals = dict(zip(keys, map(sum, histogram)))
                    exact = dict(zip(keys, (row[0] for row in histogram)))
                    counts.append(dict(zip(q['query_id'], (totals[key] for key in q['cdr3']))))
                    exact_counts.append(dict(zip(q['query_id'], (exact[key] for key in q['cdr3']))) if exclude_exact else {})
                else:
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
            [*groups, "ned_score", *PMHC],
            descending=[False]*len(groups)+[True,False,False,False,False],
            nulls_last=True,
        )
        out = out.with_columns(
            pl.int_range(1, pl.len() + 1).over(groups).cast(pl.UInt32).alias("rank")
        )
        best = out.group_by(groups).agg(
            pl.col("ned_score").first().alias("_best"),
            pl.col("ned_score")
            .get(1, null_on_oob=True)
            .fill_null(0.0)
            .alias("_second"),
        )
        out = out.join(best, on=groups).with_columns(
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

    def paired_scope_candidates(self, data, *, cdr3a='cdr3a', cdr3b='cdr3b',
                                scope=DEFAULT_SCOPE, ball=False, matrix=None,
                                va='va',vb='vb',ja='ja',jb='jb',
                                match_v=False, match_j=False, control=None, calibrate=False,
                                species=None, threads=1, progress=False, exclude_exact=False,
                                max_radius=5):
        """Separate linked and full-reference marginal scope evidence.

        Each chain uses the declared edit budget; paired radius is the bottleneck
        maximum of the two edit counts. Linked ranking retains the existing paired
        exponential mass; marginal ranking/calibration retains NED. Linked E-values
        and final linked-first/marginal-union selection calibration are unavailable.
        Empty/invalid partners retain valid marginal evidence. No detailed hit table
        is retained. Native coordinated consumption with ``skip_empty`` is required.
        """
        if isinstance(max_radius, bool) or not isinstance(max_radius, int) or not 1 <= max_radius <= 5:
            raise ValueError('max_radius must be an integer from 1 through 5')
        if not ball and max_radius != 5:
            raise ValueError('max_radius applies only to paired ball mode')
        if any(not callable(getattr(self._index.index_for(chain),'search_paired_batch_consume',None))
               for chain in self.loci if chain in {'TRA','TRB'}):
            raise RuntimeError('paired scope axes require seqtree Index.search_paired_batch_consume with skip_empty; install the compatible native build')
        _,qa = _prepare(data,cdr3a,va,ja,locus='TRA')
        _,qb = _prepare(data,cdr3b,vb,jb,locus='TRB')
        source = data if 'query_id' in data.columns else data.with_row_index('query_id')
        restored = []
        for q,seq,vcol,jcol,chain in ((qa,cdr3a,va,ja,'TRA'),(qb,cdr3b,vb,jb,'TRB')):
            # _prepare drops missing/invalid sequences. They remain original paired
            # observations here, with an explicitly empty native axis and audit flag.
            raw = source.select('query_id',
                pl.col(seq).cast(pl.String).str.strip_chars().str.to_uppercase().alias('_raw_cdr3'),
                (pl.col(vcol).cast(pl.String) if vcol in source.columns else pl.lit(None,dtype=pl.String)).alias('_raw_v'),
                (pl.col(jcol).cast(pl.String) if jcol in source.columns else pl.lit(None,dtype=pl.String)).alias('_raw_j'))
            q = raw.join(q,on='query_id',how='left',validate='1:1',maintain_order='left').with_columns(
                pl.coalesce('cdr3','_raw_cdr3').alias('cdr3'),pl.coalesce('v','_raw_v').alias('v'),
                pl.coalesce('j','_raw_j').alias('j'),pl.col('locus').fill_null(chain),
                pl.col('count').fill_null(0)).drop('_raw_cdr3','_raw_v','_raw_j')
            restored.append(q)
        qa,qb = restored
        queries = [self._matching_queries(q,match_v=match_v,match_j=match_j,species=species)
            if match_v or match_j else q for q in (qa,qb)]
        parts = {'linked':[],'TRA':[],'TRB':[]}
        gap = DEFAULT_SCALE if matrix is not None else 1
        radii = range(1,max_radius+1) if ball and matrix is not None else [None]
        for radius in radii:
            params = search_params(f'{radius},{min(radius,2)},{min(radius,2)},{radius}',
                matrix=matrix,gap_open=gap,gap_extend=gap) if radius is not None else (
                search_params(f'{max_radius},{min(max_radius,2)},{min(max_radius,2)},{max_radius}',gap_open=1,gap_extend=1) if ball else
                (search_params(scope,matrix=matrix or '',gap_open=gap,gap_extend=gap) if isinstance(scope,str) else scope))
            callbacks=[];marginals=[];references=[]
            for chain,q in zip(('TRA','TRB'),queries):
                valid_q = q.filter(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False))
                if valid_q.height and chain in self.loci:
                    callback,rows,reference = self._evidence(valid_q,params,threads=threads,
                        match_v=match_v,match_j=match_j,control=control,calibrate=calibrate,
                        species=species,score_scale=400. if matrix is not None else 1.,
                        exclude_exact=exclude_exact,radii=range(1,max_radius+1) if ball and matrix is None else None,
                        retain_hits=False,_prepare_consumer=True)
                else:
                    callback,rows,reference = lambda h:None,[],self._index.records_for(chain)
                callbacks.append(callback);marginals.append(rows);references.append(reference)
            linked = self.paired_candidates(data,cdr3a=cdr3a,cdr3b=cdr3b,scope=params,
                threads=threads,progress=progress,exclude_exact=exclude_exact,
                score_scale=400. if matrix is not None else 1.,_axis_consumers=callbacks,
                _prepared_queries=queries,_reference_views=references,
                _linked_radii=range(1,max_radius+1) if ball and matrix is None else None,
                _match_v=match_v,_match_j=match_j)
            frames = [linked,*[pl.concat(rows,how='diagonal_relaxed') if rows else
                pl.DataFrame(schema=candidate_schema(qa.schema['query_id'])) for rows in marginals]]
            for axis,frame in zip(parts,frames):
                if radius is not None:frame = frame.with_columns(pl.lit(radius,dtype=pl.UInt16).alias('radius'))
                elif 'radius' not in frame.columns:frame = frame.with_columns(pl.lit(None,dtype=pl.UInt16).alias('radius'))
                parts[axis].append(frame)
        result = {axis:pl.concat(frames,how='diagonal_relaxed').sort('radius','query_id','rank',nulls_last=True)
            for axis,frames in parts.items()}
        return result,queries

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
    gene_flags=[name for kind in ('v','j') for name in ('original_'+kind,kind+'_imputed') if name in q.columns]
    valid = q.select("query_id", pl.col("locus").alias("_locus"), *gene_flags)
    base = (
        base.with_row_index("_order")
        .join(valid, on="query_id", how="left")
        .join(top.select("query_id", *names), on="query_id", how="left")
        .sort("_order")
    )
    names += gene_flags
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


def repertoire_enrichment(observations, reference, controls, **options):
    """Reference-null repertoire test; see aggregate.reference_enrichment contracts.

    Returns summary, complete mock counts and real chain hit links. Inputs are
    normalized tables with explicit Pgen bins; the CLI owns raw-file ingestion
    and model-backed Pgen computation. Existing annotation scores are unchanged.
    """
    from .aggregate.repertoire_enrichment import reference_enrichment
    return reference_enrichment(observations, reference, controls, **options)
