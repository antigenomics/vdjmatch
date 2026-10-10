"""Fixed geometry and empirical background-mass ranking across length routes."""
from __future__ import annotations

import numpy as np
import polars as pl
from seqtree import SubstitutionMatrix, gapblock

from ..aggregate.candidates import candidate_schema, candidates
from . import search_params


def gapped_extension(queries, reference, control, *, species='human', threads=1, radius=90, temperature=12):
    """Experimental additive density on edges outside the original five-substitution ball.

    The kernel uses total TCRdist, while its empirical CDF uses weighted CDR3
    distance only. Allele-stripped equal V genes have categorical weight one;
    different genes have weight .25, independently of their native loop similarity.
    This is an uncalibrated component: no new P-value is assigned.
    Reference junctions retain their first source-order V representative. N includes
    all those representatives; unavailable genes contribute no extension edges.
    """
    from .tcrdist import distance_matrix,load_v_loops,_cdr3_options,resolve_v_alleles
    if isinstance(threads,bool) or not isinstance(threads,int) or threads<1:
        raise ValueError('threads must be positive')
    if isinstance(radius,bool) or not isinstance(radius,int) or radius<0 or not np.isfinite(temperature) or temperature<=0:
        raise ValueError('nonnegative integer radius and finite positive temperature required')
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('query_id must be unique and non-null')
    model=load_v_loops(species)
    def prepare(frame):
        return frame.with_columns(pl.col('cdr3','v').cast(pl.String)).with_columns(
            resolve_v_alleles(frame['v'],model).alias('_allele')).with_columns(
                (pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]{8,}$') &
                 pl.col('_allele').is_in(model)).fill_null(False).alias('_available'))
    ref=prepare(reference.unique('cdr3',keep='first',maintain_order=True))
    n_reference=ref.height
    ref=ref.filter(pl.col('_available'))
    loci={model[v]['locus'] for v in ref['_allele']}
    if len(loci)>1:
        raise ValueError('gapped extension reference must contain one locus')
    q=prepare(queries)
    if loci:
        q=q.with_columns((pl.col('_available') & pl.col('_allele').is_in(
            [v for v,m in model.items() if m['locus'] in loci])).alias('_available'))
    full_control=list(dict.fromkeys(control.ref_seqs()))
    full_control=[s for s in full_control if len(s)>=8]
    m_gap=len(full_control)
    if not m_gap:
        raise ValueError('gapped extension requires nonempty raw controls with junction length>=8')
    full_control_set=set(full_control)
    # Do not deduplicate these trimmed strings: distinct full junctions are observations.
    trimmed_control=[s[3:-2] for s in full_control]
    usable=q.filter(pl.col('_available')).select('cdr3','_allele').unique(maintain_order=True)
    diagnostics={'gapped_density':pl.Float64,
        'gapped_density_same_v_unweighted':pl.Float64,'gapped_density_cross_v_unweighted':pl.Float64,
        'gapped_edges_same_v':pl.UInt64,'gapped_edges_cross_v':pl.UInt64,
        'gapped_floor_density':pl.Float64,'gapped_best_total_distance':pl.Int32,
        'gapped_best_cdr3_distance':pl.Int32,'gapped_best_vloop_distance':pl.Int32}
    result=usable.with_columns(*[pl.lit(-1 if 'best_' in c else 0,dtype=dtype).alias(c)
                                for c,dtype in diagnostics.items()])
    if ref.height and usable.height:
        rs,rv=ref['cdr3'].to_list(),ref['_allele'].to_list()
        rlen=np.asarray(list(map(len,rs)))
        # Same normalization as regions.gene_family: strip allele, retain the gene.
        rgenes=ref['_allele'].str.replace(r'\*.*$', '').to_numpy()
        # Reserve component/total matrices, masks and per-length vectorized Hamming work.
        budget=32*1024**2
        if 32*ref.height>budget:
            raise ValueError('gapped reference exceeds the32MiB work budget for one query')
        batch=max(1,budget//(32*ref.height))
        values={c:[] for c in diagnostics}
        options=_cdr3_options(threads)
        for part in usable.iter_slices(batch):
            seqs=part['cdr3'].to_list()
            qgenes=part['_allele'].str.replace(r'\*.*$', '').to_numpy()
            total,cdr3=distance_matrix(seqs,rs,part['_allele'].to_list(),rv,
                species=species,threads=threads,return_cdr3=True)
            new=total<=radius
            qlen=np.asarray(list(map(len,seqs)))
            for length in np.unique(qlen):
                qi,ri=np.flatnonzero(qlen==length),np.flatnonzero(rlen==length)
                if not len(ri) or not new[np.ix_(qi,ri)].any():
                    continue
                qc=np.frombuffer(''.join(seqs[i] for i in qi).encode('ascii'),dtype=np.uint8).reshape(len(qi),int(length))
                rc=np.frombuffer(''.join(rs[i] for i in ri).encode('ascii'),dtype=np.uint8).reshape(len(ri),int(length))
                mismatches=np.zeros((len(qi),len(ri)),dtype=np.uint16)
                for pos in range(int(length)):
                    mismatches+=qc[:,None,pos]!=rc[None,:,pos]
                new[np.ix_(qi,ri)] &= mismatches>5
            active=np.flatnonzero(new.any(axis=1))
            thresholds=[np.unique(2*cdr3[i,new[i]]//3).tolist() for i in active]
            counts=gapblock.count_batch([seqs[i][3:-2] for i in active],trimmed_control,thresholds,
                                       exclude_exact=False,**options) if len(active) else []
            arrays={c:np.full(len(seqs),-1 if 'best_' in c else 0,
                    dtype=np.int32 if 'best_' in c else np.uint64 if 'edges_' in c else np.float64)
                    for c in diagnostics}
            for i,cutoffs,nc in zip(active,thresholds,counts):
                # Puncture the full identity only, preserving other trimmed-identical controls.
                nc=np.asarray(nc,dtype=np.int64)-int(seqs[i] in full_control_set)
                positions=np.searchsorted(cutoffs,2*cdr3[i,new[i]]//3)
                empirical=n_reference/m_gap*nc[positions]
                unweighted=np.exp(-total[i,new[i]]/temperature)/np.maximum(empirical,.01)
                same=qgenes[i]==rgenes[new[i]]
                weighted=unweighted*np.where(same,1.0,.25)
                arrays['gapped_density'][i]=weighted.sum()
                arrays['gapped_density_same_v_unweighted'][i]=unweighted[same].sum()
                arrays['gapped_density_cross_v_unweighted'][i]=unweighted[~same].sum()
                arrays['gapped_edges_same_v'][i]=same.sum()
                arrays['gapped_edges_cross_v'][i]=(~same).sum()
                arrays['gapped_floor_density'][i]=weighted[empirical<.01].sum()
                # argmin selects the first accepted source-order reference on ties.
                accepted=np.flatnonzero(new[i])
                best=accepted[np.argmin(total[i,accepted])]
                arrays['gapped_best_total_distance'][i]=total[i,best]
                arrays['gapped_best_cdr3_distance'][i]=cdr3[i,best]
                arrays['gapped_best_vloop_distance'][i]=total[i,best]-cdr3[i,best]
            for c in diagnostics:
                values[c].extend(arrays[c].tolist())
        result=usable.with_columns(*[pl.Series(c,values[c],dtype=dtype) for c,dtype in diagnostics.items()])
    out=q.join(result,on=['cdr3','_allele'],how='left',validate='m:1',maintain_order='left')
    return out.select('query_id',*[pl.when(pl.col(c)>=0).then(pl.col(c)).otherwise(None).alias(c)
        if 'best_' in c else pl.col(c).fill_null(0).alias(c) for c in diagnostics],
        (pl.col('_available') & pl.lit(bool(ref.height))).alias('availability'),
        pl.when(~pl.col('_available')).then(pl.lit('unavailable_query_or_gene'))
        .when(pl.lit(not bool(ref.height))).then(pl.lit('no_usable_reference'))
        .otherwise(pl.lit('available')).alias('gapped_status'),
        pl.lit(m_gap,dtype=pl.UInt64).alias('M_gap'),
        pl.lit(n_reference-ref.height,dtype=pl.UInt64).alias('n_reference_unavailable'))


def require_historical_pssm_kernel():
    """Fail before raw-control loading when either Python or native support is stale."""
    import inspect
    from seqtree import _core
    keyword='position_weights_by_length'
    if any(keyword not in inspect.signature(fn).parameters
           for fn in (gapblock.score_matrix,gapblock.count_batch)) or any(
            keyword not in (fn.__doc__ or '') for fn in (_core.gapblock_matrix,_core.gapblock_count_batch)):
        raise RuntimeError('historical-pssm geometry requires positional seqtree gapblock support; rebuild/install the reviewed native source')


def historical_pssm_extension(queries, reference, control, *, species='human', threads=1, matched_background=False, kernel_scale=400.0):
    """Opt-in original-position PSSM kernel on edges outside its substitution ball.

    Use full junctions and the native symmetric longer-sequence frame. Pin the
    gap after six matched prefix residues (Cys at index0), clamped to the
    shorter junction length. Equal-length pairs have no gap. Controls
    count the identical weighted kernel at each accepted penalty. With
    matched_background=True, include original edges under that same predicate/CDF
    instead of adding this component to the original density. Original density
    and its P-value remain separate; this component has no assigned P-value.
    """
    from .regions import significance_weights
    from .vgene import vsim
    require_historical_pssm_kernel()
    if not isinstance(matched_background,bool):
        raise ValueError('matched_background must be boolean')
    if isinstance(kernel_scale,bool) or not isinstance(kernel_scale,(int,float)) or not np.isfinite(kernel_scale) or kernel_scale<=0:
        raise ValueError('kernel_scale must be finite and positive')
    if isinstance(threads,bool) or not isinstance(threads,int) or threads<1:
        raise ValueError('threads must be positive')
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('query_id must be unique and non-null')
    def prepare(frame):
        return frame.with_columns(pl.col('cdr3','v').cast(pl.String)).with_columns(
            pl.col('v').str.replace(r'\*.*$', '').fill_null('').alias('_gene'),
            pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).alias('_available'))
    q=prepare(queries)
    ref=prepare(reference.unique('cdr3',keep='first',maintain_order=True))
    n_reference=ref.height
    ref=ref.filter(pl.col('_available'))
    full_control=list(dict.fromkeys(control.ref_seqs()))
    m_gap=len(full_control)
    if not m_gap:
        raise ValueError('historical PSSM extension requires nonempty raw controls')
    usable=q.filter(pl.col('_available')).select('cdr3','_gene').unique(maintain_order=True)
    numeric={'gapped_density':pl.Float64,'gapped_density_same_v_unweighted':pl.Float64,
        'gapped_density_cross_v_unweighted':pl.Float64,'gapped_edges_same_v':pl.UInt64,
        'gapped_edges_cross_v':pl.UInt64,'gapped_floor_density':pl.Float64}
    details={'gapped_best_kernel_penalty':pl.Int32,'gapped_best_reference_junction':pl.String,
        'gapped_best_reference_v':pl.String,'gapped_best_contribution':pl.Float64,
        'gapped_best_control_count':pl.UInt64,'gapped_best_gap_length':pl.UInt32,
        'gapped_best_block_position':pl.Int32,'gapped_best_alignment_status':pl.String,
        'gapped_best_total_distance':pl.Int32,'gapped_best_cdr3_distance':pl.Int32,
        'gapped_best_vloop_distance':pl.Int32}
    schema={**numeric,**details}
    result=usable.with_columns(*[pl.lit(0,dtype=dtype).alias(c) for c,dtype in numeric.items()],
        *[pl.lit('no_extension_edge' if c=='gapped_best_alignment_status' else None,dtype=dtype).alias(c)
          for c,dtype in details.items()])
    if ref.height and usable.height:
        matrix=SubstitutionMatrix.blosum62()
        rs,rv=ref['cdr3'].to_list(),ref['v'].to_list()
        rlen=np.asarray(list(map(len,rs)))
        rgenes=ref['_gene'].to_numpy()
        genes,codes=np.unique(rgenes,return_inverse=True)
        lengths={len(s) for s in [*usable['cdr3'].to_list(),*rs,*full_control]}
        weights={length:[max(1,round(100*w)) for w in significance_weights(length)] for length in lengths}
        options=dict(matrix=matrix,gap_open=2*matrix.scale()*100,gap_extend=matrix.scale()*100,
            gap_prior=gapblock.positions_prior((6,)),threads=threads,
            position_weights_by_length=weights)
        cutoff=5*matrix.scale()*100
        budget=32*1024**2
        if 32*ref.height>budget:
            raise ValueError('gapped reference exceeds the32MiB work budget for one query')
        batch=max(1,budget//(32*ref.height))
        values={c:[] for c in schema}
        for part in usable.iter_slices(batch):
            seqs=part['cdr3'].to_list()
            qgenes=part['_gene'].to_list()
            # V factors mirror the original component, independently of sequence geometry.
            priors={gene:np.asarray([1.0 if gene==r else .25*vsim(gene,r) for r in genes])[codes]
                    for gene in set(qgenes)}
            penalties=np.asarray(gapblock.score_matrix(seqs,rs,**options))
            new=(penalties<=cutoff) & (np.asarray(seqs)[:,None]!=np.asarray(rs)[None,:])
            if not matched_background:
                qlen=np.asarray(list(map(len,seqs)))
                for length in np.unique(qlen):
                    qi,ri=np.flatnonzero(qlen==length),np.flatnonzero(rlen==length)
                    if not len(ri) or not new[np.ix_(qi,ri)].any():
                        continue
                    qc=np.frombuffer(''.join(seqs[i] for i in qi).encode('ascii'),dtype=np.uint8).reshape(len(qi),int(length))
                    rc=np.frombuffer(''.join(rs[i] for i in ri).encode('ascii'),dtype=np.uint8).reshape(len(ri),int(length))
                    mismatches=np.zeros((len(qi),len(ri)),dtype=np.uint16)
                    for pos in range(int(length)):
                        mismatches+=qc[:,None,pos]!=rc[None,:,pos]
                    # This also excludes every exact full-junction reference.
                    new[np.ix_(qi,ri)] &= mismatches>5
            for i,gene in enumerate(qgenes):
                new[i] &= priors[gene]>0
            active=np.flatnonzero(new.any(axis=1))
            thresholds=[np.unique(penalties[i,new[i]]).tolist() for i in active]
            counts=gapblock.count_batch([seqs[i] for i in active],full_control,thresholds,
                exclude_exact=True,**options) if len(active) else []
            rows=[{**dict.fromkeys(numeric,0),**dict.fromkeys(details,None),
                   'gapped_best_alignment_status':'no_extension_edge'} for _ in seqs]
            for i,cutoffs,nc in zip(active,thresholds,counts):
                accepted=np.flatnonzero(new[i])
                positions=np.searchsorted(cutoffs,penalties[i,accepted])
                nc=np.asarray(nc,dtype=np.uint64)[positions]
                empirical=n_reference/m_gap*nc
                unweighted=np.exp(-penalties[i,accepted]/kernel_scale)/np.maximum(empirical,.01)
                same=qgenes[i]==rgenes[accepted]
                weighted=unweighted*priors[qgenes[i]][accepted]
                row=rows[i]
                row.update(gapped_density=float(weighted.sum()),
                    gapped_density_same_v_unweighted=float(unweighted[same].sum()),
                    gapped_density_cross_v_unweighted=float(unweighted[~same].sum()),
                    gapped_edges_same_v=int(same.sum()),gapped_edges_cross_v=int((~same).sum()),
                    gapped_floor_density=float(weighted[empirical<.01].sum()))
                # Highest contribution, with first source-order representative on ties.
                best_pos=int(np.argmax(weighted));best=int(accepted[best_pos])
                row.update(gapped_best_kernel_penalty=int(penalties[i,best]),
                    gapped_best_reference_junction=rs[best],gapped_best_reference_v=rv[best],
                    gapped_best_contribution=float(weighted[best_pos]),
                    gapped_best_control_count=int(nc[best_pos]),
                    gapped_best_gap_length=abs(len(seqs[i])-len(rs[best])),
                    gapped_best_alignment_status='not_returned_by_native')
            for c in schema:
                values[c].extend(row[c] for row in rows)
        result=usable.with_columns(*[pl.Series(c,values[c],dtype=dtype) for c,dtype in schema.items()])
    out=q.join(result,on=['cdr3','_gene'],how='left',validate='m:1',maintain_order='left')
    return out.select('query_id',*[pl.col(c).fill_null(0).alias(c) for c in numeric],
        pl.when(pl.col('gapped_best_reference_junction').is_not_null()).then(pl.col('v'))
        .otherwise(None).alias('gapped_best_query_v'),
        *[pl.col(c) for c in details],
        (pl.col('_available') & pl.lit(bool(ref.height))).alias('availability'),
        pl.when(~pl.col('_available')).then(pl.lit('unavailable_query_sequence'))
        .when(pl.lit(not bool(ref.height))).then(pl.lit('no_usable_reference'))
        .otherwise(pl.lit('available')).alias('gapped_status'),
        pl.lit(m_gap,dtype=pl.UInt64).alias('M_gap'),
        pl.lit(n_reference-ref.height,dtype=pl.UInt64).alias('n_reference_unavailable'))


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
        if 4 * len(refseqs) > 64 * 1024**2:
            raise ValueError('gapped reference exceeds the64MiB matrix budget for one query')
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
