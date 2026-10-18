"""Uncalibrated kernel-ladder evidence over cumulative distance-ball counts."""
from __future__ import annotations

_MAX_RETAINED_DETAIL_BYTES = 8_000_000_000

import numpy as np

KERNEL_SCORERS = ('v-bounded', 'all-bounded', 'vonly-bounded',
                  'v-rare', 'all-rare', 'vonly-rare')


def ladder_cutoffs(matrix_scale, *, radius_factor=22.5, reach_factor=100.0):
    """Recorded 90 inner/62 outer intervals, in ordinary matrix cost units."""
    if not np.isfinite(matrix_scale) or matrix_scale <= 0:
        raise ValueError('matrix scale must be finite and positive')
    if not 0 < radius_factor <= reach_factor:
        raise ValueError('positive radius factor must not exceed reach factor')
    radius, reach = round(radius_factor*matrix_scale), round(reach_factor*matrix_scale)
    return sorted({round(k*radius/90) for k in range(91)} | {
        round(radius+k*(reach-radius)/62) for k in range(1,63)})


def kernel_ladder(reference_counts, same_v_counts, control_counts, reference_sizes,
                  control_sizes, cutoffs, *, kernel_scale, cross_v=.25, floor=.01,
                  scorer='v-bounded', extras=True):
    """Recorded V-weighted bounded score plus five advanced ladder alternatives.

    Each curve counts all eligible references through each ordered cutoff. Exact
    exclusions mask hits without changing reference/control model sizes. Scores are relative evidence,
    not probabilities, E-values or calibrated specificity confidence.
    """
    if scorer not in KERNEL_SCORERS:
        raise ValueError('unknown kernel scorer')
    n, same, m = (np.asarray(x,dtype=np.float64,order='C')
                  for x in (reference_counts,same_v_counts,control_counts))
    N,M,t = (np.asarray(x,dtype=np.float64)
             for x in (reference_sizes,control_sizes,cutoffs))
    if (n.ndim!=2 or same.shape!=n.shape or m.shape!=n.shape or
            N.shape!=(n.shape[0],) or M.shape!=N.shape or t.shape!=(n.shape[1],)):
        raise ValueError('kernel curves and population axes must agree')
    if (not np.isfinite(kernel_scale) or kernel_scale<=0 or not 0<=cross_v<=1
            or not np.isfinite(floor) or floor<=0 or not len(t)
            or np.any(np.diff(t)<=0) or np.any(t<0)
            or any(not np.isfinite(x).all() or np.any(x<0) for x in (n,same,m,N,M,t))
            or any(np.any(np.diff(x,axis=1)<0) for x in (n,same,m))
            or np.any(same>n) or np.any(n>N[:,None]) or np.any(m>M[:,None])):
        raise ValueError('finite monotone counts within populations and positive kernel parameters required')
    shells=np.diff(n,axis=1,prepend=0)
    same_shells=np.diff(same,axis=1,prepend=0)
    if np.any(same_shells>shells):
        raise ValueError('same-V shells must be subsets of reference shells')
    return _kernel_scores(shells,same_shells,m,N,M,t,kernel_scale,cross_v,floor,scorer,extras)


def _kernel_scores(shells,same_shells,m,N,M,t,kernel_scale,cross_v,floor,scorer,extras, *, _reuse_buffers=False):
    """Shared arithmetic after the public or native counter's validation."""
    kernel=np.exp(-t/kernel_scale)
    usable=(N>0)&(M>0)
    if _reuse_buffers and not extras:
        # Only the compact native path supplies owned writable float buffers.
        name,mode=scorer.split('-')
        mass=shells if name=='all' else same_shells
        if name=='v':
            np.multiply(shells,cross_v,out=shells)
            np.multiply(same_shells,1-cross_v,out=same_shells)
            np.add(shells,same_shells,out=shells)
            mass=shells
        if mode=='bounded':np.add(1,m,out=m)
        else:np.maximum(m,floor,out=m)
        np.divide(kernel[None,:],m,out=m)
        score=np.einsum('ij,ij->i',mass,m)/np.where(N>0,N,1)
        score=np.where(usable,score,np.nan)
        return {'kernel_'+scorer.replace('-','_'):score,'kernel_score':score}
    out={}
    for name in ('all','v','vonly'):
        if not extras and name!=scorer.split('-')[0]:
            continue
        mass=(shells if name=='all' else same_shells if name=='vonly' else
              cross_v*shells+(1-cross_v)*same_shells)
        for mode in ('bounded','rare'):
            if not extras and name+'-'+mode!=scorer:continue
            denom=1+m if mode=='bounded' else np.maximum(m,floor)
            score=np.einsum('ij,ij->i',mass,kernel[None,:]/denom)/np.where(N>0,N,1)
            out['kernel_'+name+'_'+mode]=np.where(usable,score,np.nan)
    out['kernel_score']=out['kernel_'+scorer.replace('-','_')]
    return out


def kernel_peptides(scores):
    """Maximum finite kernel evidence and exact support per query/peptide."""
    import polars as pl
    scores=scores.with_columns(pl.when(pl.col('kernel_score').is_finite())
        .then(pl.col('kernel_score')).otherwise(None).alias('kernel_score'))
    exact = pl.col('exact_reference_count').max() if 'exact_reference_count' in scores.columns else pl.lit(0)
    return scores.group_by('query_id','epitope').agg(pl.col('kernel_score').max(),
        exact.alias('exact_reference_count'))



def _compact_kernel_selection(scores, exact, permutation, starts, peptides):
    """Maxima only: restrictions first, then lexical winner and strongest rival."""
    finite=np.where(np.isfinite(scores),scores,-np.inf)
    per_peptide=np.maximum.reduceat(finite[:,permutation],starts,axis=1)
    exact_peptides=np.maximum.reduceat(exact[:,permutation],starts,axis=1)
    winner=np.argmax(per_peptide,axis=1)
    best=per_peptide[np.arange(len(winner)),winner]
    if len(peptides)>1:
        rivals=per_peptide.copy()
        rivals[np.arange(len(winner)),winner]=-np.inf
        runner=np.argmax(rivals,axis=1)
        second=rivals[np.arange(len(winner)),runner]
    else:
        second=np.full(len(winner),-np.inf)
    labels=np.asarray(peptides,dtype=object)[winner].copy()
    labels[best<=0]=None
    return dict(selected_epitope=labels,kernel_score=np.where(best>0,best,np.nan),
        kernel_competing_score=np.where(np.isfinite(second),second,np.nan),
        tied_peptides=np.sum((per_peptide==best[:,None])&np.isfinite(per_peptide),axis=1,dtype=np.uint32),
        exact_peptides=np.sum(exact_peptides>0,axis=1,dtype=np.uint32))

def grouped_kernel_selection(queries, references, groups, *, query_keys, search_keys,
        paired, cutoffs, radius, kernel_scale, scorer, assayed=False, withheld_column=None,
        position_weighting='uniform', geometry_preset=None, _detailed=False, progress=False,
        _allow_empty_background=False, _scores_only=False, _consume=None, _return_hits=False, **geometry):
    """Compact all-label selection from one bounded native reference traversal.

    Memberships deduplicate search keys separately from pMHC annotations. Other
    peptide counts subtract only keys exclusive to that peptide; shared labels
    and unselected peptides remain in the background. No predictions persist.
    The complete-score runner retains only pMHC keys and the selected score;
    detailed diagnostics remain available to other callers. An internal consumer
    can receive each detailed observation batch directly, without retaining the
    complete table or starting another native search. Consumer batches include
    the integer execution axes ``_target`` and ``_query_order`` for keyed writes.
    Detailed mode can also return radius-hit IDs for query dispositions.
    """
    import polars as pl
    from ..match.geometry import grouped_total_distance_count_batch, distance_cutoff
    if scorer not in KERNEL_SCORERS:
        raise ValueError('unknown kernel scorer')
    if _scores_only and not _detailed:
        raise ValueError('score projection requires detailed kernel mode')
    if _consume is not None and (not _detailed or not callable(_consume)):
        raise ValueError('kernel consumer requires detailed mode and a callable')
    if not groups or not queries.height:
        return None,set()
    if progress:
        import sys,time
        progress_started=time.monotonic()
    unique=queries.select(query_keys).unique(maintain_order=True).with_row_index('_row')
    whole=references.unique(search_keys,maintain_order=True).with_row_index('_ref')
    if progress:
        print(f'Kernel batch: {unique.height} unique query keys, {whole.height} reference keys; '
              f'{unique.height*whole.height} reference-key visits',file=sys.stderr,flush=True)
    normalise=geometry_preset in ('norm','jnorm')
    same_v_totals=None
    if normalise:
        # Same-V mass and all population denominators belong to the fixed model.
        gene_columns=['_allelea','_alleleb'] if paired else ['_allele']
        reference_genes=['_allele_alpha','_allele_beta'] if paired else ['_allele']
        gene_names=['_gene'+str(i) for i in range(len(gene_columns))]
        def gene_axes(frame,columns):
            return frame.with_columns(*[pl.col(c).str.split('*').list.first().alias(g)
                for c,g in zip(columns,gene_names)])
        unique=gene_axes(unique,gene_columns)
        same_v_totals=pl.concat([gene_axes(part.unique(search_keys),reference_genes)
            .group_by(gene_names).len().with_columns(pl.lit(i).alias('_target'))
            for i,(_,_,part) in enumerate(groups)])
    targets=pl.DataFrame({'_target':range(len(groups)),
                         'epitope':[task['epitope'] for task,_,_ in groups]})
    if _detailed:
        from ..aggregate.candidates import PMHC
        targets=targets.with_columns(*[pl.Series(c,[task[c] for task,_,_ in groups],dtype=pl.String)
            for c in PMHC if c!='epitope'])
        detail_columns=(['query_id',*PMHC,'kernel_score'] if _scores_only else
            ['query_id','n_control','control_weight','control_population','n_reference',
             'reference_weight','reference_population','expected_count','enrichment_status',
             'enrichment','zero_control_hits','zero_control_weight','exact_reference_count',
             *['kernel_'+w+'_'+d for w in ('all','v','vonly') for d in ('bounded','rare')],
             'kernel_score','kernel_status',*PMHC])
    peptides=targets.select('epitope').unique(maintain_order=True).with_row_index('_peptide')
    target_peptides=targets.join(peptides,on='epitope',nulls_equal=True,maintain_order='left')['_peptide'].to_numpy()
    members=pl.concat([part.select(search_keys).unique().with_columns(pl.lit(i).alias('_target'))
                      for i,(_,_,part) in enumerate(groups)]).join(
        whole.select('_ref',*search_keys),on=search_keys,how='left',validate='m:1')
    members=members.sort('_ref','_target')
    lengths=members.group_by('_ref').len()
    lengths=whole.select('_ref').join(lengths,on='_ref',how='left',maintain_order='left')['len'].fill_null(0).to_numpy()
    offsets=np.concatenate((np.zeros(1,dtype=np.uint64),np.cumsum(lengths,dtype=np.uint64)))
    exclusive=(references.group_by(search_keys).agg(pl.col('epitope').n_unique().alias('_labels'),
        pl.col('epitope').first()).filter(pl.col('_labels')==1).join(peptides,on='epitope',how='left'))
    exclusive=whole.select('_ref',*search_keys).join(exclusive,on=search_keys,how='left',
        validate='1:1',maintain_order='left')['_peptide'].fill_null(-1).to_numpy()
    folds={}
    if withheld_column is not None:
        if withheld_column not in query_keys:
            raise ValueError('withheld peptide must be included in query deduplication keys')
        # Selected peptides keep the count-channel namespace; sparse labels follow.
        full_peptides=pl.concat([peptides.select('epitope'),references.select('epitope')]).unique(
            maintain_order=True).with_row_index('_peptide')
        labels=references.select(*search_keys,'epitope').unique().join(
            whole.select('_ref',*search_keys),on=search_keys,validate='m:1').join(
            full_peptides,on='epitope',nulls_equal=True,validate='m:1').sort('_ref','_peptide')
        lengths=whole.select('_ref').join(labels.group_by('_ref').len(),on='_ref',how='left',
            maintain_order='left')['len'].fill_null(0).to_numpy()
        withheld=unique.select(pl.col(withheld_column).alias('epitope')).join(
            full_peptides,on='epitope',how='left',nulls_equal=True,maintain_order='left')
        folds=dict(query_withheld_peptide_ids=withheld['_peptide'].fill_null(-1).to_numpy(),
            reference_peptide_offsets=np.concatenate((np.zeros(1,dtype=np.uint64),np.cumsum(lengths,dtype=np.uint64))),
            reference_peptide_ids=labels['_peptide'].to_numpy())
    def axes(frame, reference=False):
        seqs=['alpha','beta'] if reference and paired else ['cdr3a','cdr3b'] if paired else ['cdr3']
        genes=['_allele_alpha','_allele_beta'] if reference and paired else ['_allelea','_alleleb'] if paired else ['_allele']
        from ..match.geometry import J_GEOMETRY_PRESETS
        jcols=['j_alpha','j_beta'] if reference and paired else ['ja','jb'] if paired else ['j']
        return frame.select(*[pl.col(c).alias('seq'+str(i)) for i,c in enumerate(seqs)],
                            *[pl.col(c).alias('v'+str(i)) for i,c in enumerate(genes)],
                            *([pl.col(c).alias('j'+str(i)) for i,c in enumerate(jcols)]
                              if geometry_preset in J_GEOMETRY_PRESETS else []))
    native_cuts=sorted(set([int(200*t) for t in cutoffs]+[distance_cutoff(radius)]))
    ladder_indices=np.searchsorted(native_cuts,np.asarray(cutoffs)*200)
    ladder_slice=slice(None) if len(native_cuts)==len(cutoffs) else ladder_indices
    radius_index=native_cuts.index(distance_cutoff(radius))
    channels=1+2*len(groups)+peptides.height
    bytes_per_query=channels*len(native_cuts)*np.dtype(np.uint32).itemsize
    if bytes_per_query>256*1024**2:
        raise ValueError('one grouped query exceeds the 256 MiB native count-buffer budget; restrict requested targets')
    retained_bytes=0
    iterator=grouped_total_distance_count_batch(axes(unique),axes(whole,True),native_cuts,
        paired=paired,position_weighting=position_weighting,geometry_preset=geometry_preset,reference_target_offsets=offsets,
        reference_target_ids=members['_target'].to_numpy(),target_peptide_ids=target_peptides,
        reference_exclusive_peptide_ids=exclusive,
        chunk_size=min(65536,max(1,256*1024**2//bytes_per_query)),**folds,**geometry)
    mapping=queries.select('query_id',*query_keys).join(unique,on=query_keys,how='left',validate='m:1',nulls_equal=True,
        **({'maintain_order':'left'} if _detailed else {}))
    if _detailed:
        mapping=mapping.with_row_index('_query_order')
    # Native chunks span contiguous search-key rows. Sort once, then slice;
    # rescanning the complete source map per chunk is quadratic in cohort size.
    mapping=mapping.sort('_row')
    mapped_rows=mapping['_row'].to_numpy()
    allowed=None
    if assayed:
        allowed=pl.concat([part.select('query_id').with_columns(pl.lit(i).alias('_target'))
                          for i,(_,part,_) in enumerate(groups)])
    sparse=_detailed and assayed
    requested=(mapping.select('query_id','_row').join(allowed,on='query_id',how='inner')
        .select('_row','_target').unique().sort('_row','_target') if sparse else None)
    requested_rows=requested['_row'].to_numpy() if sparse else None
    states=[];hit_ids=set();T=len(groups)
    shared_single_scope=(not _detailed and allowed is not None and T==1 and allowed.height==queries.height and
        allowed['query_id'].n_unique()==queries.height and
        allowed.select('query_id').join(queries.select('query_id'),on='query_id',how='anti').is_empty())
    compact_all=not _detailed and (allowed is None or shared_single_scope)
    if compact_all:
        # The same peptide restrictions and lexical order apply to every chunk.
        names=[task['epitope'] for task,_,_ in groups]
        permutation=np.asarray(sorted(range(T),key=lambda i:(names[i] is None,names[i] or '')))
        ordered=[names[i] for i in permutation]
        starts=np.asarray([i for i,e in enumerate(ordered) if i==0 or e!=ordered[i-1]])
        peptide_names=[ordered[i] for i in starts]
        mapped_ids=mapping['query_id'].to_numpy()
    reported=-1
    for chunk in iterator:
        C=chunk.stop-chunk.start
        milestone=10*chunk.stop//unique.height
        if progress and (reported<0 or milestone>reported):
            print(f'Kernel counted {chunk.stop}/{unique.height} query keys '
                  f'in {time.monotonic()-progress_started:.1f}s',file=sys.stderr,flush=True)
            reported=milestone
        first,last=np.searchsorted(mapped_rows,[chunk.start,chunk.stop])
        chunk_mapping=mapping.slice(int(first),int(last-first))
        curves=np.asarray(chunk.counts);sizes=np.asarray(chunk.sizes)
        # Convert directly to contiguous numerical rows; advanced-cut indexing
        # otherwise transposes and recopies every query x target x cutoff array.
        if sparse:
            first,last=np.searchsorted(requested_rows,[chunk.start,chunk.stop])
            pairs=requested.slice(int(first),int(last-first))
            if not pairs.height:
                continue
            row=pairs['_row'].to_numpy();target=pairs['_target'].to_numpy();qrow=row-chunk.start
            reference=curves[qrow,1+target,:][:,ladder_slice]
            same=curves[qrow,1+T+target,:][:,ladder_slice]
            shells=np.empty(reference.shape,dtype=np.float64);same_shells=np.empty_like(shells)
            for source,out in ((reference,shells),(same,same_shells)):
                out[:,0]=source[:,0]
                out[:,1:]=np.subtract(source[:,1:],source[:,:-1],dtype=np.uint32)
            control=np.subtract(curves[qrow,0,:][:,ladder_slice],
                curves[qrow,1+2*T+target_peptides[target],:][:,ladder_slice],dtype=np.float64)
            control=np.asarray(control,order='C')
            N=sizes[qrow,1+target];M=sizes[qrow,0]-sizes[qrow,1+2*T+target_peptides[target]]
        else:
            shells=np.empty((C,T,len(cutoffs)),dtype=np.float64)
            same_shells=np.empty_like(shells)
            for source,out in ((curves[:,1:1+T,:][:,:,ladder_slice],shells),
                               (curves[:,1+T:1+2*T,:][:,:,ladder_slice],same_shells)):
                out[:,:,0]=source[:,:,0]
                # Monotone uint32 curves have exact nonnegative uint32 increments.
                out[:,:,1:]=np.subtract(source[:,:,1:],source[:,:,:-1],dtype=np.uint32)
            shells=shells.reshape(C*T,-1);same_shells=same_shells.reshape(C*T,-1)
            control=np.empty((C,T,len(native_cuts)),dtype=np.float64)
            np.subtract(curves[:,0,None,:],np.take(curves,1+2*T+target_peptides,axis=1),out=control,dtype=np.float64)
            control=np.asarray(control[:,:,ladder_slice],order='C').reshape(C*T,-1)
            N=sizes[:,1:1+T];M=sizes[:,0,None]-sizes[:,1+2*T+target_peptides]
            if not compact_all:
                row=np.repeat(np.arange(chunk.start,chunk.stop),T);target=np.tile(np.arange(T),C)
        score_N=N.reshape(-1)
        if normalise:
            rows=(pairs.join(unique.select('_row',*gene_names),on='_row',how='left') if sparse else
                unique.slice(chunk.start,C).select('_row',*gene_names).join(targets.select('_target'),how='cross')).join(same_v_totals,
                on=['_target',*gene_names],how='left',validate='m:1').sort('_row','_target')
            same_total=rows['len'].fill_null(0).to_numpy()
            weight=scorer.split('-')[0]
            score_N=(score_N if weight=='all' else same_total if weight=='vonly' else
                     .25*score_N+.75*same_total)
            score_N=np.where(score_N>0,score_N,1)
        # Immutable native curves already satisfy the validated monotone/subset
        # and population contract. Reuse the exact arithmetic without scanning
        # every generated curve again; public kernel_ladder keeps all checks.
        # The complete-score adapter can use the defined 1+m bounded density
        # at M=0; its final control randomizations calibrate this case. Reported
        # population sizes below remain actual zero, with no enrichment claim.
        scores=_kernel_scores(shells,same_shells,
            control,score_N,(np.maximum(M.reshape(-1),1) if normalise or _allow_empty_background else M.reshape(-1)),np.asarray(cutoffs,dtype=np.float64),
            kernel_scale,.25,.01,scorer,_detailed,_reuse_buffers=not _detailed)
        exact=np.zeros((C,T),dtype=np.uint32) if geometry.get('exclude_exact',False) else np.asarray(chunk.exact_counts)
        if compact_all:
            selected=_compact_kernel_selection(scores['kernel_score'].reshape(C,T),exact,
                permutation,starts,peptide_names)
            local_rows=mapped_rows[first:last]-chunk.start
            ids=mapped_ids[first:last]
            states.append(pl.DataFrame({'query_id':pl.Series(ids.tolist(),dtype=mapping.schema['query_id']),
                'selected_epitope':pl.Series(selected['selected_epitope'][local_rows].tolist(),dtype=pl.String),
                **{name:pl.Series(value[local_rows],nan_to_null=True) for name,value in selected.items()
                   if name!='selected_epitope'}}))
            supported=np.any(curves[:,1:1+T,radius_index]>0,axis=1)
            hit_ids.update(ids[supported[local_rows]].tolist())
            continue
        values=pl.DataFrame({'_row':row,'_target':target,'kernel_score':scores['kernel_score'],
            'exact_reference_count':exact[qrow,target] if sparse else exact.reshape(-1),
            '_radius_hit':(curves[qrow,1+target,radius_index]>0 if sparse else
                (curves[:,1:1+T,radius_index]>0).reshape(-1))})
        if _detailed:
            # Counts and weights in the established detailed table use the final
            # ladder cutoff; candidate radius diagnostics are a separate output.
            reference_count=(curves[qrow,1+target,ladder_indices[-1]] if sparse else
                curves[:,1:1+T,ladder_indices[-1]].reshape(-1))
            control_count=control[:,-1]
            values=values.with_columns(
                pl.Series('n_control',control_count,dtype=pl.UInt64),
                pl.Series('control_weight',control_count,dtype=pl.Float64),
                pl.Series('control_population',M.reshape(-1),dtype=pl.UInt64),
                pl.Series('n_reference',reference_count,dtype=pl.UInt64),
                pl.Series('reference_weight',reference_count,dtype=pl.Float64),
                pl.Series('reference_population',N.reshape(-1),dtype=pl.UInt64),
                *(pl.Series(k,v,nan_to_null=True) for k,v in scores.items()),
                pl.col('exact_reference_count').cast(pl.UInt64),
                pl.lit('relative_evidence_uncalibrated').alias('kernel_status'))
            available=(pl.col('reference_population')>0)&(pl.col('control_population')>0)
            values=values.with_columns(
                pl.when(available).then(pl.col('reference_population')*(pl.col('control_weight')+.5)/
                    (pl.col('control_population')+1)).alias('expected_count'),
                pl.when(available).then(pl.lit('available')).otherwise(pl.lit('empty_population')).alias('enrichment_status'))
            values=values.with_columns((pl.col('reference_weight')/pl.col('expected_count')).alias('enrichment'),
                (pl.col('n_control')==0).alias('zero_control_hits'),
                (pl.col('control_weight')==0).alias('zero_control_weight'))
            observations=chunk_mapping.select(
                'query_id','_query_order','_row')
            for part in observations.iter_slices(max(1,65536//T)):
                expanded=(part.join(allowed,on='query_id',how='inner').join(values,on=['_row','_target'],how='inner')
                    if assayed else part.join(values,on='_row',how='inner'))
                expanded=expanded.join(targets,on='_target',how='left',validate='m:1')
                if _return_hits:
                    hit_ids.update(expanded.filter(pl.col('_radius_hit'))['query_id'])
                if _consume is None:
                    retained=expanded.select('query_id',*PMHC,'kernel_score','_target','_query_order') if _scores_only else expanded.drop('_row','_radius_hit')
                    retained_bytes += retained.estimated_size()
                    if retained_bytes > _MAX_RETAINED_DETAIL_BYTES:
                        raise MemoryError('retained kernel detail exceeds8GB; use --stream-evidence or a batch consumer')
                    states.append(retained)
                else:
                    _consume(expanded.select(*detail_columns,'_target','_query_order'))
            continue
        chunk_mapping=chunk_mapping.select('query_id','_row')
        columns=['selected_epitope','kernel_score','kernel_competing_score','tied_peptides','exact_peptides']
        # Source observations can repeat one search key arbitrarily often.
        # Restrict assay memberships before joining curves; bound even the
        # case where an assay deliberately requests every target.
        for observations in chunk_mapping.iter_slices(max(1,65536//T)):
            expanded=observations.join(allowed,on='query_id',how='inner').join(
                values,on=['_row','_target'],how='inner').join(targets,on='_target',how='left',validate='m:1').filter(
                    pl.col('epitope').is_not_null())
            hit_ids.update(expanded.filter(pl.col('_radius_hit'))['query_id'])
            selected=select_kernel(expanded.select('query_id').unique(),expanded)
            states.append(selected.select('query_id',*columns))
    if _detailed and states:
        return pl.concat(states).sort('_target','_query_order').select(detail_columns),hit_ids
    return pl.concat(states) if states else None,hit_ids


def grouped_kernel_enrichment(queries,references,groups,**kwargs):
    """Detailed continuous evidence using the same bounded grouped native batch."""
    evidence,_=grouped_kernel_selection(queries,references,groups,_detailed=True,**kwargs)
    return evidence


def update_kernel_selection(state, scores):
    """Fold one completed peptide into run-local best/runner-up state.

    Restrictions must be combined before this fold, and each peptide supplied
    once. State has at most one row per query and retains all score/exact ties.
    """
    import polars as pl
    if not scores.height:
        return state
    if scores['epitope'].n_unique()!=1:
        raise ValueError('streaming kernel selection requires one completed peptide')
    peptide=kernel_peptides(scores).rename({'epitope':'_new_epitope',
        'kernel_score':'_new_score','exact_reference_count':'_new_exact'})
    if state is None:
        return peptide.select('query_id',pl.col('_new_epitope').alias('selected_epitope'),
            pl.col('_new_score').alias('kernel_score'),
            pl.lit(None,dtype=pl.Float64).alias('kernel_competing_score'),
            pl.col('_new_score').is_not_null().cast(pl.UInt32).alias('tied_peptides'),
            (pl.col('_new_exact')>0).cast(pl.UInt32).alias('exact_peptides'))
    merged=state.join(peptide,on='query_id',how='full',coalesce=True,validate='1:1')
    has_new=pl.col('_new_epitope').is_not_null()
    first=pl.col('selected_epitope').is_null()
    higher=(pl.col('_new_score').is_not_null() &
        (pl.col('kernel_score').is_null() | (pl.col('_new_score')>pl.col('kernel_score')))).fill_null(False)
    tied=(pl.col('_new_score')==pl.col('kernel_score')).fill_null(False)
    choose=has_new & (first | higher | (tied & (pl.col('_new_epitope')<pl.col('selected_epitope'))))
    return merged.select('query_id',
        pl.when(choose).then(pl.col('_new_epitope')).otherwise(pl.col('selected_epitope')).alias('selected_epitope'),
        pl.when(first | higher).then(pl.col('_new_score')).otherwise(pl.col('kernel_score')).alias('kernel_score'),
        pl.max_horizontal('kernel_competing_score',
            pl.when(first).then(None).when(higher).then(pl.col('kernel_score')).otherwise(pl.col('_new_score'))).alias('kernel_competing_score'),
        pl.when(first | higher).then(pl.col('_new_score').is_not_null().cast(pl.UInt32))
            .otherwise(pl.col('tied_peptides')+tied.cast(pl.UInt32)).alias('tied_peptides'),
        (pl.col('exact_peptides').fill_null(0)+(pl.col('_new_exact')>0).fill_null(False).cast(pl.UInt32)).alias('exact_peptides'))


def finish_kernel_selection(queries, state):
    """Format compact streaming state with every original query disposition."""
    import polars as pl
    if state is None:
        state=pl.DataFrame(schema={'query_id':queries.schema['query_id'],'selected_epitope':pl.String,
            'kernel_score':pl.Float64,'kernel_competing_score':pl.Float64,
            'tied_peptides':pl.UInt32,'exact_peptides':pl.UInt32})
    return _format_kernel_selection(queries,state)


def select_kernel(queries, scores):
    """One ranked peptide per query, with ties/no-support retained explicitly.

    Restrictions of the same peptide do not compete against one another. The
    support share is winner/(winner+strongest other peptide), a relative scale,
    not a binding probability. Original query order and metadata are retained.
    """
    import polars as pl
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('query identities must be unique and non-null')
    peptides=kernel_peptides(scores)
    top=peptides.sort(['query_id','kernel_score','epitope'],descending=[False,True,False],nulls_last=True)
    ranked=top.group_by('query_id',maintain_order=True).agg(
        pl.col('epitope').first().alias('selected_epitope'),
        pl.col('kernel_score').first(),
        pl.col('kernel_score').get(1,null_on_oob=True).alias('kernel_competing_score'),
        (pl.col('kernel_score')==pl.col('kernel_score').max()).sum().alias('tied_peptides'),
        (pl.col('exact_reference_count')>0).sum().alias('exact_peptides'))
    return _format_kernel_selection(queries,ranked)


def _format_kernel_selection(queries, ranked):
    import polars as pl
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('query identities must be unique and non-null')
    out=queries.join(ranked,on='query_id',how='left',validate='1:1',maintain_order='left').with_columns(
        pl.col('tied_peptides','exact_peptides').fill_null(0))
    supported=(pl.col('kernel_score')>0).fill_null(False)
    return out.with_columns(
        pl.when(supported).then(pl.col('kernel_score')/(pl.col('kernel_score')+
            pl.col('kernel_competing_score').fill_null(0))).alias('kernel_support_share'),
        pl.when(~supported).then(pl.lit('no_support')).when((pl.col('tied_peptides')>1)|(pl.col('exact_peptides')>1))
            .then(pl.lit('ambiguous')).otherwise(pl.lit('ranked_evidence')).alias('selection_status'),
        pl.when(supported).then(pl.col('selected_epitope')).alias('selected_epitope'),
        pl.when(supported).then(pl.col('kernel_score')).alias('kernel_score'))
