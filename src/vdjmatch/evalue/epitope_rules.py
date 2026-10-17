"""Reference-defined joint V usage and consistently V-conditioned sequence kernels.

Empirical predictive ranking, never a calibrated P/E-value or biological posterior.
Motif overlap is a kernel, not a product of supposedly independent motif hits.
"""
import polars as pl
import numpy as np

MOTIFS=('exact3','imgt3','wild4','imgt_wild4')
IMGT_CLASSES={a:group[0] for group in ('AILV','RHK','CM','ST','DE','NQ','G','F','P','W','Y') for a in group}


def _tokens(frame, axes, mode, *, full_junction=False, wildcard_positions=(1,2)):
    """Unique per-receptor motif tuples; paired tuples belong to one linked pair."""
    k=4 if 'wild' in mode else 3
    result=frame.select('_id',*[f'v{i}' for i in axes])
    for i in axes:
        seq=pl.col(f'seq{i}') if full_junction else pl.col(f'seq{i}').str.slice(3,pl.col(f'seq{i}').str.len_chars()-6)
        if 'imgt' in mode:seq=seq.str.replace_many(list(IMGT_CLASSES),list(IMGT_CLASSES.values()))
        tokens=frame.select('_id',seq.alias('_seq')).with_columns(
            pl.int_ranges(0,(pl.col('_seq').str.len_chars().cast(pl.Int64)-k+1).clip(lower_bound=0)).alias('_pos'))
        tokens=tokens.explode('_pos').filter(pl.col('_pos').is_not_null()).with_columns(
            pl.col('_seq').str.slice(pl.col('_pos'),k).alias('_token'))
        if k==4:
            tokens=pl.concat([tokens.select('_id',pl.concat_str(pl.col('_token').str.slice(0,p),pl.lit('X'),pl.col('_token').str.slice(p+1,3-p)).alias(f't{i}')) for p in wildcard_positions])
        else:tokens=tokens.select('_id',pl.col('_token').alias(f't{i}'))
        # A tuple is a single paired kernel feature; repeated/overlapping tokens
        # are averaged rather than claimed as independent observations.
        result=result.join(tokens.unique(),on='_id',how='inner',validate='m:m')
    return result


def _rule_rank_mass(q, r, b, out, *, paired, exclude_exact, species, threads, native_radius, v_loop_distances):
    """Scalar oracle geometry preparation and occupied-shell reduction."""
    from ..match.geometry import grouped_total_distance_count_batch
    axes=range(2 if paired else 1);seq=[f'seq{i}' for i in axes]
    # Exhaustive CDR3-only rank uses the same matrix-derived gaps and quantum.
    # This is a finite geometry bound, not a radius/parameter scan.
    from .local_rank import _harmonic_interval
    from scipy.special import digamma
    from ..match.geometry import distance_units
    matrix_scale,quantum=distance_units();step=3*quantum
    # One union traversal supplies both populations' CDFs; shared keys retain
    # membership in both models. Never deduplicate across model memberships.
    keys=[*seq,*['_allele'+str(i) for i in axes]]
    members=pl.concat([pop.select(keys).with_columns(pl.lit(side).alias('_side'))
                       for side,pop in enumerate((r,b))])
    union=members.select(keys).unique(maintain_order=True).with_row_index('_ref')
    members=members.join(union,on=keys,how='left',validate='m:1').sort('_ref','_side')
    offsets=np.r_[0,np.cumsum(np.bincount(members['_ref'].to_numpy(),minlength=union.height))].tolist()
    native_reference=union.select(*seq,*[pl.col('_allele'+str(i)).alias('v'+str(i)) for i in axes])
    bound=3*matrix_scale*sum(max((f[f'seq{i}'].str.len_chars().max() or 0) for f in (q,r,b)) for i in axes)
    cuts=list(range(0,200*bound+1,200*step))
    weights=np.maximum(0,native_radius-np.asarray(cuts,dtype=np.int64))
    # Seven uint32 channels: whole, two target, two same-V, two exclusive.
    # Only compact per-query masses/ranks survive each bounded native chunk.
    chunk_size=max(1,min(65536,(256*1024**2)//(7*len(cuts)*4)))
    mass=np.zeros((q.height,2));rank=np.ones(q.height)
    def ranks(part, conditional):
        if not part.height:return
        native_query=part.select(*seq,*[pl.col('_allele'+str(i)).alias('v'+str(i)) for i in axes])
        iterator=grouped_total_distance_count_batch(native_query,native_reference,cuts,
            paired=paired,species=species,threads=threads,ctrim=3,v_loop_weight=0,
            same_v_gene=conditional,exclude_exact=exclude_exact,v_loop_distances=v_loop_distances,
            chunk_size=chunk_size,reference_target_offsets=offsets,
            reference_target_ids=members['_side'].to_list(),target_peptide_ids=[0,1],
            reference_exclusive_peptide_ids=[-1]*union.height)
        for chunk in iterator:
            ids=part['_id'].slice(chunk.start,chunk.stop-chunk.start).to_numpy()
            curves=np.asarray(chunk.counts)[:,1:3,:]
            suffix='_v_count' if conditional else '_population'
            populations=[out['rule_'+side+suffix].to_numpy()[ids] for side in ('reference','control')]
            excluded_suffix='_removed_v' if conditional else '_excluded_hits'
            excluded=[out['rule_'+side+excluded_suffix].to_numpy()[ids] for side in ('reference','control')]
            if any(not np.array_equal(curves[:,i,-1],pop-removed) for i,(pop,removed) in enumerate(zip(populations,excluded))):
                raise ArithmeticError('CDR3 rank bound did not cover every eligible reference')
            if conditional:
                shells=np.diff(curves.astype(np.int64),axis=2,prepend=0)
                mass[ids]=shells@weights/native_radius
            target,control=curves[:,0,:],curves[:,1,:]
            shells=np.diff(target,axis=1,prepend=0)
            row,column=np.nonzero(shells)
            before=np.where(column>0,control[row,column-1],0).astype(np.int64)
            # Occupied shells retain their fixed distance/reduction order.
            numerator=np.bincount(row,weights=shells[row,column]*_harmonic_interval(before+1,control[row,column].astype(np.int64)+1),
                minlength=len(ids))
            N,M=populations
            result=np.ones(len(ids));usable=(N>0)&(M>0)
            result[usable]=(M[usable]+1)*numerator[usable]/(N[usable]*(digamma(M[usable]+2)+np.euler_gamma))
            rank[ids]=result
    ranks(q,True)
    fallback=~out['rule_v_supported'].to_numpy()
    if fallback.any():ranks(q.filter(pl.Series(fallback)),False)
    return mass,rank


def reference_epitope_rules(queries, reference, background, *, paired=False,
                            exclude_exact=False, species='human', threads=1, radius=90, v_loop_distances=None,
                            _rank_mass=None):
    """Frames: query_id, seq0/v0 [seq1/v1], with resolved V alleles.

    For a supported joint gene category g, shared prior p_g=(n_g+m_g+.5)/(N+M+.5K).
    Population predictive V probability is (n_g+p_g)/(N+1). Conditional
    kernel mean is (mass_g+c)/(n_g+1), c=(mass_ref+mass_bg+.5)/(n_g+m_g+1).
    The score is the ratio of their products, counting V once. K is the union
    of reference/control gene tuples plus one unseen category. Population/gene
    counts, pooled priors and motif frequencies use the fixed full models.
    Benchmark exclusion masks exact-hit numerators only; paired kernels compare
    the same linked reference pair. Excluded hits never reappear as censored ranks.
    Rank uses all CDR3-only reference distances within g; if either population
    lacks g, V is neutral and rank uses all V categories. Enriched motifs add a
    bounded correlated sequence bonus, never a calibrated independent likelihood.
    """
    from ..match.geometry import distance_cutoff,precompute_v_loop_distances
    if not 0<radius<2**28/200:raise ValueError('positive radius below native V barrier required')
    native_radius=distance_cutoff(radius)
    if native_radius==0:raise ValueError('positive radius in native units required')
    if v_loop_distances is None:v_loop_distances=precompute_v_loop_distances(species=species)
    if any(v_loop_distances.get(k)!=value for k,value in
           [('species',species),('matrix','clipped'),('alphabet','full')]):
        raise ValueError('epitope rules require matching species, clipped matrix and full alphabet')
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('unique non-null query IDs required')
    axes=range(2 if paired else 1);seq=[f'seq{i}' for i in axes];v=[f'v{i}' for i in axes]
    def prepare(frame,query=False):
        f=frame.select(*(['query_id'] if query else []),*seq,*v)
        if any(f[c].null_count() for c in seq+v):raise ValueError('complete junctions/resolved V alleles required')
        if not query:f=f.unique([*seq,*v],maintain_order=True)
        return f.with_row_index('_id').with_columns(*[pl.col(c).alias('_allele'+str(i)) for i,c in enumerate(v)]).with_columns(*[pl.col(c).str.replace(r'\*.*$','') for c in v])
    q=prepare(queries,True);r=prepare(reference);b=prepare(background)
    K=pl.concat([r.select(v),b.select(v)]).unique().height+1
    def counts(pop,side):
        out=q.select('_id','query_id',*seq,*v).join(pop.group_by(v).len().rename({'len':'g'}),on=v,how='left',validate='m:1')
        out=out.with_columns(pl.col('g').fill_null(0).cast(pl.Int64),pl.lit(pop.height,dtype=pl.Int64).alias('N'))
        removed=pl.DataFrame(schema={'_id':pl.UInt32,'_total':pl.Int64,'_g':pl.Int64})
        if exclude_exact and pop.height:
            same=q.join(pop,on=seq,how='inner',suffix='_ref',validate='m:m')
            match=pl.all_horizontal(*[pl.col(c)==pl.col(c+'_ref') for c in v])
            removed=same.group_by('_id').agg(pl.len().cast(pl.Int64).alias('_total'),match.sum().cast(pl.Int64).alias('_g'))
        out=out.join(removed,on='_id',how='left',validate='1:1').with_columns(
            pl.col('_total').fill_null(0),pl.col('_g').fill_null(0))
        return out.select('_id',pl.col('N').alias(side+'_population'),pl.col('g').alias(side+'_v_count'),
            pl.col('_total').alias(side+'_excluded_hits'),pl.col('_g').alias(side+'_removed_v'))
    out=q.select('_id','query_id').join(counts(r,'rule_reference'),on='_id',validate='1:1').join(counts(b,'rule_control'),on='_id',validate='1:1').sort('_id')
    N=pl.col('rule_reference_population');M=pl.col('rule_control_population')
    n=pl.col('rule_reference_v_count');m=pl.col('rule_control_v_count')
    pg=(n+m+.5)/(N+M+.5*K)
    out=out.with_columns(pg.alias('rule_pooled_v_probability'),((n+pg)/(N+1)).alias('rule_reference_v_probability'),((m+pg)/(M+1)).alias('rule_control_v_probability'))
    supported=(n>0)&(m>0)
    out=out.with_columns(supported.alias('rule_v_supported'),
        pl.when(supported).then(pl.col('rule_reference_v_probability')/pl.col('rule_control_v_probability')).otherwise(1.0).alias('rule_v_score'),
        pl.when((N>0)&(M>0)).then(pl.lit('empirical_reference_rule_only')).otherwise(pl.lit('empty_population')).alias('rule_status'))
    def score(kind):
        rm=pl.col('rule_reference_'+kind+'_mass');bm=pl.col('rule_control_'+kind+'_mass')
        # The pooled pseudo-observation is a model statistic, not a filtered
        # hit statistic. An exact same-V receptor has kernel one when available.
        available=pl.col('rule_'+kind+'_available') if kind in MOTIFS else pl.lit(True)
        exact=pl.when(available).then(pl.col('rule_reference_removed_v')+pl.col('rule_control_removed_v')).otherwise(0)
        c=(rm+bm+exact+.5)/(n+m+1)
        return pl.col('rule_v_score')*((rm+c)/(n+1))/((bm+c)/(m+1))
    if _rank_mass is None:
        mass,rank=_rule_rank_mass(q,r,b,out,paired=paired,exclude_exact=exclude_exact,
            species=species,threads=threads,native_radius=native_radius,v_loop_distances=v_loop_distances)
    else:
        mass,rank,conditional_totals,global_totals=_rank_mass
        if mass.shape!=(q.height,2) or rank.shape!=(q.height,):
            raise ValueError('prepared rule rank/mass axes differ from query population')
        expected=np.column_stack([out['rule_'+side+'_v_count'].to_numpy()-
            out['rule_'+side+'_removed_v'].to_numpy() for side in ('reference','control')])
        fallback=~out['rule_v_supported'].to_numpy()
        global_expected=np.column_stack([out['rule_'+side+'_population'].to_numpy()-
            out['rule_'+side+'_excluded_hits'].to_numpy() for side in ('reference','control')])
        if not np.array_equal(conditional_totals,expected) or not np.array_equal(global_totals[fallback],global_expected[fallback]):
            raise ArithmeticError('CDR3 rank bound did not cover every eligible reference')
    out=out.with_columns(*[pl.Series('rule_'+side+'_spatial_mass',mass[:,i]) for i,side in enumerate(('reference','control'))])
    out=out.with_columns(score('spatial').alias('rule_spatial_score'))
    out=out.with_columns(pl.Series('rule_cdr3_rank',rank)).with_columns(
        (pl.col('rule_v_score')*pl.col('rule_cdr3_rank')).alias('rule_rank_score'))
    for mode in MOTIFS:
        qt=_tokens(q,axes,mode);token_cols=[f't{i}' for i in axes]
        lengths=qt.group_by('_id').len().rename({'len':'_tokens'})
        frequencies=[]
        for side,pop in [('reference',r),('control',b)]:
            pt=_tokens(pop,axes,mode)
            freq=pt.group_by([*v,*token_cols]).len().rename({'len':'_count'})
            frequencies.append(freq.rename({'_count':'_'+side}))
            mass=qt.join(freq,on=[*v,*token_cols],how='left',validate='m:1').group_by('_id').agg(pl.col('_count').fill_null(0).sum().cast(pl.Float64).alias('_mass')).join(lengths,on='_id',validate='1:1')
            mass=mass.select('_id',(pl.col('_mass')/pl.col('_tokens')).alias('_mass'))
            out=out.join(mass,on='_id',how='left',validate='1:1').with_columns(
                (pl.col('_mass').fill_null(0)-pl.when(pl.col('_mass').is_not_null()).then(pl.col('rule_'+side+'_removed_v')).otherwise(0)).alias('rule_'+side+'_'+mode+'_mass')).drop('_mass')
        usable=out.select('_id').join(lengths,on='_id',how='left',validate='1:1')['_tokens'].fill_null(0)>0
        out=out.with_columns(pl.Series('rule_'+mode+'_available',usable))
        for side in ('reference','control'):
            # Fixed equal mixture of kernels: overlap is not independent evidence.
            out=out.with_columns(pl.when(pl.col('rule_'+mode+'_available')).then((pl.col('rule_'+side+'_spatial_mass')+pl.col('rule_'+side+'_'+mode+'_mass'))/2)
                .otherwise(pl.col('rule_'+side+'_spatial_mass')).alias('rule_'+side+'_mixed_'+mode+'_mass'))
        out=out.with_columns(score(mode).alias('rule_'+mode+'_score'),score('mixed_'+mode).alias('rule_mixed_'+mode+'_score'))
        # Reference enrichment is a fixed V-conditioned model statistic. Query
        # tokens only look up its frequencies; hit exclusion never refits it.
        if (int(out['rule_reference_v_count'].max() or 0) * int(out['rule_control_v_count'].max() or 0) *
                int(lengths['_tokens'].max() or 0) > 2**53-1):
            raise ValueError('conditional motif counts exceed the exactly representable integer budget')
        enriched=qt.join(frequencies[0],on=[*v,*token_cols],how='left',validate='m:1').join(frequencies[1],on=[*v,*token_cols],how='left',validate='m:1')
        enriched=enriched.join(out.select('_id','rule_reference_v_count','rule_control_v_count'),on='_id',validate='m:1')
        # Sum the exact positive cross-count differences before one division.
        # A floating group mean can change its last bit with the query batch.
        delta=(pl.col('_reference').fill_null(0).cast(pl.Int64)*pl.col('rule_control_v_count')-
            pl.col('_control').fill_null(0).cast(pl.Int64)*pl.col('rule_reference_v_count')).clip(lower_bound=0)
        enriched=enriched.group_by('_id').agg(delta.sum().alias('_sum'),pl.len().alias('_terms'),
            pl.col('rule_reference_v_count').first(),pl.col('rule_control_v_count').first()).select('_id',
                (pl.col('_sum')/(pl.col('_terms').cast(pl.Int64)*pl.col('rule_reference_v_count').clip(lower_bound=1)*
                    pl.col('rule_control_v_count').clip(lower_bound=1))).alias('_enrichment'))
        out=out.join(enriched,on='_id',how='left',validate='1:1').with_columns(
            pl.when(pl.col('rule_v_supported')).then(pl.col('_enrichment').fill_null(0)).otherwise(0.0).alias('rule_'+mode+'_enrichment')).drop('_enrichment')
        # Bounded correlated sequence bonus, not independent likelihood evidence.
        out=out.with_columns((pl.col('rule_rank_score')*(1+pl.col('rule_'+mode+'_enrichment'))).alias('rule_enriched_'+mode+'_score'))
    scores=[c for c in out.columns if c.endswith('_score')]
    return out.with_columns(*[pl.when((N>0)&(M>0)).then(pl.col(c)).otherwise(None).alias(c) for c in scores]).drop('_id')


def grouped_reference_epitope_rules(reference, groups, *, paired=False, exclude_exact=False,
                                    species='human', threads=1, radius=90, v_loop_distances=None):
    """Share one conditional geometry traversal across declared task models.

    Reference uses seq0/v0 [seq1/v1] and epitope; groups contain task, queries,
    selected reference. Only compact per-observation masses/ranks survive native
    chunks. Scalar priors, motif counts and score reductions remain the oracle.
    """
    from ..match.geometry import grouped_total_distance_count_batch, distance_units, distance_cutoff, precompute_v_loop_distances
    from .local_rank import _harmonic_interval
    from scipy.special import digamma
    if not groups:return []
    if sum(q.height for _,q,_ in groups)>5_000_000:
        raise ValueError('grouped rule observations exceed the five-million detailed-case guard')
    axes=range(2 if paired else 1)
    seq=[f'seq{i}' for i in axes];genes=[f'v{i}' for i in axes];keys=[*seq,*genes]
    if any(frame[c].null_count() for frame in [reference,*[q for _,q,_ in groups]] for c in keys):
        raise ValueError('complete junctions/resolved V alleles required')
    options=dict(paired=paired,exclude_exact=exclude_exact,species=species,threads=threads,
        radius=radius,v_loop_distances=v_loop_distances)
    if not any(q.height for _,q,_ in groups):
        return [reference_epitope_rules(q,r,reference.filter(pl.col('epitope')!=task['epitope']),**options)
            for task,q,r in groups]
    if not 0<radius<2**28/200 or not distance_cutoff(radius):
        raise ValueError('positive radius below native V barrier required')
    context=v_loop_distances if v_loop_distances is not None else precompute_v_loop_distances(species=species)
    if any(context.get(k)!=value for k,value in [('species',species),('matrix','clipped'),('alphabet','full')]):
        raise ValueError('epitope rules require matching species, clipped matrix and full alphabet')
    options['v_loop_distances']=context
    for _,q,_ in groups:
        if q['query_id'].null_count() or q['query_id'].n_unique()!=q.height:
            raise ValueError('unique non-null query IDs required')
    whole=reference.select(keys).unique(maintain_order=True).with_row_index('_ref')
    query=pl.concat([q.select(keys) for _,q,_ in groups]).unique(maintain_order=True).with_row_index('_row')
    peptides=list(dict.fromkeys(task['epitope'] for task,_,_ in groups))
    T=len(groups);models=T+1+len(peptides)
    # Explicit whole/exclusive models also expose fixed same-V denominators in
    # native target-size channels, without per-query model decrements.
    members=[r.select(keys).unique().with_columns(pl.lit(i).alias('_target')) for i,(_,_,r) in enumerate(groups)]
    members.append(whole.select(keys).with_columns(pl.lit(T).alias('_target')))
    exclusive=(reference.group_by(keys).agg(pl.col('epitope').n_unique().alias('_labels'),
        pl.col('epitope').first()).filter(pl.col('_labels')==1))
    for i,peptide in enumerate(peptides):
        members.append(exclusive.filter(pl.col('epitope')==peptide).select(keys)
            .with_columns(pl.lit(T+1+i).alias('_target')))
    members=pl.concat(members).join(whole,on=keys,how='left',validate='m:1').sort('_ref','_target')
    if members['_ref'].null_count():raise ValueError('selected rule references must belong to the whole population')
    offsets=np.r_[0,np.cumsum(np.bincount(members['_ref'].to_numpy(),minlength=whole.height))]
    requested=pl.concat([q.select(keys).with_row_index('_observation').join(query,on=keys,
        how='left',validate='m:1',maintain_order='left').select('_row','_observation')
        .with_columns(pl.lit(i).alias('_target')) for i,(_,q,_) in enumerate(groups)]).sort('_row','_target')
    matrix_scale,quantum=distance_units()
    bound=3*matrix_scale*sum(max((f[c].str.len_chars().max() or 0) for f in (query,whole)) for c in seq)
    cuts=list(range(0,200*bound+1,600*quantum))
    weights=np.maximum(0,distance_cutoff(radius)-np.asarray(cuts,dtype=np.int64))
    channels=1+2*models+1
    # Include native counts/sizes/exact counts AND gathered curves, shells and
    # rank scratch for the largest possible requested observation multiplicity.
    max_requests=int(requested.group_by('_row').len()['len'].max())
    bytes_per_query=4*(channels*(len(cuts)+1)+models)+max_requests*(len(cuts)*80+256)
    if bytes_per_query>256*1024**2:
        raise ValueError('one grouped rule query exceeds the 256 MiB working-buffer budget')
    chunk_size=min(65536,max(1,256*1024**2//bytes_per_query))
    masses=[np.zeros((q.height,2)) for _,q,_ in groups]
    ranks=[np.ones(q.height) for _,q,_ in groups]
    conditional_totals=[np.zeros((q.height,2),dtype=np.int64) for _,q,_ in groups]
    global_totals=[np.zeros((q.height,2),dtype=np.int64) for _,q,_ in groups]
    fallback=[]
    peptide_index=np.asarray([peptides.index(task['epitope']) for task,_,_ in groups])
    def traverse(part,requests,conditional):
        if not part.height or not requests.height:return
        iterator=grouped_total_distance_count_batch(part.select(keys),whole.select(keys),cuts,
            paired=paired,species=species,threads=threads,ctrim=3,v_loop_weight=0,
            matrix='clipped',same_v_gene=conditional,exclude_exact=exclude_exact,v_loop_distances=context,
            chunk_size=chunk_size,reference_target_offsets=offsets,
            reference_target_ids=members['_target'].to_numpy(),target_peptide_ids=[0]*models,
            reference_exclusive_peptide_ids=[-1]*whole.height)
        requested_rows=requests['_row'].to_numpy()
        for chunk in iterator:
            lo,hi=np.searchsorted(requested_rows,[chunk.start,chunk.stop])
            selected=requests.slice(int(lo),int(hi-lo))
            if not selected.height:continue
            row=selected['_row'].to_numpy()-chunk.start;target=selected['_target'].to_numpy()
            removed=T+1+peptide_index[target]
            counts=np.asarray(chunk.counts);sizes=np.asarray(chunk.sizes)
            curves=np.stack([counts[row,1+target,:],
                counts[row,1+T,:]-counts[row,1+removed,:]],axis=1)
            size_offset=1+models if conditional else 1
            N=sizes[row,size_offset+target]
            M=sizes[row,size_offset+T]-sizes[row,size_offset+removed]
            if np.any(curves[:,:,-1]>np.stack([N,M],axis=1)):
                raise ArithmeticError('grouped rule counts exceed fixed model populations')
            if conditional:
                mass=np.diff(curves.astype(np.int64),axis=2,prepend=0)@weights/distance_cutoff(radius)
                unsupported=(N==0)|(M==0)
                if unsupported.any():fallback.append(selected.filter(pl.Series(unsupported)))
            target_curve,control=curves[:,0,:],curves[:,1,:]
            shells=np.diff(target_curve,axis=1,prepend=0)
            occupied,columns=np.nonzero(shells)
            before=np.where(columns>0,control[occupied,columns-1],0).astype(np.int64)
            numerator=np.bincount(occupied,weights=shells[occupied,columns]*_harmonic_interval(
                before+1,control[occupied,columns].astype(np.int64)+1),minlength=len(row))
            result=np.ones(len(row));usable=(N>0)&(M>0)
            result[usable]=(M[usable]+1)*numerator[usable]/(N[usable]*(digamma(M[usable]+2)+np.euler_gamma))
            for (target_id,),items in selected.with_row_index('_local').group_by('_target',maintain_order=True):
                ids=items['_observation'].to_numpy();local=items['_local'].to_numpy()
                ranks[target_id][ids]=result[local]
                if conditional:masses[target_id][ids]=mass[local]
                (conditional_totals if conditional else global_totals)[target_id][ids]=curves[local,:,-1]
    traverse(query,requested,True)
    if fallback:
        fallback=pl.concat(fallback)
        original_query=query
        query=query.join(fallback.select('_row').unique(),on='_row',how='semi').rename({'_row':'_original_row'}).with_row_index('_row')
        requests=fallback.rename({'_row':'_original_row'}).join(query.select('_original_row','_row'),
            on='_original_row',how='left',validate='m:1').sort('_row','_target')
        traverse(query,requests,False)
        query=original_query
    return [reference_epitope_rules(q,r,reference.filter(pl.col('epitope')!=task['epitope']).select(keys),
        _rank_mass=(mass,rank,conditional,global_),**options)
        for (task,q,r),mass,rank,conditional,global_ in zip(groups,masses,ranks,conditional_totals,global_totals)]
