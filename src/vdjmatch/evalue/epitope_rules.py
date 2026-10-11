"""Reference-defined joint V usage and consistently V-conditioned sequence kernels.

Empirical predictive ranking, never a calibrated P/E-value or biological posterior.
Motif overlap is a kernel, not a product of supposedly independent motif hits.
"""
import polars as pl

MOTIFS=('exact3','imgt3','wild4','imgt_wild4')
IMGT_CLASSES={a:group[0] for group in ('AILV','RHK','CM','ST','DE','NQ','G','F','P','W','Y') for a in group}


def _tokens(frame, axes, mode):
    """Unique per-receptor motif tuples; paired tuples belong to one linked pair."""
    k=4 if 'wild' in mode else 3
    result=frame.select('_id',*[f'v{i}' for i in axes])
    for i in axes:
        seq=pl.col(f'seq{i}').str.slice(3,pl.col(f'seq{i}').str.len_chars()-6)
        if 'imgt' in mode:seq=seq.str.replace_many(list(IMGT_CLASSES),list(IMGT_CLASSES.values()))
        tokens=frame.select('_id',seq.alias('_seq')).with_columns(
            pl.int_ranges(0,(pl.col('_seq').str.len_chars().cast(pl.Int64)-k+1).clip(lower_bound=0)).alias('_pos'))
        tokens=tokens.explode('_pos').filter(pl.col('_pos').is_not_null()).with_columns(
            pl.col('_seq').str.slice(pl.col('_pos'),k).alias('_token'))
        if k==4:
            tokens=pl.concat([tokens.select('_id',pl.concat_str(pl.col('_token').str.slice(0,2),pl.lit('X'),pl.col('_token').str.slice(3,1)).alias(f't{i}')),
                tokens.select('_id',pl.concat_str(pl.col('_token').str.slice(0,1),pl.lit('X'),pl.col('_token').str.slice(2,2)).alias(f't{i}'))])
        else:tokens=tokens.select('_id',pl.col('_token').alias(f't{i}'))
        # A tuple is a single paired kernel feature; repeated/overlapping tokens
        # are averaged rather than claimed as independent observations.
        result=result.join(tokens.unique(),on='_id',how='inner',validate='m:m')
    return result


def reference_epitope_rules(queries, reference, background, *, paired=False,
                            exclude_exact=False, species='human', threads=1, radius=90):
    """Frames: query_id, seq0/v0 [seq1/v1], with resolved V alleles.

    For a joint gene category g, shared prior p_g=(n_g+m_g+.5)/(N+M+.5K).
    Population predictive V probability is (n_g+p_g)/(N+1). Conditional
    kernel mean is (mass_g+c)/(n_g+1), c=(mass_ref+mass_bg+.5)/(n_g+m_g+1).
    The score is the ratio of their products, counting V once. K is the union
    of reference/control gene tuples plus one unseen category. All counts are
    query-punctured; paired kernels compare the same linked reference pair.
    """
    from ..match.tcrdist import total_distance_count_batch,paired_total_distance_mass_batch
    if not 0<radius<2**28/200:raise ValueError('positive radius below native V barrier required')
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
            (pl.col('N')-pl.col('_total').fill_null(0)).alias('N'),
            (pl.col('g')-pl.col('_g').fill_null(0)).alias('g'),pl.col('_g').fill_null(0))
        return out.select('_id',pl.col('N').alias(side+'_population'),pl.col('g').alias(side+'_v_count'),pl.col('_g').alias(side+'_removed_v'))
    out=q.select('_id','query_id').join(counts(r,'rule_reference'),on='_id',validate='1:1').join(counts(b,'rule_control'),on='_id',validate='1:1').sort('_id')
    N=pl.col('rule_reference_population');M=pl.col('rule_control_population')
    n=pl.col('rule_reference_v_count');m=pl.col('rule_control_v_count')
    pg=(n+m+.5)/(N+M+.5*K)
    out=out.with_columns(pg.alias('rule_pooled_v_probability'),((n+pg)/(N+1)).alias('rule_reference_v_probability'),((m+pg)/(M+1)).alias('rule_control_v_probability'))
    out=out.with_columns((pl.col('rule_reference_v_probability')/pl.col('rule_control_v_probability')).alias('rule_v_score'),
        pl.when((N>0)&(M>0)).then(pl.lit('empirical_reference_rule_only')).otherwise(pl.lit('empty_population')).alias('rule_status'))
    def spatial(pop):
        if not pop.height or not q.height:return [0.0]*q.height
        options=dict(species=species,threads=threads,ctrim=3,exclude_exact=exclude_exact,v_loop_weight=0,same_v_gene=True)
        cut=[[200*radius]]*q.height
        if paired:
            _,mass=paired_total_distance_mass_batch(q['seq0'],q['seq1'],pop['seq0'],pop['seq1'],q['_allele0'],q['_allele1'],pop['_allele0'],pop['_allele1'],cut,**options)
        else:
            _,mass=total_distance_count_batch(q['seq0'],pop['seq0'],q['_allele0'],pop['_allele0'],cut,linear_mass=True,**options)
        return [row[0]/(200*radius) for row in mass]
    for side,pop in [('reference',r),('control',b)]:out=out.with_columns(pl.Series('rule_'+side+'_spatial_mass',spatial(pop),dtype=pl.Float64))
    def score(kind):
        rm=pl.col('rule_reference_'+kind+'_mass');bm=pl.col('rule_control_'+kind+'_mass')
        c=(rm+bm+.5)/(n+m+1)
        return pl.col('rule_v_score')*((rm+c)/(n+1))/((bm+c)/(m+1))
    out=out.with_columns(score('spatial').alias('rule_spatial_score'))
    for mode in MOTIFS:
        qt=_tokens(q,axes,mode);token_cols=[f't{i}' for i in axes]
        lengths=qt.group_by('_id').len().rename({'len':'_tokens'})
        for side,pop in [('reference',r),('control',b)]:
            pt=_tokens(pop,axes,mode)
            freq=pt.group_by([*v,*token_cols]).len().rename({'len':'_count'})
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
    scores=[c for c in out.columns if c.endswith('_score')]
    return out.with_columns(*[pl.when((N>0)&(M>0)).then(pl.col(c)).otherwise(None).alias(c) for c in scores]).drop('_id')
