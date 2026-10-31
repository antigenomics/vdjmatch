"""Strict Pgen-bin/count-matched repertoire enrichment from batched memberships."""
from __future__ import annotations

import numpy as np
import polars as pl
from scipy import sparse

from ..evalue.mock_reference import _sample_event_indices
from ..io.columns import _integer_count


def exact_reference_counts(observations: pl.DataFrame, panels: pl.DataFrame, *,
                           matcher='exact', mismatches=1, threads=1, radius_factor=None,
                           locus='TRB', species='human') -> pl.DataFrame:
    """Count source rows and reads against all exact V/J/junction panels at once.

    Observations contain junction_aa/v/j/count; panels additionally contain
    group_id/replicate. Source rows sharing a key remain separate observations.
    Repeated sampled panel keys do not multiply observed counts. Zero-hit panels are kept.
    This is the original bulk-repertoire counting unit, not paired-cell union.
    Optional Hamming matching changes junction membership only; exact
    V/J and one count per source observation per panel remain mandatory.
    """
    keys=['junction_aa','v','j'];groups=['group_id','replicate']
    if not {*keys,'count'}<=set(observations.columns) or not {*keys,*groups}<=set(panels.columns) or not panels.height:
        raise ValueError('source keys/counts and nonempty reference panels required')
    if any(panels[c].null_count() for c in [*groups,*keys]):
        raise ValueError('reference panel identities and keys must be non-null')
    panels=panels.unique(subset=[*groups,*keys],maintain_order=True)
    obs=observations.select(*keys,_integer_count('count').alias('count'))
    if obs['count'].null_count() or (obs['count']<0).any():
        raise ValueError('nonnegative integer source counts required')
    if obs.select(pl.col('count').cast(pl.Int128).sum()).item()>np.iinfo(np.int64).max:
        raise ValueError('source read total exceeds int64')
    counted=obs.group_by(keys).agg(pl.len().cast(pl.Int64).alias('matched_records'),pl.col('count').sum().alias('matched_reads'))
    weighted=matcher in {'blosum62','positional-blosum62'}
    from ..match.scope import GRADED_FAMILIES, graded_family_params
    if matcher not in {'exact','hamming','blosum62','positional-blosum62',*GRADED_FAMILIES}:
        raise ValueError('unknown panel matcher')
    from ..match.scope import blosum_radius_params
    if weighted:blosum_radius_params(matcher,radius_factor,1)
    elif matcher in GRADED_FAMILIES:
        params,scale=graded_family_params(matcher,mismatches,radius_factor)
    elif radius_factor is not None:raise ValueError('radius factor requires a BLOSUM matcher')
    if type(threads) is not int or threads<1 or type(mismatches) is not int or not 0<=mismatches<=3:
        raise ValueError('positive threads and zero to three mismatches required')
    if matcher=='exact' or not counted.height:
        members=panels.join(counted,on=keys,how='inner',validate='m:1')
    else:
        from seqtree import Index
        from ..match.scope import search_params
        observed=counted.select('junction_aa').unique().sort('junction_aa').with_row_index('_ri')
        queries=panels.select('junction_aa').unique().sort('junction_aa').with_row_index('_qi')
        index=Index.build(observed['junction_aa'].to_list(),alphabet='aa')
        buffers=[]
        def consume(start,stop,buffer):
            pairs=np.asarray(buffer)
            if matcher in GRADED_FAMILIES:
                pairs=pairs[pairs[:,2]/scale<=radius_factor]
            if len(pairs):
                pairs=pairs[:,:2].copy();pairs[:,0]=query_ids[pairs[:,0]]
                buffers.append(pairs)
        batches=(queries.with_columns(pl.col('junction_aa').str.len_chars().alias('_length')).partition_by('_length',as_dict=True).values() if matcher=='positional-blosum62' else [queries])
        for batch in batches:
            query_ids=batch['_qi'].to_numpy()
            if matcher not in GRADED_FAMILIES:
                params=(blosum_radius_params(matcher,radius_factor,len(batch['junction_aa'][0])) if weighted else search_params(f'{mismatches},0,0,{mismatches}'))
            index.search_batch_consume(batch['junction_aa'].to_list(),params,consume,
                threads=threads,chunk_size=2048,max_hit_bytes=64*1024**2)
        pairs=pl.DataFrame(np.concatenate(buffers) if buffers else np.empty((0,2),dtype=np.uint32),
            schema=[('_qi',pl.UInt32),('_ri',pl.UInt32)],orient='row')
        if matcher in GRADED_FAMILIES:
            from ..evalue.control import impute_gene_calls
            kinds=('v','j') if matcher.startswith('v-') else ('j',)
            sampled=counted.with_row_index('_observation')
            resolved=impute_gene_calls(sampled,species=species,locus=locus,kinds=kinds).join(observed,on='junction_aa',validate='m:1')
            roles=impute_gene_calls(panels,species=species,locus=locus,kinds=kinds).join(queries,on='junction_aa',validate='m:1')
            predicate=['_ri',*kinds]
            linked=pairs.join(roles,on='_qi',validate='m:m').join(
                resolved.select('_ri',*kinds,'_observation'),on=predicate,validate='m:m')
            members=linked.select(*groups,'_observation').unique().join(sampled,on='_observation',validate='m:1')
        else:
            roles=panels.join(queries,on='junction_aa',validate='m:1').drop('junction_aa')
            memberships=pairs.join(roles,on='_qi',validate='m:m').join(observed,on='_ri',validate='m:1').select(*groups,*keys).unique()
            members=memberships.join(counted,on=keys,validate='m:1')
    matched=members.group_by(groups).agg(pl.col('matched_records','matched_reads').sum())
    return panels.select(groups).unique().join(matched,on=groups,how='left',validate='1:1').with_columns(
        pl.col('matched_records','matched_reads').fill_null(0)).sort(groups)


def matched_enrichment(observations, control, matches, targets, *, replicates=200, seed=0, measure=None):
    """Return real/null counts, empirical upper-tail P-values and exact-bin audit.

    Observations have sample_id/query_id/count/pgen_bin. Controls have query_id,
    event_weight and pgen_bin. Memberships have query_id plus the columns of the
    explicit target inventory (normally full pMHC). One supplied source record is
    one receptor observation, including repeated amino-acid keys; counts are reads.
    Every mock preserves N, each Pgen bin and the source count multiset in that bin.
    Control events are sampled without replacement within each bin; source counts
    are assigned in uniform random order. No absent-bin remapping is permitted.
    P-values are (1 + number of mocks >= real)/(replicates + 1). Intervals describe
    the central95percent mock distribution, not uncertainty over biological donors.
    No automatic detection threshold or between-sample comparison is applied.
    Optional measure(query_ids, counts) returns two target-ordered integer vectors
    by matching only that real/mock draw; it replaces stored memberships entirely.
    """
    if type(replicates) is not int or replicates<200 or type(seed) is not int or seed<0:
        raise ValueError('at least200mocks and a nonnegative integer seed required')
    if not observations.height or not control.height or not targets.height or not targets.columns:
        raise ValueError('nonempty source, control and explicit target inventory required')
    if not {'sample_id','query_id','count','pgen_bin'}<=set(observations.columns) or not {'query_id','event_weight','pgen_bin'}<=set(control.columns):
        raise ValueError('source/control identities, integer counts and Pgen bins required')
    columns=targets.columns
    if 'epitope' not in columns or (measure is None and not {'query_id',*columns}<=set(matches.columns)):
        raise ValueError('explicit epitope target memberships required')
    obs=observations.select('sample_id','query_id','pgen_bin',_integer_count('count').alias('count'))
    ctrl=control.select('query_id','pgen_bin',_integer_count('event_weight').alias('event_weight'))
    for frame,name in [(obs,'observations'),(ctrl,'control')]:
        if any(frame[c].null_count() for c in frame.columns) or not frame.schema['query_id'].is_integer() or not frame.schema['pgen_bin'].is_integer():
            raise ValueError(name+' requires non-null integer identities/bins/counts')
        if (frame['query_id']<0).any():raise ValueError('query identities must be nonnegative')
    if (obs['count']<0).any() or (ctrl['event_weight']<=0).any() or obs['sample_id'].cast(pl.String).str.strip_chars().eq('').any():
        raise ValueError('source counts nonnegative, control weights positive and sample identities nonempty required')
    # Event-slot sampler has the established fewer-than-one-billion/bin contract.
    if ctrl['event_weight'].max()>=10**9:raise ValueError('control event bin exceeds supported slot budget')
    obs=obs.with_columns(pl.col('query_id').cast(pl.Int64),pl.col('pgen_bin').cast(pl.Int64))
    ctrl=ctrl.with_columns(pl.col('query_id','pgen_bin').cast(pl.Int64)).sort('pgen_bin','query_id')
    inventory=targets.unique(maintain_order=True)
    if measure is None:
        inventory=inventory.sort(columns,nulls_last=True)
    inventory=inventory.with_row_index('_target')
    if measure is None:
        q=pl.concat([obs.select('query_id'),ctrl.select('query_id')]).unique().sort('query_id').with_row_index('_row')
        obs=obs.join(q,on='query_id',validate='m:1')
        ctrl=ctrl.join(q,on='query_id',validate='m:1')
        member=matches.select(pl.col('query_id').cast(pl.Int64),*columns).unique().join(q,on='query_id',how='left',validate='m:1').join(inventory,on=columns,how='left',validate='m:1',nulls_equal=True)
        if member['_row'].null_count() or member['_target'].null_count():raise ValueError('membership contains unknown query or target')
        matrix=sparse.coo_matrix((np.ones(member.height,dtype=np.int64),(member['_row'].to_numpy(),member['_target'].to_numpy())),shape=(q.height,inventory.height)).tocsr()
        def measure(rows,counts):
            matched=matrix[rows]
            return np.asarray(matched.sum(axis=0)).ravel(),np.asarray(matched.T@counts).ravel()
    else:
        obs=obs.with_columns(pl.col('query_id').alias('_row'))
        ctrl=ctrl.with_columns(pl.col('query_id').alias('_row'))
    obs=obs.sort('sample_id','pgen_bin','query_id')
    pools={}
    for (bin_id,),part in ctrl.partition_by('pgen_bin',as_dict=True).items():
        cumulative=np.cumsum(part['event_weight'].to_numpy(),dtype=np.int64)
        if cumulative[-1]>=10**9:raise ValueError('control event bin exceeds supported slot budget')
        pools[bin_id]=(part['_row'].to_numpy(),cumulative)
    rng=np.random.default_rng(seed);summaries=[];mock_tables=[];audits=[]
    for (sample,),part in obs.partition_by('sample_id',as_dict=True).items():
        counts=part['count'].to_numpy();total=sum(map(int,counts))
        if total>np.iinfo(np.int64).max:raise ValueError('sample read total exceeds int64 count budget')
        rows=part['_row'].to_numpy();plan=[]
        for (bin_id,),group in part.with_row_index('_position').partition_by('pgen_bin',as_dict=True).items():
            available=int(pools[bin_id][1][-1]) if bin_id in pools else 0
            if available<group.height:raise ValueError(f'Pgen bin {bin_id}: need{group.height}control events, have{available}; exact matching cannot remap or discard source records')
            plan.append((group['_position'].to_numpy(),pools[bin_id]))
            audits.append(dict(sample_id=sample,pgen_bin=bin_id,requested_records=group.height,available_control_events=available,mock_replicates=replicates,exact_histogram_and_counts=True))
        actual_records,actual_reads=measure(rows,counts)
        mock_records=np.empty((replicates,inventory.height),dtype=np.int64);mock_reads=np.empty_like(mock_records)
        draw=np.empty(len(rows),dtype=rows.dtype)
        for replicate in range(replicates):
            for positions,(events,cumulative) in plan:
                selected=events[_sample_event_indices(cumulative,len(positions),rng)];rng.shuffle(selected)
                draw[positions]=selected
            mock_records[replicate],mock_reads[replicate]=measure(draw,counts)
        fields=dict(sample_id=pl.repeat(sample,inventory.height,eager=True),n_records=pl.repeat(len(rows),inventory.height,eager=True),n_reads=pl.repeat(total,inventory.height,eager=True),matched_records=actual_records,matched_reads=actual_reads)
        for endpoint,values,observed in [('records',mock_records,actual_records),('reads',mock_reads,actual_reads)]:
            mean=values.mean(axis=0)
            fields['excess_'+endpoint]=observed-mean
            fields['enrichment_'+endpoint]=np.divide(observed,mean,out=np.full(inventory.height,np.nan),where=mean>0)
            fields.update({'mock_mean_'+endpoint:mean,'mock_lower_'+endpoint:np.quantile(values,.025,axis=0),'mock_upper_'+endpoint:np.quantile(values,.975,axis=0),'p_'+endpoint:(1+np.count_nonzero(values>=observed,axis=0))/(replicates+1)})
        summaries.append(pl.concat([inventory.drop('_target'),pl.DataFrame(fields)],how='horizontal'))
        mock=pl.DataFrame(dict(_target=np.tile(np.arange(inventory.height),replicates),replicate=np.repeat(np.arange(replicates),inventory.height),mock_records=mock_records.ravel(),mock_reads=mock_reads.ravel())).with_columns(pl.col('_target').cast(inventory.schema['_target']),pl.lit(sample).alias('sample_id'))
        mock_tables.append(mock.join(inventory,on='_target',validate='m:1').drop('_target'))
    return pl.concat(summaries),pl.concat(mock_tables),pl.DataFrame(audits)


def reference_enrichment(observations: pl.DataFrame, reference: pl.DataFrame,
                         controls: dict[str, pl.DataFrame], *, replicates: int = 200,
                         seed: int = 0, policy: str = 'pgen', threads: int = 1,
                         mismatches: int = 1, match_j: bool = False, chain_mode: str = "union",
                         tra_matcher: str = "hamming", trb_matcher: str = "hamming", tra_min_score: float = 0., trb_min_score: float = 0., species: str = "human", target_groups=None, progress=None, kernel_scorer="v-bounded", tra_radius_factor=None, trb_radius_factor=None):
    """Test a fixed repertoire against fresh Pgen-matched mock target corpora.

    Observations: receptor_id/locus/junction_aa/v/j/count. Multiple chains with
    one receptor_id represent one receptor and must carry the same abundance.
    Reference: group_id/locus/junction_aa/v/j/pgen_bin. Controls are normalized
    unique junction/V/J keys with pgen_bin per locus. chain_mode selects separate
    tra/trb, union or intersection support; paired receptors count once. Linked and unlinked
    references both contribute. This calibrates the selected rule, not a paired score.
    Counts are unique receptor identities and their summed abundance separately.
    Returns summary, complete aggregate mock endpoints and real chain hit links.
    Full-junction Hamming matching uses one native batch per locus/population,
    with exact V and optional J gene matching. Each chain can instead use stem4
    or gapblock, accepting dimensionless scores strictly above its nonnegative min_score
    (zero by default), identically in every real/mock panel. Its
    chain_score is retained alongside strict-chain indicators, not a probability.
    chain_mode=all reuses searches and mock draws for all four reducers; outputs
    then include chain_mode. No query Pgen is required.
    Optional sample_id batches an in-memory cohort with identical seeded mock
    panels and shared native searches; receptor IDs remain local to each sample.
    target_groups restricts reported targets, never reference models or mock draws.
    Optional epitope/mhc_a/mhc_b/mhc_class columns preserve peptide identity
    across restrictions in the fuzzy background; otherwise group IDs are labels.
    Fuzzy models use production V geometry imputation, retaining raw V keys.
    Graded unit/v-unit/weighted/v-weighted families instead threshold the nearest
    individual alignment cost (unit scale1, weighted scale400) by radius_factor.
    They reuse the graded edit/indel ceiling and never sum weak hits for acceptance.
    """
    from seqtree import Index
    from ..match.scope import search_params, blosum_radius_params, graded_family_params, GRADED_FAMILIES
    from ..evalue.mock_reference import iter_epitope_mocks
    from .candidates import PMHC

    if type(replicates) is not int or replicates < 1 or type(seed) is not int or seed < 0:
        raise ValueError('positive replicates and nonnegative seed required')
    from vdjtools.sc import chain_support
    if chain_mode not in {"tra", "trb", "union", "intersection", "all"}:
        raise ValueError("chain_mode must be tra, trb, union, intersection or all")
    matchers={"TRA":tra_matcher,"TRB":trb_matcher}
    weighted={'blosum62','positional-blosum62'}
    if any(m not in {'hamming','stem4','gapblock',*weighted,*GRADED_FAMILIES} for m in matchers.values()):raise ValueError('unknown chain matcher')
    factors={'TRA':tra_radius_factor,'TRB':trb_radius_factor}
    family_params={}
    for locus,matcher in matchers.items():
        if matcher in weighted:blosum_radius_params(matcher,factors[locus],1)
        elif matcher in GRADED_FAMILIES:
            family_params[locus]=graded_family_params(matcher,mismatches,factors[locus])
        elif factors[locus] is not None:raise ValueError('radius factor requires a BLOSUM matcher')
    cutoffs={'TRA':tra_min_score,'TRB':trb_min_score}
    if any(not np.isfinite(c) or c<0 for c in cutoffs.values()):raise ValueError('finite nonnegative fuzzy score cutoffs required')
    if any(cutoffs[locus] and matcher not in {'stem4','gapblock'} for locus,matcher in matchers.items()):raise ValueError('score cutoff requires a fuzzy matcher')
    from ..evalue.kernel import KERNEL_SCORERS
    if kernel_scorer not in KERNEL_SCORERS:raise ValueError("unknown kernel scorer")
    fuzzy_loci=[locus for locus,matcher in matchers.items() if matcher in {'stem4','gapblock'}]
    family_loci=[locus for locus,matcher in matchers.items() if matcher in GRADED_FAMILIES]
    anchored=chain_mode=="intersection" and len(fuzzy_loci)==1
    required={'receptor_id','locus','junction_aa','v','j','count'}
    if not required <= set(observations.columns) or not observations.height:
        raise ValueError('nonempty observations with receptor identities, chain keys and counts required')
    if any(observations[c].null_count() for c in required):
        raise ValueError('non-null receptor keys and nonnegative integer counts required')
    cohort = 'sample_id' in observations.columns
    if cohort:
        if observations['sample_id'].null_count() or observations['sample_id'].cast(pl.String).str.strip_chars().eq('').any():
            raise ValueError('sample_id cannot be missing')
        identities = observations.select('sample_id','receptor_id').unique(maintain_order=True).with_row_index('_cohort_id')
        observations = observations.join(identities,on=['sample_id','receptor_id'],validate='m:1').drop('receptor_id').rename({'_cohort_id':'receptor_id'})
    obs=observations.select(sorted(required)).with_columns(_integer_count('count').alias('count'))
    if any(obs[c].null_count() for c in required) or (obs['count']<0).any():
        raise ValueError('non-null receptor keys and nonnegative integer counts required')
    totals=obs.group_by('receptor_id').agg(pl.col('count').n_unique().alias('_n'),pl.col('count').first())
    if (totals['_n']!=1).any():raise ValueError('all chains of a receptor must have the same count')
    total_count=sum(totals['count'].to_list())
    if total_count>np.iinfo(np.int64).max:raise ValueError('total count exceeds int64')
    if type(threads) is not int or threads<1 or type(mismatches) is not int or mismatches not in (0,1,2,3):
        raise ValueError('positive threads and Hamming mismatches from 0 to 3 required')
    if not {'group_id','locus','junction_aa','v','j','pgen_bin'}<=set(reference.columns) or not reference.height:
        raise ValueError('nonempty reference with target, locus, receptor keys and Pgen bins required')
    if set(PMHC)&set(reference.columns):
        if not set(PMHC)<=set(reference.columns):raise ValueError('complete pMHC metadata required')
        labels=reference.select('group_id',pl.col(PMHC).cast(pl.String)).unique()
        if labels.height!=reference['group_id'].n_unique() or labels.select(PMHC).n_unique()!=labels.height:
            raise ValueError('group IDs and full pMHC targets must map one-to-one')
    else:
        labels=reference.select('group_id').unique().with_columns(pl.col('group_id').cast(pl.String).alias('epitope'),
            pl.lit(None,dtype=pl.String).alias('mhc_a'),pl.lit(None,dtype=pl.String).alias('mhc_b'),pl.lit('MHCI').alias('mhc_class'))
    ref=reference.select('group_id','locus','junction_aa','v','j','pgen_bin').unique().sort('group_id','locus','junction_aa','v','j')
    for frame in [obs,ref]:
        if not set(frame['locus'])<={'TRA','TRB'} or frame['locus'].null_count():raise ValueError('TRA/TRB loci required')
        if not frame['junction_aa'].str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).all():
            raise ValueError('valid full-junction amino acid sequences required')
    modes=['tra','trb','union','intersection'] if chain_mode=='all' else [chain_mode]
    inventory=ref.select('group_id').unique().sort('group_id')
    if target_groups is not None:
        if not target_groups or not set(target_groups)<=set(inventory['group_id']):
            raise ValueError('target_groups must be a nonempty subset of reference groups')
        inventory=inventory.filter(pl.col('group_id').is_in(target_groups))
    if chain_mode=='all':inventory=inventory.join(pl.DataFrame({'chain_mode':modes}),how='cross')
    if cohort:
        inventory=inventory.join(identities.select('sample_id').unique().sort('sample_id'),how='cross')
        totals=totals.join(identities.select(pl.col('_cohort_id').alias('receptor_id'),'sample_id'),on='receptor_id',validate='1:1')
    inventory=inventory.with_row_index('_target')
    indexes={};samplers={};unavailable=[]
    for (locus,),part in obs.partition_by('locus',as_dict=True).items():
        keys=part.select('junction_aa').unique().sort('junction_aa').with_row_index('_ri')
        indexed=part.join(keys,on='junction_aa',validate='m:1').rename({'v':'_sample_v','j':'_sample_j'}).with_columns(
            pl.col('_sample_v').str.replace(r'\*.*$','').alias('_match_v'),pl.col('_sample_j').str.replace(r'\*.*$','').alias('_match_j'))
        indexes[locus]=(None if locus in family_loci else Index.build(keys['junction_aa'].to_list(),alphabet='aa'),indexed)
    for (locus,),part in ref.partition_by('locus',as_dict=True).items():
        if locus not in controls:raise ValueError(f'missing {locus} control pool')
        samplers[locus]=iter_epitope_mocks(part,controls[locus],replicates=replicates,seed=seed+(0 if locus=='TRA' else 1),policy=policy,unavailable=unavailable)
    empty=pl.DataFrame(schema={'receptor_id':obs.schema['receptor_id'],'group_id':ref.schema['group_id'],'locus':pl.String,'junction_aa':pl.String,'v':pl.String,'j':pl.String})
    fuzzy_by_locus={locus:obs.filter(pl.col('locus')==locus).with_row_index('query_id').rename({'junction_aa':'cdr3'}) for locus in fuzzy_loci+family_loci}
    weighted_params={(locus,length):blosum_radius_params(matcher,factors[locus],length)
        for locus,matcher in matchers.items() if matcher in weighted
        for length in obs.filter(pl.col('locus')==locus)['junction_aa'].str.len_chars().unique()}
    contexts = None
    if any(fuzzy_by_locus[locus].height for locus in fuzzy_loci):
        from ..match.geometry import precompute_v_loop_distances
        matrices={'blosum62' if matchers[locus]=='stem4' else 'clipped' for locus in fuzzy_loci}
        contexts={matrix:precompute_v_loop_distances(species=species,alphabet='full',matrix=matrix,v_regions='cdr12plus') for matrix in sorted(matrices)}
    if fuzzy_loci:empty=empty.with_columns(pl.lit(None,dtype=pl.Float64).alias('chain_score'))
    def measure(population, *, detailed=False):
        parts=[]
        if anchored:
            population=population.with_columns((pl.col('locus')==fuzzy_loci[0]).alias('_fuzzy_order')).sort('_fuzzy_order',maintain_order=True).drop('_fuzzy_order')
        for (locus,),targets in population.partition_by('locus',as_dict=True).items():
            if locus not in indexes:continue
            index,samples=indexes[locus]
            if locus in family_loci:
                from ..api import Annotator
                queries=fuzzy_by_locus[locus]
                params,scale=family_params[locus]
                model=targets.rename({'junction_aa':'cdr3','locus':'gene'}).join(labels,on='group_id',validate='m:1')
                scores=Annotator.from_frame(model).target_memberships(queries,
                    scope=params,threads=threads,max_score=factors[locus],
                    match_v=matchers[locus].startswith('v-'),match_j=match_j,
                    species=species,score_scale=scale)
                # One qualifying alignment suffices; weak alignments never add
                # together to cross the distance threshold.
                supported=scores.join(
                    labels,on=PMHC,nulls_equal=True,validate='m:1').join(queries,on='query_id',validate='m:1')
                if target_groups is not None:
                    supported=supported.filter(pl.col('group_id').is_in(target_groups))
                columns=['receptor_id','group_id','locus']
                if detailed:
                    supported=supported.rename({'cdr3':'junction_aa'})
                    columns+=['junction_aa','v','j']
                    if fuzzy_loci:
                        supported=supported.with_columns(pl.lit(1.).alias('chain_score'))
                        columns+=['chain_score']
                parts.append(supported.select(columns).unique())
                continue
            if locus in fuzzy_by_locus:
                fuzzy_queries=fuzzy_by_locus[locus]
                kernel=matchers[locus]
                from ..runner.combined import _kernel
                scored=[]
                kernel_reference=targets.rename({'junction_aa':'cdr3'}).join(labels,on='group_id',validate='m:1')
                def collect(name,frame):
                    scored.append(frame.filter(pl.col('kernel_score')>cutoffs[locus]).join(labels,on=PMHC,nulls_equal=True,validate='m:1').select('query_id','group_id',pl.col('kernel_score').alias('chain_score')))
                cases=None
                if anchored:
                    if not parts:continue
                    cases=pl.concat(parts).select('receptor_id','group_id').unique().join(
                        fuzzy_queries.select('query_id','receptor_id'),on='receptor_id',validate='m:m').join(
                        labels,on='group_id',validate='m:1').select('query_id',*PMHC)
                elif target_groups is not None:
                    cases=fuzzy_queries.select('query_id').join(
                        inventory.select('group_id').unique().join(labels,on='group_id',validate='1:1').select(PMHC),how='cross')
                if cases is not None:
                    if not cases.height:continue
                    queries=fuzzy_queries.join(cases.select('query_id').unique(),on='query_id',how='semi')
                    result=_kernel(queries,kernel_reference,cases,species=species,threads=threads,paired=False,
                        exclude_exact=False,contexts=contexts,kernels=(kernel,),scorer=kernel_scorer,_impute_v=True)
                    collect(kernel,result.rename({kernel+'_score':'kernel_score'}))
                else:
                    _kernel(fuzzy_queries,kernel_reference,None,species=species,threads=threads,paired=False,
                        exclude_exact=False,contexts=contexts,_consume=collect,kernels=(kernel,),scorer=kernel_scorer,_impute_v=True)
                if scored:
                    supported=pl.concat(scored).join(fuzzy_queries,on='query_id',validate='m:1').rename({'cdr3':'junction_aa'})
                    columns=['receptor_id','group_id','locus']+(['junction_aa','v','j','chain_score'] if detailed else [])
                    parts.append(supported.select(columns).unique())
                continue
            if target_groups is not None:
                targets=targets.filter(pl.col('group_id').is_in(target_groups))
                if not targets.height:continue
            keys=targets.select('junction_aa').unique().sort('junction_aa').with_row_index('_qi')
            roles=targets.join(keys,on='junction_aa',validate='m:1').select('_qi','group_id',
                pl.col('v').str.replace(r'\*.*$','').alias('_match_v'),pl.col('j').str.replace(r'\*.*$','').alias('_match_j'))
            buffers=[]
            def consume(start,stop,buffer):
                pairs=np.asarray(buffer)
                if len(pairs):
                    pairs=pairs[:,:2].copy()
                    pairs[:,0]=query_ids[pairs[:,0]]
                    buffers.append(pairs)
            batches=(keys.with_columns(pl.col('junction_aa').str.len_chars().alias('_length')).partition_by('_length',as_dict=True).values() if matchers[locus] in weighted else [keys])
            for batch in batches:
                query_ids=batch['_qi'].to_numpy()
                params=(weighted_params.get((locus,len(batch['junction_aa'][0]))) if matchers[locus] in weighted else search_params(f'{mismatches},0,0,{mismatches}'))
                if params is None:continue  # No observed junction has this length.
                index.search_batch_consume(batch['junction_aa'].to_list(),params,consume,
                    threads=threads,chunk_size=2048,max_hit_bytes=64*1024**2)
            # Join sparse hits once: native memory bounds can make callbacks tiny.
            if buffers:
                pairs=pl.DataFrame(np.concatenate(buffers),schema=[('_qi',pl.UInt32),('_ri',pl.UInt32)],orient='row')
                predicate=['_ri','_match_v']+(['_match_j'] if match_j else [])
                hits=pairs.join(roles,on='_qi',validate='m:m').join(samples,on=predicate,validate='m:m')
                if detailed:
                    detail=hits.select('receptor_id','group_id',pl.lit(locus).alias('locus'),'junction_aa',pl.col('_sample_v').alias('v'),pl.col('_sample_j').alias('j')).unique()
                    if fuzzy_loci:detail=detail.with_columns(pl.lit(1.).alias('chain_score'))
                    parts.append(detail)
                else:
                    parts.append(hits.select('receptor_id','group_id').with_columns(pl.lit(locus).alias('locus')).unique())
        members=pl.concat(parts).unique() if parts else empty
        support=chain_support(members,by=['receptor_id','group_id'])
        aggregates=[];details=[]
        for mode in modes:
            selected=support.filter(pl.col(mode+'_supported')).select('receptor_id','group_id')
            aggregate=selected.join(totals.select('receptor_id','count',*(['sample_id'] if cohort else [])),on='receptor_id',validate='m:1').group_by('group_id',*(['sample_id'] if cohort else [])).agg(pl.len().cast(pl.Int64).alias('receptors'),pl.col('count').sum())
            if chain_mode=='all':aggregate=aggregate.with_columns(pl.lit(mode).alias('chain_mode'))
            aggregates.append(aggregate)
            if detailed:
                detail=members.join(selected,on=['receptor_id','group_id'],how='inner')
                if mode in {'tra','trb'}:detail=detail.filter(pl.col('locus')==mode.upper())
                if chain_mode=='all':detail=detail.with_columns(pl.lit(mode).alias('chain_mode'))
                details.append(detail)
        keys=['group_id','chain_mode'] if chain_mode=='all' else ['group_id']
        if cohort:keys.append('sample_id')
        hits=pl.concat(details).sort('receptor_id','group_id','locus') if detailed else empty
        values=inventory.join(pl.concat(aggregates),on=keys,how='left',validate='1:1').fill_null(0).sort('_target')
        return values['receptors'].to_numpy(),values['count'].to_numpy(),hits
    actual,abundance,hits=measure(ref,detailed=True)
    if progress is not None:progress('real',0)
    null=np.empty((replicates,inventory.height,2),dtype=np.int64)
    for b in range(replicates):
        population=pl.concat([next(sampler) for sampler in samplers.values()])
        if b==0 and progress is not None:progress('first_mock_draw',1)
        null[b,:,0],null[b,:,1],_=measure(population)
        if progress is not None and ((b+1)%20==0 or b+1==replicates):progress('mock',b+1)
    if cohort:
        denominators=totals.group_by('sample_id').agg(pl.len().alias('n_receptors'),pl.col('count').sum().alias('total_count'))
        summary=inventory.drop('_target').join(denominators,on='sample_id',validate='m:1',maintain_order='left')
        hits=hits.rename({'receptor_id':'_cohort_id'}).join(identities,on='_cohort_id',validate='m:1').drop('_cohort_id')
    else:
        summary=inventory.drop('_target').with_columns(pl.lit(totals.height).alias('n_receptors'),pl.lit(total_count).alias('total_count'))
    for name,observed,values in [('receptors',actual,null[:,:,0]),('count',abundance,null[:,:,1])]:
        mean=values.mean(axis=0)
        ratio=np.divide(observed,mean,out=np.full(len(mean),np.nan),where=mean>0)
        summary=summary.with_columns(pl.Series('matched_'+name,observed),pl.Series('mock_mean_'+name,mean),
            pl.Series('mock_lower_'+name,np.quantile(values,.025,axis=0)),pl.Series('mock_upper_'+name,np.quantile(values,.975,axis=0)),
            pl.Series('excess_'+name,observed-mean),pl.Series('enrichment_'+name,ratio).fill_nan(None),
            pl.Series('p_'+name,(1+(values>=observed).sum(axis=0))/(replicates+1)))
    mocks=pl.DataFrame(dict(_target=np.tile(np.arange(inventory.height,dtype=np.uint32),replicates),
        replicate=np.repeat(np.arange(replicates),inventory.height),mock_receptors=null[:,:,0].ravel(),mock_count=null[:,:,1].ravel()))
    mocks=mocks.join(inventory,on='_target',validate='m:1').drop('_target')
    reasons={}
    for row in unavailable:reasons.setdefault(row['group_id'],[]).append(row['reason'])
    audit=inventory.with_columns(
        pl.Series('status',['untestable_null' if g in reasons else 'ok' for g in inventory['group_id']]),
        pl.Series('null_unavailable_reason',['; '.join(dict.fromkeys(reasons[g])) if g in reasons else None for g in inventory['group_id']],dtype=pl.String))
    audit=audit.select('group_id','status','null_unavailable_reason').unique()
    summary=summary.join(audit,on='group_id',validate='m:1',maintain_order='left')
    inferred=[c for c in summary.columns if c.startswith(('mock_','excess_','enrichment_','p_'))]
    summary=summary.with_columns(pl.when(pl.col('status')=='ok').then(pl.col(c)).otherwise(None).alias(c) for c in inferred)
    mocks=mocks.with_columns(pl.when(~pl.col('group_id').is_in(list(reasons))).then(pl.col(c)).otherwise(None).alias(c) for c in ['mock_receptors','mock_count'])
    return summary,mocks,hits


def kernel_reference_counts(observations, panels, inventory, full_reference, *,
                            locus, species, min_score, threads=1, progress=None):
    """Count original real/mock panels with the fixed stem4 dimensionless score.

    One model per panel scores all cohort queries together. Other-peptide
    reference records stay fixed. No query is removed from any model.
    """
    from ..aggregate.candidates import PMHC
    from ..runner.combined import _kernel, _prepare_kernel_chain
    from ..match.geometry import (precompute_v_loop_distances, grouped_total_distance_count_batch,
                                 distance_units, distance_cutoff, _substitution_costs)
    if not np.isfinite(min_score) or min_score<=0:
        raise ValueError('stem4 requires an explicit positive dimensionless score threshold')
    keys=['junction_aa','v','j']
    counted=observations.group_by('sample_id',*keys).agg(
        pl.len().cast(pl.Int64).alias('matched_records'),pl.col('count').sum().alias('matched_reads'))
    queries=counted.select(keys).unique().sort(keys).with_row_index('query_id')
    counted=counted.join(queries,on=keys,validate='m:1')
    queries=queries.rename({'junction_aa':'cdr3'}).with_columns(pl.lit(locus).alias('locus'))
    targets=inventory.select('group_id',*PMHC)
    background=full_reference.filter(~pl.col('epitope').is_in(targets['epitope'].to_list()))
    contexts={'blosum62':precompute_v_loop_distances(species=species,alphabet='full',matrix='blosum62',v_regions='cdr12plus')}
    eligible=None
    if min_score>np.exp(-1.):
        # S <= c/N + (1-c/N)*exp(-1), with c target keys inside the fixed
        # radius (= kernel scale). Model-union c and smallest target N loosen
        # this bound. Only queries unable to pass are skipped, never references.
        search_keys=['cdr3','_allele','_v_key']
        models=panels.select('replicate').unique().sort('replicate').with_row_index('_model')
        model_axis=dict(models.select('replicate','_model').iter_rows())
        prepared=_prepare_kernel_chain(panels.rename({'junction_aa':'cdr3'}),contexts['blosum62']['model'],locus)
        q=_prepare_kernel_chain(queries,contexts['blosum62']['model'],locus)
        eligible=np.zeros((models.height,queries.height),dtype=bool)
        if prepared.height and q.height:
            sizes=prepared.select('replicate','group_id',*search_keys).unique().group_by('replicate','group_id').len().group_by('replicate').agg(pl.col('len').min())
            minimum=models.join(sizes,on='replicate',how='left',maintain_order='left')['len'].cast(pl.Float64).fill_null(float('inf')).to_numpy()
            whole=prepared.select(search_keys).unique().with_row_index('_ref')
            members=prepared.select('replicate',*search_keys).unique().join(models,on='replicate',validate='m:1').join(whole,on=search_keys,validate='m:1').sort('_ref','_model')
            lengths=whole.select('_ref').join(members.group_by('_ref').len(),on='_ref',how='left',maintain_order='left')['len'].to_numpy()
            offsets=np.concatenate((np.zeros(1,dtype=np.uint64),np.cumsum(lengths,dtype=np.uint64)))
            scale,_=distance_units(_substitution_costs('full','blosum62'))
            axes=lambda frame:frame.select(pl.col('cdr3').alias('seq0'),pl.col('_allele').alias('v0'))
            chunks=grouped_total_distance_count_batch(axes(q),axes(whole),[distance_cutoff(22.5*scale)],
                species=species,threads=threads,ctrim=3,v_loop_weight=1,alphabet='full',stem_depth=4,
                stem_weight=1.,matrix='blosum62',v_loop_distances=contexts['blosum62'],
                reference_target_offsets=offsets,reference_target_ids=members['_model'].to_numpy(),
                target_peptide_ids=np.zeros(models.height,dtype=np.int32),
                reference_exclusive_peptide_ids=np.zeros(whole.height,dtype=np.int32),
                chunk_size=min(65536,max(1,256*1024**2//((2*models.height+2)*4))))
            for chunk in chunks:
                close=np.asarray(chunk.counts)[:,1:1+models.height,0]
                upper=np.exp(-1.)+(1-np.exp(-1.))*close/minimum[None,:]
                eligible[:,q['query_id'].slice(chunk.start,chunk.stop-chunk.start).to_numpy()]=(upper>=min_score).T
        del prepared,q
    results=[]
    for (replicate,),panel in panels.partition_by('replicate',as_dict=True).items():
        replacement=panel.join(targets,on='group_id',validate='m:1').select(
            pl.col('junction_aa').alias('cdr3'),'v','j',pl.lit(locus).alias('locus'),*PMHC)
        model=pl.concat([background.select(replacement.columns),replacement])
        members=[]
        def consume(name,scores):
            members.append(scores.filter(pl.col('kernel_score')>min_score).select('query_id',*PMHC))
        active=queries if eligible is None else queries.filter(pl.Series(eligible[model_axis[replicate]]))
        if active.height:
            _kernel(active,model,None,species=species,threads=threads,paired=False,
                exclude_exact=False,contexts=contexts,kernels=('stem4',),_consume=consume,_targets=targets.select(PMHC),_impute_v=True)
        supported=pl.concat(members) if members else queries.head(0).select('query_id').join(targets.select(PMHC),how='cross')
        support=supported.join(targets,on=PMHC,nulls_equal=True,validate='m:1').select('query_id','group_id').unique()
        matched=support.join(counted,on='query_id',validate='m:m').group_by('sample_id','group_id').agg(pl.col('matched_records','matched_reads').sum())
        complete=targets.select('group_id').join(observations.select('sample_id').unique(),how='cross').join(
            matched,on=['sample_id','group_id'],how='left',validate='1:1').with_columns(
            pl.col('matched_records','matched_reads').fill_null(0),pl.lit(replicate).alias('replicate'))
        results.append(complete)
        if progress is not None and (replicate==-1 or (replicate+1)%20==0):progress(replicate)
    return pl.concat(results).sort('sample_id','group_id','replicate')
