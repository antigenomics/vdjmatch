"""Strict Pgen-bin/count-matched repertoire enrichment from batched memberships."""
from __future__ import annotations

import numpy as np
import polars as pl
from scipy import sparse

from ..evalue.mock_reference import _sample_event_indices
from ..io.columns import _integer_count


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
                         mismatches: int = 1, match_j: bool = False, progress=None):
    """Test a fixed repertoire against fresh Pgen-matched mock target corpora.

    Observations: receptor_id/locus/junction_aa/v/j/count. Multiple chains with
    one receptor_id represent one receptor and must carry the same abundance.
    Reference: group_id/locus/junction_aa/v/j/pgen_bin. Controls are normalized
    unique junction/V/J keys with pgen_bin per locus. Target support is the union
    of available chain matches; paired receptors count once. Linked and unlinked
    references both contribute. This calibrates that union, not a paired score.
    Counts are unique receptor identities and their summed abundance separately.
    Returns summary, complete aggregate mock endpoints and real chain hit links.
    Full-junction Hamming matching uses one native batch per locus/population,
    with exact V and optional J gene matching. No query Pgen is required.
    """
    from seqtree import Index
    from ..match.scope import search_params
    from ..evalue.mock_reference import iter_epitope_mocks

    if type(replicates) is not int or replicates < 1 or type(seed) is not int or seed < 0:
        raise ValueError('positive replicates and nonnegative seed required')
    required={'receptor_id','locus','junction_aa','v','j','count'}
    if not required <= set(observations.columns) or not observations.height:
        raise ValueError('nonempty observations with receptor identities, chain keys and counts required')
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
    ref=reference.select('group_id','locus','junction_aa','v','j','pgen_bin').unique().sort('group_id','locus','junction_aa','v','j')
    for frame in [obs,ref]:
        if not set(frame['locus'])<={'TRA','TRB'} or frame['locus'].null_count():raise ValueError('TRA/TRB loci required')
        if not frame['junction_aa'].str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).all():
            raise ValueError('valid full-junction amino acid sequences required')
    inventory=ref.select('group_id').unique().sort('group_id').with_row_index('_target')
    indexes={};samplers={};unavailable=[]
    for (locus,),part in obs.partition_by('locus',as_dict=True).items():
        keys=part.select('junction_aa').unique().sort('junction_aa').with_row_index('_ri')
        indexed=part.join(keys,on='junction_aa',validate='m:1').rename({'v':'_sample_v','j':'_sample_j'}).with_columns(
            pl.col('_sample_v').str.replace(r'\*.*$','').alias('_match_v'),pl.col('_sample_j').str.replace(r'\*.*$','').alias('_match_j'))
        indexes[locus]=(Index.build(keys['junction_aa'].to_list(),alphabet='aa'),indexed)
    for (locus,),part in ref.partition_by('locus',as_dict=True).items():
        if locus not in controls:raise ValueError(f'missing {locus} control pool')
        samplers[locus]=iter_epitope_mocks(part,controls[locus],replicates=replicates,seed=seed+(0 if locus=='TRA' else 1),policy=policy,unavailable=unavailable)
    empty=pl.DataFrame(schema={'receptor_id':obs.schema['receptor_id'],'group_id':ref.schema['group_id'],'locus':pl.String,'junction_aa':pl.String,'v':pl.String,'j':pl.String})
    def measure(population, *, detailed=False):
        parts=[]
        for (locus,),targets in population.partition_by('locus',as_dict=True).items():
            if locus not in indexes:continue
            index,samples=indexes[locus]
            keys=targets.select('junction_aa').unique().sort('junction_aa').with_row_index('_qi')
            roles=targets.join(keys,on='junction_aa',validate='m:1').select('_qi','group_id',
                pl.col('v').str.replace(r'\*.*$','').alias('_match_v'),pl.col('j').str.replace(r'\*.*$','').alias('_match_j'))
            buffers=[]
            def consume(start,stop,buffer):
                pairs=np.asarray(buffer)
                if len(pairs):buffers.append(pairs[:,:2].copy())
            index.search_batch_consume(keys['junction_aa'].to_list(),search_params(f'{mismatches},0,0,{mismatches}'),consume,
                threads=threads,chunk_size=2048,max_hit_bytes=64*1024**2)
            # Join sparse hits once: native memory bounds can make callbacks tiny.
            if buffers:
                pairs=pl.DataFrame(np.concatenate(buffers),schema=[('_qi',pl.UInt32),('_ri',pl.UInt32)],orient='row')
                predicate=['_ri','_match_v']+(['_match_j'] if match_j else [])
                hits=pairs.join(roles,on='_qi',validate='m:m').join(samples,on=predicate,validate='m:m')
                if detailed:
                    parts.append(hits.select('receptor_id','group_id',pl.lit(locus).alias('locus'),'junction_aa',pl.col('_sample_v').alias('v'),pl.col('_sample_j').alias('j')).unique())
                else:
                    parts.append(hits.select('receptor_id','group_id').unique())
        members=pl.concat(parts).unique() if parts else empty
        hits=members.sort('receptor_id','group_id','locus') if detailed else empty
        aggregate=members.select('receptor_id','group_id').unique().join(totals.select('receptor_id','count'),on='receptor_id',validate='m:1').group_by('group_id').agg(pl.len().cast(pl.Int64).alias('receptors'),pl.col('count').sum())
        values=inventory.join(aggregate,on='group_id',how='left',validate='1:1').fill_null(0).sort('_target')
        return values['receptors'].to_numpy(),values['count'].to_numpy(),hits
    actual,abundance,hits=measure(ref,detailed=True)
    if progress is not None:progress('real',0)
    null=np.empty((replicates,inventory.height,2),dtype=np.int64)
    for b in range(replicates):
        population=pl.concat([next(sampler) for sampler in samplers.values()])
        if b==0 and progress is not None:progress('first_mock_draw',1)
        null[b,:,0],null[b,:,1],_=measure(population)
        if progress is not None and ((b+1)%20==0 or b+1==replicates):progress('mock',b+1)
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
    audit=inventory.select('group_id').with_columns(
        pl.Series('status',['untestable_null' if g in reasons else 'ok' for g in inventory['group_id']]),
        pl.Series('null_unavailable_reason',['; '.join(dict.fromkeys(reasons[g])) if g in reasons else None for g in inventory['group_id']],dtype=pl.String))
    summary=summary.join(audit,on='group_id',validate='1:1')
    inferred=[c for c in summary.columns if c.startswith(('mock_','excess_','enrichment_','p_'))]
    summary=summary.with_columns(pl.when(pl.col('status')=='ok').then(pl.col(c)).otherwise(None).alias(c) for c in inferred)
    mocks=mocks.with_columns(pl.when(~pl.col('group_id').is_in(list(reasons))).then(pl.col(c)).otherwise(None).alias(c) for c in ['mock_receptors','mock_count'])
    return summary,mocks,hits
