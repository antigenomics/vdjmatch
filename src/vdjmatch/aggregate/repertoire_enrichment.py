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
    inventory=targets.unique().sort(columns,nulls_last=True).with_row_index('_target')
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
