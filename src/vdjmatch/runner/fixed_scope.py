"""Naive bounded unit-edit neighbours, with source and pairing units explicit."""
from __future__ import annotations

import polars as pl
from seqtree import Index

from ..aggregate.candidates import PMHC
from ..match.scope import search_params


def fixed_scope_scores(queries, reference, *, species='human', paired=False,
                       scope='1,0,0,1', exclude_exact=False, threads=1):
    """Score every declared query/pMHC using one native batch per present chain.

    Queries are canonical long junction/V/J/locus rows with full pMHC metadata;
    paired mode uses declared ``pair_id`` and retains incomplete pairs. No V/J,
    abundance, background or learned weight changes a unit-edit neighbour.
    Compatible pMHC unions are formed before counting unique junctions or joint
    junction/V/J keys. Linked references are genuine complete complexes, with
    their two chain budgets applied separately. Linked nearest radius is the
    minimum bottleneck ``max(d_alpha,d_beta)``; marginal distances remain separate.
    """
    from ..cli.historical import assayed_targets
    from ..evalue.paired import build_paired_ref

    if threads < 1:
        raise ValueError('positive threads required')
    required={'query_id','cdr3','v','j','locus',*PMHC}
    if missing:=required-set(queries.columns):
        raise ValueError('missing canonical query fields: '+', '.join(sorted(missing)))
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('source query IDs must be unique and non-null')
    q=queries.with_columns(pl.col('query_id').cast(pl.String))
    identity='pair_id' if paired else 'query_id'
    if identity not in q.columns or q[identity].null_count():
        raise ValueError('paired rows require explicit non-null pair IDs')
    q=q.with_columns(pl.col(identity).cast(pl.String).alias('_case'))
    if q.group_by('_case','locus').len().filter(pl.col('len')>1).height:
        raise ValueError('query case has multiple source rows at one locus')
    if q.group_by('_case').agg(pl.col(PMHC).n_unique()).filter(pl.any_horizontal(pl.col(PMHC)>1)).height:
        raise ValueError('query case has conflicting target metadata')
    cases=q.group_by('_case',maintain_order=True).agg(pl.col(PMHC).first(),pl.len().alias('source_chain_rows'),
        pl.col('locus').n_unique().alias('_lanes'),pl.col('locus').first().alias('locus')).rename({'_case':'query_id'}).with_columns(
        pl.lit(species).alias('species'),pl.when(pl.col('_lanes')>1).then(pl.lit('paired')).otherwise(pl.col('locus')).alias('locus')).drop('_lanes')
    records=[];groups=[]
    for group_id,(task,requested,selected) in enumerate(assayed_targets(cases.select('query_id'),cases,reference,species)):
        groups.append(requested.select('query_id').with_columns(pl.lit(group_id).alias('_group')))
        records.append(selected.with_columns(pl.lit(group_id).alias('_group')))
    assigned=pl.concat(groups)
    membership=pl.concat(records)
    cases=cases.join(assigned,on='query_id',validate='1:1',maintain_order='left')
    q=q.join(assigned,left_on='_case',right_on='query_id',validate='m:1')
    canonical=pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
    params=search_params(scope,engine='seqtm',matrix='')
    hits={};report={'native_batches':0,'native_queries':{},'native_hits':{},'scope':scope}
    result=cases
    for chain,name in [('TRA','alpha'),('TRB','beta')]:
        lane=q.filter(pl.col('locus')==chain).select('_case','_group',pl.col('cdr3').alias('_query'))
        refs=membership.filter(pl.col('gene')==chain)
        nodes=reference.filter(pl.col('gene')==chain).filter(canonical).select('cdr3').unique().sort('cdr3')
        queries_unique=lane.filter(pl.col('_query').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False))['_query'].unique().sort().to_list()
        rows=[]
        if nodes.height and queries_unique:
            index=Index.build(nodes['cdr3'].to_list(),'aa')
            found=index.search_batch(queries_unique,params,threads)
            count=sum(map(len,found))
            if count>5_000_000:
                raise MemoryError('fixed-scope native hit result exceeds5000000 rows; no dense expansion performed')
            seqs=nodes['cdr3'].to_list()
            rows=[(sequence,seqs[h.ref_id],h.n_subs+h.n_ins+h.n_dels) for sequence,group in zip(queries_unique,found) for h in group]
            report['native_batches']+=1
        h=pl.DataFrame(rows,schema={'_query':pl.String,'cdr3':pl.String,'edits':pl.Int64},orient='row')
        report['native_queries'][chain]=len(queries_unique);report['native_hits'][chain]=h.height
        hits[name]=h
        selected=h.filter(pl.col('_query')!=pl.col('cdr3')) if exclude_exact else h
        matched=selected.join(refs.select('_group','cdr3','v','j'),on='cdr3',how='inner').join(lane,on=['_group','_query'],how='inner')
        summary=matched.group_by('_case').agg(pl.col('edits').min().alias(name+'_nearest_edits'),
            pl.col('cdr3').n_unique().alias(name+'_neighbor_junctions'),
            pl.struct('cdr3','v','j').n_unique().alias(name+'_neighbor_clonotypes'),pl.len().alias(name+'_matched_source_records'))
        sizes=refs.group_by('_group').agg(pl.len().alias(name+'_reference_records'),
            pl.col('cdr3').n_unique().alias(name+'_reference_junctions'),pl.struct('cdr3','v','j').n_unique().alias(name+'_reference_clonotypes'))
        result=result.join(summary,left_on='query_id',right_on='_case',how='left',validate='1:1').join(sizes,on='_group',how='left',validate='m:1')
        exact=refs.group_by('_group','cdr3').agg(pl.len().alias('_remove_records'),pl.struct('cdr3','v','j').n_unique().alias('_remove_clonotypes')).rename({'cdr3':'_query'})
        removed=lane.join(exact,on=['_group','_query'],how='left').select(pl.col('_case').alias('query_id'),
            pl.col('_remove_records').fill_null(0).alias(name+'_excluded_reference_records'),
            pl.col('_remove_clonotypes').fill_null(0).alias('_remove_clonotypes'))
        result=result.join(removed,on='query_id',how='left',validate='1:1').with_columns(
            pl.col(name+'_excluded_reference_records').fill_null(0),pl.col('_remove_clonotypes').fill_null(0))
        if not exclude_exact:
            result=result.with_columns(pl.lit(0).alias(name+'_excluded_reference_records'),pl.lit(0).alias('_remove_clonotypes'))
        result=result.with_columns(*[pl.col(name+'_reference_'+unit).alias(name+'_full_reference_'+unit) for unit in ('records','junctions','clonotypes')],
            (pl.col(name+'_reference_records')-pl.col(name+'_excluded_reference_records')).alias(name+'_reference_records'),
            (pl.col(name+'_reference_junctions')-(pl.col(name+'_excluded_reference_records')>0).cast(pl.UInt32)).alias(name+'_reference_junctions'),
            (pl.col(name+'_reference_clonotypes')-pl.col('_remove_clonotypes')).alias(name+'_reference_clonotypes')).drop('_remove_clonotypes')
    linked=build_paired_ref(reference)
    pair_members=[]
    for group_id,part in membership.partition_by('_group',as_dict=True).items():
        pair_members.append(linked.join(part.select('complex_id').unique(),on='complex_id',how='semi').with_columns(pl.lit(group_id[0]).alias('_group')))
    pair_refs=pl.concat(pair_members) if pair_members else linked.head(0).with_columns(pl.lit(0).alias('_group'))
    alpha=q.filter(pl.col('locus')=='TRA').select('_case','_group',pl.col('cdr3').alias('_qa'))
    beta=q.filter(pl.col('locus')=='TRB').select('_case','_group',pl.col('cdr3').alias('_qb'))
    paired_queries=alpha.join(beta,on=['_case','_group'],validate='1:1')
    a=hits['alpha'].rename({'_query':'_qa','cdr3':'alpha','edits':'_da'})
    b=hits['beta'].rename({'_query':'_qb','cdr3':'beta','edits':'_db'})
    matched=a.join(pair_refs,on='alpha').join(paired_queries,on=['_group','_qa']).join(b,on=['_qb','beta'])
    if exclude_exact:
        matched=matched.filter(~((pl.col('_qa')==pl.col('alpha'))&(pl.col('_qb')==pl.col('beta'))))
    pair_key=['alpha','beta',*[c for c in ('v_alpha','j_alpha','v_beta','j_beta') if c in matched.columns]]
    summary=matched.group_by('_case').agg(pl.max_horizontal('_da','_db').min().alias('linked_nearest_bottleneck_edits'),
        pl.struct('alpha','beta').n_unique().alias('linked_neighbor_junction_pairs'),
        pl.struct(pair_key).n_unique().alias('linked_neighbor_clonotypes'),pl.col('complex_id').n_unique().alias('linked_matched_complexes'))
    sizes=pair_refs.group_by('_group').agg(pl.col('complex_id').n_unique().alias('linked_reference_complexes'),
        pl.struct('alpha','beta').n_unique().alias('linked_reference_junction_pairs'))
    result=result.join(summary,left_on='query_id',right_on='_case',how='left',validate='1:1').join(sizes,on='_group',how='left',validate='m:1')
    exact_pairs=pair_refs.group_by('_group','alpha','beta').agg(pl.col('complex_id').n_unique().alias('_removed_complexes'))
    removed=paired_queries.join(exact_pairs,left_on=['_group','_qa','_qb'],right_on=['_group','alpha','beta'],how='left').select(
        pl.col('_case').alias('query_id'),pl.col('_removed_complexes').fill_null(0).alias('linked_excluded_reference_complexes'))
    result=result.join(removed,on='query_id',how='left',validate='1:1').with_columns(pl.col('linked_excluded_reference_complexes').fill_null(0))
    if not exclude_exact:result=result.with_columns(pl.lit(0).alias('linked_excluded_reference_complexes'))
    result=result.with_columns(pl.col('linked_reference_complexes').alias('linked_full_reference_complexes'),
        pl.col('linked_reference_junction_pairs').alias('linked_full_reference_junction_pairs'),
        (pl.col('linked_reference_complexes')-pl.col('linked_excluded_reference_complexes')).alias('linked_reference_complexes'),
        (pl.col('linked_reference_junction_pairs')-(pl.col('linked_excluded_reference_complexes')>0).cast(pl.UInt32)).alias('linked_reference_junction_pairs'))
    counts=[c for c in result.columns if 'neighbor_' in c or '_reference_' in c or c.endswith('_matched_source_records') or c=='linked_matched_complexes']
    result=result.with_columns(pl.col(counts).fill_null(0).cast(pl.UInt64),
        pl.lit(scope).alias('fixed_scope'),pl.lit(exclude_exact).alias('exclude_exact'))
    result=result.with_columns((-pl.col('alpha_nearest_edits')).alias('alpha_nearest_score'),
        (-pl.col('beta_nearest_edits')).alias('beta_nearest_score'),
        (-pl.col('linked_nearest_bottleneck_edits')).alias('linked_nearest_score'),
        pl.min_horizontal('alpha_nearest_edits','beta_nearest_edits').alias('marginal_any_nearest_edits'),
        *[pl.when(pl.col('locus')=='TRA').then(pl.col('alpha_'+source)).when(pl.col('locus')=='TRB').then(pl.col('beta_'+source)).otherwise(None).alias('naive_'+dest) for source,dest in [('nearest_edits','nearest_edits'),('neighbor_junctions','neighbor_junctions'),('neighbor_clonotypes','neighbor_clonotypes')]],
        ((pl.col('alpha_neighbor_junctions')>0)|(pl.col('beta_neighbor_junctions')>0)|(pl.col('linked_neighbor_junction_pairs')>0)).alias('any_hit'))
    result=result.with_columns((-pl.col('naive_nearest_edits')).alias('naive_nearest_score'))
    return result.drop('_group'),report
