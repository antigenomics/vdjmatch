"""Reproduce the manuscript's single-chain PSSM density component.

This is a baseline, not the full manuscript scorer: sparse-reference germline priors
and cohort rank fusion are separate calculations. Reference representatives follow
source order within each pMHC/junction, as in the original per-epitope producer.
"""
from __future__ import annotations

from collections import defaultdict, deque, namedtuple
from array import array
from itertools import accumulate
import math
import sys

import polars as pl
from seqtree import Index, SearchParams
from seqtree.evalue import evalue_result

from ..aggregate.candidates import PMHC
from ..evalue import first_hit
from .regions import significance_pssm, gene_family
from .vgene import precompute_v_similarities


_PREPARED_LIMIT_BYTES = 2 * 1024**3
# Union-only transient construction budget; retained prepared state stays at 2 GiB.
_UNION_PEAK_LIMIT_BYTES = 8 * 1024**3
_PreparedDensity = namedtuple('_PreparedDensity',
    'queries references query_indices offsets refs substitutions penalties unit_counts radii')


class PreparedBudget:
    """One invocation's conservative retained graph allocation budget."""
    def __init__(self):
        self.retained = 0
        self.query_axes = {}

    def check(self, allocation, temporary=0, *, peak_limit_bytes=None):
        """Always cap owned state at 2 GiB; default transient-inclusive cap is 2 GiB."""
        owned=self.retained+allocation
        peak=_PREPARED_LIMIT_BYTES if peak_limit_bytes is None else peak_limit_bytes
        if owned > _PREPARED_LIMIT_BYTES:
            raise ValueError('historical prepared neighbours exceed the 2 GiB owned-state allocation budget '
                             f'({peak/1024**3:g} GiB peak ceiling); '
                             'restrict the requested targets or query population')
        if owned+temporary > peak:
            raise ValueError(f'historical prepared neighbours exceed the {peak/1024**3:g} GiB peak '
                             'allocation budget (2 GiB owned-state ceiling); '
                             'restrict the requested targets or query population')

    def query_axis(self, seqs):
        axis = tuple(seqs)
        shared = self.query_axes.get(axis)
        if shared is None:
            size = _query_axis_size(axis)
            self.check(size)
            self.query_axes[axis] = axis
            self.retained += size
            return axis
        return shared


def _query_axis_size(seqs):
    return 128 + sys.getsizeof(seqs) + sum(sys.getsizeof(s) for s in seqs)


def _prepared_size(seqs, rs, edges, *, active_queries=None, include_query_axis=True):
    # One shared raw axis; sparse lookup/CSR/count/radius storage per target.
    active = len(seqs) if active_queries is None else active_queries
    return (1024 + (_query_axis_size(seqs) if include_query_axis else 0) + sys.getsizeof(rs)
            + sum(sys.getsizeof(s) for s in rs)
            + active * (128 + 4 + 8 + 8 + 1) + 8 + edges * 13)


def _pack_prepared(seqs, rs, hits, unit, radii):
    active = [(i,s) for i,s in enumerate(seqs) if hits[s] or unit[i][1]]
    n_edges=sum(len(row) for row in hits.values())
    indices=array('I',(i for i,_ in active))
    offsets=array('Q',[0])*(len(active)+1)
    refs=array('I',[0])*n_edges;subs=array('B',[0])*n_edges;scores=array('q',[0])*n_edges
    counts=array('Q',[0])*len(active);needed=array('B',[0])*len(active)
    k=0
    for j,(i,s) in enumerate(active):
        for h in hits[s]:
            refs[k]=h.ref_id;subs[k]=h.n_subs;scores[k]=h.score;k+=1
        offsets[j+1]=k;counts[j]=unit[i][1];needed[j]=radii[s]
    return _PreparedDensity(seqs,rs,indices,offsets,refs,subs,scores,counts,needed)


def control_histograms(control, sequences, *, threads=1, radii=None):
    """Transient cumulative unit-edit counts; no neighbour objects cross into Python."""
    if not hasattr(control, 'edit_histogram_batch'):
        raise RuntimeError('historical scoring requires seqtree edit_histogram_batch; install the repaired source')
    sequences=list(sequences)
    if radii is None:
        radii={s:5 for s in sequences}
    groups=defaultdict(list)
    for s in sequences:
        radius=radii[s]
        if isinstance(radius,bool) or not isinstance(radius,int) or not 0<=radius<=5:
            raise ValueError('control radius must be an integer from zero to five')
        if radius:groups[radius].append(s)
    out={s:[0]*6 for s in sequences}
    for radius,batch in groups.items():
        counts=control.edit_histogram_batch(batch,first_hit.scope(radius),threads,True)
        for s,row in zip(batch,counts):
            cumulative=list(accumulate(row))
            out[s]=cumulative+[cumulative[-1]]*(5-radius)
    return out


def _density_neighbours(seqs, rs, threads, *, unit_reference_ids=None):
    idx = Index.build(rs, 'aa')
    by_length = defaultdict(list)
    for s in seqs:
        by_length[len(s)].append(s)
    hits = {}
    for length, batch in by_length.items():
        p = SearchParams(max_subs=5, max_ins=0, max_dels=0, max_total_edits=5, engine='seqtm')
        p.pos_matrix = significance_pssm(length)
        hits.update(zip(batch, idx.search_batch(batch, p, threads)))
    if unit_reference_ids is None:
        unit = idx.edit_histogram_batch(seqs, first_hit.scope(1,1,1), threads, True)
    else:
        # Histogram bins cannot be partitioned among targets. Keep radius-one
        # reference identities from the same Index, including indel-only hits.
        for row in idx.search_batch(seqs, first_hit.scope(1,1,1), threads):
            unit_reference_ids.append([h.ref_id for h in row if h.n_subs+h.n_ins+h.n_dels == 1])
        unit = [(0,len(row)) for row in unit_reference_ids]
    radii = {s:max([1 if hs[1] else 0]+[h.n_subs for h in hits[s] if rs[h.ref_id]!=s])
             for s,hs in zip(seqs,unit)}
    return hits, unit, radii


def _union_control_histograms(control, batches, *, threads, prepared_batches, budget):
    """One fixed union search; 2 GiB owned state and 8 GiB transient-inclusive estimate.

    These allocation estimates are not whole-process RSS. The caller must still
    measure total process-tree memory before accepting a representative workload.
    """
    seqs=None
    first_queries=None
    static_batches=type(batches) in (list,tuple)
    references=[]
    for queries, reference in batches:
        if not queries.height or not reference.height:
            continue
        if seqs is None:
            first_queries=queries
            seqs=queries['cdr3'].unique(maintain_order=True).to_list()
        elif not static_batches or queries is not first_queries:
            current=queries['cdr3'].unique(maintain_order=True).to_list()
            if current != seqs:
                raise ValueError('identical_queries requires identical ordered unique junction axes')
        references.append(reference.unique('cdr3',keep='first',maintain_order=True)['cdr3'].to_list())
    if seqs is None:
        return {}
    axis=budget.query_axis(seqs)
    base=sum(_prepared_size(axis,rs,0,active_queries=0,include_query_axis=False) for rs in references)
    budget.check(base,peak_limit_bytes=_UNION_PEAK_LIMIT_BYTES)
    members=pl.concat([pl.DataFrame({'cdr3':rs}).with_row_index('_local_ref').with_columns(
        pl.lit(target,dtype=pl.UInt32).alias('_target')) for target,rs in enumerate(references)])
    union=members.select('cdr3').unique(maintain_order=True).with_row_index('_union_ref')
    members=members.join(union,on='cdr3',validate='m:1').drop('cdr3')
    rs=union['cdr3'].to_list()
    unit_ids=[]
    hits,unit,radii=_density_neighbours(seqs,rs,threads,unit_reference_ids=unit_ids)
    temporary=(sys.getsizeof(hits)+sys.getsizeof(radii)+_query_axis_size(seqs)+_query_axis_size(rs)
        +sum(sys.getsizeof(row)+sum(max(128,sys.getsizeof(h)) for h in row) for row in hits.values())
        +sys.getsizeof(unit_ids)+sum(sys.getsizeof(row)+len(row)*32 for row in unit_ids)
        +sys.getsizeof(unit)+len(unit)*128+members.estimated_size()*4+union.estimated_size()*2)
    budget.check(base,temporary,peak_limit_bytes=_UNION_PEAK_LIMIT_BYTES)
    if prepared_batches is None:
        return control_histograms(control,radii,threads=threads,radii=radii)
    # Generator conversion can retain Python tuples/lists and both new frames
    # while the native Hit lists and unit-reference lists are still live.
    conversion=256*(sum(len(row) for row in hits.values())+sum(map(len,unit_ids)))
    budget.check(base,temporary+conversion,peak_limit_bytes=_UNION_PEAK_LIMIT_BYTES)
    schema={'_query':pl.UInt32,'_union_ref':pl.UInt32,'_order':pl.UInt32,
            '_subs':pl.UInt8,'_penalty':pl.Int64,'_unit':pl.Boolean,'_exact':pl.Boolean}
    edges=pl.DataFrame(((i,h.ref_id,k,h.n_subs,h.score,False,rs[h.ref_id]==s)
        for i,s in enumerate(seqs) for k,h in enumerate(hits[s])),schema=schema,orient='row')
    unit_edges=pl.DataFrame(((i,r,0,1,0,True,False) for i,row in enumerate(unit_ids) for r in row),
        schema=schema,orient='row')
    del hits,unit_ids,unit
    edges=pl.concat([edges,unit_edges])
    del unit_edges
    # Count expansion before allocating the target-membership join or its sort.
    fanout=members.group_by('_union_ref').len()
    expanded=edges.join(fanout,on='_union_ref',validate='m:1').select(
        pl.col('len').cast(pl.UInt64).sum()).item() or 0
    active_bound=min(expanded,len(seqs)*len(references))
    # Includes sparse pack metadata, joined rows, sort/copy scratch and raw axes.
    temporary=(_query_axis_size(seqs)+_query_axis_size(rs)+sys.getsizeof(radii)+len(radii)*64
        +members.estimated_size()*4+edges.estimated_size()*4+expanded*256)
    budget.check(base+expanded*13+active_bound*149,temporary,peak_limit_bytes=_UNION_PEAK_LIMIT_BYTES)
    # All-hit seqtm traverses codec-ordered trie leaves, not score/ref-ID order.
    # Filtering that union order reproduces each subset trie exactly; sorting by
    # local IDs or penalties would change the historical floating reduction.
    joined=edges.join(members,on='_union_ref',validate='m:m').sort('_target','_query','_unit','_order')
    del edges,members,union,fanout
    temporary=joined.estimated_size()*4+_query_axis_size(rs)+sys.getsizeof(radii)+len(radii)*64
    if prepared_batches is not None:
        groups=iter(joined.group_by('_target',maintain_order=True))
        current=next(groups,None)
        part=stats=pssm=None
        def packed(code,values):
            out=array(code)
            out.frombytes(values.to_numpy().tobytes())
            return out
        for target,rs in enumerate(references):
            if current is None or current[0][0] != target:
                size=_prepared_size(axis,rs,0,active_queries=0,include_query_axis=False)
                budget.check(size,temporary,peak_limit_bytes=_UNION_PEAK_LIMIT_BYTES)
                result=_PreparedDensity(axis,rs,array('I'),array('Q',[0]),array('I'),array('B'),
                    array('q'),array('Q'),array('B'))
            else:
                part=current[1]
                stats=part.group_by('_query').agg(
                    pl.col('_unit').sum().cast(pl.UInt64).alias('_counts'),
                    (~pl.col('_unit')).sum().cast(pl.UInt64).alias('_edges'),
                    pl.when(pl.col('_unit')).then(pl.lit(1,dtype=pl.UInt8))
                        .when(~pl.col('_exact')).then(pl.col('_subs')).otherwise(pl.lit(0,dtype=pl.UInt8))
                        .max().alias('_radius')).sort('_query')
                pssm=part.filter(~pl.col('_unit'))
                size=_prepared_size(axis,rs,pssm.height,active_queries=stats.height,include_query_axis=False)
                budget.check(size,temporary+part.estimated_size()*4+stats.estimated_size()*4,
                    peak_limit_bytes=_UNION_PEAK_LIMIT_BYTES)
                offsets=array('Q',[0]);offsets.extend(packed('Q',stats['_edges'].cum_sum()))
                result=_PreparedDensity(axis,rs,packed('I',stats['_query']),offsets,
                    packed('I',pssm['_local_ref']),packed('B',pssm['_subs']),packed('q',pssm['_penalty']),
                    packed('Q',stats['_counts']),packed('B',stats['_radius']))
                current=next(groups,None)
            prepared_batches.append(result)
            budget.retained+=size
        del groups,current,part,stats,pssm
    del joined
    return control_histograms(control,radii,threads=threads,radii=radii)


def shared_control_histograms(control, batches, *, threads=1, prepared_batches=None,
                              prepared_budget=None, identical_queries=False):
    """One transient control CDF per junction at its largest needed target radius.

    Optional prepared batches retain flat CSR integer edges for consumption once
    in this invocation. A shared budget bounds all retained batches to 2 GiB,
    including shared raw axes, sparse query metadata and temporary Python hits.
    The explicit identical-query union path keeps the same 2 GiB owned-state
    ceiling and allows an 8 GiB transient-inclusive allocation estimate.
    Counts below the largest radius are exactly unchanged.
    ``identical_queries`` explicitly validates one ordered unique junction axis
    and searches its union reference once; source query multiplicities remain in
    the downstream reducer. Ordinary per-target preparation remains the default.
    """
    if threads < 1 or not len(control):
        raise ValueError('positive threads and a nonempty control are required')
    radii = {}
    budget = prepared_budget if prepared_budget is not None else PreparedBudget()
    if prepared_budget is None and prepared_batches:
        for p in prepared_batches:
            if budget.query_axis(p.queries) is not p.queries:
                # Independently supplied batches may still own equal raw axes.
                budget.retained += _query_axis_size(p.queries)
            budget.retained += _prepared_size(p.queries,p.references,len(p.refs),
                active_queries=len(p.query_indices),include_query_axis=False)
    if identical_queries:
        return _union_control_histograms(control,batches,threads=threads,
            prepared_batches=prepared_batches,budget=budget)
    for queries, reference in batches:
        if not queries.height or not reference.height:
            continue
        seqs = queries['cdr3'].unique(maintain_order=True).to_list()
        rs = reference.unique('cdr3', keep='first', maintain_order=True)['cdr3'].to_list()
        if prepared_batches is not None:
            axis = budget.query_axis(seqs)
            budget.check(_prepared_size(seqs,rs,0,active_queries=0,include_query_axis=False),
                _query_axis_size(seqs))
        hits, unit, needed = _density_neighbours(seqs, rs, threads)
        if prepared_batches is not None:
            # Preserve the native edge order and scalar reduction, without retaining
            # Python Neighbor objects or a computed prediction table.
            active = sum(bool(hits[s]) or bool(unit[i][1]) for i,s in enumerate(seqs))
            size = _prepared_size(seqs,rs,sum(len(row) for row in hits.values()),
                active_queries=active,include_query_axis=False)
            temporary = (sys.getsizeof(hits) + sys.getsizeof(needed)
                         + sum(sys.getsizeof(row)+sum(max(128,sys.getsizeof(h)) for h in row)
                               for row in hits.values())
                         + sys.getsizeof(unit)+sum(sys.getsizeof(row)+8*len(row) for row in unit)
                         + _query_axis_size(seqs) + active*128)
            budget.check(size, temporary)
            prepared_batches.append(_pack_prepared(axis,rs,hits,unit,needed))
            budget.retained += size
        del hits
        for s, radius in needed.items():
            radii[s] = max(radii.get(s, 0), radius)
        del unit, needed
    return control_histograms(control, radii, threads=threads, radii=radii)


def density(queries, reference, control, *, threads=1, pool_reference=False, control_counts=None,
            puncture_populations=False, v_similarities=None, prepared_neighbours=None,
            export_punctured=False):
    """Return keyed PSSM density and original radius<=1 significance component.

    The original exp(-P/400) kernel is a tempered, self-normalised BLOSUM
    likelihood ratio: P=100*(A(q,q)+A(r,r)-2*A(q,r)), where A is the
    weighted raw similarity sum. Keep its aggregation and count P-value distinct.

    Query IDs and contradictory labels remain distinct. Exact junctions are excluded
    symmetrically. Neither significance nor density is posterior specificity.
    ``pool_reference`` reproduces the original task-level reference union; the caller
    selects compatible MHC records first. Heterogeneous MHC output fields are null.
    Candidate ``puncture_populations`` also removes exact matches from population
    sizes. An empty remaining population yields missing score/significance.
    ``export_punctured`` retains the fixed-population outputs and accumulates a
    second punctured density from the same edges and control CDF.
    """
    if threads < 1 or not len(control):
        raise ValueError("positive threads and a nonempty control are required")
    if export_punctured and puncture_populations:
        raise ValueError('export_punctured requires the fixed-population primary reduction')
    if v_similarities is None:
        v_similarities=precompute_v_similarities(queries['v'],reference['v'])
    seqs = queries['cdr3'].unique(maintain_order=True).to_list()
    control_sizes = dict.fromkeys(seqs, len(control))
    punctured_control_sizes = None
    if puncture_populations or export_punctured:
        # Bin zero retains exact draw multiplicity. Native histograms require
        # radius >=1; only bin zero is used for population puncturing.
        exact = control.edit_histogram_batch(seqs, first_hit.scope(1,0,0), threads, False)
        punctured_control_sizes = {s:len(control)-int(row[0]) for s,row in zip(seqs,exact)}
        if puncture_populations:
            control_sizes = punctured_control_sizes
    rows = []
    prepared_batches = deque()
    preparation_requested = False
    if pool_reference:
        for column in ('epitope', 'species', 'gene'):
            if column in reference.columns and reference[column].n_unique() > 1:
                raise ValueError('pooled reference requires one '+column)
        key = tuple(reference[c][0] if reference.height and reference[c].n_unique() == 1
                    else None for c in PMHC)
        groups = [(key, reference)] if reference.height else []
    else:
        groups = reference.group_by(PMHC, maintain_order=True)
        if control_counts is None and reference.select(PMHC).unique().height > 1:
            control_counts = shared_control_histograms(control,
                ((queries, obs) for _, obs in groups), threads=threads, prepared_batches=prepared_batches)
            preparation_requested = control_counts is not None
    for key, obs in groups:
        # Reproduce the original per-epitope reference representative rule explicitly.
        ref = obs.unique('cdr3', keep='first', maintain_order=True)
        rs, rv = ref['cdr3'].to_list(), ref['v'].to_list()
        rv=[gene_family(v) for v in rv]
        prepared = prepared_neighbours
        if prepared is None and preparation_requested and seqs and rs:
            if not prepared_batches:
                raise ValueError('missing prepared density neighbour population')
            prepared = prepared_batches.popleft()
        if prepared is None:
            hits, unit, radii = _density_neighbours(seqs, rs, threads)
        else:
            prepared_seqs, prepared_rs, query_indices, offsets, refs, subs, penalties, unit_counts, prepared_radii = prepared
            if prepared_seqs != tuple(seqs) or prepared_rs != rs:
                raise ValueError('prepared density neighbour axes differ from query/reference population')
            seq_index = {prepared_seqs[q]:i for i,q in enumerate(query_indices)}
            radii = dict.fromkeys(seqs,0)
            radii.update((prepared_seqs[q],r) for q,r in zip(query_indices,prepared_radii))
            sparse_counts = dict(zip(query_indices,unit_counts))
            unit = ((0,sparse_counts.get(i,0)) for i in range(len(seqs)))
        columnar = (pool_reference and prepared is not None and
                    not puncture_populations and not export_punctured)
        # Original significance used the nearest nonexact target radius, capped at one.
        # Only target-supported radii enter the original formula. Narrower control
        # balls have the identical CDF at those radii; unsupported queries need none.
        counts = control_histograms(control,seqs,threads=threads,radii=radii) if control_counts is None else control_counts
        if columnar:
            import numpy as np
            scores=np.zeros(queries.height,dtype=np.float64)
            probabilities=np.ones(queries.height,dtype=np.float64)
            significance={prepared_seqs[q]: (evalue_result(n,counts[prepared_seqs[q]][1],len(rs),len(control))
                ['p_enrichment'] if n else 1.0) for q,n in zip(query_indices,unit_counts)}
            # Source appearances, including repeated IDs, are indexed by position.
            # Only sparse supported junctions enter the unchanged scalar edge loop.
            selected=queries.select('query_id','cdr3','v').with_row_index('_density_row')
            selected=(selected.join(pl.DataFrame({'cdr3':[prepared_seqs[q] for q in query_indices]}),
                on='cdr3',how='semi',maintain_order='left') if query_indices else selected.head(0))
        else:
            reference_keys = set(rs)
            reference_sizes = {s:len(rs)-int(puncture_populations and s in reference_keys) for s in seqs}
            significance = {s:(evalue_result(hs[1],counts[s][1],reference_sizes[s],control_sizes[s])['p_enrichment']
                            if hs[1] else 1.0) if reference_sizes[s] and control_sizes[s] else None
                            for s,hs in zip(seqs,unit)}
            selected=queries
        for query in (() if columnar and not selected.height else selected.iter_rows(named=True)):
            s, v = query['cdr3'], gene_family(query['v'])
            nr, nc_population = (len(rs),len(control)) if columnar else (reference_sizes[s],control_sizes[s])
            total = 0.0 if nr and nc_population else None
            if export_punctured:
                punctured_nr = len(rs)-int(s in reference_keys)
                punctured_nc = punctured_control_sizes[s]
                punctured_total = 0.0 if punctured_nr and punctured_nc else None
            i = seq_index.get(s) if prepared is not None else None
            edges = (((refs[k],subs[k],penalties[k]) for k in
                      (range(offsets[i],offsets[i+1]) if i is not None else ())) if prepared is not None else
                     ((hit.ref_id, hit.n_subs, hit.score) for hit in hits[s]))
            for ref_id, n_subs, penalty in edges if total is not None else []:
                if rs[ref_id] == s:
                    continue
                ref_v = rv[ref_id]
                if v != ref_v and (v,ref_v) not in v_similarities:
                    raise ValueError('precomputed V similarities do not cover the query/reference families')
                weight = 1.0 if v == ref_v else .25*v_similarities[(v,ref_v)]
                nc = counts[s][n_subs]
                weighted = weight*math.exp(-penalty/400.0)
                total += weighted/max(nr/nc_population*nc,.01)
                if export_punctured and punctured_total is not None:
                    punctured_total += weighted/max(punctured_nr/punctured_nc*nc,.01)
            if columnar:
                scores[query['_density_row']]=total
                probabilities[query['_density_row']]=significance[s]
                continue
            rows.append({'query_id':query['query_id'], **dict(zip(PMHC,key)),
                         'score':total, 'p_enrichment':significance[s],
                         'n_reference':nr,'control_size':nc_population,
                         'estimator':'punctured-pssm-density-v1' if puncture_populations else 'historical-pssm-density-v1',
                         **({'punctured_score':punctured_total,'punctured_n_reference':punctured_nr,
                             'punctured_control_size':punctured_nc} if export_punctured else {})})
        if columnar:
            return queries.select('query_id').with_columns(*[
                pl.lit(value,dtype=pl.String).alias(name) for name,value in zip(PMHC,key)],
                pl.Series('score',scores),pl.Series('p_enrichment',probabilities),
                pl.lit(len(rs),dtype=pl.UInt64).alias('n_reference'),
                pl.lit(len(control),dtype=pl.UInt64).alias('control_size'),
                pl.lit('historical-pssm-density-v1').alias('estimator'))
    schema={'query_id':queries.schema['query_id'], **{c:pl.String for c in PMHC},
            'score':pl.Float64,'p_enrichment':pl.Float64,'n_reference':pl.UInt64,
            'control_size':pl.UInt64,'estimator':pl.String}
    if export_punctured:
        schema.update(punctured_score=pl.Float64,punctured_n_reference=pl.UInt64,punctured_control_size=pl.UInt64)
    return pl.DataFrame(rows,schema=schema)


def legacy_components(queries, reference, control, background, *, threads=1,
                      v_similarities=None):
    """Owning held-out NED, separate V/J/length odds and global first-hit gate.

    Reference junctions take their first union-level source representative. The
    caller supplies the owning raw query projection and sampled prior background.
    Exact junction hits are excluded, while population/model counts remain fixed.
    Empirical negative-tail/Sidak fusion and accuracy summaries are external.
    """
    import numpy as np
    from collections import Counter
    if threads < 1 or not len(control) or not background.height:
        raise ValueError('legacy components require positive threads and nonempty control/background')
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('query identities must be unique and non-null')
    reference=reference.unique('cdr3',keep='first',maintain_order=True)
    if not reference.height:
        raise ValueError('legacy components require a nonempty reference union')
    seqs=queries['cdr3'].unique(maintain_order=True).to_list()
    rs=reference['cdr3'].to_list()
    hits,unit,radii=_density_neighbours(seqs,rs,threads)
    counts=control_histograms(control,seqs,threads=threads,radii=radii)
    rv=[gene_family(v) for v in reference['v']]
    if v_similarities is None:
        v_similarities=precompute_v_similarities(queries['v'],reference['v'])
    epitopes=reference['epitope'].unique(maintain_order=True).to_list()
    ei={e:i for i,e in enumerate(epitopes)}
    labels=reference['epitope'].to_list()
    # Native sparse edges are reused for every raw observation of a junction.
    edges={s:[h for h in hits[s] if rs[h.ref_id]!=s] for s in seqs}
    penalty={s:np.fromiter((math.exp(-h.score/400.)
                          for h in edges[s]),dtype=np.float64) for s in seqs}
    ned=np.zeros((queries.height,len(epitopes)),dtype=np.float64)
    for i,(seq,v) in enumerate(zip(queries['cdr3'],queries['v'])):
        family=gene_family(v)
        weights=[]
        for h in edges[seq]:
            other=rv[h.ref_id]
            if family!=other and (family,other) not in v_similarities:
                raise ValueError('precomputed V similarities do not cover the query/reference families')
            weights.append(1. if family==other else .25*v_similarities[(family,other)])
        # add.at retains the native edge addition order of the owning scalar sum.
        np.add.at(ned[i],[ei[labels[h.ref_id]] for h in edges[seq]],
                  np.asarray(weights)*penalty[seq]/np.asarray([1+counts[seq][h.n_subs] for h in edges[seq]]))
    gates={}
    for seq,curve in zip(seqs,unit):
        n=int(curve[1]);c=int(counts[seq][1])
        result=(evalue_result(n,c,reference.height,len(control)) if n else
                first_hit.pvalue([],[],reference.height,len(control)))
        gates[seq]={'first_hit_radius':1 if n else None,'first_hit_n_target':n,
                    'first_hit_n_control':c if n else 0,'first_hit_E':result['E'],
                    'first_hit_p_enrichment':result['p_enrichment'],
                    'first_hit_gate':result['p_enrichment']<.001}
    qfeatures={'v_logodds':[gene_family(v) for v in queries['v']],
               'j_logodds':[gene_family(j) for j in queries['j']],
               'length_logodds':[len(s) for s in queries['cdr3']]}
    bfeatures={'v_logodds':[gene_family(v) for v in background['v']],
               'j_logodds':[gene_family(j) for j in background['j']],
               'length_logodds':[len(s) for s in background['cdr3']]}
    bc={name:Counter(values) for name,values in bfeatures.items()}
    pieces=[]
    for epitope in epitopes:
        ref=reference.filter(pl.col('epitope')==epitope)
        available=ref.height>=20  # exact owning build_models eligibility
        rfeatures={'v_logodds':[gene_family(v) for v in ref['v']],
                   'j_logodds':[gene_family(j) for j in ref['j']],
                   'length_logodds':[len(s) for s in ref['cdr3']]}
        columns=[]
        for name,values in rfeatures.items():
            rc=Counter(values);bb=bc[name]
            nr,nb=ref.height+len(rc),background.height+len(bb)
            odds={value:(math.log if name=='length_logodds' else np.log)(((rc[value]+1)/nr)/((bb[value]+1)/nb))
                  for value in rc.keys()|bb.keys()}
            columns.append(pl.Series(name,[odds.get(value,0.) if available else None
                for value in qfeatures[name]],dtype=pl.Float64))
        frame=queries.select('query_id').with_columns(pl.lit(epitope,dtype=pl.String).alias('epitope'),
            pl.Series('ned_score',ned[:,ei[epitope]]),*columns,
            pl.lit(ref.height,dtype=pl.UInt64).alias('n_reference'),
            pl.lit(reference.height,dtype=pl.UInt64).alias('reference_union_size'),
            pl.lit(len(control),dtype=pl.UInt64).alias('control_size'),
            pl.lit(background.height,dtype=pl.UInt64).alias('prior_background_size'),
            pl.lit(available).alias('feature_model_available'))
        gate_fields={'first_hit_radius':pl.UInt32,'first_hit_n_target':pl.UInt64,
                     'first_hit_n_control':pl.UInt64,'first_hit_E':pl.Float64,
                     'first_hit_p_enrichment':pl.Float64,'first_hit_gate':pl.Boolean}
        frame=frame.with_columns(*[pl.Series(name,[gates[s][name] for s in queries['cdr3']],dtype=dtype)
                                  for name,dtype in gate_fields.items()])
        pieces.append(frame)
    return pl.concat(pieces).with_columns(pl.lit('legacy-heldout-components-v1').alias('estimator'))


def _prior_background(background):
    weighted = 'event_weight' in background.columns
    if weighted:
        from ..io.columns import _integer_count
        background = background.with_columns(_integer_count('event_weight').alias('event_weight'))
        if background['event_weight'].null_count() or (background['event_weight'] <= 0).any():
            raise ValueError('prior event weights must be positive integers')
    return background, weighted, int(background['event_weight'].sum()) if weighted else background.height


def _prior_features(frame):
    return frame.with_columns(pl.col('v', 'j').str.replace(r'\*.*$', ''),
        pl.col('cdr3').str.len_chars().alias('length'))


def _prepare_prior_background(background):
    """Operation-local sufficient counts for the fixed V/J/length control prior."""
    background, weighted, nb = _prior_background(background)
    if not nb:
        raise ValueError('germline prior requires nonempty reference and background')
    background = _prior_features(background)
    return nb, {column: (background.group_by(column).agg(pl.col('event_weight').sum().alias('_bg_count'))
        if weighted else background.group_by(column).len().rename({'len': '_bg_count'}))
        for column in ('v', 'j', 'length')}


def germline_prior(queries, reference, background, *, include_length=True, exclude_exact=False,
                   _background_counts=None):
    """Original smoothed V/J[/junction-length] log likelihood ratio, keyed by query.

    Frames use ``cdr3,v,j``. Reference junctions use their first source-order V/J;
    background rows retain raw multiplicity. Alpha reproduction deliberately passes
    the original beta background with ``include_length=False``; no correction is
    inferred here. Unknown/missing-input disposition belongs to the caller.
    ``exclude_exact`` removes every full-junction match from both populations and
    feature counts per query. This candidate option preserves historical defaults;
    an emptied population yields missing evidence, not a prior-only assignment.
    A canonical background ``event_weight`` column represents positive integer
    rearrangement multiplicities without materializing repeated rows; ordinary
    read/count columns remain ignored. Reference observations retain their usual
    first-junction unit irrespective of background weights.
    """
    ref = reference.unique('cdr3', keep='first', maintain_order=True)
    if _background_counts is not None and exclude_exact:
        raise ValueError('prepared prior background requires a fixed population')
    if _background_counts is None:
        background, weighted, nb = _prior_background(background)
    else:
        nb, background_counts = _background_counts
    nr = ref.height
    if not nr or not nb:
        raise ValueError('germline prior requires nonempty reference and background')
    columns = ['v', 'j'] + (['length'] if include_length else [])

    def counts(population, keys, name, *, background=False):
        if weighted and background:
            return population.group_by(keys).agg(pl.col('event_weight').sum().alias(name))
        return population.group_by(keys).len().rename({'len': name})

    if not exclude_exact:
        # Preserve the accepted scalar division and inline null-fill operation order.
        out = _prior_features(queries).select('query_id', *columns)
        ref = _prior_features(ref)
        if _background_counts is None:
            bg = _prior_features(background)
        out = out.with_columns(pl.lit(0.0).alias('germline_lr'))
        for column in columns:
            rc = ref.group_by(column).len().rename({'len': '_ref_count'})
            bc = (counts(bg, column, '_bg_count', background=True) if _background_counts is None
                  else background_counts[column])
            out = (out.join(rc, on=column, how='left', maintain_order='left')
                   .join(bc, on=column, how='left', maintain_order='left')
                   .with_columns((pl.col('germline_lr') +
                       (((pl.col('_ref_count').fill_null(0) + .5) / (nr + 30)) /
                        ((pl.col('_bg_count').fill_null(0) + .5) / (nb + 30))).log())
                       .alias('germline_lr'))
                   .drop('_ref_count', '_bg_count'))
        return out.select('query_id', 'germline_lr').with_columns(
            pl.lit(nr, dtype=pl.UInt64).alias('n_reference'),
            pl.lit(nb, dtype=pl.UInt64).alias('background_size'))

    out = _prior_features(queries).select('query_id', *(['cdr3'] if exclude_exact else []), *columns)
    ref, bg = _prior_features(ref), _prior_features(background)
    out = out.with_columns(pl.lit(0.0).alias('germline_lr'),
        pl.lit(nr, dtype=pl.UInt64).alias('n_reference'),
        pl.lit(nb, dtype=pl.UInt64).alias('background_size'))
    if exclude_exact:
        for population, size in [(ref, 'n_reference'), (bg, 'background_size')]:
            removed = counts(population, 'cdr3', '_removed', background=size == 'background_size')
            out = (out.join(removed, on='cdr3', how='left', validate='m:1', maintain_order='left')
                .with_columns((pl.col(size) - pl.col('_removed').fill_null(0)).alias(size))
                .drop('_removed'))
    for column in columns:
        rc = ref.group_by(column).len().rename({'len': '_ref_count'})
        bc = counts(bg, column, '_bg_count', background=True)
        out = (out.join(rc, on=column, how='left', maintain_order='left')
               .join(bc, on=column, how='left', maintain_order='left')
               .with_columns(pl.col('_ref_count', '_bg_count').fill_null(0)))
        if exclude_exact:
            keys = ['cdr3', column]
            for population, count in [(ref, '_ref_count'), (bg, '_bg_count')]:
                removed = counts(population, keys, '_removed', background=count == '_bg_count')
                out = (out.join(removed, on=keys, how='left', validate='m:1', maintain_order='left')
                    .with_columns((pl.col(count) - pl.col('_removed').fill_null(0)).alias(count))
                    .drop('_removed'))
        out = (out
               .with_columns((pl.col('germline_lr') +
                   (((pl.col('_ref_count') + .5) / (pl.col('n_reference') + 30)) /
                    ((pl.col('_bg_count') + .5) / (pl.col('background_size') + 30))).log())
                   .alias('germline_lr'))
               .drop('_ref_count', '_bg_count'))
    if exclude_exact:
        out = out.with_columns(pl.when((pl.col('n_reference') > 0) & (pl.col('background_size') > 0))
            .then(pl.col('germline_lr')).otherwise(None).alias('germline_lr'))
    return out.select('query_id', 'germline_lr', 'n_reference', 'background_size')


def _score_series(values):
    values = pl.Series(values, dtype=pl.Float64)
    if values.null_count() or not values.is_finite().all():
        raise ValueError('original fusion requires finite, nonmissing component scores')
    return values


def unpaired_score(densities, priors, reference_size, *, log_scale=False):
    """Original score; log_scale puts dense targets on the same log-density units.

    Germline terms retain the original sparse-reference activation. Historical
    defaults and sparse scores are unchanged; dense ordering is unchanged.
    """
    density_values, prior_values = _score_series(densities), _score_series(priors)
    if isinstance(reference_size, pl.Series):
        sizes = _score_series(reference_size)
        if len(sizes) != len(density_values) or (sizes < 1).any() or (sizes != sizes.floor()).any():
            raise ValueError('aligned components and a positive integer reference size are required')
        if len(density_values) != len(prior_values) or (density_values < 0).any():
            raise ValueError('aligned components and nonnegative density are required')
        return pl.DataFrame({'density': density_values, 'prior': prior_values, 'size': sizes}).select(
            pl.when(pl.col('size') < 500).then((pl.col('density') + 1e-6).log() + pl.col('prior'))
            .otherwise((pl.col('density')+1e-6).log() if log_scale else pl.col('density'))
            .alias('score'))['score']
    if len(density_values) != len(prior_values) or reference_size < 1:
        raise ValueError('aligned components and a positive reference size are required')
    if (density_values < 0).any():
        raise ValueError('density must be nonnegative')
    score = (density_values+1e-6).log() if log_scale or reference_size<500 else density_values
    if reference_size<500:
        score=score+prior_values
    return score.rename('score')


def paired_score(beta_densities, alpha_densities, combined_priors, *, log_scale=False):
    """Original equal-weight cohort rank sum; ties receive mean zero-based rank.

    Inputs are aligned observations, retaining repeated junctions and pair context.
    ``combined_priors`` is beta V/J/length plus alpha V/J. Optional log_scale
    adds original chain log-densities and this same prior, without cohort ranks.
    It is factorised original evidence, not calibrated joint binding probability.
    """
    beta, alpha, prior = map(_score_series, (beta_densities, alpha_densities, combined_priors))
    if len(beta) != len(alpha) or len(beta) != len(prior):
        raise ValueError('paired components must describe the same aligned observations')
    if (beta < 0).any() or (alpha < 0).any():
        raise ValueError('density must be nonnegative')
    if log_scale:
        return ((beta+1e-6).log()+(alpha+1e-6).log()+prior).rename('score')
    return (beta.rank('average') + alpha.rank('average') + prior.rank('average') - 3).rename('score')
