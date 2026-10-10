"""Reproduce the manuscript's single-chain PSSM density component.

This is a baseline, not the full manuscript scorer: sparse-reference germline priors
and cohort rank fusion are separate calculations. Reference representatives follow
source order within each pMHC/junction, as in the original per-epitope producer.
"""
from __future__ import annotations

import bisect
from collections import defaultdict
import math

import polars as pl
from seqtree import Index, SearchParams

from ..aggregate.candidates import PMHC
from ..evalue import first_hit
from .regions import significance_pssm, gene_family
from .vgene import vsim


def density(queries, reference, control, *, threads=1):
    """Return keyed PSSM density and original radius<=1 significance component.

    Query IDs and contradictory labels remain distinct. Exact junctions are excluded
    symmetrically. Neither significance nor density is posterior specificity.
    """
    if threads < 1 or not len(control):
        raise ValueError("positive threads and a nonempty control are required")
    seqs = queries['cdr3'].unique(maintain_order=True).to_list()
    costs = first_hit._cost_lists(control, seqs, first_hit.scope(), threads, True,
                                  max(1, len(seqs)), 'control', False)
    counts = {s: [d for d, _ in row] for s, row in zip(seqs, costs)}
    by_length = defaultdict(list)
    for s in seqs:
        by_length[len(s)].append(s)
    rows = []
    for key, obs in reference.group_by(PMHC, maintain_order=True):
        # Reproduce the original per-epitope reference representative rule explicitly.
        ref = obs.unique('cdr3', keep='first', maintain_order=True)
        rs, rv = ref['cdr3'].to_list(), ref['v'].to_list()
        idx = Index.build(rs, 'aa')
        hits = {}
        for length, batch in by_length.items():
            p = SearchParams(max_subs=5, max_ins=0, max_dels=0, max_total_edits=5,
                             engine='seqtm')
            p.pos_matrix = significance_pssm(length)
            hits.update(zip(batch, idx.search_batch(batch, p, threads)))
        # Original significance used the nearest nonexact target radius, capped at one.
        unit = first_hit._cost_lists(idx, seqs, first_hit.scope(1,1,1), threads, True,
                                     max(1,len(seqs)), 'target', False)
        significance = {s: first_hit.pvalue([(d,'target') for d,_ in hs], counts[s],len(rs),len(control))
                        for s,hs in zip(seqs,unit)}
        for query in queries.iter_rows(named=True):
            s, v = query['cdr3'], query['v']
            total = 0.0
            for hit in hits[s]:
                if rs[hit.ref_id] == s:
                    continue
                ref_v = rv[hit.ref_id]
                weight = 1.0 if gene_family(v) == gene_family(ref_v) else .25*vsim(v,ref_v)
                nc = bisect.bisect_right(counts[s],hit.n_subs)
                total += weight*math.exp(-hit.score/400.0)/max(len(rs)/len(control)*nc,.01)
            rows.append({'query_id':query['query_id'], **dict(zip(PMHC,key)),
                         'score':total, 'p_enrichment':significance[s]['p_enrichment'],
                         'n_reference':len(rs),'control_size':len(control),
                         'estimator':'historical-pssm-density-v1'})
    schema={'query_id':queries.schema['query_id'], **{c:pl.String for c in PMHC},
            'score':pl.Float64,'p_enrichment':pl.Float64,'n_reference':pl.UInt64,
            'control_size':pl.UInt64,'estimator':pl.String}
    return pl.DataFrame(rows,schema=schema)
