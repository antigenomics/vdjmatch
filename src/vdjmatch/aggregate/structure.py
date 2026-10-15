"""Reference graph summaries with exact leave-one-junction-out connectivity.

A single iterative depth-first traversal gives the components remaining after
removing each vertex, including articulation vertices. No per-query graph search.
"""
from __future__ import annotations
import math
import operator


def graph_profiles(n, edges, *, include_punctures=True):
    """Return full and vertex-punctured graph counts in O(vertices + edges).

    Vertices are unique full junctions; edges join one-substitution neighbours.
    Removing a vertex removes every associated V annotation from the caller's
    composition table. This function supplies connectivity, not recognition odds.
    Normal scoring sets include_punctures=False and returns no deletion diagnostics.
    """
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError('nonnegative integer vertex count required')
    unique = set()
    for a, b in edges:
        try:
            a, b = operator.index(a), operator.index(b)
        except TypeError as exc:
            raise ValueError('integer vertex IDs required') from exc
        if not (0 <= a < n and 0 <= b < n) or a == b:
            raise ValueError('edges must join distinct valid vertex IDs')
        unique.add((min(a,b), max(a,b)))
    neighbours = [[] for _ in range(n)]
    for a,b in unique:
        neighbours[a].append(b); neighbours[b].append(a)
    seen = [-1]*n; low = [0]*n; parent = [-1]*n; sizes = [1]*n
    component = [-1]*n; separated = [[] for _ in range(n)]; component_sizes = []
    clock = 0
    for root in range(n):
        if seen[root] >= 0:
            continue
        cid = len(component_sizes)
        seen[root] = low[root] = clock; clock += 1; component[root] = cid
        stack = [(root, iter(neighbours[root]))]
        while stack:
            u, iterator = stack[-1]
            v = next(iterator, None)
            if v is not None:
                if v == parent[u]:
                    continue
                if seen[v] < 0:
                    parent[v] = u; component[v] = cid
                    seen[v] = low[v] = clock; clock += 1
                    stack.append((v, iter(neighbours[v])))
                else:
                    low[u] = min(low[u], seen[v])
            else:
                stack.pop()
                p = parent[u]
                if p >= 0:
                    sizes[p] += sizes[u]; low[p] = min(low[p], low[u])
                    if low[u] >= seen[p]:
                        separated[p].append(sizes[u])
        component_sizes.append(sizes[root])
    largest = [0,0]
    for size in component_sizes:
        if size >= largest[0]:
            largest = [size,largest[0]]
        else:
            largest[1] = max(largest[1],size)
    size_log = lambda s: s*math.log(s) if s else 0.
    total_log = math.fsum(size_log(s) for s in component_sizes)
    total_square = sum(s*s for s in component_sizes)
    singletons = sum(s == 1 for s in component_sizes)

    def summary(nodes, count, isolated, largest_size, edges_count, log_sum, square_sum):
        return dict(junctions=nodes, components=count, isolated_junctions=isolated,
            graph_edges=edges_count, largest_component=largest_size,
            largest_component_fraction=largest_size/nodes if nodes else None,
            clustered_fraction=1-isolated/nodes if nodes else None,
            component_entropy=max(0.,math.log(nodes)-log_sum/nodes) if nodes else None,
            size_weighted_component_mean=square_sum/nodes if nodes else None)

    full = summary(n,len(component_sizes),singletons,largest[0] if largest else 0,
                   len(unique),total_log,total_square)
    if not include_punctures:
        return full, []
    punctured = []
    for u in range(n):
        old = component_sizes[component[u]]
        fragments = separated[u]
        remainder = old-1-sum(fragments)
        if remainder:
            fragments = [*fragments, remainder]
        unaffected = largest[0] if old != largest[0] else largest[1] if len(largest)>1 else 0
        punctured.append(summary(n-1,len(component_sizes)-1+len(fragments),
            singletons-(old == 1)+sum(s == 1 for s in fragments),
            max(unaffected,max(fragments,default=0)),len(unique)-len(neighbours[u]),
            total_log-size_log(old)+math.fsum(size_log(s) for s in fragments),
            total_square-old*old+sum(s*s for s in fragments)))
    return full, punctured


def composition(records, junctions, *, include_punctures=True):
    """Full and exact-junction-punctured length/V summaries; no query labels."""
    import math
    import numpy as np
    import polars as pl

    seqs = junctions['cdr3'].to_list()
    lengths = np.array([len(s) for s in seqs], dtype=np.float64)
    keys = records.select('cdr3','_allele').unique().with_columns(
        pl.col('_allele').str.split('*').list.first().alias('_gene'))
    counts = keys.group_by('_gene').len().sort('_gene')
    sizes = counts['len'].to_numpy().astype(float)
    total = keys.height
    flog = lambda x: x*math.log(x) if x else 0.
    v_log = math.fsum(flog(x) for x in sizes)
    length_counts = dict(zip(*np.unique(lengths, return_counts=True)))
    length_log = math.fsum(flog(x) for x in length_counts.values())
    sum_length = lengths.sum(); sum_square = (lengths*lengths).sum()
    records_per_seq = dict(records.group_by('cdr3').len().iter_rows())
    removed_keys = dict(keys.group_by('cdr3').len().iter_rows())
    # Each sequence removes its full contribution from every V category.
    v_changes = {}
    if include_punctures:
        changes = keys.group_by('cdr3','_gene').len().rename({'len':'_removed'}).join(
            counts,on='_gene',validate='m:1').sort('cdr3','_gene')
        for sequence, part in changes.group_by('cdr3',maintain_order=True):
            before = part['len'].to_list(); after = (part['len']-part['_removed']).to_list()
            v_changes[sequence[0]] = (math.fsum(flog(new)-flog(old) for old,new in zip(before,after)),
                                      sum(x==0 for x in after))

    def row(removed=None):
        n = len(seqs)-(removed is not None)
        length = len(removed) if removed is not None else 0
        k = removed_keys[removed] if removed is not None else 0
        remain = total-k
        delta, lost_genes = v_changes[removed] if removed is not None else (0.,0)
        v_entropy = max(0., math.log(remain)-(v_log+delta)/remain) if remain else None
        f_length = length_log
        if removed is not None:
            count = length_counts[float(length)]
            f_length += flog(count-1)-flog(count)
        mean = (sum_length-length)/n if n else None
        variance = max(0.,(sum_square-length*length)/n-mean*mean) if n else None
        return dict(clonotype_keys=remain, v_genes=len(sizes)-lost_genes,
            v_entropy=v_entropy, effective_v_genes=math.exp(v_entropy) if v_entropy is not None else None,
            length_mean=mean, length_sd=math.sqrt(variance) if variance is not None else None,
            length_entropy=max(0.,math.log(n)-f_length/n) if n else None,
            reference_records=records.height-(records_per_seq[removed] if removed is not None else 0),
            excluded_reference_records=records_per_seq[removed] if removed is not None else 0,
            excluded_clonotype_keys=k)
    return row(), [row(sequence) for sequence in seqs] if include_punctures else []
