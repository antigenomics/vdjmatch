"""Reference graph summaries with exact leave-one-junction-out connectivity.

A single iterative depth-first traversal gives the components remaining after
removing each vertex, including articulation vertices. No per-query graph search.
"""
from __future__ import annotations
import math
import operator


def graph_profiles(n, edges):
    """Return full and vertex-punctured graph counts in O(vertices + edges).

    Vertices are unique full junctions; edges join one-substitution neighbours.
    Removing a vertex removes every associated V annotation from the caller's
    composition table. This function supplies connectivity, not recognition odds.
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
