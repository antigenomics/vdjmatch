"""Occurrence-weighted reference communities and vdjtools Chao/iNEXT adapters."""
from itertools import islice

import numpy as np
import polars as pl


def abundance_coverage(counts, *, nboot=50, seed=0, threads=1):
    """Chao1 richness, unseen mass, entropy and effort/coverage curves.

    Counts are database occurrences, not read/cell abundance. The iNEXT
    finite-sample Chao1 and zero-doubleton handling, augmented-assemblage
    bootstrap and coverage function belong to vdjtools. Only the coverage
    percentile envelope is assembled here from that owner's bootstrap law.
    """
    from vdjtools.stats import asymptotic_diversity, rarefaction_batch, sample_coverage
    from vdjtools.stats import bootstrap_abundances

    x = np.asarray(counts)
    if (any(isinstance(v, bool) or not isinstance(v, int) for v in (nboot, seed, threads)) or
            x.ndim != 1 or not x.size or x.dtype.kind not in 'iu' or
            np.any(x <= 0) or nboot < 2 or seed < 0 or threads < 1):
        raise ValueError('positive integer occurrence counts, nboot>=2, seed>=0 and threads>=1 required')
    total = sum(map(int, x))
    if total >= 2**31:
        raise ValueError('occurrence total exceeds vdjtools native bootstrap int32 range')
    estimate = float(asymptotic_diversity(x, q=0, se=False)['estimator'].item())
    working_bytes = (len(x) + total) * 32 + nboot * 41 * 16
    # Include retained draws, float coercion, Python lists, native copies and sample scratch.
    batch_size = min(nboot, (256 * 1024**2 - working_bytes) // (96 * total))
    if batch_size < 1:
        raise ValueError('abundance bootstrap working memory exceeds 256 MiB')
    curve = rarefaction_batch([x], q=0, nboot=nboot, seed=seed, threads=threads).drop('sample')
    depths = curve['m'].to_list()
    coverage_draws = np.empty((nboot, len(depths)))
    draws = bootstrap_abundances(x, nboot=nboot, seed=seed)
    for start in range(0, nboot, batch_size):
        end = min(start + batch_size, nboot)
        block = rarefaction_batch(list(islice(draws, end - start)), q=0, sizes=depths,
            se=False, threads=threads).sort('sample', 'm')
        coverage_draws[start:end] = block['sample_coverage'].to_numpy().reshape(end - start, len(depths))
    curve = curve.with_columns(pl.Series('coverage_lo', np.quantile(coverage_draws, .025, axis=0)),
        pl.Series('coverage_hi', np.quantile(coverage_draws, .975, axis=0)))
    coverage = float(sample_coverage(x))
    p = x.astype(float) / total
    summary = dict(observations=total, observed_taxa=len(x), singletons=int(np.count_nonzero(x == 1)),
        doubletons=int(np.count_nonzero(x == 2)), estimated_taxa=estimate, coverage=coverage,
        unseen_mass=1-coverage, unseen_type_fraction=(estimate-len(x))/estimate,
        entropy_nats=float(-np.sum(p*np.log(p))),
        coverage_status='zero_doubleton_formula' if not np.any(x == 2) and np.any(x == 1) else 'available')
    return summary, curve


def _sw_geometry(search, *, species, locus, cdr_loops, gap_open, gap_extend, threads):
    """Prepare the existing full-loop geometry; unavailable keys stay explicit."""
    from ..match.geometry import load_v_loops, resolve_v_alleles, require_geometry_capabilities
    from seqtree import SubstitutionMatrix
    from seqtree.pairwise import dist_matrix
    require_geometry_capabilities('sw')
    matrix = SubstitutionMatrix.blosum62()
    if not cdr_loops:
        return search.with_columns(pl.lit(0).alias('_vgroup')), [[0]], matrix
    model = load_v_loops(species)
    alleles = resolve_v_alleles(search['v'], model, locus=locus)
    usable = {a for a in alleles.drop_nulls().unique() if all(
        model[a][loop].replace('.', '') and
        set(model[a][loop].replace('.', '')) <= set('ACDEFGHIKLMNPQRSTVWY') for loop in cdr_loops)}
    unique = sorted(usable)
    if len(unique)**2 * 32 > 64 * 1024**2:
        raise ValueError('germline loop distance tables exceed 64 MiB')
    table = np.zeros((len(unique), len(unique)), dtype=np.int64)
    for loop in cdr_loops:
        loops = [model[a][loop].replace('.', '') for a in unique]
        if loops:
            table += np.asarray(dist_matrix(loops, loops, matrix, mode='local',
                gap_open=gap_open, gap_extend=gap_extend, threads=threads))
    groups = alleles.replace_strict(dict(zip(unique, range(len(unique)))), default=None).cast(pl.Int64)
    return search.with_columns(pl.Series('_vgroup', groups, dtype=pl.Int64)), table.tolist(), matrix


def _target_csr(mapping, n):
    """Fixed key-to-population memberships, without a dense incidence matrix."""
    ordered = mapping.sort('_search', 'group_id')
    offsets = np.r_[0, np.cumsum(np.bincount(ordered['_search'].to_numpy(), minlength=n))]
    return offsets.tolist(), ordered['group_id'].to_list()


def sw_neighbor_knees(records, *, species, locus, cdr_loops, gap_open=11, gap_extend=1, threads=1,
            pair_budget=2_000_000, neighbor_orders=(1,)):
    """Kth-other-clonotype SW distances and pooled knees from one count batch.

    Memberships and full clonotype keys stay fixed; only the query's own search
    index is omitted from its own target's sparse native hits. Ties at
    zero distance remain other clonotypes. Unmodeled/undersized groups stay null.
    Knees use the existing unsmoothed convex Kneedle S=1 rule, separately per k.
    """
    orders = tuple(neighbor_orders)
    if not orders or any(isinstance(k, bool) or not isinstance(k, (int, np.integer)) or k < 1 for k in orders):
        raise ValueError('neighbor_orders must contain positive integers')
    orders = tuple(sorted(set(map(int, orders))))
    from seqtree import gapblock
    from ..match.geometry import kneedle_radius
    nodes = records.select('group_id', 'cdr3', 'v', 'j').unique().sort(
        'group_id', 'cdr3', 'v', 'j', nulls_last=True).with_row_index('node')
    search = nodes.select('cdr3', 'v', 'j').unique().sort('cdr3', 'v', 'j', nulls_last=True)
    search, distances, matrix = _sw_geometry(search, species=species, locus=locus,
        cdr_loops=cdr_loops, gap_open=gap_open, gap_extend=gap_extend, threads=threads)
    search = search.filter(pl.col('_vgroup').is_not_null()).with_row_index('_search')
    modeled = nodes.join(search, on=['cdr3', 'v', 'j'], nulls_equal=True, validate='m:1').sort('node')
    pairs = modeled.height * search.height
    if pairs > pair_budget:
        raise ValueError(f'knee pair budget exceeded: {pairs} > {pair_budget}')
    nearest = {k: np.full(modeled.height, -1, dtype=np.int64) for k in orders}
    if modeled.height:
        groups = int(nodes['group_id'].max())+1
        offsets, targets = _target_csr(modeled, search.height)
        # BLOSUM62 diagonal <=11, local similarity >=0: this encloses every score.
        bound = 22 * int(search['cdr3'].str.len_chars().max()) + int(np.max(distances))
        largest = int(modeled.group_by('group_id').len()['len'].max())
        # Bound native hit output and chunk-local ragged rank work before allocation.
        chunk_size = min(16, (64*1024**2)//(12*largest))
        if chunk_size < 1:
            raise ValueError('one target exceeds bounded sparse nearest-distance working memory')
        self_refs = modeled['_search'].to_numpy()
        for chunk in gapblock.grouped_count_batch(modeled['cdr3'].to_list(), search['cdr3'].to_list(),
                [bound], reference_target_offsets=offsets, reference_target_ids=targets,
                target_peptide_ids=[0]*groups, reference_exclusive_peptide_ids=[-1]*search.height,
                query_gene_ids=[0]*modeled.height, reference_gene_ids=[0]*search.height,
                query_group_ids=modeled['_vgroup'].to_list(), reference_group_ids=search['_vgroup'].to_list(),
                group_distances=distances, matrix=matrix, alignment_mode='local', gap_open=gap_open,
                gap_extend=gap_extend, threads=threads, chunk_size=chunk_size, return_hits=True,
                query_target_ids=modeled['group_id'].to_list(), max_hit_bytes=64*1024**2):
            hits = np.asarray(chunk.hits)
            other = hits[hits[:, 1] != self_refs[hits[:, 0]]]
            rows = other[:, 0].astype(np.int64)-chunk.start
            counts = np.bincount(rows, minlength=chunk.stop-chunk.start)
            width = int(counts.max(initial=0))
            if width:
                starts = np.cumsum(counts)-counts
                columns = np.arange(len(other))-np.repeat(starts, counts)
                ranked = np.full((len(counts), width), bound+1, dtype=np.int64)
                ranked[rows, columns] = other[:, 2]
                valid = [k-1 for k in orders if k <= width]
                if valid:
                    ranked.partition(valid, axis=1)
                    for k in orders:
                        if k <= width:
                            nearest[k][chunk.start:chunk.stop] = np.where(counts >= k, ranked[:, k-1], -1)
            del hits, other, chunk
    curves, summaries = [], []
    for k in orders:
        values = modeled.select('node').with_columns(pl.Series('neighbor_distance', nearest[k]),
            pl.lit(True).alias('geometry_available'))
        curve = nodes.join(values, on='node', how='left', validate='1:1').with_columns(
            pl.col('geometry_available').fill_null(False),
            pl.when(pl.col('neighbor_distance') >= 0).then(pl.col('neighbor_distance'))
                .otherwise(None).alias('neighbor_distance'),
            pl.lit(k, dtype=pl.Int64).alias('neighbor_order'))
        curve = curve.with_columns(
            pl.when(~pl.col('geometry_available')).then(pl.lit('geometry_unavailable'))
            .when(pl.col('neighbor_distance').is_null()).then(pl.lit('insufficient_neighbors'))
            .otherwise(pl.lit('available')).alias('distance_status'))
        distances = curve['neighbor_distance'].drop_nulls().to_numpy()
        radius = kneedle_radius(distances)
        summaries.append(dict(neighbor_order=k, epsilon_knee=None if radius is None else int(radius),
            distances_count=len(distances), reference_clonotypes=nodes.height,
            modeled_clonotypes=modeled.height, unmodeled_clonotypes=nodes.height-modeled.height,
            insufficient_neighbor_clonotypes=modeled.height-len(distances),
            knee_status=('available' if radius is not None else 'no_modeled_vertices' if not modeled.height
                else 'no_neighbors' if not len(distances) else 'no_knee')))
        curves.append(curve)
    summary = pl.DataFrame(summaries, schema_overrides={'epsilon_knee': pl.Int64})
    return pl.concat(curves).sort('node', 'neighbor_order'), summary


def sw_knee(records, *, species, locus, cdr_loops, gap_open=11, gap_extend=1, threads=1,
            pair_budget=2_000_000):
    """Existing k=1 pooled radius and node-order curve; unsupported knees raise."""
    curves, summary = sw_neighbor_knees(records, species=species, locus=locus,
        cdr_loops=cdr_loops, gap_open=gap_open, gap_extend=gap_extend, threads=threads,
        pair_budget=pair_budget)
    if not summary['modeled_clonotypes'].item():
        raise ValueError('no modeled SW vertices for radius selection')
    radius = summary['epsilon_knee'].item()
    if radius is None:
        raise ValueError('reference nearest-distance curve has no confirmed Kneedle boundary; supply --sw-radius')
    curve = curves.select('node', 'group_id', 'cdr3', 'v', 'j',
        pl.col('neighbor_distance').alias('nearest_distance'))
    return int(radius), curve


def reference_graph(records, *, geometry, radius=None, species='human', locus='TRB',
                    threads=1, scope='1,1,1,1', gap_open=11, gap_extend=1,
                    cdr_loops=('cdr1', 'cdr2', 'cdr25'), pair_budget=2_000_000, edge_sink=None):
    """One native graph batch; full pMHC vertices retain source occurrences.

    An optional sink consumes bounded edge frames with their population IDs,
    allowing CLI callers to partition on disk instead of retaining every graph.
    The returned edge table is empty when a sink owns the edges.
    """
    from ..cluster import overlap
    from ..match.scope import parse_scope
    from seqtree import gapblock
    if geometry not in ('levenshtein', 'sw') or threads < 1 or pair_budget < 1:
        raise ValueError('known geometry, positive threads and pair budget required')
    if records.filter(~pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)).height:
        raise ValueError('graph requires valid full junctions')
    nodes = records.group_by('group_id', 'cdr3', 'v', 'j').len(name='occurrences').sort(
        'group_id', 'cdr3', 'v', 'j', nulls_last=True).with_row_index('node')
    sizes = nodes.group_by('group_id').len()['len']
    pairs = sum(int(n)*(int(n)-1)//2 for n in sizes)
    search = nodes.select('cdr3', 'v', 'j').unique().sort('cdr3', 'v', 'j', nulls_last=True)
    if geometry == 'levenshtein':
        unique_junctions = nodes['cdr3'].n_unique()
        pairs = unique_junctions*(unique_junctions-1)//2
    if pairs > pair_budget:
        raise ValueError(f'graph pair budget exceeded: {pairs} candidate pairs > {pair_budget}; raise --pair-budget explicitly')
    edge_schema = {'a': pl.Int64, 'b': pl.Int64, 'distance_ceiling': pl.Int64}
    if geometry == 'sw':
        if (not isinstance(radius, int) or isinstance(radius, bool) or not 0 <= radius < 2**28 or
                gap_open < 0 or gap_extend < 0 or not set(cdr_loops) <= {'cdr1', 'cdr2', 'cdr25'}):
            raise ValueError('SW requires an explicit nonnegative integer radius below 2**28, nonnegative gaps and known CDR loops')
        search, distances, matrix = _sw_geometry(search, species=species, locus=locus,
            cdr_loops=cdr_loops, gap_open=gap_open, gap_extend=gap_extend, threads=threads)
        available = search.select('cdr3', 'v', 'j', pl.col('_vgroup').is_not_null().alias('geometry_available'))
        nodes = nodes.join(available, on=['cdr3', 'v', 'j'], nulls_equal=True, validate='m:1').sort('node')
        search = search.filter(pl.col('_vgroup').is_not_null()).with_row_index('_search')
        mapping = nodes.join(search, on=['cdr3', 'v', 'j'], nulls_equal=True, validate='m:1').sort('node')
        if not mapping.height:
            return nodes, pl.DataFrame(schema=edge_schema)
        offsets, targets = _target_csr(mapping, search.height)
        groups = int(nodes['group_id'].max())+1
        edge_parts, edge_count = [], 0
        node_ids = mapping['node'].to_numpy().astype(np.int64)
        group_ids = mapping['group_id'].to_numpy().astype(np.int64)
        search_indices = mapping['_search'].to_numpy()
        group_offsets = np.r_[0, np.cumsum(np.bincount(group_ids, minlength=groups))]
        reference_nodes = np.empty(search.height, dtype=np.int64)
        last_group = None
        for chunk in gapblock.grouped_count_batch(mapping['cdr3'].to_list(), search['cdr3'].to_list(), [radius],
                reference_target_offsets=offsets, reference_target_ids=targets, target_peptide_ids=[0]*groups,
                reference_exclusive_peptide_ids=[-1]*search.height, query_gene_ids=[0]*mapping.height,
                reference_gene_ids=[0]*search.height, chunk_size=128, matrix=matrix,
                gap_open=gap_open, gap_extend=gap_extend, alignment_mode='local',
                query_group_ids=mapping['_vgroup'].to_list(), reference_group_ids=search['_vgroup'].to_list(),
                group_distances=distances, threads=threads, return_hits=True,
                query_target_ids=mapping['group_id'].to_list(), query_reference_ids=mapping['_search'].to_list()):
            hits = np.asarray(chunk.hits)
            edge_count += len(hits)
            if edge_sink is None and edge_count * 64 > 1024**3:
                raise ValueError('sparse population edges exceed 1 GiB budget')
            if len(hits):
                query = hits[:, 0]
                hit_groups = group_ids[query]
                boundaries = np.r_[0, np.flatnonzero(hit_groups[1:] != hit_groups[:-1])+1, len(hits)]
                matched_nodes = np.empty(len(hits), dtype=np.int64)
                for lo, hi in zip(boundaries[:-1], boundaries[1:]):
                    group = hit_groups[lo]
                    if group != last_group:
                        begin, end = group_offsets[group:group+2]
                        reference_nodes[search_indices[begin:end]] = node_ids[begin:end]
                        last_group = group
                    matched_nodes[lo:hi] = reference_nodes[hits[lo:hi, 1]]
                part = pl.DataFrame({'a': node_ids[query], 'b': matched_nodes,
                    'distance_ceiling': np.full(len(hits), radius, dtype=np.int64),
                    'group_id': group_ids[query]})
                if edge_sink is None:
                    edge_parts.append(part.drop('group_id'))
                else:
                    edge_sink(part.select(pl.min_horizontal('a', 'b').alias('a'),
                        pl.max_horizontal('a', 'b').alias('b'), 'distance_ceiling', 'group_id'))
        edges = pl.concat(edge_parts) if edge_parts else pl.DataFrame(schema=edge_schema)
        # Query order is node order; the native upper triangle uses search order.
        return nodes, edges.select(pl.min_horizontal('a', 'b').alias('a'),
            pl.max_horizontal('a', 'b').alias('b'), 'distance_ceiling').sort('a', 'b')
    _, insertions, deletions, radius = parse_scope(scope)
    if insertions != deletions:
        raise ValueError('undirected edit graphs require symmetric insertion/deletion budgets')
    nodes = nodes.with_columns(pl.lit(True).alias('geometry_available'))
    # Edit geometry does not depend on V/J: search each junction once, then
    # expand its fixed full-clonotype memberships in bounded table chunks.
    search = nodes.select('cdr3').unique().sort('cdr3').with_row_index('_search')
    mapping = nodes.join(search, on='cdr3', validate='m:1')
    hits = overlap(search['cdr3'].to_list(), scope=scope, threads=threads) if search.height else None
    if hits is None:
        return nodes, pl.DataFrame(schema=edge_schema)
    a, b = hits['a_idx'].to_numpy(), hits['b_idx'].to_numpy()
    keep = a < b
    global_edges = pl.concat([pl.DataFrame({'_a': a[keep], '_b': b[keep]}),
        search.select(pl.col('_search').cast(pl.Int64).alias('_a'), pl.col('_search').cast(pl.Int64).alias('_b'))])
    del hits
    left = mapping.select(pl.col('_search').cast(pl.Int64).alias('_a'), 'group_id', pl.col('node').cast(pl.Int64).alias('a'))
    right = mapping.select(pl.col('_search').cast(pl.Int64).alias('_b'), 'group_id', pl.col('node').cast(pl.Int64).alias('b'))
    multiplicity = max(1, int(mapping.group_by('_search').len()['len'].max() or 1))
    local = max(1, int(mapping.group_by('_search', 'group_id').len()['len'].max() or 1))
    if multiplicity*local*64 > 256*1024**2:
        raise ValueError('one edit edge membership expansion exceeds 256 MiB budget')
    expanded_edges = global_edges.lazy().join(left.lazy(), on='_a', validate='m:m').join(
        right.lazy(), on=['_b', 'group_id'], validate='m:m').filter(pl.col('a') < pl.col('b')).select(
            'a', 'b', pl.lit(radius, dtype=pl.Int64).alias('distance_ceiling'), 'group_id')
    parts, count = [], 0
    for expanded in expanded_edges.collect_batches(chunk_size=65536, maintain_order=False, engine='streaming'):
        count += expanded.height
        if edge_sink is None and count * 64 > 1024**3:
            raise ValueError('sparse population edges exceed 1 GiB budget')
        if edge_sink is None:
            parts.append(expanded.drop('group_id'))
        elif expanded.height:
            edge_sink(expanded)
    edges = pl.concat(parts) if parts else pl.DataFrame(schema=edge_schema)
    return nodes, edges.sort('a', 'b')


def leiden_taxa(nodes, edges, *, resolution=1., seed=0):
    """Seeded unweighted modularity Leiden within each fixed full-pMHC graph.

    Partitions are fitted to the complete observed graph and then held fixed
    during occurrence rarefaction. Isolates remain singleton communities.
    """
    import igraph as ig
    import leidenalg

    if not np.isfinite(resolution) or resolution <= 0 or seed < 0:
        raise ValueError('positive finite Leiden resolution and nonnegative seed required')
    result = []
    for (group,), selected in nodes.group_by('group_id', maintain_order=True):
        ids = selected['node'].to_numpy()
        local = edges.filter(pl.col('a').is_in(ids) & pl.col('b').is_in(ids))
        endpoints = np.searchsorted(ids, local.select('a', 'b').to_numpy())
        graph = ig.Graph(n=len(ids), edges=endpoints, directed=False)
        membership = leidenalg.find_partition(graph, leidenalg.RBConfigurationVertexPartition,
            resolution_parameter=resolution, n_iterations=-1, seed=seed % 2**31).membership if local.height else list(range(len(ids)))
        result.append(selected.with_columns(pl.Series('community', membership, dtype=pl.Int64)))
    return pl.concat(result) if result else nodes.with_columns(pl.Series('community', [], dtype=pl.Int64))
