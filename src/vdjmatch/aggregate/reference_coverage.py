"""Occurrence-weighted reference communities and vdjtools Chao/iNEXT adapters."""
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
    if (len(x) + total) * 32 + nboot * 41 * 16 > 256 * 1024**2:
        raise ValueError('abundance bootstrap working memory exceeds 256 MiB')
    curve = rarefaction_batch([x], q=0, nboot=nboot, seed=seed, threads=threads).drop('sample')
    depths = curve['m'].to_list()
    coverage_draws = np.empty((nboot, len(depths)))
    for i, draw in enumerate(bootstrap_abundances(x, nboot=nboot, seed=seed)):
        coverage_draws[i] = sample_coverage(draw, m=depths)['sample_coverage'].to_numpy()
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


def sw_knee(records, *, species, locus, cdr_loops, gap_open=11, gap_extend=1, threads=1,
            pair_budget=2_000_000):
    """One reference-only batch; nearest other clonotype within each full pMHC.

    Pool those distances once across the selected species/chain reference. The
    fixed convex Kneedle S=1 helper supplies one radius for real and mock graphs.
    Singleton and unresolved-loop vertices have no nearest distance, never zero.
    """
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
    if not modeled.height:
        raise ValueError('no modeled SW vertices for radius selection')
    groups = int(nodes['group_id'].max())+1
    offsets, targets = _target_csr(modeled, search.height)
    # BLOSUM62 diagonal <=11, local similarity >=0: this encloses every score.
    bound = 22 * int(search['cdr3'].str.len_chars().max()) + int(np.max(distances))
    row_bytes = 4 * ((2+2*groups)*(bound+1) + 2+3*groups)
    chunk_size = min(16, (256*1024**2)//row_bytes)
    if chunk_size < 1:
        raise ValueError('nearest-distance count row exceeds 256 MiB budget')
    nearest = []
    for chunk in gapblock.grouped_count_batch(modeled['cdr3'].to_list(), search['cdr3'].to_list(),
            list(range(bound+1)), reference_target_offsets=offsets, reference_target_ids=targets,
            target_peptide_ids=[0]*groups, reference_exclusive_peptide_ids=[-1]*search.height,
            query_gene_ids=[0]*modeled.height, reference_gene_ids=[0]*search.height,
            query_group_ids=modeled['_vgroup'].to_list(), reference_group_ids=search['_vgroup'].to_list(),
            group_distances=distances, matrix=matrix, alignment_mode='local', gap_open=gap_open,
            gap_extend=gap_extend, threads=threads, chunk_size=chunk_size):
        own = modeled['group_id'][chunk.start:chunk.stop].to_numpy()
        curves = np.asarray(chunk.counts)[np.arange(chunk.stop-chunk.start), 1+own, :]
        found = curves > 1  # only this vertex is excluded; the fixed model is intact
        nearest.extend(int(row.argmax()) if row.any() else None for row in found)
    curve = nodes.join(modeled.select('node').with_columns(
        pl.Series('nearest_distance', nearest, dtype=pl.Int64)), on='node', how='left', validate='1:1')
    radius = kneedle_radius(curve['nearest_distance'].drop_nulls().to_numpy())
    if radius is None:
        raise ValueError('reference nearest-distance curve has no confirmed Kneedle boundary; supply --sw-radius')
    return int(radius), curve


def reference_graph(records, *, geometry, radius=None, species='human', locus='TRB',
                    threads=1, scope='1,1,1,1', gap_open=11, gap_extend=1,
                    cdr_loops=('cdr1', 'cdr2', 'cdr25'), pair_budget=2_000_000):
    """One native graph batch; full pMHC vertices retain source occurrences."""
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
        codes = group_ids * search.height + mapping['_search'].to_numpy()
        order = np.argsort(codes)
        sorted_codes = codes[order]
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
            if edge_count * 64 > 1024**3:
                raise ValueError('sparse population edges exceed 1 GiB budget')
            if len(hits):
                query = hits[:, 0]
                positions = np.searchsorted(sorted_codes, group_ids[query]*search.height+hits[:, 1])
                edge_parts.append(pl.DataFrame({'a': node_ids[query], 'b': node_ids[order[positions]],
                    'distance_ceiling': np.full(len(hits), radius, dtype=np.int64)}))
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
    chunk_size = max(1, min(65536, (256*1024**2)//(multiplicity*local*64)))
    parts, count = [], 0
    for chunk in global_edges.iter_slices(chunk_size):
        expanded = chunk.join(left, on='_a', validate='m:m').join(right,
            on=['_b', 'group_id'], validate='m:m').filter(pl.col('a') < pl.col('b')).select(
                'a', 'b', pl.lit(radius, dtype=pl.Int64).alias('distance_ceiling'))
        count += expanded.height
        if count * 64 > 1024**3:
            raise ValueError('sparse population edges exceed 1 GiB budget')
        parts.append(expanded)
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
        graph = ig.Graph(n=len(ids), edges=endpoints.tolist(), directed=False)
        membership = leidenalg.find_partition(graph, leidenalg.RBConfigurationVertexPartition,
            resolution_parameter=resolution, n_iterations=-1, seed=seed % 2**31).membership if local.height else list(range(len(ids)))
        result.append(selected.with_columns(pl.Series('community', membership, dtype=pl.Int64)))
    return pl.concat(result) if result else nodes.with_columns(pl.Series('community', [], dtype=pl.Int64))
