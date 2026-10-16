"""Occurrence-weighted reference communities and vdjtools Chao/iNEXT adapters."""
import math

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
    from vdjtools.stats.inext import _bootstrap_probs

    x = np.asarray(counts)
    if (any(isinstance(v, bool) or not isinstance(v, int) for v in (nboot, seed, threads)) or
            x.ndim != 1 or not x.size or x.dtype.kind not in 'iu' or
            np.any(x <= 0) or nboot < 2 or seed < 0 or threads < 1):
        raise ValueError('positive integer occurrence counts, nboot>=2, seed>=0 and threads>=1 required')
    total = sum(map(int, x))
    if total >= 2**31:
        raise ValueError('occurrence total exceeds vdjtools native bootstrap int32 range')
    estimate = float(asymptotic_diversity(x, q=0, se=False)['estimator'].item())
    # The owning augmented law introduces ceil(estimated unseen richness) taxa.
    # Bound its probabilities/counts before either native or Python bootstrap.
    if 6 * math.ceil(estimate) * 8 + nboot * 41 * 16 > 256 * 1024**2:
        raise ValueError('augmented Chao bootstrap exceeds 256 MiB; reduce the selected population')
    curve = rarefaction_batch([x], q=0, nboot=nboot, seed=seed, threads=threads).drop('sample')
    depths = curve['m'].to_list()
    probabilities = _bootstrap_probs(x.astype(float))
    probabilities /= probabilities.sum()
    rng = np.random.default_rng(seed)
    coverage_draws = np.empty((nboot, len(depths)))
    for i in range(nboot):
        draw = rng.multinomial(total, probabilities)
        coverage_draws[i] = sample_coverage(draw[draw > 0], m=depths)['sample_coverage'].to_numpy()
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


def reference_graph(records, *, geometry, radius=None, species='human', locus='TRB',
                    threads=1, scope='1,1,1,1', gap_open=11, gap_extend=1,
                    cdr_loops=('cdr1', 'cdr2', 'cdr25'), pair_budget=2_000_000):
    """One population batch; sparse graph on distinct full-junction/V/J calls.

    ``group_id`` identifies full pMHC strata; cross-stratum edges are discarded.
    Full SW uses the sum of raw BLOSUM62 self-score distances over the complete
    junction and requested germline loops, without trimming or kernel surrogates.
    Edge ``distance_ceiling`` is the inclusive graph radius, not an exact score.
    """
    from ..cluster import overlap
    from ..match.geometry import load_v_loops, resolve_v_alleles, require_geometry_capabilities
    from ..match.scope import parse_scope
    from seqtree import SubstitutionMatrix, gapblock
    from seqtree.pairwise import dist_matrix

    if geometry not in ('levenshtein', 'sw') or threads < 1 or pair_budget < 1:
        raise ValueError('known geometry, positive threads and pair budget required')
    if records.filter(~pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)).height:
        raise ValueError('graph requires valid full junctions')
    nodes = records.group_by('group_id', 'cdr3', 'v', 'j').len(name='occurrences').sort(
        'group_id', 'cdr3', 'v', 'j', nulls_last=True).with_row_index('node')
    search = nodes.select('cdr3', 'v', 'j').unique().sort('cdr3', 'v', 'j', nulls_last=True).with_row_index('_search')
    mapping = nodes.join(search, on=['cdr3', 'v', 'j'], nulls_equal=True, validate='m:1')
    n = search.height
    if n*n > pair_budget:
        raise ValueError(f'graph pair budget exceeded: {n*n} comparisons > {pair_budget}; select fewer pMHC targets or raise --pair-budget explicitly')
    edge_schema = {'a': pl.Int64, 'b': pl.Int64, 'distance_ceiling': pl.Int64}
    sequences = search['cdr3'].to_list()
    if geometry == 'levenshtein':
        _, insertions, deletions, radius = parse_scope(scope)
        if insertions != deletions:
            raise ValueError('undirected edit graphs require symmetric insertion/deletion budgets')
        hits = overlap(sequences, scope=scope, threads=threads) if n else None
        a = hits['a_idx'].to_numpy() if n else np.array([], dtype=int)
        b = hits['b_idx'].to_numpy() if n else np.array([], dtype=int)
    else:
        require_geometry_capabilities('sw')
        if (not isinstance(radius, int) or isinstance(radius, bool) or not 0 <= radius < 2**28 or
                gap_open < 0 or gap_extend < 0 or not set(cdr_loops) <= {'cdr1', 'cdr2', 'cdr25'}):
            raise ValueError('SW requires an explicit nonnegative integer radius below 2**28, nonnegative gaps and known CDR loops')
        if not n:
            return nodes, pl.DataFrame(schema=edge_schema)
        matrix = SubstitutionMatrix.blosum62()
        group_ids, distances = [0]*n, [[0]]
        if cdr_loops:
            model = load_v_loops(species)
            alleles = resolve_v_alleles(search['v'], model, locus=locus)
            if alleles.null_count():
                raise ValueError(f'full SW germline loops unavailable for {alleles.null_count()} clonotypes; supply resolved V calls or explicitly select --cdr-loops junction')
            unique = alleles.unique().sort().to_list()
            if len(unique)**2 * 32 > 64 * 1024**2:
                raise ValueError('germline loop distance tables exceed 64 MiB')
            table = np.zeros((len(unique), len(unique)), dtype=np.int64)
            for loop in cdr_loops:
                loops = [model[a][loop].replace('.', '') for a in unique]
                if any(not s or not set(s) <= set('ACDEFGHIKLMNPQRSTVWY') for s in loops):
                    raise ValueError('full SW requires complete canonical germline-loop templates')
                table += np.asarray(dist_matrix(loops, loops, matrix, mode='local',
                    gap_open=gap_open, gap_extend=gap_extend, threads=threads))
            group_ids = alleles.replace_strict(unique, list(range(len(unique)))).cast(pl.Int64).to_list()
            distances = table.tolist()
        # One native full-alignment traversal. Each reference key is one target,
        # so bounded binary target counts expose edges without an N x N matrix.
        rows_a, rows_b = [], []
        chunks = gapblock.grouped_count_batch(sequences, sequences, [radius],
            reference_target_offsets=list(range(n+1)), reference_target_ids=list(range(n)),
            target_peptide_ids=[0]*n, reference_exclusive_peptide_ids=[-1]*n,
            query_gene_ids=[0]*n, reference_gene_ids=[0]*n, chunk_size=32,
            matrix=matrix, gap_open=gap_open, gap_extend=gap_extend, alignment_mode='local',
            query_group_ids=group_ids, reference_group_ids=group_ids, group_distances=distances,
            threads=threads)
        for chunk in chunks:
            row, col = np.nonzero(np.asarray(chunk.counts)[:, 1:1+n, 0])
            row = row + chunk.start
            keep = row < col
            rows_a.append(row[keep]); rows_b.append(col[keep])
        a, b = np.concatenate(rows_a), np.concatenate(rows_b)
    keep = a < b
    global_edges = pl.DataFrame({'_a': a[keep], '_b': b[keep]})
    # Reuse the same searched key in every real/mock graph, retaining each
    # population's own vertex and occurrence count. No graph-search reruns.
    left = mapping.select(pl.col('_search').cast(pl.Int64).alias('_a'), 'group_id', pl.col('node').cast(pl.Int64).alias('a'))
    right = mapping.select(pl.col('_search').cast(pl.Int64).alias('_b'), 'group_id', pl.col('node').cast(pl.Int64).alias('b'))
    worst_edges = sum(count*(count-1)//2 for count in nodes.group_by('group_id').len()['len'])
    expanded = int(np.dot(np.bincount(a[keep], minlength=n),
        np.bincount(mapping['_search'].to_numpy(), minlength=n)))
    if max(worst_edges, expanded) * 64 > 256 * 1024**2:
        raise ValueError('population edge expansion exceeds 256 MiB; select fewer targets/replicates')
    edges = global_edges.join(left, on='_a', validate='m:m').join(right,
        on=['_b', 'group_id'], validate='m:1').select('a', 'b', pl.lit(radius, dtype=pl.Int64).alias('distance_ceiling'))
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
