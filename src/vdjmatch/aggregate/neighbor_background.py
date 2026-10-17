"""Uniform same-size clonotype draws, with one native Hamming graph batch."""
import numpy as np
import polars as pl

DRAW_SCHEMA = {name: pl.Int64 for name in ('replicate', 'clonotypes', 'population_clonotypes', 'graph_edges', 'pair_opportunities')}
DRAW_SCHEMA.update({name: pl.Float64 for name in ('pair_density', 'connected_fraction', 'mean_degree')})


def uniform_neighbor_draws(junctions, sizes, *, replicates, seed=0, threads=1):
    """Sample distinct population rows without replacement, ignoring abundance.

    Rows are distinct full-junction/V-gene/J-gene clonotypes; repeated junctions
    retain exact links between different clonotypes. Prefixes of each independent
    maximum-size draw supply all requested N; these nested draws are correlated
    across N. No knee selection, Pgen matching or significance test is implied.
    """
    from .. import cluster

    pool, sizes = list(junctions), np.asarray(sizes)
    if (sizes.ndim != 1 or not sizes.size or sizes.dtype.kind not in 'iu'
            or (sizes < 1).any() or (sizes > len(pool)).any()
            or type(replicates) is not int or replicates < 1
            or type(seed) is not int or seed < 0
            or type(threads) is not int or threads < 1):
        raise ValueError('integer draw sizes within the nonempty population, positive replicates/threads and nonnegative seed required')
    sequences = pl.Series(pool)
    if sequences.dtype != pl.String or not sequences.str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).all():
        raise ValueError('uniform graph draws require valid full amino-acid junctions')
    sizes = np.unique(sizes)
    rng = np.random.default_rng(seed)
    draws = np.stack([rng.choice(len(pool), size=int(sizes[-1]), replace=False) for _ in range(replicates)])
    selected, inverse = np.unique(draws, return_inverse=True)
    union, node_ids = np.unique(np.asarray(pool, dtype=object)[selected], return_inverse=True)
    draw_nodes = node_ids[inverse.reshape(draws.shape)]
    pairs = cluster.overlap(union.tolist(), scope='1,0,0,1', threads=threads)
    edges = pairs.filter(pl.col('a_idx') < pl.col('b_idx')).select('a_idx', 'b_idx').to_numpy()
    rows = []
    for replicate, nodes in enumerate(draw_nodes):
        unique, local_ids = np.unique(nodes, return_inverse=True)
        lookup = np.full(len(union), -1, dtype=np.int64)
        lookup[unique] = np.arange(len(unique))
        local_edges = lookup[edges]
        local_edges = local_edges[(local_edges >= 0).all(axis=1)]
        for n in sizes:
            counts = np.bincount(local_ids[:n], minlength=len(unique))
            inclusive = counts.copy()
            np.add.at(inclusive, local_edges[:, 0], counts[local_edges[:, 1]])
            np.add.at(inclusive, local_edges[:, 1], counts[local_edges[:, 0]])
            degree_sum = int(np.dot(counts, inclusive-1))
            rows.append(dict(replicate=replicate, clonotypes=int(n), population_clonotypes=len(pool),
                graph_edges=degree_sum//2, pair_opportunities=int(n)*(int(n)-1)//2,
                pair_density=degree_sum/(int(n)*(int(n)-1)) if n > 1 else None,
                connected_fraction=float(counts[inclusive > 1].sum()/n),
                mean_degree=degree_sum/int(n)))
    return pl.DataFrame(rows, schema=DRAW_SCHEMA)


def compare_neighbor_draws(summary, draws, *, population_clonotypes):
    """Keep target identities, observed values and empirical same-N denominators."""
    context = [c for c in ('species', 'locus', 'epitope', 'mhc_a', 'mhc_b', 'mhc_class') if c in summary.columns]
    if 'degree_singletons' not in summary.columns:
        summary = summary.with_columns(pl.lit(None, dtype=pl.Int64).alias('degree_singletons'))
    observed = summary.select(*context, pl.col('degree_clonotypes').alias('observed_clonotypes'),
        pl.col('degree_clonotype_edges').alias('observed_graph_edges'),
        pl.when(pl.col('degree_clonotypes') > 0).then(2*pl.col('degree_clonotype_edges')/pl.col('degree_clonotypes')).alias('observed_mean_degree'),
        pl.when(pl.col('degree_clonotypes') > 0).then(1-pl.col('degree_singletons')/pl.col('degree_clonotypes')).alias('observed_connected_fraction'),
        pl.when(pl.col('degree_clonotypes') > 1).then(2*pl.col('degree_clonotype_edges')/
            (pl.col('degree_clonotypes')*(pl.col('degree_clonotypes')-1))).alias('observed_pair_density'))
    metrics = ('pair_density', 'connected_fraction', 'mean_degree')
    baseline = draws.group_by('clonotypes').agg(pl.len().alias('draws'), *[
        expression.alias(stat+'_'+metric) for metric in metrics for stat, expression in
        [('mean', pl.col(metric).mean()), ('median', pl.col(metric).median()),
         ('q10', pl.col(metric).quantile(.1, interpolation='linear')),
         ('q90', pl.col(metric).quantile(.9, interpolation='linear'))]])
    result = observed.join(baseline, left_on='observed_clonotypes', right_on='clonotypes', how='left', validate='m:1')
    return result.with_columns(pl.lit(population_clonotypes).alias('population_clonotypes'),
        pl.col('draws').fill_null(0),
        pl.when(pl.col('observed_clonotypes') == 0).then(pl.lit('empty_reference'))
        .when(pl.col('observed_clonotypes') > population_clonotypes).then(pl.lit('insufficient_background'))
        .otherwise(pl.lit('observed')).alias('status'),
        *[pl.when(pl.col('mean_'+metric) > 0).then(pl.col('observed_'+metric)/pl.col('mean_'+metric))
          .when((pl.col('mean_'+metric) == 0) & (pl.col('observed_'+metric) > 0)).then(float('inf'))
          .alias(metric+'_ratio') for metric in metrics])
