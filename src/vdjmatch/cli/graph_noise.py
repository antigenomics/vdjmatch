"""Rank nodes of a supplied graph by triangle-supported degree."""
from pathlib import Path

from . import cpu_budget


def register(subparsers):
    p = subparsers.add_parser('graph-noise', help=__doc__)
    p.add_argument('nodes', type=Path, help='TSV with a unique node_id column')
    p.add_argument('--edges', type=Path, required=True, help='undirected edge TSV with a and b node IDs')
    p.add_argument('--output', type=Path, required=True, help='new TSV with degrees and noise rank')
    p.add_argument('--retain', type=int, help='node budget, with fractional weights for the cutoff tie')
    p.add_argument('--threads', type=cpu_budget, default='0')
    p.set_defaults(func=main)


def main(args):
    import polars as pl
    from ..aggregate.graph_acquisition import triangle_supported_degrees

    nodes = pl.read_csv(args.nodes, separator='\t', infer_schema=False)
    edges = pl.read_csv(args.edges, separator='\t', infer_schema=False)
    reserved = {'degree', 'triangle_degree', 'noise_rank'}
    if args.retain is not None:
        if not 0 <= args.retain <= nodes.height:
            raise ValueError('retain must be between zero and the number of nodes')
        reserved.add('retention_weight')
    if 'node_id' not in nodes or reserved.intersection(nodes.columns) or not {'a', 'b'} <= set(edges.columns):
        raise ValueError('nodes require node_id and unreserved result fields; edges require a and b')
    ids = nodes['node_id']
    if ids.null_count() or ids.n_unique() != nodes.height or (ids.str.strip_chars() == '').any():
        raise ValueError('node IDs must be unique and nonblank')
    index = nodes.select('node_id').with_row_index('position')
    endpoints = edges.select('a', 'b')
    for key in ('a', 'b'):
        endpoints = endpoints.join(index.rename({'node_id': key, 'position': key + '_index'}),
            on=key, how='left', validate='m:1')
        if endpoints[key + '_index'].null_count():
            raise ValueError('edge endpoint is absent from the node table')
    degree, triangle = triangle_supported_degrees(nodes.height, endpoints.select('a_index', 'b_index').to_numpy())
    scores = nodes.with_columns(pl.Series('degree', degree), pl.Series('triangle_degree', triangle))
    order = 2 * pl.col('triangle_degree').cast(pl.Int64) + (pl.col('degree') > 0).cast(pl.Int64)
    scores = scores.with_columns(order.rank('min', descending=True).alias('noise_rank'))
    if args.retain is not None:
        remaining = args.retain + 1 - pl.col('noise_rank').cast(pl.Int64)
        scores = scores.with_columns((remaining / pl.len().over('noise_rank')).clip(0, 1).alias('retention_weight'))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('xb') as stream:
        scores.write_csv(stream, separator='\t')
    return 0
