"""Reference-only kth-neighbour SW distances and epsilon knees.

All requested neighbour orders share one fixed full-pMHC count batch. Outputs
are density diagnostics; this command does not fit graphs, communities or calls.
"""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('reference-density', help=__doc__.splitlines()[0])
    p.add_argument('--vdjdb', required=True)
    p.add_argument('--species', default='human')
    p.add_argument('--locus', choices=['TRA', 'TRB'], required=True)
    p.add_argument('--epitope', action='append')
    p.add_argument('--reference-population', choices=['curated-benchmark', 'all-valid'], default='curated-benchmark')
    p.add_argument('--neighbor-order', action='append', type=int, dest='neighbor_orders',
                   help='positive kth other clonotype; repeatable (default: 1, 2, 5, 10)')
    p.add_argument('--cdr-loops', choices=['junction', 'cdr12', 'cdr12plus'], default='cdr12plus')
    p.add_argument('--gap-open', type=int, default=11)
    p.add_argument('--gap-extend', type=int, default=1)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--pair-budget', type=int, default=2_000_000)
    p.add_argument('--output-prefix', required=True)
    from .epitope import add_options
    add_options(p)
    p.set_defaults(func=main, targets_from_sample=None)


def main(a):
    import json
    import time
    from importlib.metadata import version
    import polars as pl
    import seqtree._core as native
    from .. import db
    from ..aggregate.reference_coverage import sw_neighbor_knees
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..match.geometry import kneedle_radius, load_v_loops
    from ..runner.annotation import _validate_output_paths
    from .structure import reference_groups

    started = time.monotonic()
    orders = a.neighbor_orders if a.neighbor_orders is not None else (1, 2, 5, 10)
    if (a.threads < 1 or a.pair_budget < 1 or a.gap_open < 0 or a.gap_extend < 0 or
            not orders or any(k < 1 for k in orders)):
        raise ValueError('positive threads, pair budget and neighbour orders; nonnegative gap costs required')
    prefix = _validate_output_paths(a.output_prefix, [a.vdjdb])
    organism = _organism(a.species)
    species = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[organism]
    from .epitope import load_reference
    raw, epitope_info = load_reference(a, species)
    valid = pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
    if a.reference_population == 'curated-benchmark':
        valid &= pl.col('reference_valid').fill_null(False)
    groups = reference_groups(raw, a)
    records, contexts = [], []
    for group_id, (task, source) in enumerate(groups):
        retained = source.filter(valid).select('cdr3', 'v', 'j').with_columns(
            pl.lit(group_id, dtype=pl.Int64).alias('group_id'))
        records.append(retained)
        contexts.append(dict(group_id=group_id, **task, source_records=source.height,
            retained_records=retained.height, excluded_records=source.height-retained.height))
    records = (pl.concat(records) if records else raw.head(0).select('cdr3', 'v', 'j').with_columns(
        pl.lit(0, dtype=pl.Int64).alias('group_id')))
    loops = () if a.cdr_loops == 'junction' else ('cdr1', 'cdr2') if a.cdr_loops == 'cdr12' else ('cdr1', 'cdr2', 'cdr25')
    curves, knees = sw_neighbor_knees(records, species=organism, locus=a.locus, cdr_loops=loops,
        gap_open=a.gap_open, gap_extend=a.gap_extend, threads=a.threads,
        pair_budget=a.pair_budget, neighbor_orders=orders)
    context_schema = {'group_id': pl.Int64, **{c: pl.String for c in ('epitope', 'mhc_a', 'mhc_b', 'mhc_class')},
        **{c: pl.Int64 for c in ('source_records', 'retained_records', 'excluded_records')}}
    metadata = pl.DataFrame(contexts, schema=context_schema)
    curves = curves.join(metadata, on='group_id', validate='m:1').with_columns(
        pl.lit(species).alias('species'), pl.lit(a.locus).alias('locus'))
    knees = knees.with_columns(pl.lit(species).alias('species'), pl.lit(a.locus).alias('locus'),
        pl.lit('pooled within exact full-pMHC groups').alias('distance_population'))
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for suffix, frame in [('.neighbor-distances.tsv', curves), ('.neighbor-knees.tsv', knees)]:
        path = Path(str(prefix)+suffix)
        frame.write_csv(path, separator='\t')
        outputs[path.name] = sha256(path)
    from .epitope import write_tables
    outputs.update({path.name: sha256(path) for path in write_tables(prefix, epitope_info)})
    manifest = dict(reference=db.provenance(a.vdjdb), species=species, locus=a.locus, epitope=a.epitope,
        reference_population=a.reference_population, targets_sha256=None, **epitope_info, neighbor_orders=sorted(set(orders)),
        geometry='sw', cdr_loops=a.cdr_loops, gap_open=a.gap_open, gap_extend=a.gap_extend,
        distance='full junction + selected V-loop local BLOSUM62 diagonal-sum distances; raw integer units',
        neighbors='other distinct full junction/V/J clonotype vertices within the same exact full-pMHC; self count included once, kth-other is count>k; tied zero distances retained',
        knee='pooled within-full-pMHC kth distances; existing convex Kneedle S=1; no smoothing/flooring/fallback radius',
        unavailable='null distance for unresolved geometry or fewer than k other modeled vertices; null epsilon/status when no knee',
        search_batches=int(knees['modeled_clonotypes'].max() > 0), threads=a.threads, pair_budget=a.pair_budget,
        input_chain_records=raw.height, retained_chain_records=records.height,
        source_sha256=sha256(Path(__file__)), reducer_sha256=sha256(Path(sw_neighbor_knees.__code__.co_filename)),
        native_sha256=sha256(Path(native.__file__)),
        knee_source_sha256=sha256(Path(kneedle_radius.__code__.co_filename)),
        v_model_sha256=(sha256(Path(load_v_loops.__code__.co_filename).parents[1]/'resources/vloops'/f'{organism}_v_loops.tsv')
            if loops else None), tables=outputs,
        software={name: version(name) for name in ('vdjmatch', 'seqtree', 'polars', 'numpy')},
        seconds=time.monotonic()-started)
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0
