"""Reference-only nearest-distance knees, with explicit unsupported dispositions."""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('reference-radius', help='Kneedle distance radii from nonexact reference neighbours')
    p.add_argument('--vdjdb', required=True)
    p.add_argument('--species', default='HomoSapiens')
    p.add_argument('--locus', choices=['TRA', 'TRB', 'paired'], required=True)
    p.add_argument('--epitope', action='append')
    p.add_argument('--targets-from-sample', help='profile compatible pMHC unions declared by this AIRR sample')
    p.add_argument('--junction-ends', choices=['trim32', 'trim3'], default='trim3')
    p.add_argument('--matrix', choices=['clipped', 'blosum62', 'pam250', 'mj', 'vdjam'],
                   default='clipped', help='shared substitution costs for junctions, V loops and distance units')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--knee-floor-factor', type=float, default=.5)
    p.add_argument('--coverage-radius-factor', type=float, default=22.5,
                   help='coverage radius in shared substitution-matrix scale units')
    p.add_argument('--output-prefix', required=True)
    p.set_defaults(func=main)


def main(a):
    import json
    import time
    import numpy as np
    import polars as pl
    from .. import db
    from ..aggregate.candidates import PMHC
    from ..evalue.control import _organism
    from ..evalue.paired import build_paired_ref
    from ..match.geometry import (load_v_loops, resolve_v_alleles, reference_nearest_distances,
                                kneedle_radius, regularize_reference_radii, reference_structure,
                                distance_units, _substitution_costs)

    if a.threads < 1:
        raise ValueError('positive threads required')
    if not np.isfinite(a.knee_floor_factor) or a.knee_floor_factor <= 0:
        raise ValueError('knee floor factor must be finite and positive')
    start = time.perf_counter()
    paired = a.locus == 'paired'
    organism = _organism(a.species)
    reference_species = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[organism]
    raw = db.load(a.vdjdb, species=reference_species, gene=None if paired else a.locus, epitope=a.epitope)
    source = build_paired_ref(raw.filter(pl.col('reference_valid'))) if paired else raw
    model = load_v_loops(organism)
    costs = _substitution_costs(matrix=a.matrix)
    matrix_scale = distance_units(costs)[0]
    radii, nearest = [], []
    if a.targets_from_sample:
        from ..io.airr import _read_table
        from .historical import assayed_targets
        targets = _read_table(a.targets_from_sample)
        if 'query_id' not in targets.columns:
            targets = targets.with_row_index('query_id')
        groups = ((tuple(task[c] for c in PMHC), selected) for task, _, selected in
            assayed_targets(targets.select('query_id'), targets, source, a.species))
    else:
        groups = source.group_by(PMHC, maintain_order=True)
    for key, observations in groups:
        fields = []
        for i in range(2 if paired else 1):
            suffix = ('alpha', 'beta')[i]
            fields.extend([pl.col(suffix if paired else 'cdr3').alias(f'seq{i}'),
                pl.col('v_'+suffix if paired else 'v').alias(f'v{i}')])
        frame = observations.select(fields)
        for i in range(2 if paired else 1):
            frame = frame.with_columns(resolve_v_alleles(frame[f'v{i}'], model,
                locus=('TRA', 'TRB')[i] if paired else a.locus).alias(f'v{i}'))
        available = pl.all_horizontal([pl.col(f'v{i}').is_not_null() &
            pl.col(f'seq{i}').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]{8,}$').fill_null(False)
            for i in range(2 if paired else 1)])
        if not paired:
            frame = frame.with_columns(observations['reference_valid'].alias('_valid'))
            available &= pl.col('_valid')
        retained = frame.filter(available).drop('_valid', strict=False)
        frame = retained.unique(maintain_order=True)
        distances = reference_nearest_distances(frame, paired=paired, species=organism,
            threads=a.threads, ctrim=2 if a.junction_ends == 'trim32' else 3, matrix=a.matrix)
        finite = np.isfinite(distances)
        radius = kneedle_radius(distances[finite])
        coverage, shape = reference_structure(distances, radius_factor=a.coverage_radius_factor, costs=costs)
        context = {'species': reference_species, 'locus': a.locus, **dict(zip(PMHC, key))}
        radii.append({**context, 'source_observations': observations.height,
            'source_unit': 'linked_pair_observation' if paired else 'chain_record',
            'unavailable_observations': observations.height-retained.height,
            'reference_keys': frame.height, 'nearest_keys': int(finite.sum()),
            'no_nonexact_neighbour_keys': int((~finite).sum()), 'raw_radius': radius,
            'reference_coverage': coverage, 'knee_shape': shape,
            'radius_status': 'knee' if radius is not None else 'no_knee'})
        nearest.append(frame.with_row_index('reference_key_id').with_columns(
            pl.Series('nearest_distance', distances, nan_to_null=True),
            *(pl.lit(v, dtype=pl.String).alias(k) for k, v in context.items())))
    adjusted, floor, floor_source = regularize_reference_radii(
        [r['raw_radius'] if r['raw_radius'] is not None else np.nan for r in radii],
        np.concatenate([n['nearest_distance'].to_numpy() for n in nearest]) if nearest else [],
        factor=a.knee_floor_factor)
    for row, radius in zip(radii, adjusted):
        row['radius'] = float(radius) if np.isfinite(radius) else None
        row['knee_radius_ratio'] = (radius/(radius+a.coverage_radius_factor*matrix_scale)
                                   if np.isfinite(radius) else 0.)
        row['radius_floor'] = floor
        row['floor_source'] = floor_source
        if row['raw_radius'] is not None:
            row['radius_status'] = ('degenerate_no_scale' if row['radius'] is None else
                'knee_floored' if radius > row['raw_radius'] else 'knee')
    prefix = Path(a.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    context_schema = {c: pl.String for c in ['species', 'locus', *PMHC]}
    radius_schema = {**context_schema, 'source_observations': pl.UInt64, 'source_unit': pl.String,
        **{c: pl.UInt64 for c in ['unavailable_observations', 'reference_keys', 'nearest_keys',
                                  'no_nonexact_neighbour_keys']},
        **{c: pl.Float64 for c in ['raw_radius', 'radius', 'radius_floor',
                                  'reference_coverage', 'knee_shape', 'knee_radius_ratio']},
        'floor_source': pl.String, 'radius_status': pl.String}
    pl.DataFrame(radii, schema=radius_schema).write_csv(str(prefix)+'.radii.tsv', separator='\t')
    if nearest:
        pl.concat(nearest).write_csv(str(prefix)+'.nearest.tsv', separator='\t')
    manifest = {'reference': db.provenance(a.vdjdb), 'options': vars(a) | {'func': 'reference-radius'},
        'source_chain_records': raw.height, 'valid_source_chain_records': int(raw['reference_valid'].sum()),
        'radius_selection': 'reference-only nonexact nearest distances; strongest confirmed offline Kneedle S=1 bend, convex sorted curve, no smoothing',
        'structure_features': {'coverage_radius_factor': a.coverage_radius_factor,
            'matrix_scale': matrix_scale,
            'coverage_denominator': 'all distinct full-junction/resolved-V reference keys',
            'knee_shape': 'maximum normalized convex bend; zero for flat/short curves',
            'knee_radius_ratio': 'regularized epsilon/(epsilon+matrix-relative coverage radius); zero when no knee, with radius_status retained'},
        'degenerate_knee_guard': {'factor': a.knee_floor_factor, 'floor': floor,
            'source': floor_source, 'scope': 'same species, locus and junction geometry',
            'raised_knees': sum(r['radius_status'] == 'knee_floored' for r in radii)},
        'reference_unit': 'distinct full junctions and resolved V alleles; same linked pair for paired mode',
        'geometry': f'uniform {a.matrix} CDR3 plus V loops; paired alpha+beta sum',
        'no_knee': 'no radius is selected; reference keys remain available',
        'probability_interpretation': 'none; selected radii carry no transferred calibration',
        'seconds': time.perf_counter()-start}
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0
