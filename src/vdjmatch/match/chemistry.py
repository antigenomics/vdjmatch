"""Experimental V-conditioned local KF4 mixture; no transferred P/E-value.

Each receptor has unit mixture weight over its ordered four-residue core
windows. Paired windows are the Cartesian product within one linked receptor.
Query windows are averaged, never multiplied as independent observations.
"""
import math

import numpy as np
import polars as pl
from scipy.spatial.distance import cdist


def _windows(frame, axes):
    from vdjtools.features.kmer import _explode_kmers
    from vdjtools.features.physchem import load_property_table

    table = load_property_table().select('amino_acid', 'kf4')
    scale = float(table['kf4'].std(ddof=0))
    result = frame.select('_id')
    for axis in axes:
        t = _explode_kmers(frame.select('_id', pl.col(f'seq{axis}').alias('junction_aa')), 4, flank=3)
        for pos in range(4):
            name = f'x{axis}{pos}'
            t = t.with_columns(pl.col('kmer').str.slice(pos, 1).alias('amino_acid')).join(
                table.rename({'kf4': name}), on='amino_acid', how='left', validate='m:1').drop('amino_acid')
            if t[name].null_count():
                raise ValueError('local chemistry requires canonical amino acids')
            t = t.with_columns(pl.col(name) / scale)
        # Repeated windows are retained: the mixture chooses a core position,
        # not a unique word. Per-receptor weights sum to one at any length.
        result = result.join(t.select('_id', *[f'x{axis}{p}' for p in range(4)]),
                             on='_id', how='inner', validate='m:m')
    return result.with_columns((1 / pl.len().over('_id')).alias('_weight'))


def local_chemistry_evidence(queries, reference, background, *, paired=False,
                             exclude_exact=False, bandwidth_factor=1.0):
    """Frames: unique query_id, seq0/v0 [seq1/v1], full anchored junctions.

    Gaussian KDE uses bandwidth_factor times the AA-table KF4 standard
    deviation on each coordinate. Both populations use exactly the same
    geometry. Exact full-junction identity punctures all V categories, target
    and background. Rows without core windows/categories remain explicit.
    The conditional log-density ratio and once-only V usage ratio are separate
    columns. Their sum is experimental ranking, not a biological likelihood
    ratio or confidence. Requires the existing ``bench`` extra's scikit-learn.
    """
    try:
        from sklearn.neighbors import KernelDensity
    except ImportError as error:
        raise RuntimeError('local chemistry requires the vdjmatch[bench] extra') from error

    if isinstance(bandwidth_factor, bool) or not np.isfinite(bandwidth_factor) or bandwidth_factor <= 0:
        raise ValueError('bandwidth factor must be finite and positive')
    if queries['query_id'].null_count() or queries['query_id'].n_unique() != queries.height:
        raise ValueError('unique non-null query IDs required')
    axes = range(2 if paired else 1)
    seq, vs = [f'seq{i}' for i in axes], [f'v{i}' for i in axes]
    coordinates = [f'x{i}{p}' for i in axes for p in range(4)]
    def prepare(frame, query=False):
        f = frame.select(*(['query_id'] if query else []), *seq, *vs)
        if any(f[c].null_count() for c in seq + vs):
            raise ValueError('complete junctions and V genes required')
        if f.filter(pl.any_horizontal(*[~pl.col(c).str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$') for c in seq],
                                      *[pl.col(c).str.strip_chars()=='' for c in vs])).height:
            raise ValueError('canonical junctions and nonempty V genes required')
        f = f.with_columns(*[pl.col(c).str.replace(r'\*.*$', '') for c in vs])
        if not query:
            f = f.unique([*seq, *vs], maintain_order=True)
        return f.with_row_index('_id')
    query_rows, r, b = prepare(queries, True), prepare(reference), prepare(background)
    # Deduplicate search keys within this batch, then restore every query ID.
    q = query_rows.unique([*seq, *vs], maintain_order=True).drop('_id').with_row_index('_id')
    query_map = query_rows.select('query_id', *seq, *vs).join(q.select('_id', *seq, *vs),
        on=[*seq, *vs], how='left', validate='m:1', maintain_order='left').select('query_id', '_id')
    qt, rt, bt = (_windows(f, axes) for f in (q, r, b))
    # Only keys with irreducibly available local windows define this component;
    # the principal query population is retained, including unavailable rows.
    populations = []
    for pop, tokens in ((r, rt), (b, bt)):
        available = pop.join(tokens.select('_id').unique(), on='_id', how='semi')
        populations.append((available, tokens.join(available.select('_id', *vs), on='_id', validate='m:1')))
    K = pl.concat([pop.select(vs) for pop, _ in populations]).unique().height + 1
    out = q.select('_id', 'query_id')
    for side, (pop, _) in zip(('reference', 'control'), populations):
        counts = q.select('_id', *seq, *vs).join(pop.group_by(vs).len().rename({'len': '_n'}),
                                               on=vs, how='left', validate='m:1')
        removed = pl.DataFrame(schema={'_id': pl.UInt32, '_total': pl.Int64, '_v': pl.Int64})
        if exclude_exact and pop.height:
            same = q.join(pop, on=seq, how='inner', suffix='_ref', validate='m:m')
            match = pl.all_horizontal(*[pl.col(v) == pl.col(v + '_ref') for v in vs])
            removed = same.group_by('_id').agg(pl.len().cast(pl.Int64).alias('_total'),
                                               match.sum().cast(pl.Int64).alias('_v'))
        counts = counts.join(removed, on='_id', how='left', validate='1:1').with_columns(
            pl.col('_n', '_total', '_v').fill_null(0).cast(pl.Int64))
        counts = counts.select('_id', (pl.col('_n')-pl.col('_v')).alias('chemistry_'+side+'_v_keys'),
                               (pop.height-pl.col('_total')).alias('chemistry_'+side+'_population'),
                               pl.col('_v').alias('_removed_'+side))
        out = out.join(counts, on='_id', validate='1:1').with_columns(
            pl.lit((r if side=='reference' else b).height-pop.height).alias('chemistry_'+side+'_no_window_keys'))
    self_rows = []
    if exclude_exact:
        for (key,), part in qt.group_by('_id', maintain_order=True):
            points = part.select(coordinates).to_numpy()
            # Bound the transient distance matrix by actual point count/bytes.
            batch = max(1, 32*1024**2 // (8*len(points)))
            total = 0.0
            for start in range(0, len(points), batch):
                total += np.exp(-cdist(points[start:start+batch], points, 'sqeuclidean') /
                                (2*bandwidth_factor**2)).sum()
            self_rows.append((key, total/len(points)**2))
    self_values = pl.DataFrame(self_rows, schema={'_id': pl.UInt32, '_self': pl.Float64}, orient='row')
    out = out.join(self_values, on='_id', how='left', validate='1:1').with_columns(pl.col('_self').fill_null(0))
    norm = len(coordinates) * math.log(math.sqrt(2 * math.pi) * bandwidth_factor)
    mass_rows = []
    for category, part in q.group_by(vs, maintain_order=True):
        mask = pl.all_horizontal(*[pl.col(v) == value for v, value in zip(vs, category)])
        windows = qt.join(part.select('_id'), on='_id', how='semi')
        if not windows.height:
            continue
        x, inverse = np.unique(windows.select(coordinates).to_numpy(), axis=0, return_inverse=True)
        masses = windows.select('_id')
        for side, (pop, tokens) in zip(('reference', 'control'), populations):
            selected, points = pop.filter(mask), tokens.filter(mask)
            values = np.zeros(windows.height)
            if points.height:
                # Identical coordinates have additive mixture weight; retain
                # every receptor/window contribution without repeating KDE work.
                points = points.group_by(coordinates).agg(pl.col('_weight').sum()).sort(coordinates)
                model = KernelDensity(bandwidth=bandwidth_factor, kernel='gaussian',
                                      algorithm='ball_tree', atol=0, rtol=0)
                model.fit(points.select(coordinates).to_numpy(), sample_weight=points['_weight'].to_numpy())
                values = np.exp(model.score_samples(x)[inverse] + norm) * selected.height
            masses = masses.with_columns(pl.Series('_mass_'+side, values))
        mass_rows.append(masses.group_by('_id').agg(pl.col('_mass_reference', '_mass_control').mean(),
                                                    pl.len().alias('chemistry_query_windows')))
    schema = {'_id': pl.UInt32, '_mass_reference': pl.Float64, '_mass_control': pl.Float64,
              'chemistry_query_windows': pl.UInt32}
    masses = pl.concat(mass_rows) if mass_rows else pl.DataFrame(schema=schema)
    out = out.join(masses, on='_id', how='left', validate='1:1')
    for side in ('reference', 'control'):
        mass = pl.col('_mass_'+side)-pl.col('_removed_'+side)*pl.col('_self')
        if out.filter(mass < -64*np.finfo(float).eps*pl.col('chemistry_'+side+'_v_keys').clip(lower_bound=1)).height:
            raise ArithmeticError('negative punctured chemistry mass')
        out = out.with_columns(pl.when(pl.col('chemistry_'+side+'_v_keys')>0).then(
            mass.clip(lower_bound=0)/pl.col('chemistry_'+side+'_v_keys')).otherwise(None).alias('chemistry_'+side+'_density'))
    n, m = pl.col('chemistry_reference_v_keys'), pl.col('chemistry_control_v_keys')
    N, M = pl.col('chemistry_reference_population'), pl.col('chemistry_control_population')
    out = out.with_columns(pl.when(pl.col('chemistry_query_windows').is_null()).then(pl.lit('no_core_window'))
        .when(n==0).then(pl.lit('no_reference_v')).when(m==0).then(pl.lit('no_control_v'))
        .when((pl.col('chemistry_reference_density')<=0)|(pl.col('chemistry_control_density')<=0)).then(pl.lit('numerical_tail'))
        .otherwise(pl.lit('scored')).alias('chemistry_status'))
    pg = (n+m+.5)/(N+M+.5*K)
    out = out.with_columns(pl.when(pl.col('chemistry_status')=='scored').then(
        pl.col('chemistry_reference_density').log()-pl.col('chemistry_control_density').log()).alias('chemistry_log_conditional_ratio'),
        pl.when(pl.col('chemistry_status')=='scored').then(((n+pg)/(N+1)).log()-((m+pg)/(M+1)).log()).alias('chemistry_log_v_ratio'))
    out = out.with_columns((pl.col('chemistry_log_conditional_ratio')+pl.col('chemistry_log_v_ratio')).alias('chemistry_log_score')).drop(
        'query_id', '_self', '_mass_reference', '_mass_control', '_removed_reference', '_removed_control')
    return query_map.join(out, on='_id', how='left', validate='m:1', maintain_order='left').drop('_id')
