"""Declared-target ranking by two empirical tails over fixed other-peptide receptor keys."""

import numpy as np
import polars as pl

from ..aggregate.candidates import PMHC
from ..cli.historical import assayed_targets
from ..db.schema import mhc_compatible
from ..match.geometry import (
    distance_matrix, impute_reference_v_alleles, paired_distance_matrix,
    precompute_v_loop_distances, resolve_v_alleles,
)
from .kernel import grouped_kernel_enrichment, ladder_cutoffs
from .paired import build_paired_ref

_NO_DISTANCE = np.iinfo(np.int32).max
_NEAREST_BLOCK_CELLS = 32 * 1024**2 // 16
DEFAULT_KERNEL_SCORE_BUDGET_BYTES = 128 * 1024**2


def _prepare(frame, chains, model):
    for c in chains:
        locus = 'TRA' if c == 'a' else 'TRB'
        raw = frame['v' + c]
        allele, imputed = impute_reference_v_alleles(raw, model, locus=locus)
        # A fallback template supplies geometry without merging distinct unknown V calls.
        raw_key = pl.struct(pl.col('v' + c).alias('v')).struct.json_encode()
        frame = frame.with_columns(
            resolve_v_alleles(raw, model, locus=locus).alias('native' + c),
            allele.alias('allele' + c), imputed.alias('imputed' + c),
        ).with_columns(pl.when('imputed' + c).then(raw_key).otherwise('allele' + c).alias('key' + c))
    return frame


def _eligible(chains, allele, minimum):
    return pl.all_horizontal([
        pl.col(allele + c).is_not_null()
        & pl.col('cdr3' + c).str.contains(f'^[ACDEFGHIKLMNPQRSTVWY]{{{minimum},}}$').fill_null(False)
        for c in chains
    ])


def _nearest(queries, targets, chains, context, species, threads, exclude_exact):
    """Reduce bounded native query/reference blocks with a full-junction hit mask."""
    result = np.full((len(targets), len(queries)), _NO_DISTANCE, dtype=np.int32)
    eligible = _eligible(chains, 'native', 5)
    queries = queries.with_row_index('_i').filter(eligible)
    # Contiguous membership rows share query validation and native alignment across targets.
    # A reference belonging to two targets contributes to each target's own minimum.
    columns = ['_junction', *['cdr3' + c for c in chains], *['native' + c for c in chains]]
    if len(targets) == 1:
        refs = targets[0].filter(eligible)
    else:
        refs = pl.concat([target.select(columns).with_columns(pl.lit(i).alias('_target'))
                          for i, target in enumerate(targets)]).filter(eligible)
    distance = paired_distance_matrix if len(chains) == 2 else distance_matrix
    for ref in refs.iter_slices(65536):
        if len(targets) > 1:
            target_ids = ref['_target'].to_numpy()
            starts = np.r_[0, np.flatnonzero(target_ids[1:] != target_ids[:-1]) + 1]
            groups = target_ids[starts]
        refseq = [ref['cdr3' + c].to_list() for c in chains]
        refv = [ref['native' + c].to_list() for c in chains]
        refid = ref['_junction'].to_numpy()
        for part in queries.iter_slices(max(1, _NEAREST_BLOCK_CELLS // len(ref))):
            matrix = distance(*[part['cdr3' + c].to_list() for c in chains], *refseq,
                              *[part['native' + c].to_list() for c in chains], *refv,
                              species=species, threads=threads, allow_short=True,
                              matrix=context['matrix'], v_loop_distances=context)
            if exclude_exact:
                matrix[part['_junction'].to_numpy()[:, None] == refid[None, :]] = _NO_DISTANCE
            idx = part['_i'].to_numpy()
            if len(targets) == 1:
                result[0, idx] = np.minimum(result[0, idx], matrix.min(axis=1))
            else:
                reduced = np.minimum.reduceat(matrix, starts, axis=1).T
                block = np.ix_(groups, idx)
                result[block] = np.minimum(result[block], reduced)
    return result


def _sparse_minima(queries, targets, chains, context, species, threads, exclude_exact):
    # Dense targets already fill native rectangles. Leave them to the dedicated
    # reference/query blocks, which also omit background rows unused by that target.
    small = [i for i, target in enumerate(targets)
             if len(target) * len(queries) <= _NEAREST_BLOCK_CELLS]
    if len(small) < 2:
        return {}
    tile = _nearest(queries, [targets[i] for i in small], chains, context, species, threads, exclude_exact)
    return dict(zip(small, tile))


def equal_tail_scores(queries, reference, *, locus, species, threads, exclude_exact=False,
                      background_scope='all-hla', distance_geometry='native',
                      kernel_score_budget_bytes=DEFAULT_KERNEL_SCORE_BUDGET_BYTES, _consume=None):
    """Score prepared query rows against one curated, species-specific reference.

    Queries carry query_id, species, full candidate pMHC and cdr3a/va and/or cdr3b/vb.
    Linked queries and references use the same pair in both native geometry terms.
    The native distance uses clipped BLOSUM62. ``blosum62-double-cdr12`` uses full
    BLOSUM62 and doubles CDR1/CDR2 in both query and background distances.
    The kernel uses the same fixed geometry for either distance choice. Its score
    slab batches targets within kernel_score_budget_bytes, with a minimum of one target.
    Counts retain exact hits in the model. Missing distance/kernel support contributes
    a unit tail; an absent target or empty calibration population has no combined score.
    """
    if type(kernel_score_budget_bytes) is not int or kernel_score_budget_bytes < 1:
        raise ValueError('kernel_score_budget_bytes must be a positive integer')
    if locus not in ('TRA', 'TRB', 'paired') or background_scope not in ('all-hla', 'compatible-hla'):
        raise ValueError('unsupported equal-tail locus or background scope')
    chains = {'TRA': 'a', 'TRB': 'b', 'paired': 'ab'}[locus]
    paired = len(chains) == 2
    if distance_geometry == 'native':
        dctx = precompute_v_loop_distances(species=species)
    elif distance_geometry == 'blosum62-double-cdr12':
        dctx = precompute_v_loop_distances(species=species, matrix='blosum62', v_regions='double-cdr12')
    else:
        raise ValueError('unsupported equal-tail distance geometry')
    kctx = precompute_v_loop_distances(species=species, matrix='blosum62')
    if paired:
        ref = build_paired_ref(reference).rename({'alpha': 'cdr3a', 'beta': 'cdr3b',
                                                   'v_alpha': 'va', 'v_beta': 'vb'})
    else:
        ref = reference.filter(pl.col('gene') == locus).rename({'cdr3': 'cdr3' + chains, 'v': 'v' + chains})
    ref = _prepare(ref, chains, kctx['model'])
    q = _prepare(queries, chains, kctx['model'])
    seq = ['cdr3' + c for c in chains]
    kernel_keys = [p + c for p in ('cdr3', 'allele', 'key') for c in chains]
    identity = kernel_keys + ['native' + c for c in chains]
    ref = ref.unique(identity + PMHC, maintain_order=True)
    cohort = ref.select(identity).unique(maintain_order=True).with_row_index('_sid')
    # The integer full-junction ID matches both chains jointly, independently of V.
    ids = pl.concat([cohort.select(seq), q.select(seq)]).unique().with_row_index('_junction')
    fixed = cohort.join(ids, on=seq, validate='m:1', nulls_equal=True, maintain_order='left')
    q = q.with_row_index('_sid', offset=len(fixed)).join(
        ids, on=seq, validate='m:1', nulls_equal=True, maintain_order='left')
    ref = ref.join(fixed, on=identity, validate='m:1', nulls_equal=True)
    if _consume is None:
        declared = q.select(pl.col('_sid').alias('query_id'), 'species', *PMHC)
    else:
        # Enumerate the complete species catalogue once, without expanding
        # a query-by-target table. Each candidate uses the ordinary compatible union.
        declared = (reference.filter(pl.col('gene').is_in(['TRA','TRB'])).select(PMHC).unique().drop_nulls()
                    .filter(pl.all_horizontal(pl.col(*PMHC).str.strip_chars() != ''))
                    .sort(PMHC).with_row_index('query_id').with_columns(pl.lit(species).alias('species')))
    groups = list(assayed_targets(declared.select('query_id'), declared, ref, species))
    kernel_ok = _eligible(chains, 'allele', 8)
    query_names = ['cdr3a','cdr3b','_allelea','_alleleb','_v_keya','_v_keyb'] if paired else ['cdr3','_allele','_v_key']
    ref_names = ['alpha','beta','_allele_alpha','_allele_beta','_v_key_alpha','_v_key_beta'] if paired else query_names
    qr = dict(zip(kernel_keys, query_names))
    rr = dict(zip(kernel_keys, ref_names))
    kernel_ref = ref.filter(kernel_ok).select(*kernel_keys, *PMHC).rename(rr)
    allkeys = pl.concat([fixed.select('_sid', *kernel_keys), q.select('_sid', *kernel_keys)])
    kernel_q = allkeys.filter(kernel_ok).rename(qr | {'_sid':'query_id'})
    scale = kctx['distance_units'][0]
    rows = []
    cols = ['_sid', '_junction', *identity]
    complete = pl.concat([fixed.select(cols), q.select(cols)]) if _consume is not None else None
    # Limit the retained int32 minima independently of the requested kernel-score budget.
    nearest_batch = max(1, 128 * 1024**2 // max(1, len(allkeys) * np.dtype(np.int32).itemsize))
    # Cap retained kernel values per target batch; native search itself streams query blocks.
    batch = max(1, kernel_score_budget_bytes // max(1, len(allkeys) * np.dtype(np.float64).itemsize))
    for offset in range(0, len(groups), batch):
        tasks = groups[offset:offset + batch]
        values = np.full((len(tasks), len(allkeys)), np.nan)
        kernel_groups = [(task, part, target.filter(kernel_ok).select(*kernel_keys, *PMHC).rename(rr))
                         for task, part, target in tasks]
        def consume(frame):
            values[frame['_target'].to_numpy(), frame['query_id'].to_numpy()] = frame['kernel_score'].to_numpy()
        if kernel_ref.height:
            grouped_kernel_enrichment(
                kernel_q, kernel_ref, kernel_groups, query_keys=query_names, search_keys=ref_names,
                paired=paired, cutoffs=ladder_cutoffs(scale), radius=22.5*scale, kernel_scale=22.5*scale,
                scorer='v-bounded', _scores_only=True, _consume=consume, v_loop_distances=kctx,
                species=species, threads=threads, ctrim=3, v_loop_weight=1, alphabet='full',
                stem_depth=4, stem_weight=1., matrix='blosum62', exclude_exact=exclude_exact,
            )
        targets = [target.unique('_sid', maintain_order=True) for _, _, target in tasks]
        distances = {}
        for i, (task, part, _) in enumerate(tasks):
            if complete is not None and i % nearest_batch == 0:
                distances.clear()
                distances = _sparse_minima(complete, targets[i:i + nearest_batch], chains,
                                           dctx, species, threads, exclude_exact)
            observed = q if _consume is not None else q.join(part.select(pl.col('query_id').alias('_sid')), on='_sid', how='semi')
            other = ref.filter(pl.col('epitope') != task['epitope'])
            if background_scope == 'compatible-hla':
                other = other.filter((pl.col('mhc_class') == task['mhc_class'])
                                     & mhc_compatible('mhc_a', task['mhc_a'])
                                     & mhc_compatible('mhc_b', task['mhc_b']))
            bg = other.unique('_sid', maintain_order=True)
            target = targets[i]
            combined = pl.concat([bg.select(cols), observed.select(cols)])
            m, n = len(bg), len(target)
            position = i % nearest_batch
            if position not in distances:
                dist = _nearest(combined, [target], chains, dctx, species, threads, exclude_exact)[0]
            else:
                dist = distances[position][combined['_sid'].to_numpy()]
            kernel = values[i, combined['_sid'].to_numpy()]
            missing = np.isnan(kernel[m:])
            kernel = np.nan_to_num(kernel, nan=0.)
            d, k = dist[m:], kernel[m:]
            dn = 1 + np.searchsorted(np.sort(dist[:m]), d, side='right')
            kn = 1 + m - np.searchsorted(np.sort(kernel[:m]), k, side='left')
            dn, kn = dn.astype(np.uint64), kn.astype(np.uint64)
            product = dn * kn
            has_evidence = (d != _NO_DISTANCE) | (k > 0)
            status = 'no_target_reference' if not n else 'empty_background' if not m else 'available'
            result = observed.select('query_id', '_sid').with_columns(
                pl.lit(n, dtype=pl.UInt64).alias('N'), pl.lit(m, dtype=pl.UInt64).alias('M'),
                pl.Series('d', d), pl.Series('k', k), pl.Series('distance_missing', d == _NO_DISTANCE),
                pl.Series('kernel_missing', missing), pl.Series('has_evidence', has_evidence),
                pl.lit(bool(n and m)).alias('calibration_available'),
                pl.when(pl.Series(~has_evidence) & pl.lit(bool(n and m))).then(pl.lit('no_evidence'))
                .otherwise(pl.lit(status)).alias('status'),
            )
            score = -np.log(product.astype(float) / (m+1)**2)
            for name, v in [('D_tail_n', dn), ('K_tail_n', kn), ('product_numerator', product),
                            ('tail_denominator', np.full(len(d), m+1, dtype=np.uint64)),
                            ('equal_tail_score', score)]:
                result = result.with_columns(pl.Series(name, v) if n and m else
                                             pl.lit(None, dtype=pl.Float64 if name=='equal_tail_score' else pl.UInt64).alias(name))
            if _consume is None:
                rows.append(result)
            else:
                _consume(task, result)
    if _consume is not None:
        return
    if not rows:
        raise ValueError('no queries for the selected species and locus')
    return queries.join(pl.concat(rows, how='vertical_relaxed').drop('_sid'), on='query_id',
                        validate='1:1', maintain_order='left')
