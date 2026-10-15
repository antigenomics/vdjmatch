"""Source-preserving frozen manuscript marginal scoring in one native batch/chain."""
from __future__ import annotations

import numpy as np
import polars as pl
from seqtree import Index, SearchParams

from ..aggregate.candidates import PMHC
from ..match.classical import frozen_model, gene_key
from ..match.scope import parse_scope

_AA = r'^[ACDEFGHIKLMNPQRSTVWY]+$'
_OUTPUT_SCHEMA = {
    **{name: pl.Float64 for name in ('classical_score', 'classical_pfit', 'classical_informativeness',
        'classical_raw_similarity', 'classical_normalized_similarity')},
    **{name: pl.UInt64 for name in ('classical_reference_records', 'classical_reference_junctions',
        'classical_full_reference_records', 'classical_full_reference_junctions',
        'classical_reference_epitope_roles', 'classical_hit_records', 'classical_modeled_hit_records',
        'classical_hit_junctions', 'classical_excluded_exact_records', 'classical_candidate_junctions', 'classical_nw_pairs')},
    'classical_gap_residues': pl.Int64, 'classical_query_supported': pl.Boolean,
    **{name: pl.String for name in ('classical_best_reference_junction', 'classical_best_reference_v',
        'classical_best_reference_j', 'classical_status')},
}


def _csr(frame, key, values, size):
    """Pack already-sorted numeric rows without materializing per-key lists."""
    counts = np.bincount(frame[key].to_numpy().astype(np.int64), minlength=size)
    return np.r_[0, np.cumsum(counts, dtype=np.int64)], frame[values].to_numpy().astype(np.int64)


def _other_counts(index, roles, lanes, params, max_indels, threads):
    """Distinct other-epitope AA neighbours: global minus epitope-exclusive AA.

    Counts are inclusive of exact AA here; the fused native reducer punctures the
    query from every AA/epitope role and adjusts these counts consistently.
    Histogram calls are one batch per global/exclusive population, not per item.
    """
    needed = lanes.select('cdr3', 'epitope').unique().sort('cdr3', 'epitope')
    sequences = needed['cdr3'].unique().sort().to_list()
    if not sequences:
        return np.empty(0, dtype=np.int64), 0
    calls = 0

    def count(idx, batch):
        nonlocal calls
        try:
            histogram = idx.edit_histogram_batch(batch, params, threads, False, max_indels=max_indels)
        except TypeError as exc:
            raise RuntimeError('trained scoring requires native edit_histogram_batch max_indels support') from exc
        calls += 1
        return np.asarray(histogram, dtype=np.int64).sum(axis=1)

    # Exact-only scope needs no native histogram (whose API requires radius>0).
    if not params.max_total_edits:
        global_counts = np.ones(len(sequences), dtype=np.int64)
    else:
        global_counts = count(index, sequences)
    totals = pl.DataFrame({'cdr3': sequences, '_global': global_counts})
    exclusive = roles.group_by('cdr3').agg(pl.len().alias('_roles'), pl.col('epitope').first()).filter(pl.col('_roles') == 1)
    parts = []
    for (epitope,), part in needed.group_by('epitope', maintain_order=True):
        only = exclusive.filter(pl.col('epitope') == epitope)['cdr3'].sort().to_list()
        batch = part['cdr3'].to_list()
        if only and params.max_total_edits:
            removed = count(Index.build(only, 'aa'), batch)
        elif only:
            present = set(only)
            removed = np.asarray([sequence in present for sequence in batch], dtype=np.int64)
        else:
            removed = np.zeros(len(batch), dtype=np.int64)
        parts.append(part.with_columns(pl.Series('_exclusive', removed)))
    counts = pl.concat(parts).join(totals, on='cdr3', validate='m:1').with_columns(
        (pl.col('_global') - pl.col('_exclusive')).alias('_other'))
    result = lanes.select('cdr3', 'epitope').join(counts.select('cdr3', 'epitope', '_other'),
        on=['cdr3', 'epitope'], validate='m:1', maintain_order='left')['_other'].to_numpy()
    if np.any(result < 0):
        raise ValueError('other-epitope neighbour counts cannot be negative')
    return result, calls


def classical_scores(queries, reference, *, species='human', paired=False,
                     scope='4,2,2,6', max_indels=2, exclude_exact=False, threads=1):
    """Return all source chain rows and flat population provenance for frozen Pfit*I.

    Canonical queries have unique ``query_id,cdr3,v,j,locus`` and full pMHC fields;
    paired mode additionally requires a ``pair_id`` field but retains orphan,
    missing and ambiguous pair assignments. Only trained marginal scores exist.
    A full query AA puncture removes all V/J/epitope roles when requested.

    ``I=-log10((k_other_remaining+1)/N_remaining)`` uses distinct AA neighbours
    and AA/epitope roles in the global species/chain population. Compatible
    assayed pMHC unions restrict requested target lanes before alignment. Counts
    preserve source multiplicities separately. Unknown genes remain unscored;
    they are not replaced by fabricated gene annotations or average offsets.
    """
    from ..cli.historical import assayed_targets
    from ..evalue.control import _organism

    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError('positive integer threads required')
    if isinstance(max_indels, bool) or not isinstance(max_indels, int) or not 0 <= max_indels <= 65535:
        raise ValueError('max_indels must be an integer from zero to 65535')
    s, i, d, t = parse_scope(scope)
    if max(s, i, d, t) > 65535:
        raise ValueError('scope operation budgets cannot exceed 65535')
    if exclude_exact and i != d:
        raise ValueError('exact model punctures require symmetric insertion/deletion caps')
    params = SearchParams(max_subs=s, max_ins=i, max_dels=d, max_total_edits=t, engine='seqtm', matrix='')
    if missing := {'query_id', 'cdr3', 'v', 'j', 'locus', *PMHC} - set(queries.columns):
        raise ValueError('missing canonical query fields: ' + ', '.join(sorted(missing)))
    if missing := {'cdr3', 'v', 'j', 'gene', *PMHC} - set(reference.columns):
        raise ValueError('missing canonical reference fields: ' + ', '.join(sorted(missing)))
    if queries['query_id'].null_count() or queries['query_id'].n_unique() != queries.height:
        raise ValueError('source query IDs must be unique and non-null')
    if paired and 'pair_id' not in queries.columns:
        raise ValueError('paired marginal scoring requires a pair_id field')
    if any(name.startswith('_classical_') or name.startswith('classical_') for name in queries.columns):
        raise ValueError('query fields collide with reserved classical scorer outputs')
    if not hasattr(Index, 'alignment_max_batch'):
        raise RuntimeError('trained scoring requires native Index.alignment_max_batch; install a capable seqtree build')
    organism = _organism(species)
    aliases = ['human', 'homosapiens', 'homo sapiens'] if organism == 'human' else ['mouse', 'musmusculus', 'mus musculus']
    q = queries.with_row_index('_classical_row')
    if 'species' not in q.columns:
        q = q.with_columns(pl.lit(organism).alias('species'))
    selected_species = pl.col('species').cast(pl.String).str.strip_chars().str.to_lowercase().is_in(aliases).fill_null(False)
    q = q.with_columns(selected_species.alias('_classical_species'),
        pl.col('cdr3').cast(pl.String).str.contains(_AA).fill_null(False).alias('_classical_valid'),
        pl.all_horizontal(*[pl.col(name).cast(pl.String).str.strip_chars().ne('').fill_null(False) for name in PMHC]).alias('_classical_target'))
    ref = reference.with_columns(pl.col(PMHC).cast(pl.String))
    if 'species' in ref.columns:
        ref = ref.filter(selected_species)
    ref = ref.filter(pl.col('cdr3').cast(pl.String).str.contains(_AA).fill_null(False) &
        pl.col('epitope').cast(pl.String).str.strip_chars().ne('').fill_null(False))
    if 'reference_valid' in ref.columns:
        ref = ref.filter(pl.col('reference_valid').fill_null(False))
    groups = []; memberships = []; targets = []
    metadata = q.select('query_id', pl.col('species', *PMHC).cast(pl.String))
    for target, (task, requested, selected) in enumerate(assayed_targets(q.select('query_id'), metadata, ref, species)):
        groups.append(requested.select('query_id').with_columns(pl.lit(target, dtype=pl.UInt32).alias('_classical_target_id')))
        memberships.append(selected.with_columns(pl.lit(target, dtype=pl.UInt32).alias('_classical_target_id')))
        targets.append(dict(_classical_target_id=target, **task))
    assigned = pl.concat(groups) if groups else q.head(0).select('query_id').with_columns(pl.lit(None, dtype=pl.UInt32).alias('_classical_target_id'))
    q = q.join(assigned, on='query_id', how='left', validate='1:1', maintain_order='left')
    membership = pl.concat(memberships) if memberships else ref.head(0).with_columns(pl.lit(None, dtype=pl.UInt32).alias('_classical_target_id'))
    target_frame = pl.DataFrame(targets, infer_schema_length=None).with_columns(pl.col('_classical_target_id').cast(pl.UInt32)) if targets else None
    q = q.with_columns(gene_key('v').alias('_classical_v'), gene_key('j').alias('_classical_j'))
    populations = []; outputs = []
    for chain in ('TRA', 'TRB'):
        lane = q.filter((pl.col('locus') == chain) & pl.col('_classical_species'))
        if not lane.height:
            continue
        model = frozen_model(species, chain)
        genes = [{gene: i for i, gene in enumerate(names)} for names in model['group_names']]
        lane = lane.with_columns(pl.col('_classical_v').replace_strict(genes[0], default=-1).alias('_classical_v_id'),
            pl.col('_classical_j').replace_strict(genes[1], default=-1).alias('_classical_j_id'))
        used = lane.filter(pl.col('_classical_valid') & pl.col('_classical_target'))
        keys = ['cdr3', '_classical_v_id', '_classical_j_id']
        unique = used.select(keys).unique().sort(keys).with_row_index('_classical_q')
        source = lane.join(unique, on=keys, how='left', validate='m:1', nulls_equal=True, maintain_order='left')
        requests = source.filter(pl.col('_classical_q').is_not_null() & pl.col('_classical_valid') & pl.col('_classical_target'))\
            .select('_classical_q', '_classical_target_id').unique().sort('_classical_q', '_classical_target_id')
        refs = ref.filter(pl.col('gene') == chain)
        nodes = refs.select('cdr3').unique().sort('cdr3').with_row_index('_classical_r')
        role = refs.select('cdr3', 'epitope').unique()
        selected = membership.filter(pl.col('gene') == chain)
        annotations = selected.group_by('cdr3', 'v', 'j', '_classical_target_id').len(name='_multiplicity').join(
            nodes, on='cdr3', validate='m:1')
        annotations = annotations.join(target_frame.select('_classical_target_id', 'epitope'), on='_classical_target_id', validate='m:1') if target_frame is not None else annotations.with_columns(pl.lit(None, dtype=pl.String).alias('epitope'))
        annotations = annotations.sort('_classical_r', '_classical_target_id', 'v', 'j').with_row_index('_classical_lane')
        annotations = annotations.with_columns(gene_key('v').replace_strict(genes[0], default=-1).alias('_classical_v_id'),
            gene_key('j').replace_strict(genes[1], default=-1).alias('_classical_j_id'))
        sizes = selected.group_by('_classical_target_id').agg(pl.len().alias('classical_reference_records'),
            pl.col('cdr3').n_unique().alias('classical_reference_junctions'))
        source = source.join(sizes, on='_classical_target_id', how='left', validate='m:1', maintain_order='left')
        source = source.join(role.group_by('cdr3').len().rename({'len': '_classical_own_roles'}),
            on='cdr3', how='left', validate='m:1', maintain_order='left').with_columns(
                (pl.lit(role.height) - (pl.col('_classical_own_roles').fill_null(0) if exclude_exact else 0))
                .cast(pl.UInt64).alias('classical_reference_epitope_roles'))
        hist_calls = 0; align_calls = 0; nw_pairs = 0; candidate_count = 0
        if unique.height and nodes.height and requests.height:
            index = Index.build(nodes['cdr3'].to_list(), 'aa')
            other, hist_calls = _other_counts(index, role, annotations, params, max_indels, threads)
            classes = pl.concat([role.select('epitope'), target_frame.select('epitope')]).unique().sort('epitope').with_row_index('_classical_class')
            class_role = role.join(nodes, on='cdr3', validate='m:1').join(classes, on='epitope', validate='m:1', nulls_equal=True).sort('_classical_r', '_classical_class')
            class_targets = target_frame.join(classes, on='epitope', validate='m:1', nulls_equal=True).sort('_classical_target_id')
            q_offsets, q_targets = _csr(requests, '_classical_q', '_classical_target_id', unique.height)
            r_offsets, r_targets = _csr(annotations, '_classical_r', '_classical_target_id', nodes.height)
            c_offsets, c_ids = _csr(class_role, '_classical_r', '_classical_class', nodes.height)
            result = index.alignment_max_batch(unique['cdr3'].to_list(), params, model['similarity'], model['similarity_symbols'], max_indels,
                query_target_offsets=q_offsets, query_target_ids=q_targets,
                reference_target_offsets=r_offsets, reference_target_ids=r_targets,
                target_class_ids=class_targets['_classical_class'].to_numpy(),
                reference_class_offsets=c_offsets, reference_class_ids=c_ids,
                query_group_ids=[unique[name].to_numpy() for name in ('_classical_v_id', '_classical_j_id')],
                reference_group_ids=[annotations[name].to_numpy() for name in ('_classical_v_id', '_classical_j_id')],
                group_offsets=model['group_offsets'], reference_other_counts=other,
                reference_multiplicities=annotations['_multiplicity'].to_numpy(), coefficients=model['coefficients'],
                threads=threads, exclude_exact=exclude_exact)
            align_calls = 1
            arrays = result.arrays()
            nw_pairs = int(np.asarray(arrays['nw_pairs']).sum())
            candidate_count = int(np.asarray(arrays['candidate_counts']).sum())
            names = {'scores': 'classical_score', 'counts': 'classical_hit_records', 'modeled_counts': 'classical_modeled_hit_records',
                'aa_counts': 'classical_hit_junctions', 'exact_counts': 'classical_excluded_exact_records',
                'raw_similarity': 'classical_raw_similarity', 'normalized_similarity': 'classical_normalized_similarity',
                'gap_counts': 'classical_gap_residues', 'pfit': 'classical_pfit', 'informativeness': 'classical_informativeness'}
            evidence = requests.with_columns(*[pl.Series(name, np.asarray(arrays[key]), nan_to_null=True) for key, name in names.items()],
                pl.Series('_classical_best_lane', np.asarray(arrays['best_lane_ids'])))
            evidence = evidence.join(annotations.select(pl.col('_classical_lane').cast(pl.Int64).alias('_classical_best_lane'),
                pl.col('cdr3').alias('classical_best_reference_junction'), pl.col('v').alias('classical_best_reference_v'),
                pl.col('j').alias('classical_best_reference_j')), on='_classical_best_lane', how='left', validate='m:1', maintain_order='left').drop('_classical_best_lane')
            query_counts = unique.select('_classical_q').with_columns(
                pl.Series('classical_candidate_junctions', np.asarray(arrays['candidate_counts'])),
                pl.Series('classical_nw_pairs', np.asarray(arrays['nw_pairs'])),
                pl.Series('classical_query_supported', np.asarray(arrays['query_supported'])))
            source = source.join(query_counts, on='_classical_q', how='left', validate='m:1', maintain_order='left')
            source = source.join(evidence, on=['_classical_q', '_classical_target_id'], how='left', validate='m:1', maintain_order='left')
        else:
            source = source.with_columns(*[pl.lit(None, dtype=pl.Float64).alias(name) for name in
                ('classical_score', 'classical_pfit', 'classical_informativeness', 'classical_raw_similarity', 'classical_normalized_similarity')],
                *[pl.lit(0, dtype=pl.UInt64).alias(name) for name in ('classical_hit_records', 'classical_modeled_hit_records', 'classical_hit_junctions', 'classical_excluded_exact_records', 'classical_gap_residues')],
                *[pl.lit(None, dtype=pl.String).alias(name) for name in ('classical_best_reference_junction', 'classical_best_reference_v', 'classical_best_reference_j')])
        source = source.with_columns(
            pl.col('classical_reference_records').fill_null(0).alias('classical_full_reference_records'),
            pl.col('classical_reference_junctions').fill_null(0).alias('classical_full_reference_junctions'),
            (pl.col('classical_excluded_exact_records').fill_null(0) if exclude_exact else pl.lit(0))
                .alias('classical_excluded_exact_records')).with_columns(
            (pl.col('classical_full_reference_records') - pl.col('classical_excluded_exact_records')).alias('classical_reference_records'),
            (pl.col('classical_full_reference_junctions') - (pl.col('classical_excluded_exact_records') > 0).cast(pl.UInt32))
                .alias('classical_reference_junctions'))
        status = pl.when(~pl.col('_classical_valid')).then(pl.lit('invalid_junction')).when(~pl.col('_classical_target')).then(pl.lit('invalid_target'))\
            .when((pl.col('_classical_v_id') < 0) & (pl.col('_classical_j_id') < 0)).then(pl.lit('unknown_vj_genes'))\
            .when(pl.col('_classical_v_id') < 0).then(pl.lit('unknown_v_gene')).when(pl.col('_classical_j_id') < 0).then(pl.lit('unknown_j_gene'))\
            .when(pl.col('classical_reference_records').fill_null(0) == 0).then(pl.lit('no_reference'))\
            .when((pl.col('classical_hit_records') > 0) & (pl.col('classical_modeled_hit_records') == 0)).then(pl.lit('unknown_reference_genes'))\
            .when(pl.col('classical_hit_records').fill_null(0) == 0).then(pl.lit('no_target_hit')).otherwise(pl.lit('scored'))
        source = source.with_columns(status.alias('classical_status')).with_columns(
            pl.when(pl.col('classical_status').is_in(['no_reference', 'no_target_hit'])).then(0.)
                .when(pl.col('classical_status') == 'scored').then(pl.col('classical_score')).otherwise(None).alias('classical_score'))
        if 'classical_query_supported' not in source.columns:
            source = source.with_columns((pl.col('_classical_valid') & pl.col('_classical_target') &
                (pl.col('_classical_v_id') >= 0) & (pl.col('_classical_j_id') >= 0)).alias('classical_query_supported'))
        source = source.with_columns(pl.col('classical_reference_records', 'classical_reference_junctions').fill_null(0),
            *[pl.when(pl.col('classical_status') == 'scored').then(pl.col(name)).otherwise(None).alias(name) for name in
                ('classical_pfit', 'classical_informativeness', 'classical_raw_similarity', 'classical_normalized_similarity', 'classical_gap_residues')])
        outputs.append(source.select('_classical_row', *[name for name in source.columns if name.startswith('classical_')]))
        populations.append(dict(locus=chain, frozen_species=model['species'], source_query_rows=lane.height, native_query_keys=unique.height,
            source_reference_records=refs.height, reference_junctions=nodes.height, reference_epitope_roles=role.height,
            annotation_lanes=annotations.height, alignment_batches=align_calls, histogram_batches=hist_calls,
            native_nw_pairs=nw_pairs, native_candidate_junctions=candidate_count,
            scope=scope, max_indels=max_indels, selection_gap=model['selection_gap'], exclude_exact=exclude_exact,
            unknown_query_gene_rows=lane.filter((pl.col('_classical_v_id') < 0) | (pl.col('_classical_j_id') < 0)).height,
            unknown_reference_gene_lanes=annotations.filter((pl.col('_classical_v_id') < 0) | (pl.col('_classical_j_id') < 0)).height,
            **{key: value for key, value in model.items() if key.startswith('sha256_')}))
    if outputs:
        evidence = pl.concat(outputs, how='diagonal_relaxed')
        output = q.join(evidence, on='_classical_row', how='left', validate='1:1', maintain_order='left')
    else:
        output = q.with_columns(pl.lit(None, dtype=pl.Float64).alias('classical_score'), pl.lit(None, dtype=pl.String).alias('classical_status'))
    output = output.with_columns(*[pl.lit(None, dtype=dtype).alias(name)
        for name, dtype in _OUTPUT_SCHEMA.items() if name not in output.columns])
    output = output.with_columns(*[pl.col(name).cast(dtype) for name, dtype in _OUTPUT_SCHEMA.items()])
    output = output.with_columns(pl.when(~pl.col('_classical_species')).then(pl.lit('unselected_species'))
        .when(~pl.col('locus').is_in(['TRA', 'TRB']).fill_null(False)).then(pl.lit('unknown_locus'))
        .otherwise(pl.col('classical_status')).alias('classical_status'))
    joint = pl.lit('not_paired')
    if paired:
        pair = q.filter(pl.col('pair_id').is_not_null()).group_by('pair_id').agg(
            pl.col('locus').eq('TRA').sum().alias('_classical_alpha'), pl.col('locus').eq('TRB').sum().alias('_classical_beta'),
            pl.struct(PMHC).n_unique().alias('_classical_pair_targets'))
        output = output.join(pair, on='pair_id', how='left', validate='m:1', maintain_order='left')
        joint = pl.when(pl.col('pair_id').cast(pl.String).str.strip_chars().eq('').fill_null(True)).then(pl.lit('unavailable_missing_pair_id'))\
            .when(pl.col('_classical_pair_targets') > 1).then(pl.lit('unavailable_conflicting_pair_metadata'))\
            .when((pl.col('_classical_alpha') > 1) | (pl.col('_classical_beta') > 1)).then(pl.lit('unavailable_ambiguous_pair'))\
            .when((pl.col('_classical_alpha') == 0) | (pl.col('_classical_beta') == 0)).then(pl.lit('unavailable_incomplete_pair'))\
            .otherwise(pl.lit('unavailable_no_trained_joint_rule'))
    output = output.with_columns(joint.alias('classical_joint_status')).sort('_classical_row')
    return output.drop(*[name for name in output.columns if name.startswith('_classical_')]), {'populations': populations}
