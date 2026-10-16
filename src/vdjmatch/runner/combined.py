"""In-memory retained original/gap-block/stem4 component scoring.

The caller supplies and identifies its historical or one-control model context.
No intermediate files or subprocess CLI calls are used here.
"""
from __future__ import annotations

import polars as pl

from ..aggregate.candidates import PMHC

_CASES = ['query_id', *PMHC]
_CASE_LIMIT = 5_000_000


def _check_case_budget(cases):
    if cases > _CASE_LIMIT:
        raise ValueError('combined detailed case budget exceeded')


def _cases(queries, reference, masks):
    groups = reference.select(PMHC).unique(maintain_order=True)
    if masks is None:
        _check_case_budget(queries.height * groups.height)
        return queries.select('query_id').join(groups, how='cross')
    if not set(_CASES).issubset(masks.columns):
        raise ValueError('target masks require query_id and full pMHC')
    masks = masks.select(_CASES).unique(maintain_order=True)
    if masks.select('query_id').unique().join(queries.select('query_id'), on='query_id', how='anti').height:
        raise ValueError('target masks contain unknown query identities')
    _check_case_budget(masks.height)
    return masks


def _valid_single(queries):
    return queries.filter(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False))


def _original(queries, reference, cases, index, prior_background, *, locus, threads, exclude_exact):
    from collections import deque
    from ..match.historical import (density, germline_prior, unpaired_score,
        shared_control_histograms, PreparedBudget)
    from ..match.vgene import precompute_v_similarities

    rows = []
    query = _valid_single(queries)
    similarities = precompute_v_similarities(query['v'], reference['v'])
    work = []
    for key, requested in cases.partition_by(PMHC, as_dict=True).items():
        ref = reference.filter(pl.all_horizontal(*[
            pl.col(name).eq_missing(value) for name, value in zip(PMHC, key)]))
        selected = query.join(requested.select('query_id'), on='query_id', how='semi')
        if ref.height and selected.height:
            work.append((selected, ref))
    prepared = deque()
    counts = shared_control_histograms(index, work, threads=threads,
        prepared_batches=prepared, prepared_budget=PreparedBudget()) if work else None
    exact_cdf = {}
    similarity_rows = None
    if not exclude_exact and work:
        from ..evalue import first_hit
        # Accepted density control CDFs exclude exact identities, so bin zero
        # cannot supply an ordinary exact-edge denominator. One all-query batch
        # includes those identities without changing the retained CDF/model.
        sequences = query['cdr3'].unique(maintain_order=True).to_list()
        histograms = index.edit_histogram_batch(sequences, first_hit.scope(1, 0, 0), threads, False)
        exact_cdf = {sequence: int(row[0]) for sequence, row in zip(sequences, histograms)}
        similarity_rows = pl.DataFrame([(q, r, 1.0 if q == r else 0.25 * value)
            for (q, r), value in similarities.items()], orient='row',
            schema={'_query_family': pl.String, '_reference_family': pl.String, '_exact_weight': pl.Float64})
    for selected, ref in work:
        scored = density(selected, ref, index, threads=threads, pool_reference=True,
            puncture_populations=False, v_similarities=similarities,
            control_counts=counts, prepared_neighbours=prepared.popleft())
        scored = scored.with_columns(pl.col('score').alias('historical_original_density'))
        if not exclude_exact:
            exact = selected.select('query_id', 'cdr3', pl.col('v').str.replace(r'\*.*$', '').fill_null('').alias('_query_family')).join(
                ref.unique('cdr3', keep='first', maintain_order=True).select('cdr3',
                    pl.col('v').str.replace(r'\*.*$', '').fill_null('').alias('_reference_family')),
                on='cdr3', how='inner', validate='m:1').join(similarity_rows,
                    on=['_query_family', '_reference_family'], validate='m:1').join(
                        scored.select('query_id', 'n_reference', 'control_size'), on='query_id', validate='1:1')
            if exact.height:
                exact = exact.with_columns(pl.col('cdr3').replace_strict(exact_cdf, return_dtype=pl.Int64).alias('_exact_cdf'))
                exact = exact.select('query_id', (pl.col('_exact_weight') / pl.max_horizontal(
                    pl.col('n_reference') / pl.col('control_size') * pl.col('_exact_cdf'), pl.lit(0.01)))
                    .alias('original_exact_contribution'))
            else:
                exact = exact.select('query_id').with_columns(pl.lit(None, dtype=pl.Float64).alias('original_exact_contribution'))
            scored = scored.join(exact, on='query_id', how='left', validate='1:1').with_columns(
                pl.col('original_exact_contribution').fill_null(0.0)).with_columns(
                    (pl.col('score') + pl.col('original_exact_contribution')).alias('score'))
        else:
            scored = scored.with_columns(pl.lit(0.0).alias('original_exact_contribution'))
        with_genes = selected.filter(pl.all_horizontal(*[
            pl.col(c).cast(pl.String).str.strip_chars().ne('').fill_null(False) for c in ('v', 'j')]))
        priors = germline_prior(with_genes, ref, prior_background, include_length=locus == 'TRB', exclude_exact=False)
        scored = scored.join(priors.select('query_id', 'germline_lr'), on='query_id', how='left', validate='1:1')
        available = scored.filter(pl.col('score').is_not_null() & pl.col('germline_lr').is_not_null() & (pl.col('n_reference') > 0))
        if available.height:
            available = available.with_columns(unpaired_score(available['score'], available['germline_lr'],
                available['n_reference'], log_scale=True).alias('original_score'))
        else:
            available = available.with_columns(pl.lit(None, dtype=pl.Float64).alias('original_score'))
        scored = scored.rename({'score': 'original_density'}).join(available.select('query_id', 'original_score'),
            on='query_id', how='left', validate='1:1')
        rows.append(scored.select(*_CASES, 'original_score', 'original_density', 'historical_original_density',
            'original_exact_contribution', 'germline_lr').with_columns(
                (pl.col('original_density') > 0).fill_null(False).alias('original_supported')))
    schema = {**{c: cases.schema[c] for c in _CASES}, 'original_score': pl.Float64,
              'original_density': pl.Float64, 'historical_original_density': pl.Float64,
              'original_exact_contribution': pl.Float64, 'germline_lr': pl.Float64, 'original_supported': pl.Boolean}
    scores = pl.concat(rows) if rows else pl.DataFrame(schema=schema)
    return cases.join(scores, on=_CASES, how='left', nulls_equal=True, validate='1:1', maintain_order='left').with_columns(
        pl.col('original_supported').fill_null(False))


def _kernel(queries, reference, cases, *, species, threads, paired, exclude_exact, contexts):
    from ..evalue.kernel import grouped_kernel_enrichment, ladder_cutoffs
    from ..match.geometry import (resolve_v_alleles, distance_units, _substitution_costs)

    model = contexts['clipped']['model']
    if paired:
        q = queries
        r = reference
        for suffix, chain, ref_suffix in [('a', 'TRA', 'alpha'), ('b', 'TRB', 'beta')]:
            q = q.with_columns(resolve_v_alleles(q['v' + suffix], model, locus=chain).alias('_allele' + suffix))
            r = r.with_columns(resolve_v_alleles(r['v_' + ref_suffix], model, locus=chain).alias('_allele_' + ref_suffix))
            q = q.filter(pl.col('_allele' + suffix).is_not_null() & pl.col('cdr3' + suffix).str.contains(
                r'^[ACDEFGHIKLMNPQRSTVWY]{8,}$').fill_null(False))
            r = r.filter(pl.col('_allele_' + ref_suffix).is_not_null() & pl.col(ref_suffix).str.contains(
                r'^[ACDEFGHIKLMNPQRSTVWY]{8,}$').fill_null(False))
        query_keys = ['cdr3a', 'cdr3b', '_allelea', '_alleleb']
        search_keys = ['alpha', 'beta', '_allele_alpha', '_allele_beta']
    else:
        chain = queries['locus'][0] if queries.height else 'TRB'
        q = queries.with_columns(resolve_v_alleles(queries['v'], model, locus=chain).alias('_allele')).filter(
            pl.col('_allele').is_not_null() & pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]{8,}$').fill_null(False))
        r = reference.with_columns(resolve_v_alleles(reference['v'], model, locus=chain).alias('_allele')).filter(
            pl.col('_allele').is_not_null() & pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]{8,}$').fill_null(False))
        query_keys, search_keys = ['cdr3', '_allele'], ['cdr3', '_allele']
    groups = []
    for key, requested in cases.partition_by(PMHC, as_dict=True).items():
        task = dict(zip(PMHC, key))
        part = r.filter(pl.all_horizontal(*[pl.col(name).eq_missing(value) for name, value in task.items()]))
        subset = q.join(requested.select('query_id'), on='query_id', how='semi')
        if task['epitope'] is not None and part.height and subset.height:
            groups.append((task, subset, part))
    result = cases
    for name, matrix, depth in [('gapblock', 'clipped', 0), ('stem4', 'blosum62', 4)]:
        scale, _ = distance_units(_substitution_costs('full', matrix))
        scores = grouped_kernel_enrichment(q, r, groups, query_keys=query_keys, search_keys=search_keys,
            paired=paired, cutoffs=ladder_cutoffs(scale), radius=22.5 * scale, kernel_scale=22.5 * scale,
            scorer='v-bounded', assayed=True, species=species, threads=threads, ctrim=3, v_loop_weight=1,
            alphabet='full', stem_depth=depth, stem_weight=1.0, matrix=matrix,
            v_loop_distances=contexts[matrix], exclude_exact=exclude_exact, _allow_empty_background=True,
            _scores_only=True)
        if scores is None:
            scores = cases.head(0).with_columns(pl.lit(None, dtype=pl.Float64).alias('kernel_score'))
        result = result.join(scores.select(*_CASES, pl.when(pl.col('kernel_score') > 0).then(
            pl.col('kernel_score')).otherwise(None).alias(name + '_score')),
            on=_CASES, how='left', nulls_equal=True, validate='1:1', maintain_order='left')
    return result


def _reference_population(reference, organism):
    if 'species' in reference.columns:
        reference = reference.filter(pl.col('species').str.to_lowercase().str.replace_all(' ', '').is_in(
            ['human', 'homosapiens'] if organism == 'human' else ['mouse', 'musmusculus']))
    if 'reference_valid' in reference.columns:
        reference = reference.filter(pl.col('reference_valid').fill_null(False))
    reference = reference.filter(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False))
    reference = reference.with_columns(*[pl.lit(None, dtype=pl.String).alias(c) for c in PMHC if c not in reference.columns])
    return reference


def score_components(queries, reference, control, *, species, locus, threads=1,
                     exclude_exact=False, target_masks=None):
    """Score retained components in memory, preserving unsupported source cases.

    Queries use canonical single ``query_id,cdr3,v,j,locus`` or paired
    ``query_id,cdr3a,cdr3b,va,vb,ja,jb`` fields. Reference is canonical long VDJdb.
    ``control`` is INTERNAL model context: ``density_index`` is an Index
    for singles or a TRA/TRB Index mapping for pairs; ``prior_background`` is the
    declared V/J prior population. The caller supplies either pinned historical
    assets or its one-control context and records that choice; this function
    does not infer a background. Optional ``v_loop_distances`` contains existing
    clipped/blosum62 model contexts for run-local parameter reuse.

    ``target_masks`` can restrict requested query/full-pMHC cases while retaining
    the whole reference complement. Pair output keeps linked kernels, both full
    same-target unlinked marginals and factorized original log/prior evidence
    separate. No axes are maximized here. Historical density excludes exact
    edges; ordinary scoring adds their explicitly exported contribution once.
    ``exclude_exact`` removes exact junction hits only. Model populations, priors,
    graph summaries and mixture weights remain fixed for all queries.
    Original zero-density log(epsilon)+prior scores remain in continuous outputs,
    with an explicit support flag for call/no-hit handling. No synthetic
    availability flag substitutes for a zero-weight exact edge. Scores are uncalibrated.
    """
    from ..evalue.control import _organism
    from ..match.geometry import precompute_v_loop_distances
    from ..match.historical import paired_score
    from ..evalue.paired import build_paired_ref

    if threads < 1 or locus not in {'TRA', 'TRB', 'paired'}:
        raise ValueError('component scoring requires positive threads and TRA/TRB/paired locus')
    if not isinstance(control, dict) or not {'density_index', 'prior_background'}.issubset(control):
        raise ValueError('accepted internal density_index/prior_background model context required')
    if queries['query_id'].null_count() or queries['query_id'].n_unique() != queries.height:
        raise ValueError('component query identities must be unique and non-null')
    organism = _organism(species)
    reference = _reference_population(reference, organism)
    if locus != 'paired':
        reference = reference.filter(pl.col('gene') == locus)
    cases = _cases(queries, reference, target_masks)
    contexts = control.get('v_loop_distances')
    if contexts is None:
        contexts = {matrix: precompute_v_loop_distances(species=organism, alphabet='full', matrix=matrix,
            v_regions='cdr12plus') for matrix in ('clipped', 'blosum62')}
    def metadata(result):
        return result.join(queries.select('query_id', 'pgen_bin'), on='query_id', how='left', validate='m:1',
            maintain_order='left') if 'pgen_bin' in queries.columns else result
    if locus != 'paired':
        reference = reference.filter(pl.col('gene') == locus)
        query = queries.with_columns(pl.lit(locus).alias('locus'))
        original = _original(query, reference, cases, control['density_index'], control['prior_background'],
            locus=locus, threads=threads, exclude_exact=exclude_exact)
        kernels = _kernel(query, reference, cases, species=organism, threads=threads, paired=False,
            exclude_exact=exclude_exact, contexts=contexts)
        return metadata(original.join(kernels, on=_CASES, how='left', nulls_equal=True, validate='1:1', maintain_order='left'))
    if not isinstance(control['density_index'], dict) or not {'TRA', 'TRB'}.issubset(control['density_index']):
        raise ValueError('paired accepted density context requires TRA and TRB indexes')
    result = cases
    for chain, suffix, name in [('TRA', 'a', 'alpha'), ('TRB', 'b', 'beta')]:
        query = queries.select('query_id', pl.col('cdr3' + suffix).alias('cdr3'),
            pl.col('v' + suffix).alias('v'), pl.col('j' + suffix).alias('j')).with_columns(pl.lit(chain).alias('locus'))
        ref = reference.filter(pl.col('gene') == chain)
        original = _original(query, ref, cases, control['density_index'][chain], control['prior_background'],
            locus=chain, threads=threads, exclude_exact=exclude_exact)
        kernels = _kernel(query, ref, cases, species=organism, threads=threads, paired=False,
            exclude_exact=exclude_exact, contexts=contexts)
        marginal = original.join(kernels, on=_CASES, how='left', nulls_equal=True, validate='1:1')
        result = result.join(marginal.rename({c: c.removesuffix('_score') + '_' + name + '_score'
            if c.endswith('_score') else c + '_' + name for c in marginal.columns if c not in _CASES}),
            on=_CASES, how='left', nulls_equal=True, validate='1:1', maintain_order='left')
    complete = result.filter(pl.all_horizontal(pl.col('original_density_alpha', 'original_density_beta',
        'germline_lr_alpha', 'germline_lr_beta').is_not_null()))
    complete = complete.with_columns(paired_score(complete['original_density_beta'], complete['original_density_alpha'],
        complete['germline_lr_alpha'] + complete['germline_lr_beta'], log_scale=True).alias('original_score'))
    result = result.join(complete.select(*_CASES, 'original_score'), on=_CASES, how='left', nulls_equal=True,
        validate='1:1', maintain_order='left').with_columns(
            (pl.col('original_supported_alpha') & pl.col('original_supported_beta')).alias('original_supported'))
    if 'complex_id' not in reference.columns:
        reference = reference.with_columns(pl.lit('0').alias('complex_id'))
    linked = _kernel(queries, build_paired_ref(reference), cases, species=organism, threads=threads, paired=True,
        exclude_exact=exclude_exact, contexts=contexts).rename({'gapblock_score': 'gapblock_linked_score', 'stem4_score': 'stem4_linked_score'})
    return metadata(result.join(linked, on=_CASES, how='left', nulls_equal=True, validate='1:1', maintain_order='left'))


def _graph_weights(reference, *, chain, threads, hamming_edges=None):
    """One fixed graph summary per pMHC; query expansion belongs to its axis."""
    from ..aggregate.reference_combination import graph_mixture_weights
    from ..aggregate.combination import THREE_WEIGHTS

    groups = reference.select(PMHC).unique().sort(PMHC).with_row_index('group_id')
    records = reference.filter(pl.col('gene') == chain).join(groups, on=PMHC, nulls_equal=True, validate='m:1')
    fields = [*THREE_WEIGHTS, 'weight_status', 'junctions', 'components', 'isolated_junctions',
              'largest_component_fraction', 'component_entropy', 'length_mean', 'length_sd',
              'length_entropy', 'v_genes', 'v_entropy', 'effective_v_genes', 'clonotype_keys', 'source_records']
    if not records.height:
        return groups.with_columns(*[pl.lit(None, dtype=pl.Float64).alias(c) for c in fields]).drop('group_id')
    full, _ = graph_mixture_weights(records.select('group_id', 'cdr3', 'v'), threads=threads, hamming_edges=hamming_edges, include_punctures=False)
    result = groups.join(full.select('group_id', *fields), on='group_id', how='left', validate='1:1', maintain_order='left')
    return result.drop('group_id')


def _select_axis_frames(frames):
    """Fold complete single/paired axes without retaining the stacked wide table."""
    import numpy as np

    selected = keys = None
    names, score_columns = [], []
    for axis, frame in frames:
        if axis not in {'single', 'linked', 'alpha', 'beta'} or axis in names or len(names) >= 3:
            raise ValueError('expected distinct single or paired scoring axes')
        frame = frame.with_columns(pl.lit(axis).alias('axis')).sort(_CASES, nulls_last=True)
        if keys is None:
            keys = frame.select(_CASES)
            best = [np.full(frame.height, np.nan), np.full(frame.height, np.nan)]
            ties = [np.zeros(frame.height, dtype=np.uint8), np.zeros(frame.height, dtype=np.uint8)]
            winner = np.zeros(frame.height, dtype=np.uint8)
        elif not keys.equals(frame.select(_CASES)):
            raise RuntimeError('axis selection did not preserve complete case identities')
        values = [frame[c].cast(pl.Float64).to_numpy() for c in ('score', 'call_score')]
        if any(not frame[c].drop_nulls().is_finite().all() for c in ('score', 'call_score')):
            raise ValueError('axis scores must be finite or null')
        index = len(names)
        if selected is None:
            selected = frame
        else:
            higher = ~np.isnan(values[0]) & (np.isnan(best[0]) | (values[0] > best[0]))
            equal = (values[0] == best[0]) | (np.isnan(values[0]) & np.isnan(best[0]))
            earlier = np.array([axis < name for name in names])[winner]
            take = higher | (equal & earlier)
            selected = pl.concat([selected.filter(~take), frame.filter(take)],
                how='diagonal_relaxed').sort(_CASES, nulls_last=True)
            winner[take] = index
        bit = np.uint8(1 << index)
        for channel, value in enumerate(values):
            higher = ~np.isnan(value) & (np.isnan(best[channel]) | (value > best[channel]))
            equal = ~np.isnan(value) & (value == best[channel])
            ties[channel] = np.where(higher, bit, ties[channel] | np.where(equal, bit, np.uint8(0)))
            best[channel][higher] = value[higher]
        names.append(axis)
        score_columns.append(frame['score'].rename('score_' + axis))
    if selected is None:
        raise ValueError('axis selection requires at least one complete axis')
    masks = {name: [(ties[c] & (1 << i)) != 0 for c in (0, 1)] for i, name in enumerate(names)}
    # Only the selected wide row and scalar per-axis scores survive the fold.
    detail = keys.with_columns(*[
        pl.Series('_tie_' + str(c) + '_' + name, masks[name][c])
        for c in (0, 1) for name in names]).select(
            pl.concat_str([pl.when(pl.col('_tie_0_' + name)).then(pl.lit(name))
                for name in sorted(names)], separator=',', ignore_nulls=True).alias('tied_axes'),
            pl.sum_horizontal([pl.col('_tie_0_' + name).cast(pl.UInt32) for name in names]).alias('tied_axis_count'),
            pl.Series('call_score', best[1], nan_to_null=True),
            pl.concat_str([pl.when(pl.col('_tie_1_' + name)).then(pl.lit(name))
                for name in sorted(names)], separator=',', ignore_nulls=True).alias('tied_call_axes'),
            *score_columns)
    return selected.drop('call_score').hstack(detail).with_columns(
        pl.when(pl.col('call_score').is_not_null()).then(pl.lit('supported'))
        .otherwise(pl.lit('unsupported')).alias('combined_status'))


def combine_components(components, normalizer_components, reference, queries, *, species, locus,
                       threads=1, exclude_exact=False, hamming_edges=None):
    """Normalize matching components, graph-weight each axis, then select its max.

    The normalizer is an independent internally drawn bank, scored on the SAME
    reference/model/exclusion policy. Optional ``pgen_bin`` on both component
    frames conditions normalization on a declared single or paired stratum.
    Component tails use ties and all bank no-hits; they are not final P-values.
    An available original zero-density score remains in component output, but
    cannot alone turn a no-support case into a call. Missing components add zero
    without redistributing weights. Continuous scores retain finite prior evidence;
    call_score masks unsupported axes and has its own maximum and tied_call_axes.

    Paired alpha/beta axes use full same-target marginal models. Linked-axis
    composition weights are the equal mean of those two marginal weight vectors;
    marginal composition models remain fixed, while the linked native hit list
    excludes complete paired identities only. This rule
    never deletes linked reference pairs merely for sharing one query chain.
    Linked/original-factorized and unlinked marginal evidence stay separate until
    the maximum of their complete three-channel sums. Tied axis names survive.
    Final selected-label/axis calibration belongs to the complete caller event.
    Graph summaries stay at pMHC level until consumed. Axes are reduced one at
    a time, retaining the selected metadata, every axis score and both tie sets.
    The detailed query-by-pMHC case budget remains enforced.

    ``hamming_edges`` optionally supplies a transient full-population edge table,
    or a TRA/TRB mapping of such tables, to reuse one union graph across models.
    """
    from ..aggregate.combination import THREE_WEIGHTS, null_normalized_evidence, weighted_component_support
    from ..evalue.control import _organism

    if locus not in {'TRA', 'TRB', 'paired'} or threads < 1:
        raise ValueError('combined evidence requires TRA/TRB/paired and positive threads')
    required = set(_CASES)
    if not required.issubset(components.columns) or not required.issubset(normalizer_components.columns):
        raise ValueError('component evidence requires query identities and full pMHC')
    if components.select(_CASES).is_duplicated().any() or normalizer_components.select(_CASES).is_duplicated().any():
        raise ValueError('component query/group cases must be unique')
    conditioned = 'pgen_bin' in components.columns
    if conditioned != ('pgen_bin' in normalizer_components.columns):
        raise ValueError('observed and normalizer must declare the same Pgen strata')
    species_name = _organism(species)
    reference = _reference_population(reference, species_name)
    if locus == 'paired':
        alpha = _graph_weights(reference, chain='TRA', threads=threads,
            hamming_edges=hamming_edges.get('TRA') if isinstance(hamming_edges, dict) else hamming_edges)
        beta = _graph_weights(reference, chain='TRB', threads=threads,
            hamming_edges=hamming_edges.get('TRB') if isinstance(hamming_edges, dict) else hamming_edges)
        linked = alpha.select(*PMHC, *THREE_WEIGHTS).join(beta.select(*PMHC, *THREE_WEIGHTS),
            on=PMHC, nulls_equal=True, suffix='_beta', validate='1:1').with_columns(*[
                ((pl.col(c) + pl.col(c + '_beta')) / 2).alias(c) for c in THREE_WEIGHTS]).drop(
                    *[c + '_beta' for c in THREE_WEIGHTS])
        definitions = [('linked', 'original_score', 'gapblock_linked_score', 'stem4_linked_score', 'original_supported', linked),
            ('alpha', 'original_alpha_score', 'gapblock_alpha_score', 'stem4_alpha_score', 'original_supported_alpha', alpha),
            ('beta', 'original_beta_score', 'gapblock_beta_score', 'stem4_beta_score', 'original_supported_beta', beta)]
    else:
        weights = _graph_weights(reference, chain=locus, threads=threads,
            hamming_edges=hamming_edges)
        definitions = [('single', 'original_score', 'gapblock_score', 'stem4_score', 'original_supported', weights)]
    def axis_frames():
        columns = ['original_score', 'gapblock_score', 'stem4_score']
        for axis, original, gapblock, stem4, support, weights in definitions:
            metadata = [*_CASES, *(['pgen_bin'] if conditioned else [])]
            raw_columns = [original, gapblock, stem4]
            observed = components.select(*metadata, *[pl.col(c).alias(name) for c, name in zip(raw_columns, columns)],
                pl.col(support).fill_null(False).alias('_original_supported')).join(weights,
                    on=PMHC, how='left', nulls_equal=True, validate='m:1', maintain_order='left')
            null = normalizer_components.select(*metadata, *[pl.col(c).alias(name) for c, name in zip(raw_columns, columns)])
            eligible = pl.all_horizontal(pl.col(THREE_WEIGHTS).is_not_null())
            parts = []
            by = observed.partition_by('pgen_bin', as_dict=True) if conditioned else {(): observed}
            for bin_key, part in by.items():
                bank = null.filter(pl.col('pgen_bin').eq_missing(bin_key[0])) if conditioned else null
                if not bank.height:
                    raise ValueError('normalizer bank is missing a represented Pgen stratum')
                selected = part.filter(eligible)
                if selected.height:
                    if conditioned:
                        selected = selected.with_columns(pl.col('pgen_bin').alias('normalization_stratum'))
                        bank = bank.with_columns(pl.col('pgen_bin').alias('normalization_stratum'))
                    scored = null_normalized_evidence(selected.with_columns(pl.lit(species_name).alias('species'),
                        pl.lit(locus).alias('locus')), bank.with_columns(pl.lit(species_name).alias('species'),
                            pl.lit(locus).alias('locus'), pl.lit('null').alias('outcome'),
                            pl.lit('independent-normalizer-bank').alias('null_source')), three_channel=True)
                    scored = scored.with_columns(pl.col('_original_supported').alias('original_supported'))
                    scored = scored.with_columns(pl.col('normalized_score').alias('score'),
                        pl.when(weighted_component_support(scored)).then(
                            pl.col('normalized_score')).otherwise(None).alias('call_score'))
                    parts.append(scored)
                unavailable = part.filter(~eligible).with_columns(*[pl.lit(None, dtype=pl.Float64).alias(c) for c in ('score', 'call_score')])
                if unavailable.height:
                    parts.append(unavailable)
            axis_scores = pl.concat(parts, how='diagonal_relaxed') if parts else observed.with_columns(*[pl.lit(None, dtype=pl.Float64).alias(c) for c in ('score', 'call_score')])
            yield axis, axis_scores.drop('_original_supported')
    return _select_axis_frames(axis_frames())
