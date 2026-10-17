"""Fresh calibration of the complete reference-weighted annotation statistic."""
from __future__ import annotations

import numpy as np
import polars as pl

from ..aggregate.candidates import PMHC
from .mock_reference import control_events, compute_pgen, pgen_strata, iter_mock_models

_KEYS = ['query_id', *PMHC]


def _query_banks(queries, controls, probabilities, *, locus, replicates, seed):
    """Independent normalization bank and Pgen-matched test-query bank."""
    paired = locus == 'paired'
    lanes = [('TRA', 'a'), ('TRB', 'b')] if paired else [(locus, '')]
    tagged = queries
    pools = {}
    for chain, suffix in lanes:
        bins = pgen_strata(probabilities[chain]).rename({'pgen_bin': '_bin' + suffix})
        tagged = tagged.join(bins, left_on='cdr3' + suffix, right_on='junction_aa',
                             how='left', validate='m:1')
        if tagged['_bin' + suffix].null_count():
            raise ValueError('Pgen map does not cover every query')
        pool = controls[chain].join(bins.rename({'_bin' + suffix: '_bin'}),
                                   on='junction_aa', validate='m:1').sort('junction_aa', 'v', 'j')
        pools[chain] = pool.partition_by('_bin', as_dict=True)
    columns = ['_bin' + suffix for _, suffix in lanes]
    tagged = tagged.with_columns(pl.concat_str([pl.col(c).cast(pl.String) for c in columns],
                                              separator=':').alias('pgen_bin'))
    streams = np.random.SeedSequence(seed).spawn(2)
    banks = []
    for name, stream in zip(('normalize', 'test'), streams):
        rng = np.random.default_rng(stream)
        parts = []
        if name == 'normalize':
            fields = {'query_id': [f'__calibration_normalize_{i}' for i in range(replicates)],
                      'replicate': np.arange(replicates)}
            for chain, suffix in lanes:
                pool = controls[chain].sort('junction_aa', 'v', 'j')
                weight = pool['event_weight'].to_numpy().astype(float)
                selected = pool[rng.choice(pool.height, replicates, replace=True, p=weight / weight.sum())]
                for source, dest in [('junction_aa', 'cdr3'), ('v', 'v'), ('j', 'j')]:
                    fields[dest + suffix] = selected[source]
            if not paired:
                fields['locus'] = [locus] * replicates
            banks.append(pl.DataFrame(fields))
            continue
        for number, (key, part) in enumerate(tagged.sort('pgen_bin').partition_by('pgen_bin', as_dict=True).items()):
            fields = {'query_id': [f'__calibration_{name}_{number}_{i}' for i in range(replicates)],
                      'pgen_bin': [key[0]] * replicates, 'replicate': np.arange(replicates)}
            for chain, suffix in lanes:
                pool = pools[chain].get((part['_bin' + suffix][0],))
                if pool is None:
                    raise ValueError(f'control has no {chain} events in query Pgen stratum')
                weights = pool['event_weight'].to_numpy().astype(float)
                selected = pool[rng.choice(pool.height, replicates, replace=True, p=weights / weights.sum())]
                for source, dest in [('junction_aa', 'cdr3'), ('v', 'v'), ('j', 'j')]:
                    fields[dest + suffix] = selected[source]
            if not paired:
                fields['locus'] = [locus] * replicates
            parts.append(pl.DataFrame(fields))
        banks.append(pl.concat(parts))
    return tagged.drop(columns), *banks


def _overall(scores):
    """Repeat full-label selection for every observed/randomized receptor."""
    return scores.group_by('query_id').agg(pl.col('score').max()).with_columns(
        *[pl.lit(None, dtype=scores.schema[c]).alias(c) for c in PMHC])


def _with_overall(scores):
    # Directional tails and match rates use the complete support-gated event.
    # The unmasked continuous score is retained separately in observed_detail.
    scores = scores.select(*_KEYS, pl.col('call_score').alias('score'))
    return pl.concat([scores.select(*_KEYS, 'score').with_columns(pl.lit('epitope').alias('test_scope')),
        _overall(scores).select(*_KEYS, 'score').with_columns(pl.lit('database').alias('test_scope'))])


def calibrate_combined(queries, reference, raw_control, *, species, locus, replicates=200,
                       seed=0, threads=1, exclude_exact=False, model_context=None,
                       detail_targets=None,report_p=.01):
    """Calibrate the actual sum and axis/label selection from one raw control.

    One independent unconditional internal normalization bank is held fixed in both directions.
    Mock references retain all source roles, full pMHC and random chain pairing.
    Test clonotypes have the observed chain-wise Pgen bins and independently
    drawn TRA/TRB events. No-hit draws remain in every denominator. Returned
    database tails test the maximum over all searched labels; epitope tails test
    a fixed full-pMHC group. The maximum of the two directional tails tests their
    union null (no independence assumption); it is not a binding posterior.

    Default model context explicitly uses the supplied control for original
    density and germline priors (beta priors for paired inputs). This candidate recipe is distinct from the
    historical benchmark's separately pinned original prior. An internal context
    can reproduce that comparator without another end-user calibration input.

    ``report_p`` is the configurable any-call reporting threshold (default 0.01).
    It filters compact reported specificities only; declared targets remain.
    """
    from seqtree import Index
    from .control import _organism
    from ..runner.combined import score_components, combine_components, _reference_population, _check_case_budget
    from ..runner.component_bank import ComponentBank

    if not np.isfinite(report_p) or not 0 < report_p < 1:
        raise ValueError('report_p must be finite and between zero and one')
    if replicates < 1 or threads < 1 or not queries.height:
        raise ValueError('combined calibration requires queries and positive replicates/threads')
    if queries['query_id'].n_unique() != queries.height or queries['query_id'].null_count():
        raise ValueError('combined calibration requires unique non-null query IDs')
    if queries['query_id'].cast(pl.String).str.starts_with('__calibration_').any():
        raise ValueError('query IDs use the reserved calibration prefix')
    queries = queries.with_columns(pl.col('query_id').cast(pl.String))
    chains = ['TRA', 'TRB'] if locus == 'paired' else [locus]
    reference = _reference_population(reference, _organism(species))
    reference = reference.filter(pl.col('gene').is_in(chains))
    groups = reference.select(PMHC).n_unique()
    if detail_targets is None:
        _check_case_budget((queries.height + replicates) * groups)
    controls = {chain: control_events(raw_control, species=species, locus=chain) for chain in chains}
    probabilities = {}
    for chain in chains:
        suffix = ('a' if chain == 'TRA' else 'b') if locus == 'paired' else ''
        sequences = pl.concat([reference.filter(pl.col('gene') == chain).select(
            pl.col('cdr3').alias('junction_aa')), queries.select(pl.col('cdr3' + suffix).alias('junction_aa'))])
        probabilities[chain] = compute_pgen(sequences, controls[chain], species=species, locus=chain, threads=threads)
    query_seed, reference_seed = np.random.SeedSequence(seed).spawn(2)
    tagged, normalizers, test_queries = _query_banks(queries, controls, probabilities, locus=locus,
        replicates=replicates, seed=int(query_seed.generate_state(1)[0]))
    _check_case_budget((tagged.height + normalizers.height + test_queries.height) * groups
        if detail_targets is None else normalizers.height*groups)
    pgen_columns = []
    for chain in chains:
        suffix = ('a' if chain == 'TRA' else 'b') if locus == 'paired' else ''
        name = 'pgen_' + ('alpha' if chain == 'TRA' else 'beta') if locus == 'paired' else 'pgen'
        tagged = tagged.join(probabilities[chain].rename({'junction_aa': 'cdr3' + suffix, 'pgen': name}),
            on='cdr3' + suffix, how='left', validate='m:1')
        pgen_columns.append(name)
    if model_context is None:
        indexes = {chain: Index.build(controls[chain]['junction_aa'].unique().sort().to_list(), 'aa') for chain in chains}
        beta = controls['TRB' if locus == 'paired' else locus]
        # Preserve event multiplicity as integer weights, without expanding rows.
        prior = beta.select(pl.col('junction_aa').alias('cdr3'), 'v', 'j', 'event_weight')
        model_context = {'density_index': indexes if locus == 'paired' else indexes[locus], 'prior_background': prior}
    from ..match.geometry import precompute_v_loop_distances
    from ..match.historical import _prepare_prior_background
    model_context = dict(model_context)
    model_context['_prior_counts'] = _prepare_prior_background(model_context['prior_background'])
    if 'v_loop_distances' not in model_context:
        model_context['v_loop_distances'] = {matrix: precompute_v_loop_distances(species=species,
            alphabet='full', matrix=matrix, v_regions='cdr12plus') for matrix in ('clipped', 'blosum62')}
    def components(q, model):
        return score_components(q, model, model_context, species=species, locus=locus,
            threads=threads, exclude_exact=exclude_exact,_bank=True)
    def subset(raw,q):
        return (pl.concat(raw.frames(q)) if isinstance(raw,ComponentBank) else
                raw.join(q.select('query_id'),on='query_id',how='semi'))
    def unstratified(raw):
        return raw if isinstance(raw,ComponentBank) else raw.drop('pgen_bin')
    if detail_targets is not None:
        from .compact_combined import calibrate_compact
        return calibrate_compact(tagged,normalizers,test_queries,reference,components=components,
            combine=combine_components,models=iter_mock_models(reference,controls,pgen=probabilities,
                replicates=replicates,seed=int(reference_seed.generate_state(1)[0])),
            detail_targets=detail_targets,replicates=replicates,report_p=report_p,pgen_columns=pgen_columns,
            combine_options=dict(species=species,locus=locus,threads=threads,exclude_exact=exclude_exact))
    bank = pl.concat([normalizers.drop('replicate'), test_queries.drop('replicate')], how='diagonal_relaxed')
    all_queries = pl.concat([tagged, bank], how='diagonal_relaxed')
    real = components(all_queries, reference)
    norm = subset(real,normalizers)
    combined = combine_components(unstratified(real), unstratified(norm) if 'pgen_bin' in norm.columns else norm, reference, all_queries, species=species, locus=locus,
        threads=threads, exclude_exact=exclude_exact)
    observed_detail = combined.join(tagged.select('query_id', 'pgen_bin', *pgen_columns),
        on='query_id', how='inner', validate='m:1')
    observed = _with_overall(observed_detail).with_row_index('_case')
    exceed = np.zeros(observed.height, dtype=np.int64)
    reference_hit_draws = np.zeros(observed.height, dtype=np.int64)
    observed_queries = pl.concat([tagged, normalizers.drop('replicate')], how='diagonal_relaxed')
    test = _with_overall(combined.join(test_queries.select('query_id'), on='query_id', how='semi')).join(
        test_queries.select('query_id', 'pgen_bin'), on='query_id', validate='m:1')
    lookup_keys = ['pgen_bin', *PMHC, 'test_scope']
    tagged_observed = observed.join(tagged.select('query_id', 'pgen_bin'), on='query_id', validate='m:1')
    if (test.group_by(lookup_keys).len().filter(pl.col('len') != replicates).height or
            tagged_observed.select(lookup_keys).unique().join(test.select(lookup_keys).unique(),
                on=lookup_keys, how='anti', nulls_equal=True).height):
        raise RuntimeError('query-null scoring did not preserve every draw in every testing group')
    distributions = {key: np.sort(part['score'].drop_nulls().to_numpy()) for key, part in
        test.partition_by(lookup_keys, as_dict=True).items()}
    query_exceed = np.zeros(observed.height, dtype=np.int64)
    query_hit_draws = np.zeros(observed.height, dtype=np.int64)
    for key, part in tagged_observed.partition_by(lookup_keys, as_dict=True).items():
        values = distributions[key]
        scores = part['score'].to_numpy()
        counts = len(values) - np.searchsorted(values, scores, side='left')
        counts[np.isnan(scores)] = replicates
        query_exceed[part['_case'].to_numpy()] = counts
        query_hit_draws[part['_case'].to_numpy()] = len(values)
    # Retain requested observed detail and compact tail counts only. Wide
    # real/bank frames must not overlap every subsequent mock population.
    del real, norm, combined, bank, all_queries, test, tagged_observed, test_queries, distributions
    epitope_observed = observed.filter(pl.col('test_scope') == 'epitope').select(*_KEYS, '_case')
    database_observed = observed.filter(pl.col('test_scope') == 'database').with_row_index('_query')
    observed_values = observed['score'].to_numpy()
    for model in iter_mock_models(reference, {chain: controls[chain] for chain in chains},
            pgen=probabilities, replicates=replicates, seed=int(reference_seed.generate_state(1)[0])):
        raw = components(observed_queries, model)
        mock_norm = subset(raw,normalizers)
        seen = np.zeros(observed.height, dtype=bool)
        maxima = np.full(database_observed.height, np.nan)

        def count(indices, values):
            reference_hit_draws[indices] += ~np.isnan(values)
            exceed[indices] += np.isnan(observed_values[indices]) | (values >= observed_values[indices])

        def consume(frame):
            scores = frame.join(tagged.select('query_id'), on='query_id', how='semi').select(
                *_KEYS, pl.col('call_score').alias('_null_score'))
            if scores.select(_KEYS).is_duplicated().any():
                raise RuntimeError('randomized scoring did not preserve the complete observed case grid')
            matched = scores.join(epitope_observed, on=_KEYS, how='left', nulls_equal=True, validate='1:1')
            if matched['_case'].null_count():
                raise RuntimeError('randomized scoring did not preserve the complete observed case grid')
            indices = matched['_case'].to_numpy()
            if seen[indices].any():
                raise RuntimeError('randomized scoring did not preserve the complete observed case grid')
            seen[indices] = True
            count(indices, matched['_null_score'].to_numpy())
            by_query = scores.join(database_observed.select('query_id', '_query'),
                on='query_id', how='left', validate='m:1')
            np.fmax.at(maxima, by_query['_query'].to_numpy(), by_query['_null_score'].to_numpy())

        mock = combine_components(unstratified(raw), unstratified(mock_norm) if 'pgen_bin' in mock_norm.columns else mock_norm, model, observed_queries,
            species=species, locus=locus, threads=threads, exclude_exact=exclude_exact, _consume=consume)
        # Internal adapters may still return one materialized frame.
        if mock is not None:
            consume(mock)
        indices = database_observed['_case'].to_numpy()
        count(indices, maxima)
        seen[indices] = True
        if not seen.all():
            raise RuntimeError('randomized scoring did not preserve the complete observed case grid')
        del raw, mock_norm, mock, consume, count, seen, maxima
    tails = observed.drop('_case').with_columns(pl.Series('reference_exceedances', exceed),
        pl.Series('query_exceedances', query_exceed), pl.lit(replicates).alias('reference_draws'),
        pl.lit(replicates).alias('query_draws'), pl.Series('reference_hit_draws', reference_hit_draws),
        pl.Series('query_hit_draws', query_hit_draws)).with_columns(
        (pl.col('reference_hit_draws') / replicates).alias('null_reference_match_rate'),
        (pl.col('query_hit_draws') / replicates).alias('null_query_match_rate'),
        ((pl.col('reference_exceedances') + 1) / (replicates + 1)).alias('p_reference'),
        ((pl.col('query_exceedances') + 1) / (replicates + 1)).alias('p_query')).with_columns(
        pl.max_horizontal('p_reference', 'p_query').alias('p_both'))
    return observed_detail, tails
