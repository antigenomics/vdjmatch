"""Fresh, scope-specific directional junction tails from one raw control repertoire."""
from pathlib import Path

import numpy as np
import polars as pl

from ..evalue.mock_reference import normalize_control, compute_pgen, pgen_strata, iter_mock_references, EVENT_SAMPLER

_PMHC = ['species', 'gene', 'epitope', 'mhc_a', 'mhc_b', 'mhc_class']
_LIMIT = 5_000_000  # Bound sparse hit/membership expansions, not a dense query/control matrix.


def _scores(hits, membership, statistic):
    membership = membership.unique()
    expansion = hits.group_by('junction_aa').len().join(
        membership.group_by('junction_aa').len(name='groups'), on='junction_aa').select(
            (pl.col('len').cast(pl.Int64) * pl.col('groups')).sum()).item()
    if expansion > _LIMIT:
        raise ValueError('scope produces too many hit/group memberships; reduce scope or reference population')
    joined = hits.join(membership, on='junction_aa', how='inner')
    return joined.group_by('qindex', 'group_id').agg(
        pl.col('distance').min().alias('nearest_edits'),
        pl.col('junction_aa').n_unique().alias('hit_count')).with_columns(
            (-pl.col('nearest_edits').cast(pl.Float64) if statistic == 'nearest' else
             pl.col('hit_count').cast(pl.Float64)).alias('score'))


def main(a):
    import hashlib
    from ..precursor.model import bundled_model_path
    import time
    from importlib.metadata import version
    import seqtree._core as seqtree_native
    from vdjtools import _core as pgen_native
    from seqtree import Index
    from .. import db
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..io.airr import _read_table
    from ..match.scope import search_params
    from ..runner.annotation import _validate_output_paths
    from .database import prepare, exact_clonotypes
    from .__main__ import _flat_table

    start = time.monotonic()
    if a.threads < 1 or a.mock_replicates < 1 or a.max_edits < 0:
        raise ValueError('threads/replicates must be positive and max-edits nonnegative')
    if a.statistic not in {'nearest', 'count'}:
        raise ValueError('statistic must be nearest or count')
    sample_path, control_path = Path(a.sample), Path(a.control)
    if sample_path.resolve() == control_path.resolve() or (
            sample_path.exists() and control_path.exists() and sample_path.samefile(control_path)):
        raise ValueError('query sample cannot also define its control population')
    prefix = _validate_output_paths(a.output_prefix, [a.sample, a.vdjdb, a.control])
    species = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[_organism(a.species)]
    model_source = 'olga' if _organism(a.species) == 'human' else 'arda'
    model_root = bundled_model_path(a.locus, model_source, _organism(a.species))
    model_digest = hashlib.sha256()
    for path in sorted(model_root.iterdir(), key=lambda p: p.name):
        if path.is_file():
            model_digest.update(path.name.encode() + b'\0' + hashlib.sha256(path.read_bytes()).digest())
    reference = db.load(a.vdjdb, species=species, gene=a.locus,
        epitope=getattr(a, 'epitope', None), evidence_type=getattr(a, 'evidence_type', None),
        min_score=getattr(a, 'min_score', 0))
    reference = reference.filter(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False))
    if not reference.height:
        raise ValueError('selected reference has no valid junctions')
    query, _ = prepare(a.sample, a)
    reserved = {'p_reference', 'p_query', 'p_both', 'reference_exceedances', 'query_exceedances',
        'reference_hit_draws', 'query_hit_draws', 'null_reference_match_rate', 'null_query_match_rate',
        'reference_draws', 'query_draws'}
    if reserved.intersection(query.columns) or {'pgen', 'source_pgen'}.issubset(query.columns):
        raise ValueError('sample contains reserved calibration output columns')
    if 'pgen' in query.columns:
        query = query.rename({'pgen': 'source_pgen'})
    control_raw = _read_table(a.control)
    control = normalize_control(control_raw, species=a.species, locus=a.locus)
    for column in _PMHC:
        if column not in reference.columns:
            reference = reference.with_columns(pl.lit(None, dtype=pl.String).alias(column))
    groups = reference.select(_PMHC).unique().sort(_PMHC).with_row_index('group_id', offset=1)
    real_membership = reference.join(groups, on=_PMHC, nulls_equal=True, validate='m:1').select(
        'group_id', pl.col('cdr3').alias('junction_aa'))
    membership = pl.concat([reference.select(pl.lit(0, dtype=pl.UInt32).alias('group_id'),
        pl.col('cdr3').alias('junction_aa')), real_membership])
    if membership.height > _LIMIT:
        raise ValueError('reference membership budget exceeded')
    groups = pl.concat([pl.DataFrame({'group_id': [0]}, schema={'group_id': pl.UInt32}).with_columns(
        *[pl.lit(None, dtype=reference.schema[c]).alias(c) for c in _PMHC]), groups])
    valid = query.filter(pl.col('database_status') == 'ready').select(
        pl.col('database_junction').alias('junction_aa')).unique().sort('junction_aa')
    # One fresh Pgen batch includes query-only junctions as well as reference/control.
    pgen = compute_pgen(pl.concat([membership.select('junction_aa'), valid]), control,
        species=a.species, locus=a.locus, threads=a.threads)
    strata = pgen_strata(pgen)
    valid = valid.join(strata, on='junction_aa', validate='1:1')
    pool = control.join(strata, on='junction_aa', validate='1:1').sort('pgen_bin', 'junction_aa')
    rng = np.random.default_rng(np.random.SeedSequence(a.seed).spawn(2)[1])
    draws = []
    for key, part in pool.partition_by('pgen_bin', as_dict=True).items():
        if key[0] not in valid['pgen_bin']:
            continue
        weight = part['event_weight'].to_numpy().astype(float)
        selected = rng.choice(part['junction_aa'].to_numpy(), size=a.mock_replicates,
                              replace=True, p=weight / weight.sum())
        draws.append(pl.DataFrame({'pgen_bin': key[0], 'draw': np.arange(a.mock_replicates), 'junction_aa': selected}))
    if set(valid['pgen_bin']) - set(pool['pgen_bin']):
        raise ValueError('control has no junctions in a query Pgen stratum')
    draw_table = pl.concat(draws) if draws else pl.DataFrame(schema={
        'pgen_bin': pl.Int64, 'draw': pl.Int64, 'junction_aa': pl.String})
    searched = pl.concat([valid.select('junction_aa'), draw_table.select('junction_aa')]).unique().sort('junction_aa').with_row_index('qindex')
    targets = pl.concat([membership.select('junction_aa'), control.select('junction_aa')]).unique().sort('junction_aa')['junction_aa'].to_list()
    scope = getattr(a, 'scope', None) or ','.join([str(a.max_edits)] * 4)
    # seqtrie enforces total edits but ignores individual operation budgets.
    engine = 'seqtm' if getattr(a, 'scope', None) else 'seqtrie'
    params = search_params(scope, engine=engine, mode='all')
    load_seconds = time.monotonic() - start
    search_start = time.monotonic()
    native_hits = Index.build(targets, 'aa').search_batch(searched['junction_aa'].to_list(), params, a.threads)
    if sum(map(len, native_hits)) > _LIMIT:
        raise ValueError('scope produces too many sparse hits; reduce scope or population')
    hits = pl.DataFrame(((qi, targets[h.ref_id], h.score) for qi, rows in enumerate(native_hits) for h in rows),
        schema={'qindex': pl.UInt32, 'junction_aa': pl.String, 'distance': pl.Int64}, orient='row')
    del native_hits
    search_seconds = time.monotonic() - search_start
    score_start = time.monotonic()
    real_scores = _scores(hits, membership, a.statistic)
    valid = valid.join(searched, on='junction_aa', validate='1:1')
    observed = real_scores.join(valid.select('qindex', 'pgen_bin'), on='qindex', how='inner')
    # Only database no-hits are expanded. Missing epitope rows mean no-hit, P=1.
    missing_db = valid.select('qindex', 'pgen_bin').join(
        observed.filter(pl.col('group_id') == 0).select('qindex'), on='qindex', how='anti').with_columns(
            pl.lit(0, dtype=pl.UInt32).alias('group_id'),
            pl.lit(None, dtype=pl.Int64).alias('nearest_edits'),
            pl.lit(0, dtype=pl.UInt32).alias('hit_count'), pl.lit(None, dtype=pl.Float64).alias('score'))
    observed = pl.concat([observed, missing_db.select(observed.columns)]).sort('qindex', 'group_id').with_row_index('_case')
    if observed.height > _LIMIT:
        raise ValueError('observed candidate budget exceeded')
    observed_score_seconds = time.monotonic() - score_start
    ref_start = time.monotonic()
    exceed = np.zeros(observed.height, dtype=np.int64)
    reference_hit_draws = np.zeros(observed.height, dtype=np.int64)
    ref_seed = int(np.random.SeedSequence(a.seed).spawn(2)[0].generate_state(1)[0])
    # Iterator preflights exact source-row strata, even for a sample with no valid rows.
    observed_hits = hits.join(valid.select('qindex'), on='qindex', how='semi')
    for mock in iter_mock_references(membership, control, replicates=a.mock_replicates,
                                    seed=ref_seed, pgen=pgen, species=a.species, locus=a.locus, threads=a.threads):
        scored = _scores(observed_hits, mock.select('group_id', 'junction_aa'), a.statistic)
        comparison = observed.select('qindex', 'group_id', 'score').join(
            scored.select('qindex', 'group_id', pl.col('score').alias('_null_score')),
            on=['qindex', 'group_id'], how='left', validate='1:1', maintain_order='left')
        reference_hit_draws += comparison['_null_score'].is_not_null().cast(pl.Int64).to_numpy()
        exceed += comparison.select((pl.col('score').is_null() |
            (pl.col('_null_score') >= pl.col('score')).fill_null(False)).cast(pl.Int64)).to_series().to_numpy()
    ref_seconds = time.monotonic() - ref_start
    query_start = time.monotonic()
    indexed_draws = draw_table.join(searched, on='junction_aa', validate='m:1')
    expansion = indexed_draws.group_by('qindex').len().join(real_scores.group_by('qindex').len(name='groups'),
        on='qindex').select((pl.col('len').cast(pl.Int64) * pl.col('groups')).sum()).item()
    if expansion > _LIMIT:
        raise ValueError('query null score membership budget exceeded')
    null_query_scores = indexed_draws.join(
        real_scores.select('qindex', 'group_id', 'score'), on='qindex', how='inner')
    query_exceed = np.zeros(observed.height, dtype=np.int64)
    query_hit_draws = np.zeros(observed.height, dtype=np.int64)
    distributions = {key: np.sort(part['score'].to_numpy()) for key, part in
                     null_query_scores.partition_by('pgen_bin', 'group_id', as_dict=True).items()}
    for key, part in observed.partition_by('pgen_bin', 'group_id', as_dict=True).items():
        values = distributions.get(key, np.empty(0))
        scores = part['score'].to_numpy()
        count = len(values) - np.searchsorted(values, scores, side='left')
        count[np.isnan(scores)] = a.mock_replicates
        query_exceed[part['_case'].to_numpy()] = count
        query_hit_draws[part['_case'].to_numpy()] = len(values)
    observed = observed.with_columns(pl.Series('reference_hit_draws', reference_hit_draws),
        pl.Series('query_hit_draws', query_hit_draws),
        pl.Series('reference_exceedances', exceed),
        pl.Series('query_exceedances', query_exceed), pl.lit(a.mock_replicates).alias('reference_draws'),
        pl.lit(a.mock_replicates).alias('query_draws')).with_columns(
            (pl.col('reference_hit_draws') / a.mock_replicates).alias('null_reference_match_rate'),
            (pl.col('query_hit_draws') / a.mock_replicates).alias('null_query_match_rate'),
            ((pl.col('reference_exceedances') + 1) / (a.mock_replicates + 1)).alias('p_reference'),
            ((pl.col('query_exceedances') + 1) / (a.mock_replicates + 1)).alias('p_query')).with_columns(
                pl.max_horizontal('p_reference', 'p_query').alias('p_both'))
    query_seconds = time.monotonic() - query_start
    lookup = observed.filter(pl.col('group_id') == 0).drop('group_id', 'pgen_bin', '_case').join(
        searched, on='qindex', validate='1:1').drop('qindex').rename({
            'junction_aa': 'database_junction', 'nearest_edits': 'database_nearest_edits',
            'hit_count': 'database_hit_count', 'score': 'database_score'})
    query = query.join(lookup, on='database_junction', how='left', validate='m:1', maintain_order='left').join(
        pgen.rename({'junction_aa': 'database_junction'}), on='database_junction', how='left', validate='m:1', maintain_order='left')
    output_columns = ['database_nearest_edits', 'database_hit_count', 'database_score', 'pgen',
        'reference_hit_draws', 'query_hit_draws', 'null_reference_match_rate', 'null_query_match_rate',
        'reference_exceedances', 'query_exceedances', 'reference_draws', 'query_draws', 'p_reference', 'p_query', 'p_both']
    query = query.with_columns(*[pl.when(pl.col('database_status') == 'ready').then(pl.col(c)).otherwise(None).alias(c)
        for c in output_columns], pl.when(pl.col('database_status') == 'ready').then(
            pl.col('database_junction').is_in(reference['cdr3'].unique().to_list())).otherwise(None).alias('database_exact_junction'))
    query = query.with_columns(pl.when(pl.col('database_status') == 'ready').then(
        pl.when(pl.col('database_nearest_edits').is_not_null()).then(pl.lit('matched')).otherwise(pl.lit('no_hit')))
        .otherwise(pl.col('database_status')).alias('database_status'))
    query = exact_clonotypes(query, reference, a.locus)
    candidates = observed.filter(pl.col('group_id') != 0).drop('_case', 'pgen_bin').join(
        searched, on='qindex', validate='m:1').join(groups, on='group_id', validate='m:1').join(
            query.filter(pl.col('database_status') == 'matched').select('_database_row', 'database_junction'),
            left_on='junction_aa', right_on='database_junction', how='inner').drop('qindex').rename({'_database_row': 'database_input_row'})
    bin_counts = membership.join(strata, on='junction_aa', validate='m:1').group_by('group_id', 'pgen_bin').agg(
        pl.len().alias('pgen_bin_reference_source_records'))
    control_bin_counts = pool.group_by('pgen_bin').agg(
        pl.len().alias('pgen_bin_control_unique_junctions'),
        pl.col('event_weight').sum().alias('pgen_bin_control_events'))
    calibration = groups.join(membership.group_by('group_id').agg(pl.len().alias('reference_source_records'),
        pl.col('junction_aa').n_unique().alias('reference_unique_junctions')), on='group_id', validate='1:1').with_columns(
        pl.lit(a.mock_replicates).alias('mock_replicates'), pl.lit(a.seed).alias('seed'),
        pl.lit(EVENT_SAMPLER).alias('reference_event_sampler'),
        pl.lit(a.species).alias('species_requested'), pl.lit(a.locus).alias('locus_requested'),
        pl.lit(getattr(a, 'min_score', 0)).alias('min_score'),
        pl.lit(str(getattr(a, 'evidence_type', None))).alias('evidence_filter'),
        pl.lit(str(getattr(a, 'epitope', None))).alias('epitope_filter'),
        pl.lit(sha256(Path(a.vdjdb)) if Path(a.vdjdb).is_file() else None, dtype=pl.String).alias('reference_sha256'),
        pl.lit(a.statistic).alias('statistic'), pl.lit(scope).alias('scope'), pl.lit(a.threads).alias('native_threads'),
        pl.lit(pl.thread_pool_size()).alias('table_threads'), pl.lit(1).alias('native_batches'),
        pl.lit(control_raw.height).alias('control_source_rows'), pl.lit(control.height).alias('control_unique_junctions'),
        pl.lit(control['event_weight'].sum()).alias('control_events'), pl.lit(searched.height).alias('searched_unique_junctions'),
        pl.lit(hits.height).alias('sparse_hits'), pl.lit(load_seconds).alias('load_pgen_seconds'),
        pl.lit(search_seconds).alias('search_seconds'), pl.lit(observed_score_seconds).alias('observed_score_seconds'),
        pl.lit(ref_seconds).alias('reference_tail_seconds'), pl.lit(query_seconds).alias('query_tail_seconds'),
        pl.lit(time.monotonic() - start).alias('calculation_seconds'),
        pl.lit(sha256(Path(a.sample))).alias('sample_sha256'), pl.lit(sha256(Path(a.control))).alias('control_sha256'),
        pl.lit('supported epitope candidates only; omitted valid query/group = no-hit, both P=1').alias('epitope_row_policy'),
        pl.lit('fixed group tails; selected epitope requires its own selection event; match rates k/B are Monte Carlo estimates').alias('tail_scope'),
        pl.lit('full-junction scope/statistic only; no V-loop, combined-score, binding or posterior calibration').alias('interpretation'))
    calibration = calibration.join(bin_counts, on='group_id', validate='1:m').join(
        control_bin_counts, on='pgen_bin', validate='m:1').with_columns(
            pl.lit(engine).alias('search_engine'),
            pl.lit(model_source).alias('pgen_model_source'),
            pl.lit(model_digest.hexdigest()).alias('pgen_model_sha256'),
            pl.lit(False).alias('pgen_model_collapse_alleles'),
            pl.lit('one fresh native batch over reference/control/query union').alias('pgen_policy'),
            pl.lit(version('vdjtools')).alias('vdjtools_version'),
            pl.lit(version('seqtree')).alias('seqtree_version'),
            pl.lit(version('vdjmatch')).alias('vdjmatch_version'),
            pl.lit(sha256(Path(seqtree_native.__file__))).alias('seqtree_native_sha256'),
            pl.lit(sha256(Path(pgen_native.__file__))).alias('pgen_native_sha256'),
            pl.lit('post-search hit count and preflight sparse-join checks; native hit allocation is not bounded').alias('allocation_guard'))
    prefix.parent.mkdir(parents=True, exist_ok=True)
    _flat_table(query.rename({'_database_row': 'database_input_row'})).write_csv(str(prefix) + '.database.tsv', separator='\t')
    _flat_table(candidates.sort('database_input_row', 'group_id')).write_csv(str(prefix) + '.epitopes.tsv', separator='\t')
    _flat_table(calibration.sort('group_id', 'pgen_bin')).write_csv(str(prefix) + '.calibration.tsv', separator='\t')
    return 0
