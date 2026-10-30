"""Complete-score annotation and two directional tails from one raw control."""
from pathlib import Path


def _validate(a):
    import math
    if not math.isfinite(a.fpr) or not 0 < a.fpr < 1:
        raise ValueError('fpr must be finite and between zero and one')
    expected = dict(distance_matrix='blosum62', junction_ends='trim3', stem_depth=4,
        stem_weight=1., v_loop_weight=1, v_regions='cdr12plus', alphabet='full',
        position_weighting='uniform', neighbour_weighting='kernel', kernel_scorer='v-bounded',
        background='vdjdb-other', rank_signal_prior=.5, knee_floor_factor=.5,
        chemistry_bandwidth_factor=1.)
    unsupported = [name.replace('_', '-') for name, value in expected.items() if getattr(a, name) != value]
    unsupported += [name.replace('_', '-') for name in ('geometry_extra', 'targets_from_sample',
        'withhold_epitope', 'reference_radius', 'unbounded_nearest', 'conditional_j',
        'epitope_rules', 'chemistry', 'nearest_rank', 'shortlist_audit') if getattr(a, name)]
    if a.radius_factor not in (None, 22.5):
        unsupported.append('radius-factor')
    if a.radius is not None:
        unsupported.append('radius')
    if unsupported:
        raise ValueError('one-control combined annotation does not support: ' + ', '.join(unsupported))
    if a.threads < 1 or a.mock_replicates < 1 or a.seed < 0 or a.min_score < 0:
        raise ValueError('threads/replicates must be positive; seed/min-score nonnegative')


def _queries(rows, a, organism):
    import polars as pl
    from ..io.airr import normalize_annotation_input
    from ..evalue.control import _organism

    parts, paired, offset = [], [], 0
    for sample_id, path, info in rows:
        raw, report = normalize_annotation_input(path)
        if any(c.startswith(('calibration_', '_combined_')) for c in raw.columns):
            raise ValueError('input contains reserved calibration identity/status columns')
        if ('sample_id' in raw.columns or set(info) & set(raw.columns) or
                any(c.startswith(('calibration_', '_combined_')) for c in info)):
            raise ValueError('sample metadata conflicts with repertoire columns')
        raw = raw.with_row_index('_combined_row', offset=offset).with_columns(
            pl.lit(sample_id).alias('sample_id'),
            *[pl.lit(v).alias(c) for c, v in info.items()])
        if 'species' not in raw.columns:
            raw = raw.with_columns(pl.lit(organism).alias('species'))
        offset += raw.height
        species_ok = pl.lit(True)
        for column in ('species', 'organism'):
            if column in raw.columns:
                aliases = {}
                for value in raw[column].drop_nulls().unique():
                    try:
                        aliases[value] = _organism(value) == organism
                    except ValueError:
                        aliases[value] = False
                species_ok &= pl.col(column).replace_strict(aliases, default=False, return_dtype=pl.Boolean)
        locus = pl.col('locus').str.strip_chars().str.to_uppercase().replace({'ALPHA': 'TRA', 'BETA': 'TRB'})
        if a.locus in ('TRA', 'TRB'):
            locus = locus.fill_null(a.locus)
        raw = raw.with_columns(locus.alias('calibration_locus'),
            pl.lit(organism).alias('calibration_species'),
            pl.col('junction_aa').str.strip_chars().str.to_uppercase().alias('_combined_sequence'),
            pl.concat_str([pl.lit('row:'), pl.col('_combined_row')]).alias('calibration_query_id'),
            pl.lit('single').alias('calibration_pair_status'))
        sequence = pl.col('_combined_sequence')
        locus_ok = pl.col('calibration_locus').is_in(['TRA', 'TRB'])
        if a.locus in ('TRA', 'TRB'):
            locus_ok &= pl.col('calibration_locus') == a.locus
        raw = raw.with_columns(pl.when(sequence.is_null() | sequence.eq('')).then(pl.lit('missing_sequence'))
            .when(~sequence.str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$')).then(pl.lit('invalid_sequence'))
            .when(~species_ok.fill_null(False)).then(pl.lit('species_mismatch'))
            .when(~locus_ok.fill_null(False)).then(pl.lit('locus_mismatch'))
            .otherwise(pl.lit('ready')).alias('calibration_status'))
        linked = a.locus == 'paired' or report['wide_input'] or ('pair_id' in raw.columns and
            raw.filter(pl.col('pair_id').is_not_null()).group_by('pair_id').agg(
                pl.col('calibration_locus').n_unique().alias('_loci')).filter(pl.col('_loci') > 1).height > 0)
        if linked and a.locus in ('auto', 'paired'):
            if 'pair_id' not in raw.columns:
                raise ValueError('paired annotation requires supplied pair/cell linkage')
            linked_rows = raw.filter(pl.col('pair_id').is_not_null() & pl.col('pair_id').str.strip_chars().ne(''))
            if linked_rows.group_by('pair_id', 'calibration_locus').len().filter(pl.col('len') > 1).height:
                raise ValueError('ambiguous paired input: multiple rows per pair/locus')
            ready = linked_rows.filter(pl.col('calibration_status') == 'ready')
            complete = ready.group_by('pair_id').agg(pl.col('calibration_locus').n_unique().alias('_loci'))
            complete = ready.join(complete.filter(pl.col('_loci') == 2).select('pair_id'), on='pair_id', how='semi')
            pairs = complete.group_by('pair_id', maintain_order=True).agg(
                pl.col('_combined_row').min().alias('_first'),
                *[pl.col(source).filter(pl.col('calibration_locus') == chain).first().alias(dest + suffix)
                  for chain, suffix in [('TRA', 'a'), ('TRB', 'b')]
                  for source, dest in [('_combined_sequence', 'cdr3'), ('v_call', 'v'), ('j_call', 'j')]])
            pairs = pairs.with_columns(pl.concat_str([pl.lit('pair:'), pl.col('_first')]).alias('query_id'))
            raw = raw.join(pairs.select('pair_id', pl.col('query_id').alias('_pair_query')), on='pair_id',
                how='left', validate='m:1', maintain_order='left').with_columns(
                    pl.coalesce('_pair_query', 'calibration_query_id').alias('calibration_query_id'),
                    pl.when(pl.col('_pair_query').is_not_null()).then(pl.lit('paired')).otherwise(
                        pl.col('calibration_locus')).alias('calibration_locus'),
                    pl.when(pl.col('_pair_query').is_not_null()).then(pl.lit('complete_pair'))
                    .when(pl.col('calibration_status') == 'ready').then(pl.lit('incomplete_pair_marginal'))
                    .otherwise(pl.lit('incomplete_pair')).alias('calibration_pair_status')).drop('_pair_query')
            if pairs.height:
                paired.append(pairs.drop('pair_id', '_first'))
        parts.append(raw)
    raw = pl.concat(parts, how='diagonal_relaxed')
    populations = [(chain, raw.filter((pl.col('calibration_locus') == chain) &
        (pl.col('calibration_status') == 'ready')).select(
            pl.col('calibration_query_id').alias('query_id'), pl.col('_combined_sequence').alias('cdr3'),
            pl.col('v_call').alias('v'), pl.col('j_call').alias('j'), pl.lit(chain).alias('locus')))
        for chain in ('TRA', 'TRB')]
    if paired:
        populations.append(('paired', pl.concat(paired)))
    return raw.drop('_combined_sequence').rename({'_combined_row': 'calibration_input_row'}), populations


def _raw_control(path, *, retain_missing_genes=False):
    import polars as pl
    from ..io.airr import _read_table
    from ..io.columns import validate_counts, _resolve

    raw = _read_table(path)
    validate_counts(raw)
    source_rows = raw.height
    productive = next((c for c in raw.columns if c.lower() == 'productive'), None)
    if productive:
        values = pl.col(productive).cast(pl.String).str.strip_chars().str.to_lowercase()
        if raw.filter(values.is_not_null() & ~values.is_in(['t', 'true', '1', 'f', 'false', '0'])).height:
            raise ValueError('control productive must contain AIRR true/false values')
        raw = raw.filter(values.is_in(['t', 'true', '1']).fill_null(False))
    productive_rows = raw.height
    lower = {c.lower(): c for c in raw.columns}
    sequence = next((lower[c] for c in ('junction_aa', 'cdr3aa', 'cdr3', 'cdr3_aa') if c in lower), None)
    if sequence is None:
        raise ValueError('control requires a full junction_aa or legacy cdr3aa column')
    if sequence.lower() == 'cdr3_aa' or (sequence.lower() == 'cdr3' and
            any(c in lower for c in ('v_call', 'j_call', 'sequence_id'))):
        raise ValueError('AIRR control requires junction_aa; bare CDR3 excludes anchors')
    raw = raw.with_columns(pl.col(sequence).cast(pl.String).str.strip_chars().str.to_uppercase())
    raw = raw.filter(pl.col(sequence).str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False))
    sequence_rows = raw.height
    fields = _resolve(raw)
    if retain_missing_genes:
        return raw, dict(control_source_rows=source_rows,control_productive_rows=productive_rows,control_excluded_nonproductive_rows=source_rows-productive_rows,control_excluded_invalid_sequence_rows=productive_rows-sequence_rows,control_excluded_invalid_gene_rows=0,control_valid_rows=raw.height,missing_gene_policy='retained for flagged control-event imputation')
    if not {'v', 'j'}.issubset(fields):
        raise ValueError('combined-score control requires V and J gene columns')
    v = pl.col(fields['v']).cast(pl.String).str.strip_chars().str.to_uppercase()
    j = pl.col(fields['j']).cast(pl.String).str.strip_chars().str.to_uppercase()
    valid_genes = (v.str.contains(r'^TR[AB]V[A-Z0-9/.-]+(?:\*[A-Z0-9]+)?$') &
        j.str.contains(r'^TR[AB]J[A-Z0-9/.-]+(?:\*[A-Z0-9]+)?$') &
        (v.str.slice(0, 3) == j.str.slice(0, 3)))
    if 'locus' in fields:
        declared = pl.col(fields['locus']).cast(pl.String).str.strip_chars().str.to_uppercase().replace(
            {'ALPHA': 'TRA', 'BETA': 'TRB'})
        valid_genes &= declared == v.str.slice(0, 3)
    raw = raw.filter(valid_genes.fill_null(False))
    return raw, dict(control_source_rows=source_rows, control_productive_rows=productive_rows,
        control_excluded_nonproductive_rows=source_rows - productive_rows,
        control_excluded_invalid_sequence_rows=productive_rows - sequence_rows,
        control_excluded_invalid_gene_rows=sequence_rows - raw.height,
        control_valid_rows=raw.height)


def _detail_targets(raw,receptors,reference):
    """Expand declarations to compatible full reference groups, without pooling models."""
    import polars as pl
    from ..aggregate.candidates import PMHC
    from ..db.schema import mhc_compatible
    schema={c:pl.String for c in ['query_id',*PMHC]}
    if 'epitope' not in raw.columns:
        return pl.DataFrame(schema=schema)
    declared=raw.rename({'calibration_query_id':'_declared_query'}).join(
        receptors.select(pl.col('query_id').alias('_declared_query')),on='_declared_query',how='semi')
    declared=declared.with_columns(*[pl.lit(None,dtype=pl.String).alias(c) for c in PMHC if c not in declared.columns])
    declared=declared.select(pl.col('_declared_query').alias('query_id'),pl.col(PMHC).cast(pl.String)).filter(
        pl.col('epitope').is_not_null()&pl.col('epitope').str.strip_chars().ne('')).unique()
    catalog=reference.select(PMHC).unique()
    parts=[]
    for key,part in declared.group_by(PMHC,maintain_order=True):
        if any(value is None or not value.strip() for value in key):
            selected=catalog.head(0)
        else:
            epitope,a,b,kind=key
            selected=catalog.filter((pl.col('epitope')==epitope)&(pl.col('mhc_class')==kind)&
                mhc_compatible('mhc_a',a)&mhc_compatible('mhc_b',b))
        parts.append(part.select('query_id').join(selected,how='cross') if selected.height else part)
    return pl.concat(parts).unique() if parts else pl.DataFrame(schema=schema)


def main(a):
    import hashlib
    from importlib.metadata import version
    from ..precursor.model import bundled_model_path
    import time
    import polars as pl
    import seqtree._core as native
    from vdjtools import _core as pgen_native
    from .. import db
    from ..aggregate.candidates import PMHC
    from ..db.cache import sha256
    from ..evalue.combined import calibrate_combined
    from ..evalue.control import _organism
    from ..evalue.mock_reference import control_events, EVENT_SAMPLER
    from ..runner.annotation import _input_rows, _validate_output_paths
    from .__main__ import _flat_table

    start = time.monotonic()
    _validate(a)
    organism = _organism(a.species)
    rows = _input_rows(a.sample, a.metadata)
    control_path = Path(a.control).expanduser().resolve()
    for _, path, _ in rows:
        if path == control_path or path.samefile(control_path):
            raise ValueError('query sample cannot also define its control population')
    reference_path = Path(a.vdjdb).expanduser() if a.vdjdb else db.fetch_latest()
    protected = [path for _, path, _ in rows] + [reference_path, control_path]
    if a.metadata:
        protected.append(Path(a.metadata))
    prefix = _validate_output_paths(a.output_prefix, protected)
    reference = db.load(reference_path, epitope=a.epitope, mhc_a=a.mhc_a, mhc_b=a.mhc_b,
        mhc_match=a.mhc_match, evidence_type=a.evidence_type, min_score=a.min_score)
    species_name = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[organism]
    species_policy = 'caller_declared' if reference['species'].null_count() == reference.height else 'source_annotated'
    if species_policy == 'caller_declared':
        reference = reference.with_columns(pl.lit(species_name).alias('species'))
    else:
        aliases = {}
        for value in reference['species'].drop_nulls().unique():
            try:
                aliases[value] = _organism(value) == organism
            except ValueError:
                aliases[value] = False
        source_reference_rows = reference.height
        reference = reference.filter(pl.col('species').replace_strict(aliases,
            default=False, return_dtype=pl.Boolean)).with_columns(pl.lit(species_name).alias('species'))
    excluded_reference_species_rows = 0 if species_policy == 'caller_declared' else source_reference_rows - reference.height
    valid = pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
    if a.reference_population == 'all-valid':
        reference = reference.with_columns(valid.alias('reference_valid'))
    reference = reference.filter(valid & pl.col('reference_valid'))
    if a.locus in ('TRA', 'TRB'):
        reference = reference.filter(pl.col('gene') == a.locus)
    if not reference.height:
        raise ValueError('selected reference has no valid matching-species/locus records')
    query, populations = _queries(rows, a, organism)
    control, control_counts = _raw_control(control_path, retain_missing_genes=True)
    chains_needed = set()
    for locus, receptors in populations:
        chains = ['TRA', 'TRB'] if locus == 'paired' else [locus]
        if receptors.height and reference.filter(pl.col('gene').is_in(chains)).height:
            chains_needed.update(chains)
    events_by_chain = {chain: control_events(control, species=organism, locus=chain) for chain in sorted(chains_needed)}
    details, tails, metadata = [], [], []
    for locus, receptors in populations:
        if not receptors.height:
            continue
        chains = ['TRA', 'TRB'] if locus == 'paired' else [locus]
        if not reference.filter(pl.col('gene').is_in(chains)).height:
            query = query.with_columns(pl.when(pl.col('calibration_query_id').is_in(receptors['query_id'].implode()))
                .then(pl.lit('no_reference')).otherwise(pl.col('calibration_status')).alias('calibration_status'))
            continue
        events = {chain: events_by_chain[chain] for chain in chains}
        detail, tail = calibrate_combined(receptors, reference, control, species=organism, locus=locus,
            replicates=a.mock_replicates, seed=a.seed, threads=a.threads, exclude_exact=a.exclude_exact,
            report_p=a.fpr,
            detail_targets=_detail_targets(query,receptors,reference.filter(pl.col('gene').is_in(chains)))
                if a.selection_only else None)
        details.append(detail.with_columns(pl.lit(organism).alias('species'), pl.lit(locus).alias('locus')))
        tails.append(tail.with_columns(pl.lit(organism).alias('species'), pl.lit(locus).alias('locus')))
        source = 'olga' if organism == 'human' else 'arda'
        digest = hashlib.sha256()
        for chain in chains:
            model = bundled_model_path(chain, source, organism)
            for path in sorted(model.iterdir(), key=lambda p: p.name):
                if path.is_file():
                    digest.update(chain.encode() + path.name.encode() + b'\0' + hashlib.sha256(path.read_bytes()).digest())
        metadata.append(dict(species=organism, locus=locus, query_receptors=receptors.height,
            reference_records=reference.filter(pl.col('gene').is_in(chains)).height,
            reference_population=a.reference_population, reference_species_policy=species_policy,
            mock_replicates=a.mock_replicates, seed=a.seed, native_threads=a.threads,
            table_threads=pl.thread_pool_size(), **control_counts,
            reference_excluded_species_rows=excluded_reference_species_rows,
            control_events=sum(int(x['event_weight'].sum()) for x in events.values()),
            control_keys=sum(x.height for x in events.values()), exclude_exact=a.exclude_exact,
            reference_sha256=sha256(reference_path) if reference_path.is_file() else None,
            control_sha256=sha256(control_path), min_score=a.min_score,
            evidence_type_filter=','.join(a.evidence_type or []), epitope_filter=','.join(a.epitope or []),
            mhc_a_filter=a.mhc_a, mhc_b_filter=a.mhc_b, mhc_match=a.mhc_match,
            original_model_context='matching-chain control prior' if locus != 'paired' else 'beta-control prior for both paired chains',
            context_policy='candidate one-control recipe; distinct from historical benchmark model assets',
            null_pairing='independently randomized TRA/TRB events' if locus == 'paired' else 'single-chain events',
            normalization='independent unconditional bank, fixed across both directions; final query nulls Pgen-matched',
            selection='database maximum over full pMHC and available paired axes',
            detail_scope='declared targets plus reported specificities' if a.selection_only else 'all query-pMHC pairs',
            report_p=a.fpr,
            reported_tail='both directional complete-label maxima at each reported label score',
            interpretation='conditional randomization tails, not binding probability; k/B match rates are separate',
            pgen_model_source=source, pgen_model_sha256=digest.hexdigest(),
            seqtree_native_sha256=sha256(Path(native.__file__)), pgen_native_sha256=sha256(Path(pgen_native.__file__)),
            vdjmatch_version=version('vdjmatch'), seqtree_version=version('seqtree'), vdjtools_version=version('vdjtools')))
    base = query.select('calibration_query_id', 'sample_id', 'calibration_status', 'calibration_pair_status',
        pl.col('calibration_locus').alias('locus')).unique(maintain_order=True).with_columns(pl.lit(organism).alias('species'))
    keys = ['query_id', *PMHC]
    if tails:
        tails = pl.concat(tails, how='diagonal_relaxed')
        database = base.join(tails.filter(pl.col('test_scope') == 'database').drop(*PMHC).rename(
            {'query_id': 'calibration_query_id'}), on=['calibration_query_id', 'species', 'locus'],
            how='left', validate='1:1', maintain_order='left')
        detail = pl.concat(details, how='diagonal_relaxed')
        epitopes = detail.join(tails.filter(pl.col('test_scope') == 'epitope').drop('score'),
            on=[*keys, 'species', 'locus'], how='left', nulls_equal=True, validate='1:1').rename(
                {'query_id': 'calibration_query_id'})
        epitopes = epitopes.join(base.select('calibration_query_id', 'sample_id'), on='calibration_query_id', validate='m:1')
    else:
        database = base.with_columns(pl.lit(None, dtype=pl.Float64).alias('score'),
            *[pl.lit(None, dtype=pl.Float64).alias(c) for c in ('p_reference', 'p_query', 'p_both')])
        epitopes = pl.DataFrame(schema={**{c: pl.String for c in ['calibration_query_id', 'species', 'locus', *PMHC]},
            'score': pl.Float64, 'call_score': pl.Float64})
    database = database.with_columns(pl.when(pl.col('calibration_status') == 'ready').then(
        pl.when(pl.col('score').is_not_null()).then(pl.lit('scored')).otherwise(pl.lit('no_support')))
        .otherwise(pl.col('calibration_status')).alias('calibration_status'))
    query = query.drop('calibration_status').join(database.select('calibration_query_id', 'calibration_status'),
        on='calibration_query_id', how='left', validate='m:1', maintain_order='left')
    calibration = pl.DataFrame(metadata) if metadata else pl.DataFrame({'species': [organism], 'locus': [a.locus],
        'query_receptors': [0], 'mock_replicates': [a.mock_replicates], **{k: [v] for k,v in control_counts.items()}})
    calibration = calibration.with_columns(pl.lit(time.monotonic() - start).alias('calculation_seconds'),
        pl.lit(a.fpr).alias('report_p'),
        pl.lit(EVENT_SAMPLER).alias('reference_event_sampler'))
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for name, frame in [('queries', query), ('database', database), ('epitopes', epitopes), ('calibration', calibration)]:
        _flat_table(frame).write_csv(str(prefix) + '.' + name + '.tsv', separator='\t')
    return 0
