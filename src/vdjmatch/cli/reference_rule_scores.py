"""Fresh fixed-model reference-rule retrieval, with complete continuous exports."""
from pathlib import Path


def register(subparsers):
    parser = subparsers.add_parser('reference-rule-scores', help=__doc__)
    parser.add_argument('--vdjdb', required=True)
    parser.add_argument('--species', choices=['human', 'mouse'], default='human')
    parser.add_argument('--locus', choices=['TRA', 'TRB', 'paired'], required=True)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--min-records', type=int, default=30)
    parser.add_argument('--reference-population', choices=['curated-benchmark', 'all-valid'], default='curated-benchmark')
    parser.add_argument('--class-policy',help='frozen reference-policy.tsv and sibling manifest; export only the experimental class-policy scorer')
    parser.add_argument('--kernel-only',action='store_true',help='export only the six current real-control-free stem4 kernel alternatives; do not rerun fixed rules')
    parser.add_argument('--v-model-policy', choices=['extended', 'resolved-first'], default='extended',
        help='kernel-only: preserve the original resolved-V model for resolved queries; extended remains the default')
    parser.add_argument('--components', action='store_true', help='also score the six retained single-chain components')
    parser.add_argument('--control', help='matching raw density control; default is bundled human TRB')
    parser.add_argument('--germline-background', help='explicit raw germline prior background; required with --components')
    parser.add_argument('--output-prefix', required=True)
    parser.set_defaults(func=main)


def main(a):
    import json
    import time
    from importlib.metadata import version
    import numpy as np
    import polars as pl
    import seqtree._core as native
    from .. import db
    from ..db.cache import sha256
    from ..evalue.epitope_rules import MOTIFS
    from ..evalue.reference_rule_bank import stream_reference_rules
    from ..match.geometry import precompute_v_loop_distances, impute_reference_v_alleles, resolve_v_alleles
    from ..runner.annotation import _validate_output_paths
    from ..runner.retrieval_winners import WinnerAccumulator
    from .epitope import pool_epitopes

    started = time.monotonic()
    if a.threads < 1 or a.min_records < 1:
        raise ValueError('positive threads and min-records required')
    class_policy = getattr(a,'class_policy',None)
    v_model_policy = getattr(a, 'v_model_policy', 'extended')
    if v_model_policy != 'extended' and not a.kernel_only:
        raise ValueError('--v-model-policy resolved-first requires --kernel-only')
    if class_policy and (a.components or a.kernel_only):
        raise ValueError('--class-policy, --kernel-only and --components are separate calculations')
    policy_paths = {}
    policy_frame = policy_info = None
    if class_policy:
        from .score_policy import FAST_ROUTES
        model = Path(class_policy).expanduser().resolve()
        companion = model.parent/'manifest.json'
        policy_info = json.loads(companion.read_text())
        if (policy_info.get('reference_model_asset') is not True or policy_info.get('no_benchmark_inputs') is not True
            or policy_info.get('outputs_sha256',{}).get(model.name)!=sha256(model)
            or policy_info.get('routes')!=FAST_ROUTES):
            raise ValueError('class policy must match a frozen reference-only FAST_ROUTES manifest')
        policy_frame = pl.read_csv(model,separator='\t')
        if not set(policy_frame['reference_class']) <= set(FAST_ROUTES):
            raise ValueError('class-policy scoring requires existing FAST_ROUTES classes')
        policy_paths = {'policy':model,'manifest':companion}
    def policy_freeze():
        return {name:dict(path=str(path),sha256=sha256(path)) for name,path in policy_paths.items()}
    policy_before = policy_freeze()
    paired=a.locus=='paired'
    if paired and a.components:
        raise ValueError('historical --components supports single chains only')
    if a.kernel_only and a.components:
        raise ValueError('--kernel-only and historical --components are separate calculations')
    if a.components and not a.germline_background:
        raise ValueError('--components requires --germline-background')
    if not a.components and (a.control or a.germline_background):
        raise ValueError('--control and --germline-background require --components')
    prefix = _validate_output_paths(a.output_prefix, [a.vdjdb, *[
        path for path in (a.control, a.germline_background,*policy_paths.values()) if path]])
    suffixes = ['queries.tsv', 'truth.tsv', 'scores.tsv', 'target-support.tsv', 'manifest.json', 'continuous'] + (['policy-audit.tsv'] if class_policy else [])
    destinations = [Path(str(prefix)+'.'+suffix) for suffix in suffixes]
    if any(path.exists() or path.is_symlink() for path in destinations):
        raise FileExistsError('fresh reference-rule output paths must not already exist')
    organism, species = a.species, {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[a.species]
    def freeze():
        source_paths = sorted(Path(__file__).resolve().parents[1].rglob('*.py')) + [
            Path(native.__file__), Path(__file__).parents[1]/'resources/vloops'/f'{organism}_v_loops.tsv']
        return {str(path): sha256(path) for path in source_paths}
    component_paths = {}
    if a.components:
        from importlib.resources import files
        if not a.control and (organism != 'human' or a.locus != 'TRB'):
            raise ValueError('this control is not bundled; supply --control for the matching species/locus')
        component_paths = {'control': (Path(a.control).expanduser().resolve() if a.control else
            Path(str(files('seqtree').joinpath('data/control_human_trb_aa.txt.gz')))),
            'germline_background': Path(a.germline_background).expanduser().resolve()}
    def component_freeze():
        return {name: dict(path=str(path), sha256=sha256(path)) for name, path in component_paths.items()}
    component_before = component_freeze()
    before, raw_before = freeze(), db.provenance(a.vdjdb)
    if class_policy:
        bound = policy_info.get('source_reference_identity',{}).get(organism)
        identity_key = 'files' if 'files' in raw_before else 'input_sha256'
        if not isinstance(bound,dict) or identity_key not in bound or bound[identity_key]!=raw_before[identity_key]:
            raise ValueError('class policy was prepared from a different raw reference identity')
    raw = db.load(a.vdjdb, species=species)
    pooled, info = pool_epitopes(raw, reference_population=a.reference_population, min_records=1)
    info.pop('min_records')
    for row in info['epitope_inventory']:
        row['classification_eligible'] = row['valid_source_records'] >= a.min_records
    valid = pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
    if a.reference_population == 'curated-benchmark':
        valid &= pl.col('reference_valid').fill_null(False)
    if paired:
        from ..evalue.paired import build_paired_ref
        chains=pooled.filter(valid)
        selected=build_paired_ref(chains.with_columns(pl.lit(True).alias('reference_valid')))
        keys=['alpha','beta','v_alpha','v_beta','j_alpha','j_beta']
        support=chains.group_by('epitope','gene').len().group_by('epitope').agg(pl.col('len').max().alias('valid_source_records')).join(selected.group_by('epitope').len(name='paired_source_records'),on='epitope',how='left',validate='1:1').with_columns(pl.col('paired_source_records').fill_null(0)).sort('epitope')
        axes=[('alpha','v_alpha','TRA'),('beta','v_beta','TRB')]
    else:
        selected=pooled.filter((pl.col('gene')==a.locus)&valid)
        keys=['cdr3','v','j']
        support=selected.group_by('epitope').len(name='valid_source_records').sort('epitope')
        axes=[('cdr3','v',a.locus)]
    if selected.is_empty():raise ValueError('no valid reference records for the requested species/locus')
    support=support.with_columns((pl.col('valid_source_records')>=a.min_records).alias('classification_eligible'))
    eligible=support.filter(pl.col('classification_eligible')).select('epitope')
    if eligible.is_empty():raise ValueError('no reference epitopes meet min-records')
    queries=selected.join(eligible,on='epitope',how='semi').select(keys).unique().sort(keys,nulls_last=True).with_row_index('query_id')
    context=precompute_v_loop_distances(species=organism)
    resolved=['resolved_v'+str(i) for i in range(len(axes))] if paired else ['resolved_v']
    def resolve(frame):
        # Retain original unresolved keys even when their modeled loops coincide.
        for i, (_, v, locus) in enumerate(axes):
            if v_model_policy == 'resolved-first':
                frame = frame.with_columns(resolve_v_alleles(frame[v], context['model'], locus=locus)
                    .alias('baseline_v'+str(i)))
            alleles, imputed = impute_reference_v_alleles(frame[v], context['model'], locus=locus)
            flag = 'v_geometry_imputed'+str(i) if paired else 'v_geometry_imputed'
            frame = frame.with_columns(alleles.alias(resolved[i]), imputed.alias(flag)).with_columns(
                pl.when(pl.col(flag)).then(pl.struct(pl.col(v).alias('v')).struct.json_encode())
                    .otherwise(pl.col(resolved[i])).alias('_v_key'+str(i)))
        return frame.with_columns(pl.all_horizontal(*[(pl.col(seq).str.len_chars()>=8) for seq,_,_ in axes]).alias('model_available'))
    reference=resolve(selected);queries=resolve(queries)
    feature_columns=[column for i,(seq,_,_) in enumerate(axes) for column in [pl.col(seq).alias('seq'+str(i)),pl.col(resolved[i]).alias('v'+str(i))]]
    identity_keys=['_v_key'+str(i) for i in range(len(axes))]
    flags=['v_geometry_imputed'+str(i) for i in range(len(axes))] if paired else ['v_geometry_imputed']
    baseline_columns = ['baseline_v'+str(i) for i in range(len(axes))] if v_model_policy == 'resolved-first' else []
    modeled_reference=reference.filter(pl.col('model_available')).select(*feature_columns,*identity_keys,*baseline_columns,'epitope')
    modeled_queries=queries.filter(pl.col('model_available')).select('query_id',*feature_columns,*identity_keys,*baseline_columns)
    modeled_counts=reference.filter(pl.col('model_available')).group_by('epitope').agg(pl.len().alias('modeled_source_records'),pl.struct([*[x[0] for x in axes],*resolved,*identity_keys]).n_unique().alias('modeled_search_keys'),*[pl.col(c).sum().alias(c+'_source_records') for c in flags])
    support=support.join(modeled_counts,on='epitope',how='left',validate='1:1').with_columns(pl.col('modeled_source_records','modeled_search_keys').fill_null(0)).with_columns((pl.col('modeled_search_keys')>0).alias('model_available'))
    truth = selected.join(queries.select('query_id', *keys), on=keys, how='inner', nulls_equal=True, validate='m:1').select(
        'query_id', 'epitope', 'mhc_a', 'mhc_b', 'mhc_class').unique().sort('query_id', 'epitope')
    def metadata(frame):
        return frame.with_columns(pl.lit(species).alias('species'), pl.lit('PAIRED' if paired else a.locus).alias('locus'))
    scorers = ['rule_v', 'rule_spatial', 'rule_rank', *[
        'rule_'+prefix+mode for mode in MOTIFS for prefix in ('', 'mixed_', 'enriched_')]]
    if class_policy:
        scorers=['reference_class_policy']
    elif a.kernel_only:
        from ..evalue.kernel import KERNEL_SCORERS
        scorers=['kernel_'+name.replace('-','_') for name in KERNEL_SCORERS]
    rule_scorers = scorers.copy()
    component_scorers = (['accepted_original', 'current_original', 'current_original_log',
        'gapblock', 'stem4', 'motif'] if a.components else [])
    scorers += component_scorers
    targets = support['epitope'].to_list()
    winners = WinnerAccumulator(queries.height, targets, scorers)
    target_index = dict(zip(targets, range(len(targets))))
    prefix.parent.mkdir(parents=True, exist_ok=True)
    continuous = destinations[5]
    continuous.mkdir()
    outputs, pending, pending_rows, continuous_rows = {}, [], 0, 0
    policy_counts=[]
    policy_columns=['reference_class','selected_score_column','policy_status']
    def write_part(frame):
        path = continuous/f'part{len(outputs):06d}.parquet'
        frame.write_parquet(path)
        outputs[str(path.relative_to(prefix.parent))] = sha256(path)
    def export(frame):
        nonlocal pending, pending_rows, continuous_rows
        frame = metadata(frame).select('species', 'locus', 'query_id', 'epitope', 'scorer', 'score', 'supported',*(policy_columns if class_policy else []))
        pending.append(frame)
        pending_rows += frame.height
        continuous_rows += frame.height
        if pending_rows >= 65536:
            combined = pl.concat(pending)
            while combined.height >= 65536:
                write_part(combined.head(65536))
                combined = combined.slice(65536)
            pending, pending_rows = ([combined] if combined.height else []), combined.height
    def consume(scorer, frame):
        if class_policy:
            policy_counts.append(frame.group_by('epitope',*policy_columns).len(name='query_target_cases'))
        winners.consume(scorer, frame['query_id'].to_numpy(),
            frame['epitope'].replace_strict(target_index, return_dtype=pl.UInt32).to_numpy(),
            frame['score'].to_numpy(), frame['supported'].to_numpy())
        retained = frame.join(eligible, on='epitope', how='semi')
        if retained.height:
            export(retained.with_columns(pl.lit(scorer).alias('scorer')))
    if class_policy:
        from ..evalue.reference_rule_bank import stream_reference_policy
        def stream_backend(q,r,**kwargs):
            return stream_reference_policy(q,r,policy_frame,locus=a.locus,**kwargs)
    elif a.kernel_only:
        from ..evalue.reference_rule_bank import stream_reference_kernels
        stream_backend=stream_reference_kernels
    else:stream_backend=stream_reference_rules
    backend = (stream_backend(modeled_queries, modeled_reference, species=organism,
        threads=a.threads, paired=paired, consume=consume, identity_keys=identity_keys,
        **({'v_model_policy': v_model_policy} if a.kernel_only else {})) if modeled_reference.height and modeled_queries.height else {})
    # Fill only coordinates the backend cannot emit. Existing continuous parts
    # are final deliverables and are never reopened to calculate predictions.
    available_targets = support.filter(pl.col('model_available')).select('epitope')
    complete_labels = support.select('epitope') if class_policy else eligible
    absent_targets = complete_labels.join(available_targets, on='epitope', how='anti')
    for query_subset, label_subset in [
            (queries.filter(~pl.col('model_available')).select('query_id'), complete_labels),
            (queries.filter(pl.col('model_available')).select('query_id'), absent_targets)]:
        if query_subset.is_empty() or label_subset.is_empty():
            continue
        for start in range(0, query_subset.height, max(1, 65536//label_subset.height)):
            missing = query_subset.slice(start, max(1, 65536//label_subset.height)).join(label_subset, how='cross').with_columns(
                pl.lit(None, dtype=pl.Float64).alias('score'), pl.lit(False).alias('supported'))
            if class_policy:
                from .score_policy import apply_reference_policy
                selected_missing = apply_reference_policy(metadata(missing).with_columns(
                    pl.lit(None,dtype=pl.Float64).alias('kernel_v_bounded'),
                    pl.lit(None,dtype=pl.Float64).alias('rule_v_score'),
                    pl.lit(None,dtype=pl.Float64).alias('rule_exact3_score'),
                    pl.lit(None,dtype=pl.Float64).alias('rule_rank_score')),policy_frame,
                    original_column='kernel_v_bounded')
                consume('reference_class_policy',selected_missing.select('query_id','epitope',
                    pl.col('policy_score').alias('score'),pl.lit(False).alias('supported'),*policy_columns))
            else:
                for scorer in rule_scorers:
                    export(missing.with_columns(pl.lit(scorer).alias('scorer')))
    control_info = None
    if a.components:
        from ..evalue.control import raw_background
        from ..io import read_rearrangement
        from ..runner.retrieval_components import stream_components
        density_index, control_info = raw_background(a.locus, organism, path=a.control)
        if control_info['sha256'] != component_before['control']['sha256']:
            raise ValueError('raw control changed while loading')
        prior = read_rearrangement(component_paths['germline_background'], source='legacy')
        eligible_indices = np.array([target_index[label] for label in eligible['epitope']])
        emitted = np.zeros((len(component_scorers), queries.height, (len(targets)+7)//8), dtype=np.uint8)
        def component_consume(scorer, frame):
            component = component_scorers.index(scorer)
            frame = frame.select('query_id', 'epitope', 'score', 'supported')
            q = frame['query_id'].to_numpy()
            t = frame['epitope'].replace_strict(target_index, return_dtype=pl.UInt32).to_numpy()
            flags = np.left_shift(np.uint8(1), (t % 8).astype(np.uint8))
            if (np.any(emitted[component, q, t//8] & flags) or
                    np.unique(q.astype(np.uint64)*len(targets)+t).size != len(q)):
                raise ValueError('duplicate component query/target score coordinate')
            np.bitwise_or.at(emitted[component], (q, t//8), flags)
            consume(scorer, frame)
        component_reference = reference.with_columns(pl.col('v').alias('original_v'),
            pl.when(pl.col('v_geometry_imputed')).then(pl.col('resolved_v')).otherwise(pl.col('v')).alias('v'))
        if a.reference_population == 'all-valid':
            component_reference = component_reference.with_columns(pl.lit(True).alias('reference_valid'))
        component_queries = queries.select('query_id', *keys, 'resolved_v', 'v_geometry_imputed').with_columns(
            pl.col('v').alias('original_v'),
            pl.when(pl.col('v_geometry_imputed')).then(pl.col('resolved_v')).otherwise(pl.col('v')).alias('v'),
            pl.lit(a.locus).alias('locus')).drop('resolved_v')
        stream_components(component_queries,
            component_reference, {'density_index': density_index, 'prior_background': prior},
            species=organism, locus=a.locus, threads=a.threads, consume=component_consume)
        chunk = max(1, 65536//len(targets))
        for component, scorer in enumerate(component_scorers):
            for start in range(0, queries.height, chunk):
                present = np.unpackbits(emitted[component, start:start+chunk], axis=1,
                    count=len(targets), bitorder='little')[:, eligible_indices]
                rows, labels = np.nonzero(present == 0)
                if len(rows):
                    export(pl.DataFrame({'query_id': queries['query_id'].gather(rows+start),
                        'epitope': eligible['epitope'].gather(labels), 'scorer': [scorer]*len(rows),
                        'score': pl.Series([None]*len(rows), dtype=pl.Float64), 'supported': False}))
    if pending:
        write_part(pl.concat(pending))
    expected_rows = queries.height*eligible.height*len(scorers)
    if continuous_rows != expected_rows:
        raise ValueError(f'incomplete continuous output: {continuous_rows} != {expected_rows}')
    winner_rows = 0
    score_path = destinations[2]
    with score_path.open('x') as handle:
        for frame in winners.frames(queries['query_id'], max_rows=max(65536, len(targets))):
            metadata(frame).write_csv(handle, separator='\t', include_header=winner_rows == 0)
            winner_rows += frame.height
    outputs[score_path.name] = sha256(score_path)
    for path, frame in [(destinations[0], queries), (destinations[1], truth), (destinations[3], support)]:
        metadata(frame).write_csv(path, separator='\t')
        outputs[path.name] = sha256(path)
    if class_policy:
        audit = pl.concat(policy_counts).group_by('epitope',*policy_columns).agg(
            pl.col('query_target_cases').sum()).sort('epitope',*policy_columns)
        metadata(audit).write_csv(destinations[6],separator='\t')
        outputs[destinations[6].name]=sha256(destinations[6])
    after, raw_after = freeze(), db.provenance(a.vdjdb)
    component_after = component_freeze()
    if before != after or raw_before != raw_after or component_before != component_after or policy_before != policy_freeze():
        raise ValueError('reference or producing source/native changed during fresh scoring')
    destinations[4].write_text(json.dumps(dict(reference=raw_before, reference_after=raw_after,
        source_before=before, source_after=after, source_frozen=True, reference_frozen=True,
        species=species, locus=a.locus, reference_population=a.reference_population,
        min_records=a.min_records, threads=a.threads, scorers=scorers,
        raw_records=raw.height, valid_source_records=selected.height, queries=queries.height,
        modeled_queries=modeled_queries.height, candidate_targets=len(targets), eligible_targets=eligible.height,
        winner_rows=winner_rows, continuous_rows=continuous_rows, expected_continuous_rows=expected_rows,
        backend=backend, components=a.components, kernel_only=a.kernel_only, paired=paired,
        v_model_policy=v_model_policy,
        kernel_model_binding=('resolved-first candidate; resolved queries use fixed original resolved-allele references; '
            'unresolved queries use fixed extended references; existing extended calibration does not transfer'
            if v_model_policy == 'resolved-first' else 'unchanged extended reference model'),
        class_policy=dict(inputs_before=policy_before,inputs_after=policy_freeze(),inputs_frozen=True,
            routes=policy_info['routes'],original_column='kernel_v_bounded',experimental=True,
            default_changed=False,labels_used=False,calibration_status='OPEN; component P/E/posteriors do not transfer',
            candidate_models='all valid fixed labels including missing-profile tails; eligibility limits only queries and continuous exports') if class_policy else None,
        paired_eligibility="maximum qualifying chain source-record count before exact-hit filtering" if paired else None, component_inputs_before=component_before,
        model_availability='length>=8 on every present scoring axis; unresolved V receives a flagged fixed first complete co-locus loop template; original keys retained separately from geometry; source units are genuine complete complexes for paired mode, chain observations for singles; component availability remains explicit',
        component_inputs_after=component_after, component_inputs_frozen=True, control_info=control_info,
        **info, exact_exclusion='hit-list only; full fixed reference models retained',
        score_meaning=('experimental frozen reference-class routing of existing scores; larger is stronger; final calibration OPEN' if class_policy else
            'uncalibrated current stem4 kernel evidence; larger is stronger' if a.kernel_only else 'uncalibrated fixed-reference rule evidence; larger is stronger'),
        calculation='fresh raw reference and scores; no intermediate reuse or resume',
        software={name: version(name) for name in ('vdjmatch', 'seqtree', 'vdjtools', 'polars', 'numpy', 'scipy')},
        tables=outputs, seconds=time.monotonic()-started), indent=2)+'\n')
    return 0
