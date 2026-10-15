"""Tie-aware empirical thresholds for one winning-label decision per query.

Inputs are final evidence, not predictions generated here. Null sources retain separate
query denominators. This reduction does not establish out-of-sample biological confidence.
"""
from pathlib import Path
import math


def register(subparsers):
    parser = subparsers.add_parser('calibrate-calls', help='empirical any-call operating points')
    parser.add_argument('cases', help='TSV: query_id, score, outcome, null_source')
    parser.add_argument('--target-fpr', type=float, required=True)
    parser.add_argument('--lower-is-better', action='store_true')
    parser.add_argument('--score-column',default='score',help='numeric selected-score column; frozen in the recipe')
    parser.add_argument('--input-manifest', help='optional score producer manifest; freeze its pipeline binding')
    parser.add_argument('--output-prefix', required=True)
    parser.set_defaults(func=main)
    apply = subparsers.add_parser('apply-call-calibration', help='reject independent query calls using a frozen empirical cutoff')
    apply.add_argument('cases', help='TSV: unique query_id and score; other columns are preserved evidence')
    apply.add_argument('--calibration', required=True, help='matching scorer/reference calibrate-calls JSON')
    apply.add_argument('--input-manifest', help='score producer manifest; required for a pipeline-bound recipe')
    apply.add_argument('--null-cases', help='original calibrate-calls input, hash-bound to the frozen recipe; add per-query null tails')
    apply.add_argument('--output', required=True)
    apply.set_defaults(func=apply_main)
    return parser



PAIRED_STATISTIC='maximum of final selected linked, alpha and beta kernel scores'


def paired_call_event():
    """The declared any-axis event; no combined label or evidence ranking."""
    return dict(score_column='kernel_any_axis_score',statistic=PAIRED_STATISTIC,
        axes=['linked','alpha','beta'],selection='each axis: max restriction then max peptide; ties retained',
        linked='alpha+beta distance to genuine linked references',
        marginals='full compatible TRA/TRB references including unlinked observations',
        exact_policy='linked: both full junctions; marginals: own full junction; V/J do not change puncture')


def paired_pipeline_sources():
    """Pin both shipped paired producers and their shared event/binding definition."""
    import hashlib
    return {name:hashlib.sha256(Path(__file__).with_name(name+'.py').read_bytes()).hexdigest()
            for name in ('reference_calls','neighbours','calibration')}


def _bound_score_columns(binding):
    return ('score','kernel_any_axis_score') if binding['context'][1]=='paired' else ('score','kernel_score')

def pipeline_binding(manifest):
    """Canonical kernel and selected-event identity from either shipped producer.

    This is a provenance check, not calibration. A single-population cohort
    manifest can wrap its core producer; mixed populations require their
    own declared selection and null contract.
    """
    import hashlib
    import json
    from importlib.resources import files
    from ..match.geometry import load_v_loops

    try:
        if 'core_manifests' in manifest:
            if len(manifest['core_manifests']) != 1:
                raise ValueError('binding requires one single-chain population')
            manifest = manifest['core_manifests'][0]
        reference = manifest['reference']['input_sha256']
        policy = manifest['reference_population']
        if policy not in ('all-valid', 'curated-benchmark'):
            raise ValueError('binding requires a known reference population')
        if 'geometry' in manifest:
            options = manifest['options']
            if options['scorer'] != 'stem4':
                raise ValueError('binding requires final kernel scoring')
            if (manifest['selection'] != 'maximum restriction per peptide, then maximum peptide; exact score ties remain ambiguous'
                    or manifest['calibration_statistic'] != (PAIRED_STATISTIC if options['locus']=='paired' else 'final selected kernel_score; no component P/E transfer')):
                raise ValueError('binding selection/statistic is unsupported')
            geometry = {k:v for k,v in manifest['geometry'].items() if k != 'threads'}
            sources = {k:manifest['source_sha256'][k] for k in ('kernel','geometry','seqtree_native','seqtree_python')}
        else:
            options = manifest['parameters']
            if (not options['all_targets'] or options['background'] != 'vdjdb-other'
                    or options['neighbour_weighting'] != 'kernel'
                    or any(options.get(k) for k in ('epitope','mhc_a','mhc_b','targets_from_sample',
                        'withhold_epitope','geometry_extra','reference_radius','nearest_rank','conditional_j','epitope_rules','chemistry'))):
                raise ValueError('binding requires the complete supported all-label kernel selection')
            g, kernel = manifest['distance_scale'], manifest['kernel']
            if kernel['score_column'] != 'kernel_score' or manifest['junction_end_trim']['n_terminal'] != 3:
                raise ValueError('binding statistic/trim is unsupported')
            geometry = dict(species={'HomoSapiens':'human','MusMusculus':'mouse'}[options['species']],
                matrix=g['matrix'], ctrim=manifest['junction_end_trim']['c_terminal'], v_loop_weight=g['v_loop_weight'],
                alphabet=g['alphabet'], stem_depth=g['stem_depth'], stem_weight=g['stem_weight'],
                exclude_exact=options['exclude_exact'], radius=g['radius'], kernel_scale=kernel['kernel_scale'],
                cutoffs=kernel['cutoffs'], cross_v=kernel['cross_v'], floor=.01, scorer=kernel['scorer'],
                matrix_scale=g['matrix_unit'], v_regions=g['v_regions'], position_weighting=options['position_weighting'])
            sources = {name:manifest['software'][key] for name,key in (
                ('kernel','kernel_source_sha256'), ('geometry','scorer_source_sha256'),
                ('seqtree_native','seqtree_native_sha256'), ('seqtree_python','seqtree_python_source_sha256'))}
        species, locus = options['species'], options['locus']
        organism = {'HomoSapiens':'human','MusMusculus':'mouse'}[species]
        if locus not in ('TRA','TRB','paired') or geometry['species'] != organism or not isinstance(geometry['exclude_exact'], bool):
            raise ValueError('binding requires one receptor population, species and exact policy')
        required_geometry = {'species','matrix','ctrim','v_loop_weight','alphabet','stem_depth','stem_weight',
            'exclude_exact','radius','kernel_scale','cutoffs','cross_v','floor','scorer','matrix_scale','v_regions','position_weighting'}
        if set(geometry) != required_geometry:
            raise ValueError('binding geometry is incomplete or unsupported')
        model_file = Path(str(files('vdjmatch.resources')/'vloops'/f'{organism}_v_loops.tsv'))
        file_hash = hashlib.sha256(model_file.read_bytes()).hexdigest()
        semantic_hash = hashlib.sha256(json.dumps(load_v_loops(organism), sort_keys=True, default=str).encode()).hexdigest()
        if (manifest['model_sha256'] != semantic_hash if 'geometry' in manifest
                else manifest['v_loop_model']['sha256'] != file_hash):
            raise ValueError('binding model differs from installed germline model')
        hashes = [reference, *sources.values()]
        if any(not isinstance(x,str) or len(x)!=64 or any(c not in '0123456789abcdef' for c in x) for x in hashes):
            raise ValueError('binding requires complete SHA-256 identities')
        if locus=='paired':
            if manifest.get('paired_call_event')!=paired_call_event():
                raise ValueError('binding requires the declared paired any-axis event')
            if manifest.get('paired_source_sha256')!=paired_pipeline_sources():
                raise ValueError('paired producer sources differ from installed sources')
            sources=dict(sources,**{'paired_'+k:v for k,v in manifest['paired_source_sha256'].items()})
        return dict(context=[species,locus], reference_population=policy, reference_sha256=reference,
            geometry=geometry, sources=sources, model_file_sha256=file_hash, model_semantic_sha256=semantic_hash,
            statistic=PAIRED_STATISTIC if locus=='paired' else 'final selected kernel_score',
            selection='independent linked/alpha/beta selections; maximum scalar only' if locus=='paired' else 'max restriction then max peptide; ties retained',
            exact_policy=(('linked: both full junctions; marginals: own full junction; no ordinary exact-retained transfer' if locus=='paired' else 'full-junction puncture across V/J; no ordinary exact-retained transfer') if geometry['exclude_exact']
                          else 'exact evidence retained across all labels; ordinary annotation'))
    except (KeyError, TypeError, AttributeError) as exc:
        raise ValueError('binding manifest is missing required producer identities') from exc


def apply_calibration(frame, calibration, *, input_manifest=None):
    """Apply a frozen cutoff without reading query labels or refitting it.

    Input scores remain ranking evidence; only call_score represents an accepted
    call. No-score queries and rejected candidates are no_supported_match, not
    proven nonbinders. Bound recipes require a matching producer manifest;
    callers of generic unbound recipes own their scorer/reference compatibility.
    """
    import polars as pl
    if 'pipeline_binding' in calibration:
        if input_manifest is None:
            raise ValueError('pipeline binding requires --input-manifest')
        if pipeline_binding(input_manifest) != calibration['pipeline_binding']:
            raise ValueError('pipeline binding does not match score producer manifest')
        if calibration.get('score_column','score') not in _bound_score_columns(calibration['pipeline_binding']) or calibration.get('comparison') != '>=':
            raise ValueError('pipeline binding requires the final higher-is-better kernel score')
    score_column=calibration.get('score_column','score')
    if not isinstance(score_column,str) or not score_column.strip():
        raise ValueError('score_column must name a numeric input column')
    if missing := {'query_id', score_column} - set(frame.columns):
        raise ValueError('missing columns: ' + ', '.join(sorted(missing)))
    if {'call_score', 'call_status'} & set(frame.columns):
        raise ValueError('call output columns already exist')
    if frame['query_id'].null_count() or frame['query_id'].n_unique() != frame.height:
        raise ValueError('query_id must be unique and non-null; reduce candidate labels first')
    scores = frame[score_column].cast(pl.Float64)
    if not scores.drop_nulls().is_finite().all():
        raise ValueError('non-null scores must be finite; use null for no evidence')
    action, cutoff = calibration.get('action'), calibration.get('cutoff')
    comparison = calibration.get('comparison')
    if comparison not in ('>=', '<='):
        raise ValueError('calibration comparison must be >= or <=')
    if action == 'abstain_all':
        if cutoff is not None:
            raise ValueError('abstain_all requires a null cutoff')
        accepted = pl.Series('accepted', [False] * frame.height, dtype=pl.Boolean)
    elif action == 'apply_cutoff':
        if isinstance(cutoff, bool) or not isinstance(cutoff, (int, float)) or not math.isfinite(cutoff):
            raise ValueError('apply_cutoff requires a finite numeric cutoff')
        accepted = (scores >= cutoff if comparison == '>=' else scores <= cutoff).fill_null(False)
    else:
        raise ValueError('calibration action must be apply_cutoff or abstain_all')
    return frame.with_columns(
        pl.when(accepted).then(scores).otherwise(None).alias('call_score'),
        pl.when(accepted).then(pl.lit('accepted')).otherwise(pl.lit('no_supported_match')).alias('call_status'))



from ..evalue.directional import null_tail_probabilities


def apply_main(args):
    import json
    from ..db.cache import sha256
    from ..io.airr import _read_table
    recipe = Path(args.calibration)
    frozen = json.loads(recipe.read_text())
    input_manifest_path = Path(args.input_manifest) if getattr(args,'input_manifest',None) else None
    input_manifest = json.loads(input_manifest_path.read_text()) if input_manifest_path else None
    out = apply_calibration(_read_table(args.cases), frozen, input_manifest=input_manifest)
    null_path = Path(args.null_cases) if getattr(args,'null_cases',None) else None
    tails = None
    if null_path is not None:
        if sha256(null_path) != frozen.get('input_sha256'):
            raise ValueError('null cases must match the frozen calibration input hash')
        nulls = _read_table(null_path)
        out,tails = null_tail_probabilities(out,nulls.filter(nulls['outcome']=='null'),
            higher_is_better=frozen['comparison']=='>=',score_column=frozen.get('score_column','score'))
    target = Path(args.output)
    manifest_path = Path(str(target)+'.manifest.json')
    tail_path=Path(str(target)+'.null-tails.tsv')
    inputs={Path(args.cases).resolve(),recipe.resolve()} | ({null_path.resolve()} if null_path else set())
    if input_manifest_path is not None:
        inputs.add(input_manifest_path.resolve())
    if {target.resolve(),manifest_path.resolve(),tail_path.resolve()} & inputs:
        raise ValueError('output must not overwrite inputs')
    target.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(target, separator='\t')
    if tails is not None:
        tails.write_csv(tail_path,separator='\t')
    accepted = int((out['call_status'] == 'accepted').sum())
    manifest = {'input_sha256': sha256(Path(args.cases)), 'calibration_sha256': sha256(recipe),
        'queries': out.height, 'accepted_calls': accepted, 'rejected_queries': out.height-accepted,
        'no_score_queries': out[frozen.get('score_column','score')].null_count(), 'query_labels_used': False,
        'score_column':frozen.get('score_column','score'),
        **({'pipeline_binding':frozen['pipeline_binding'], 'input_manifest_sha256':sha256(input_manifest_path)}
           if 'pipeline_binding' in frozen else {}),
        **({'null_cases_sha256':sha256(null_path),'null_tail_table_sha256':sha256(tail_path),
            'P_definition':'max source (k+1)/(N+1), including ties and no-hit null queries; exchangeability-conditional rank P, no-score P1',
            'E_definition':'max source k/N, observed false calls per query at the given score; no-score E0; not biological binding probability',
            'tail_selection':'same frozen score and winning-label pipeline required; source-specific tails exported separately'} if tails is not None else {}),
        'interpretation': 'Rejected evidence means no supported VDJdb match, not biological nonbinding. Input scores and metadata remain evidence, not accepted calls.',
        'limitations': 'The frozen cutoff must match the entire scorer/reference/geometry/selection pipeline. No biological probability or out-of-sample null-rate guarantee.'}
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    return 0


def operating_points(frame, target_fpr, *, higher_is_better=True):
    """Select the most permissive tied cutoff meeting every empirical null-source rate.

    ``outcome`` is correct/incorrect for labelled reference queries and null for null
    queries. Null score means abstention and stays in all denominators. Correctness
    describes the single winning label; callers must resolve ambiguity before input.
    Labelled queries are optional; null-only fitting leaves recall unavailable.
    """
    import numpy as np
    import polars as pl

    if not math.isfinite(target_fpr) or not 0 <= target_fpr < 1:
        raise ValueError('target-fpr must be finite and in [0, 1)')
    required = {'query_id', 'score', 'outcome', 'null_source'}
    if missing := required - set(frame.columns):
        raise ValueError('missing columns: ' + ', '.join(sorted(missing)))
    if frame['query_id'].null_count() or frame['query_id'].n_unique() != frame.height:
        raise ValueError('query_id must be unique and non-null; reduce candidate labels first')
    if frame['outcome'].null_count() or not frame['outcome'].is_in(['correct', 'incorrect', 'null']).all():
        raise ValueError('outcome must be correct, incorrect or null')
    frame = frame.with_columns(pl.col('score').cast(pl.Float64))
    if not frame['score'].drop_nulls().is_finite().all():
        raise ValueError('non-null scores must be finite; use null for no call')
    if frame.filter((pl.col('outcome') == 'correct') & pl.col('score').is_null()).height:
        raise ValueError('a correct winning label requires a score')
    null = frame.filter(pl.col('outcome') == 'null')
    labelled = frame.filter(pl.col('outcome') != 'null')
    if not null.height:
        raise ValueError('null queries are required')
    if null['null_source'].null_count() or null['null_source'].cast(pl.String).str.strip_chars().eq('').any():
        raise ValueError('every null query requires a nonempty null_source')
    sign = 1 if higher_is_better else -1
    cuts = np.unique(sign * frame['score'].drop_nulls().to_numpy())[::-1]

    def counts(part):
        scores = np.sort(sign * part['score'].drop_nulls().to_numpy())
        return np.r_[0, len(scores) - np.searchsorted(scores, cuts, side='left')]

    correct = counts(labelled.filter(pl.col('outcome') == 'correct'))
    incorrect = counts(labelled.filter(pl.col('outcome') == 'incorrect'))
    eligible = np.ones(len(cuts) + 1, dtype=bool)
    sources = []
    empirical_expected = np.zeros(len(cuts) + 1, dtype=float)
    for source in sorted(null['null_source'].unique().to_list()):
        part = null.filter(pl.col('null_source') == source)
        calls = counts(part)
        eligible &= calls / part.height <= target_fpr
        empirical_expected = np.maximum(empirical_expected, calls / part.height)
        sources.append((source, part.height, calls))
    selected = int(np.flatnonzero(eligible)[-1])
    cutoff = None if selected == 0 else float(sign * cuts[selected - 1])
    precision = np.divide(correct, correct + incorrect, out=np.full(len(correct), np.nan),
                          where=(correct + incorrect) > 0)
    result = {
        'status': 'empirical_calibration_only',
        'action': 'abstain_all' if selected == 0 else 'apply_cutoff',
        'cutoff': cutoff, 'comparison': '>=' if higher_is_better else '<=',
        'target_null_query_call_rate': target_fpr,
        'vdjdb_queries': labelled.height,
        'correct_calls': int(correct[selected]), 'incorrect_calls': int(incorrect[selected]),
        'vdjdb_correct_call_recall': float(correct[selected] / labelled.height) if labelled.height else None,
        'vdjdb_call_precision': float(precision[selected]) if np.isfinite(precision[selected]) else None,
        'no_score_queries': frame['score'].null_count(),
        'empirical_E_false_calls_per_query': float(empirical_expected[selected]),
        'E_definition': ('Maximum observed null-source any-call tail after winning-label selection; '
                         'expected false calls per query under the supplied null populations. '
                         'For Q queries at a common cutoff, empirical expected false calls = Q * E. '
                         'Not a finite-sample upper confidence bound, Poisson neighbour E, '
                         'test-martingale e-value or biological specificity probability.'),
        'nulls': [{'source': s, 'queries': n, 'calls': int(c[selected]),
                  'call_rate': float(c[selected] / n)} for s, n, c in sources],
    }
    # A long flat table keeps each null mechanism and denominator explicit.
    base = pl.DataFrame({
        'cutoff': pl.Series([None, *(sign * cuts).tolist()], dtype=pl.Float64),
        'action': ['abstain_all', *(['apply_cutoff'] * len(cuts))],
        'eligible': eligible, 'selected': np.arange(len(eligible)) == selected,
        'correct_calls': correct, 'incorrect_calls': incorrect,
        'vdjdb_queries': [labelled.height] * len(eligible),
        'vdjdb_correct_call_recall': (correct / labelled.height if labelled.height else
            pl.Series([None] * len(eligible), dtype=pl.Float64)),
        'vdjdb_call_precision': pl.Series(precision).fill_nan(None),
        'empirical_E_false_calls_per_query': empirical_expected,
    })
    curve = pl.concat([base.with_columns(
        pl.lit(s).alias('null_source'), pl.lit(n).alias('null_queries'),
        pl.Series('null_calls', c), pl.Series('null_query_call_rate', c / n))
        for s, n, c in sources])
    return result, curve


def main(args):
    import json
    from ..db.cache import sha256
    from ..io.airr import _read_table

    frame=_read_table(args.cases)
    score_column=getattr(args,'score_column','score')
    input_manifest_path=Path(args.input_manifest) if getattr(args,'input_manifest',None) else None
    binding=None
    if input_manifest_path is not None:
        binding=pipeline_binding(json.loads(input_manifest_path.read_text()))
        if score_column not in _bound_score_columns(binding) or args.lower_is_better:
            raise ValueError('pipeline binding requires its final higher-is-better kernel score')
    if score_column not in frame.columns:
        raise ValueError('missing score column: '+score_column)
    result, curve = operating_points(frame.with_columns(frame[score_column].alias('score')), args.target_fpr,
                                    higher_is_better=not args.lower_is_better)
    result['score_column']=score_column
    result['input_sha256'] = sha256(Path(args.cases))
    if binding is not None:
        result['pipeline_binding']=binding
        result['input_manifest_sha256']=sha256(input_manifest_path)
    result['limitations'] = ('Empirical rates on supplied query populations; not independent validation '
                            'or biological specificity probability. Bind source/scorer/reference identities '
                            'before recommending or applying this cutoff to another run.')
    prefix = Path(args.output_prefix)
    if input_manifest_path is not None:
        outputs={Path(str(prefix)+suffix).resolve() for suffix in ('.operating_points.tsv','.calibration.json')}
        if outputs & {Path(args.cases).resolve(),input_manifest_path.resolve()}:
            raise ValueError('output must not overwrite binding inputs')
    prefix.parent.mkdir(parents=True, exist_ok=True)
    curve.write_csv(str(prefix) + '.operating_points.tsv', separator='\t')
    Path(str(prefix) + '.calibration.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return 0
