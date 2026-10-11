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
    parser.add_argument('--output-prefix', required=True)
    parser.set_defaults(func=main)
    apply = subparsers.add_parser('apply-call-calibration', help='reject independent query calls using a frozen empirical cutoff')
    apply.add_argument('cases', help='TSV: unique query_id and score; other columns are preserved evidence')
    apply.add_argument('--calibration', required=True, help='matching scorer/reference calibrate-calls JSON')
    apply.add_argument('--output', required=True)
    apply.set_defaults(func=apply_main)
    composite = subparsers.add_parser('composite-evidence', help='experimental two-scorer evidence on a fixed empirical background scale')
    composite.add_argument('cases', help='TSV: unique query_id, original_score, conditional_score')
    composite.add_argument('--background', required=True, help='same scored background queries for both scorers, independent of cases')
    composite.add_argument('--paired', action='store_true', help='rank original alpha/beta density and germline components against the fixed background')
    composite.add_argument('--output', required=True)
    composite.set_defaults(func=composite_main)
    bg = subparsers.add_parser('composite-background', help='raw reference-derived queries for one composite target; no predictions')
    bg.add_argument('--vdjdb', required=True)
    bg.add_argument('--species', default='HomoSapiens')
    bg.add_argument('--locus', choices=['TRA', 'TRB', 'paired'], required=True)
    for field in ['epitope', 'mhc-a', 'mhc-b', 'mhc-class']:
        bg.add_argument('--'+field, required=True)
    bg.add_argument('--output', required=True)
    bg.set_defaults(func=background_main)
    return parser


def background_queries(reference, target, *, paired=False):
    """Other-peptide raw keys; full chain identities and observed pairing only."""
    import polars as pl
    from ..aggregate.candidates import PMHC
    from ..evalue.paired import build_paired_ref
    if set(target) != {'species', *PMHC} or any(not isinstance(v, str) or not v.strip() for v in target.values()):
        raise ValueError('one complete species/peptide/MHC target required')
    valid = reference.filter(pl.col('reference_valid').fill_null(False))
    if valid['species'].n_unique() != 1 or valid['species'][0] != target['species']:
        raise ValueError('reference must contain the declared species only')
    if paired:
        population = build_paired_ref(valid).filter(pl.col('epitope') != target['epitope'])
        keys = ['alpha', 'v_alpha', 'j_alpha', 'beta', 'v_beta', 'j_beta']
        population = population.select(keys).unique().sort(keys).with_row_index('clone_id')
        rows = [population.select('clone_id', pl.lit(locus).alias('locus'),
            pl.col(chain).alias('junction_aa'), pl.col('v_'+chain).alias('v_call'),
            pl.col('j_'+chain).alias('j_call')) for chain, locus in [('alpha', 'TRA'), ('beta', 'TRB')]]
        raw = pl.concat(rows).sort('clone_id', 'locus').with_row_index('query_id')
    else:
        if valid['gene'].n_unique() != 1 or valid['gene'][0] not in ('TRA', 'TRB'):
            raise ValueError('single-chain reference requires one TRA/TRB locus')
        raw = valid.filter(pl.col('epitope') != target['epitope']).select(
            pl.col('gene').alias('locus'), pl.col('cdr3').alias('junction_aa'),
            pl.col('v').alias('v_call'), pl.col('j').alias('j_call'))
        raw = raw.unique().sort(raw.columns).with_row_index('query_id')
    return raw.with_columns(*[pl.lit(v).alias(k) for k, v in target.items()])


def background_main(args):
    import json
    import time
    from .. import db
    from ..db.cache import sha256
    from ..evalue.control import _organism
    started = time.monotonic()
    species = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[_organism(args.species)]
    reference = db.load(args.vdjdb, species=species, gene=None if args.locus == 'paired' else args.locus)
    target = {key: getattr(args, key) for key in ['epitope', 'mhc_a', 'mhc_b', 'mhc_class']}
    target['species'] = species
    raw = background_queries(reference, target, paired=args.locus == 'paired')
    output = Path(args.output)
    manifest = Path(str(output)+'.manifest.json')
    if Path(args.vdjdb).resolve() in {output.resolve(), manifest.resolve()}:
        raise ValueError('output must not overwrite reference input')
    output.parent.mkdir(parents=True, exist_ok=True)
    raw.write_csv(output, separator='\t')
    manifest.write_text(json.dumps({'reference_sha256': sha256(Path(args.vdjdb)), 'target': target,
        'locus': args.locus, 'source_records': reference.height,
        'valid_source_records': reference.filter('reference_valid').height,
        'query_rows': raw.height, 'queries': raw['clone_id'].n_unique() if args.locus == 'paired' else raw.height,
        'query_sha256': sha256(output), 'wall_seconds': time.monotonic()-started,
        'population': 'distinct full junction/V/J keys from other peptides; genuine linked pairs only in paired mode',
        'predictions_computed': False, 'query_labels_used': False,
        'policy': 'Score against the entire valid reference with benchmark exact-clonotype exclusion. This population does not establish null exchangeability or biological nonbinding.'}, indent=2)+'\n')
    return 0


def composite_evidence(frame, background, *, paired=False):
    """Bonferroni minimum of two upper empirical tails on a fixed background.

    One target/species/locus/geometry per call. Background rows retain draw
    multiplicity and no-score queries. Only scores/IDs are read, never labels.
    Exchangeability of background and null test scores is required for a null
    bound; it is not established by these tables or by database annotations.
    """
    import numpy as np
    import polars as pl
    components = ['original_alpha_density_score', 'original_beta_density_score', 'original_germline_lr']
    fields = [*(components if paired else ['original_score']), 'conditional_score']
    outputs = ['p_original_empirical', 'p_conditional_empirical', 'p_composite_empirical', 'score', 'composite_status']
    if set(outputs) & set(frame.columns):
        raise ValueError('composite output columns already exist')
    for name, table in [('cases', frame), ('background', background)]:
        if missing := {'query_id', *fields} - set(table.columns):
            raise ValueError(name+' missing columns: '+', '.join(sorted(missing)))
        if table['query_id'].null_count() or table['query_id'].n_unique() != table.height:
            raise ValueError(name+' query_id must be unique and non-null')
        for field in fields:
            if not table[field].cast(pl.Float64).drop_nulls().is_finite().all():
                raise ValueError(name+' non-null scores must be finite; use null for no evidence')
    for field in ['species', 'locus', 'epitope', 'mhc_a', 'mhc_b', 'mhc_class', 'reference_sha256', 'scoring_sha256']:
        if field in frame.columns or field in background.columns:
            if field not in frame.columns or field not in background.columns:
                raise ValueError('context column required in both tables: '+field)
            values = background[field].unique().to_list()
            if len(values) != 1 or any(v != values[0] for v in frame[field].unique().to_list()):
                raise ValueError('one matching target/scorer context required: '+field)
    if not background.height:
        raise ValueError('nonempty background required')
    if paired:
        raw_query = frame.select(pl.col(components).cast(pl.Float64)).to_numpy()
        raw_control = background.select(pl.col(components).cast(pl.Float64)).to_numpy()
        valid_query = np.isfinite(raw_query).all(axis=1)
        valid_control = np.isfinite(raw_control).all(axis=1)
        query_rank = np.zeros(frame.height)
        control_rank = np.zeros(background.height)
        for i in range(3):
            axis = np.sort(raw_control[valid_control, i])
            query_rank += np.searchsorted(axis, raw_query[:, i], side='right')
            control_rank += np.searchsorted(axis, raw_control[:, i], side='right')
        original_query = np.where(valid_query, query_rank, np.nan)
        original_control = control_rank[valid_control]
    else:
        original_query = frame['original_score'].cast(pl.Float64).to_numpy()
        original_control = background['original_score'].cast(pl.Float64).drop_nulls().to_numpy()
    pairs = [(original_query, original_control),
        (frame['conditional_score'].cast(pl.Float64).to_numpy(),
         background['conditional_score'].cast(pl.Float64).drop_nulls().to_numpy())]
    tails, available = [], []
    for scores, values in pairs:
        controls = np.sort(values)
        present = np.isfinite(scores)
        tail = np.ones(frame.height)
        tail[present] = (1+len(controls)-np.searchsorted(controls, scores[present], side='left'))/(background.height+1)
        tails.append(tail)
        available.append(present)
    combined = np.minimum(1., 2*np.minimum(*tails))
    observed = available[0] | available[1]
    return frame.with_columns(
        pl.Series(outputs[0], tails[0]), pl.Series(outputs[1], tails[1]), pl.Series(outputs[2], combined),
        pl.Series('score', np.where(observed, -np.log10(combined), np.nan), nan_to_null=True),
        pl.Series('composite_status', np.where(~observed, 'no_evidence',
            np.where(available[0] & available[1], 'available', 'single_scorer'))))


def composite_main(args):
    import json
    from ..db.cache import sha256
    from ..io.airr import _read_table
    cases, background = _read_table(args.cases), _read_table(args.background)
    if sha256(Path(args.cases)) == sha256(Path(args.background)):
        raise ValueError('cases and background must not be the same scored table')
    out = composite_evidence(cases, background, paired=args.paired)
    target = Path(args.output)
    manifest = Path(str(target)+'.manifest.json')
    if {target.resolve(), manifest.resolve()} & {Path(args.cases).resolve(), Path(args.background).resolve()}:
        raise ValueError('output must not overwrite inputs')
    target.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(target, separator='\t')
    manifest.write_text(json.dumps({'input_sha256': sha256(Path(args.cases)),
        'background_sha256': sha256(Path(args.background)), 'queries': cases.height,
        'background_queries': background.height, 'query_labels_used': False,
        'paired_original': 'sum of inclusive background component ranks; conservative compared with pooled-query/full conformal rank' if args.paired else None,
        'formula': 'p_j=(1+number of background scores >= query score)/(background queries+1); p=min(1,2*min(p_original,p_conditional)); score=-log10(p)',
        'status': 'experimental_unvalidated_composite',
        'limitations': 'One target/species/locus/geometry and identical frozen scorer/reference/selection policies for cases/background are required but not inferred from tables. Background/null-score exchangeability is unproved. Not a binder probability or default promotion.'}, indent=2)+'\n')
    return 0


def apply_calibration(frame, calibration):
    """Apply a frozen cutoff without reading query labels or refitting it.

    Input scores remain ranking evidence; only call_score represents an accepted
    call. No-score queries and rejected candidates are no_supported_match, not
    proven nonbinders. The caller must use matching scorer/reference identities.
    """
    import polars as pl
    if missing := {'query_id', 'score'} - set(frame.columns):
        raise ValueError('missing columns: ' + ', '.join(sorted(missing)))
    if {'call_score', 'call_status'} & set(frame.columns):
        raise ValueError('call output columns already exist')
    if frame['query_id'].null_count() or frame['query_id'].n_unique() != frame.height:
        raise ValueError('query_id must be unique and non-null; reduce candidate labels first')
    scores = frame['score'].cast(pl.Float64)
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


def apply_main(args):
    import json
    from ..db.cache import sha256
    from ..io.airr import _read_table
    recipe = Path(args.calibration)
    out = apply_calibration(_read_table(args.cases), json.loads(recipe.read_text()))
    target = Path(args.output)
    manifest_path = Path(str(target)+'.manifest.json')
    if {target.resolve(), manifest_path.resolve()} & {Path(args.cases).resolve(), recipe.resolve()}:
        raise ValueError('output must not overwrite inputs')
    target.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(target, separator='\t')
    accepted = int((out['call_status'] == 'accepted').sum())
    manifest = {'input_sha256': sha256(Path(args.cases)), 'calibration_sha256': sha256(recipe),
        'queries': out.height, 'accepted_calls': accepted, 'rejected_queries': out.height-accepted,
        'no_score_queries': out['score'].null_count(), 'query_labels_used': False,
        'interpretation': 'Rejected evidence means no supported VDJdb match, not biological nonbinding. Input scores and metadata remain evidence, not accepted calls.',
        'limitations': 'The frozen cutoff must match the entire scorer/reference/geometry/selection pipeline. No biological probability or out-of-sample null-rate guarantee.'}
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    return 0


def operating_points(frame, target_fpr, *, higher_is_better=True):
    """Select the most permissive tied cutoff meeting every empirical null-source rate.

    ``outcome`` is correct/incorrect for labelled reference queries and null for null
    queries. Null score means abstention and stays in all denominators. Correctness
    describes the single winning label; callers must resolve ambiguity before input.
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
    if not labelled.height or not null.height:
        raise ValueError('both labelled VDJdb queries and null queries are required')
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
        'vdjdb_correct_call_recall': float(correct[selected] / labelled.height),
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
        'vdjdb_correct_call_recall': correct / labelled.height,
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

    result, curve = operating_points(_read_table(args.cases), args.target_fpr,
                                    higher_is_better=not args.lower_is_better)
    result['input_sha256'] = sha256(Path(args.cases))
    result['limitations'] = ('Empirical rates on supplied query populations; not independent validation '
                            'or biological specificity probability. Bind source/scorer/reference identities '
                            'before recommending or applying this cutoff to another run.')
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    curve.write_csv(str(prefix) + '.operating_points.tsv', separator='\t')
    Path(str(prefix) + '.calibration.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return 0
