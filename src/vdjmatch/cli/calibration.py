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
    parser.add_argument('--output-prefix', required=True)
    parser.set_defaults(func=main)
    apply = subparsers.add_parser('apply-call-calibration', help='reject independent query calls using a frozen empirical cutoff')
    apply.add_argument('cases', help='TSV: unique query_id and score; other columns are preserved evidence')
    apply.add_argument('--calibration', required=True, help='matching scorer/reference calibrate-calls JSON')
    apply.add_argument('--null-cases', help='original calibrate-calls input, hash-bound to the frozen recipe; add per-query null tails')
    apply.add_argument('--output', required=True)
    apply.set_defaults(func=apply_main)
    return parser


def apply_calibration(frame, calibration):
    """Apply a frozen cutoff without reading query labels or refitting it.

    Input scores remain ranking evidence; only call_score represents an accepted
    call. No-score queries and rejected candidates are no_supported_match, not
    proven nonbinders. The caller must use matching scorer/reference identities.
    """
    import polars as pl
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



def null_tail_probabilities(frame, nulls, *, higher_is_better=True, score_column="score"):
    """Tail of the existing selected score against each frozen null population.

    P=(1+number of null scores at least as extreme)/(1+all null queries),
    retaining ties and no-hits. It is a conservative rank P-value if query and
    null scores are exchangeable under the same complete pipeline. E=k/N is
    the observed expected false calls per query at that score threshold, not a
    confidence bound. Maxima retain the most conservative represented source.
    Unsupported queries have P=1/E=0 because they make no call.
    """
    import numpy as np
    import polars as pl
    if not isinstance(score_column,str) or not score_column.strip():
        raise ValueError('score_column must name a numeric input column')
    if {'p_null','E_null'} & set(frame.columns):
        raise ValueError('null-tail output columns already exist')
    for part in (frame, nulls):
        if missing := {'query_id',score_column} - set(part.columns):
            raise ValueError('missing columns: '+', '.join(sorted(missing)))
        if part['query_id'].null_count() or part['query_id'].n_unique()!=part.height:
            raise ValueError('query_id must be unique and non-null; reduce candidate labels first')
        if not part[score_column].cast(pl.Float64).drop_nulls().is_finite().all():
            raise ValueError('non-null scores must be finite')
    if not nulls.height or not {'outcome','null_source'} <= set(nulls.columns):
        raise ValueError('nonempty null cases require outcome and null_source')
    if nulls['outcome'].null_count() or not nulls['outcome'].eq('null').all():
        raise ValueError('only null outcomes define null tails')
    if nulls['null_source'].null_count() or nulls['null_source'].cast(pl.String).str.strip_chars().eq('').any():
        raise ValueError('every null query requires a nonempty null_source')
    sign=1 if higher_is_better else -1
    values=sign*frame[score_column].cast(pl.Float64).to_numpy()
    missing=np.isnan(values)
    p=np.zeros(frame.height); e=np.zeros(frame.height); rows=[]
    for source in sorted(nulls['null_source'].unique().to_list()):
        part=nulls.filter(pl.col('null_source')==source)
        scores=np.sort(sign*part[score_column].cast(pl.Float64).drop_nulls().to_numpy())
        counts=len(scores)-np.searchsorted(scores,values,side='left')
        counts[missing]=0
        tail=np.where(missing,1.,(counts+1)/(part.height+1))
        rate=counts/part.height
        p=np.maximum(p,tail); e=np.maximum(e,rate)
        rows.append(frame.select('query_id').with_columns(
            pl.lit(source).alias('null_source'), pl.lit(part.height).alias('null_queries'),
            pl.Series('null_tail_count',counts),pl.Series('p_null',tail),pl.Series('E_null',rate)))
    return frame.with_columns(pl.Series('p_null',p),pl.Series('E_null',e)),pl.concat(rows)


def apply_main(args):
    import json
    from ..db.cache import sha256
    from ..io.airr import _read_table
    recipe = Path(args.calibration)
    frozen = json.loads(recipe.read_text())
    out = apply_calibration(_read_table(args.cases), frozen)
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

    frame=_read_table(args.cases)
    score_column=getattr(args,'score_column','score')
    if score_column not in frame.columns:
        raise ValueError('missing score column: '+score_column)
    result, curve = operating_points(frame.with_columns(frame[score_column].alias('score')), args.target_fpr,
                                    higher_is_better=not args.lower_is_better)
    result['score_column']=score_column
    result['input_sha256'] = sha256(Path(args.cases))
    result['limitations'] = ('Empirical rates on supplied query populations; not independent validation '
                            'or biological specificity probability. Bind source/scorer/reference identities '
                            'before recommending or applying this cutoff to another run.')
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    curve.write_csv(str(prefix) + '.operating_points.tsv', separator='\t')
    Path(str(prefix) + '.calibration.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return 0
