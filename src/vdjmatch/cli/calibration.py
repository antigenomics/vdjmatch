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
    return parser


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
    for source in sorted(null['null_source'].unique().to_list()):
        part = null.filter(pl.col('null_source') == source)
        calls = counts(part)
        eligible &= calls / part.height <= target_fpr
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
