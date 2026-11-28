"""Fit a strict whole-query cutoff on one independent null sample."""
from pathlib import Path
import math
import sys


def register(subparsers):
    parser = subparsers.add_parser(
        'calibrate-ordered-calls', help='fit a whole-query cutoff with order-statistic false-call control')
    parser.add_argument('cases', help='TSV: query_id, score, outcome=null, null_source')
    parser.add_argument('--target-fpr', type=float, required=True)
    parser.add_argument('--confidence', type=float, required=True)
    parser.add_argument('--lower-is-better', action='store_true')
    parser.add_argument('--score-column', default='score')
    parser.add_argument('--input-manifest', help='optional score producer manifest to bind to the recipe')
    parser.add_argument('--output', required=True, help='calibration recipe JSON')
    parser.set_defaults(func=main)


def fit_ordered_calls(frame, target_fpr, confidence, *, higher_is_better=True,
                      score_column='score'):
    """Use the largest allowed binomial count, then reject the boundary tie."""
    import numpy as np
    import polars as pl
    from scipy.stats import binom
    from .calibration import apply_calibration

    if not math.isfinite(target_fpr) or not 0 <= target_fpr < 1:
        raise ValueError('target-fpr must be finite and in [0, 1)')
    if not math.isfinite(confidence) or not 0 < confidence < 1:
        raise ValueError('confidence must be finite and in (0, 1)')
    required = {'query_id', score_column, 'outcome', 'null_source'}
    if missing := required - set(frame.columns):
        raise ValueError('missing columns: ' + ', '.join(sorted(missing)))
    if not frame.height or not (frame['outcome'] == 'null').fill_null(False).all():
        raise ValueError('calibration requires at least one query and only null outcomes')
    sources = frame['null_source'].cast(pl.String).unique()
    if len(sources) != 1 or not sources[0] or not sources[0].strip():
        raise ValueError('calibration requires one declared null source')

    size = frame.height
    rejected_risk = 1 - confidence
    low, high = -1, size
    while high - low > 1:
        middle = (low + high) // 2
        # Evaluate the smaller tail to avoid rounding probabilities near one.
        if confidence <= 0.5:
            qualifies = binom.sf(middle, size, target_fpr) >= confidence
        else:
            qualifies = binom.cdf(middle, size, target_fpr) <= rejected_risk
        if qualifies:
            low = middle
        else:
            high = middle
    allowed = low if low >= 0 else None

    scores = frame[score_column].cast(pl.Float64).drop_nulls().to_numpy()
    if not np.isfinite(scores).all():
        raise ValueError('non-null scores must be finite; use null for no evidence')
    sign = 1 if higher_is_better else -1
    boundary, cutoff = None, None
    at_abstention = allowed is not None and allowed >= len(scores)
    if allowed is not None:
        if at_abstention:
            cutoff = -sign * sys.float_info.max
        else:
            index = len(scores) - allowed - 1
            boundary = sign * float(np.partition(sign * scores, index)[index])
            first_accepted = math.nextafter(boundary, sign * math.inf)
            if math.isfinite(first_accepted):
                cutoff = first_accepted
    recipe = dict(
        status='ordered_null_calibration', score_column=score_column,
        action='apply_cutoff' if cutoff is not None else 'abstain_all',
        comparison='>=' if higher_is_better else '<=', cutoff=cutoff,
        strict_boundary=boundary, boundary_at_abstention=at_abstention,
        target_null_query_call_rate=target_fpr, confidence=confidence,
        allowed_calibration_calls=allowed,
        binomial_violation_bound=(0.0 if allowed is None else
                                  float(binom.cdf(allowed, size, target_fpr))),
        binomial_confidence_bound=(1.0 if allowed is None else
                                   float(binom.sf(allowed, size, target_fpr))),
        limitations=('The bound assumes IID draws from the declared null and a complete '
                     'score, reference and winning-label rule fixed before calibration.'),
    )
    applied = apply_calibration(frame, recipe)
    calls = int((applied['call_status'] == 'accepted').sum())
    recipe['nulls'] = [dict(source=sources[0], queries=size, calls=calls, call_rate=calls / size)]
    return recipe


def main(args):
    import json
    from .calibration import _bound_score_columns, pipeline_binding
    from ..db.cache import sha256
    from ..io.airr import _read_table

    cases, output = Path(args.cases), Path(args.output)
    manifest = Path(args.input_manifest) if args.input_manifest else None
    inputs = {cases.resolve()}
    if manifest is not None:
        inputs.add(manifest.resolve())
    if output.resolve() in inputs:
        raise ValueError('output must not overwrite inputs')
    recipe = fit_ordered_calls(_read_table(cases), args.target_fpr, args.confidence,
                               higher_is_better=not args.lower_is_better,
                               score_column=args.score_column)
    recipe['input_sha256'] = sha256(cases)
    if manifest is not None:
        binding = pipeline_binding(json.loads(manifest.read_text()))
        if args.lower_is_better or args.score_column not in _bound_score_columns(binding):
            raise ValueError('pipeline binding requires its final higher-is-better kernel score')
        recipe['pipeline_binding'] = binding
        recipe['input_manifest_sha256'] = sha256(manifest)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(recipe, indent=2, allow_nan=False) + '\n')
    return 0
