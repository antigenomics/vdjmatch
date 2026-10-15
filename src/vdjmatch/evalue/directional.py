"""Two directional empirical tails for a frozen score and testing scope."""
from __future__ import annotations

import polars as pl

_KEYS = ['query_id', 'group_id']


def _validate_scores(frame: pl.DataFrame) -> pl.DataFrame:
    if not {'score', *_KEYS}.issubset(frame.columns):
        raise ValueError('tail inputs require query_id, group_id and score')
    if any(frame[key].null_count() for key in _KEYS):
        raise ValueError('tail case keys cannot be missing')
    if not frame.schema['score'].is_numeric() and frame.schema['score'] != pl.Null:
        raise ValueError('tail scores must be numeric or null for no-hit')
    frame = frame.with_columns(pl.col('score').cast(pl.Float64))
    if frame.filter(pl.col('score').is_not_null() & ~pl.col('score').is_finite()).height:
        raise ValueError('tail scores must be finite or null for no-hit')
    return frame


def _tail(observed: pl.DataFrame, null: pl.DataFrame, direction: str) -> pl.DataFrame:
    null = _validate_scores(null)
    if 'replicate' not in null.columns or not null.height or null['replicate'].null_count():
        raise ValueError(f'{direction} null requires nonempty, non-null replicate identities')
    if null.select(*_KEYS, 'replicate').is_duplicated().any():
        raise ValueError(f'{direction} null contains duplicate case/replicate rows')
    cases = observed.select(_KEYS)
    supplied = null.select(_KEYS).unique()
    if supplied.join(cases, on=_KEYS, how='anti').height or cases.join(supplied, on=_KEYS, how='anti').height:
        raise ValueError(f'{direction} null case keys must match observed cases exactly')
    if null.group_by('replicate').len().filter(pl.col('len') != observed.height).height:
        raise ValueError(f'{direction} null must include every case for every replicate')
    joined = null.select(*_KEYS, 'score').join(
        observed.select(*_KEYS, pl.col('score').alias('_observed')), on=_KEYS, validate='m:1')
    counts = joined.group_by(_KEYS).agg(
        pl.len().cast(pl.Int64).alias(f'{direction}_draws'),
        (pl.col('_observed').is_null() | (pl.col('score') >= pl.col('_observed')).fill_null(False))
        .sum().cast(pl.Int64).alias(f'{direction}_exceedances'))
    return counts.with_columns(
        ((pl.col(f'{direction}_exceedances') + 1) / (pl.col(f'{direction}_draws') + 1)).alias(f'p_{direction}'))


def directional_tails(observed: pl.DataFrame, reference_null: pl.DataFrame,
                      query_null: pl.DataFrame) -> pl.DataFrame:
    """Return reference/query randomization tails and their maximum per case.

    Frames have ``query_id,group_id,score``; null frames additionally have a
    ``replicate`` identity and a complete row for every observed case on every
    draw. Higher scores are better, ties count as exceedances, and null no-hits
    (``score=None``) stay in the denominator. An observed no-hit has tail one.
    Each tail is ``(1 + exceedances)/(1 + draws)``. The returned ``p_both`` is
    their maximum, a test against the union of the two null hypotheses whenever
    each applicable marginal test is valid; it requires no independence.

    The caller must freeze the same statistic and scope across observed and
    randomized scores and own the preselection event. A component tail cannot
    be reused for a selected winning group or another score. These tails are
    conditional randomization tests under exchangeability, not binding evidence
    probabilities or posterior probabilities.
    """
    reserved = {'p_reference', 'p_query', 'p_both', 'reference_exceedances',
                'query_exceedances', 'reference_draws', 'query_draws'}
    if reserved.intersection(observed.columns):
        raise ValueError('observed contains reserved tail output columns')
    normalized = _validate_scores(observed)
    if not normalized.height or normalized.select(_KEYS).is_duplicated().any():
        raise ValueError('observed requires nonempty unique case keys')
    reference = _tail(normalized, reference_null, 'reference')
    query = _tail(normalized, query_null, 'query')
    return observed.join(reference, on=_KEYS, validate='1:1', maintain_order='left').join(
        query, on=_KEYS, validate='1:1', maintain_order='left').with_columns(
            pl.max_horizontal('p_reference', 'p_query').alias('p_both'))
