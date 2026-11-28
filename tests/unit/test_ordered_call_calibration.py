import math
import sys

import polars as pl
import pytest

from vdjmatch.cli.calibration import apply_calibration, null_tail_probabilities
from vdjmatch.cli.ordered_calls import fit_ordered_calls


@pytest.mark.parametrize('sign', [1, -1])
def test_ordered_cutoff_excludes_boundary_ties_and_accepts_the_next_float(sign):
    nulls = pl.DataFrame({
        'query_id': range(10),
        'outcome': ['null'] * 10,
        'null_source': ['generated'] * 10,
        'score': [10, 9, 8, 8, 7, None, None, None, None, None],
    }).with_columns(pl.col('score') * sign)
    recipe = fit_ordered_calls(nulls, .5, .75, higher_is_better=sign == 1)
    assert recipe['allowed_calibration_calls'] == 3
    assert recipe['empirical_E_cutoff'] == .3
    assert recipe['binomial_violation_bound'] == pytest.approx(.171875)
    assert recipe['strict_boundary'] == sign * 8
    assert recipe['nulls'] == [{'source': 'generated', 'queries': 10, 'calls': 2, 'call_rate': .2}]
    adjacent = math.nextafter(sign * 8, sign * math.inf)
    queries = pl.DataFrame({'query_id': [0, 1, 2], 'score': [sign * 8.0, adjacent, None]})
    applied = apply_calibration(queries, recipe)
    assert applied['call_score'].to_list() == [None, adjacent, None]
    tails, _ = null_tail_probabilities(queries, nulls, higher_is_better=sign == 1)
    accepted_on_e = queries['score'].is_not_null() & (tails['E_null'] <= recipe['empirical_E_cutoff'])
    assert accepted_on_e.to_list() == applied['call_score'].is_not_null().to_list()


def test_ordered_abstention_boundaries_and_input_contract():
    nulls = pl.DataFrame({
        'query_id': range(5), 'outcome': ['null'] * 5, 'null_source': ['generated'] * 5,
        'score': [None] * 5,
    })
    recipes = [fit_ordered_calls(nulls, .01, .95), fit_ordered_calls(nulls, .5, .95)]
    assert recipes[0]['allowed_calibration_calls'] is None
    assert recipes[0]['empirical_E_cutoff'] is None
    assert recipes[0]['action'] == 'abstain_all'
    assert recipes[1]['allowed_calibration_calls'] == 0
    assert recipes[1]['boundary_at_abstention']
    assert recipes[1]['cutoff'] == -sys.float_info.max

    extreme = nulls.with_columns(pl.lit(sys.float_info.max).alias('score'))
    frozen = fit_ordered_calls(extreme, .5, .95)
    assert frozen['action'] == 'abstain_all'
    assert frozen['strict_boundary'] == sys.float_info.max
    assert frozen['binomial_violation_bound'] == pytest.approx(.03125)
    # Subtracting this confidence from one rounds to one in binary64.
    assert fit_ordered_calls(nulls.head(1), 1e-30, 1e-20)['allowed_calibration_calls'] is None
    with pytest.raises(ValueError, match='one declared null source'):
        fit_ordered_calls(nulls.with_columns(pl.col('query_id').cast(pl.String).alias('null_source')), .1, .95)
    with pytest.raises(ValueError, match='query_id must be unique'):
        fit_ordered_calls(nulls.with_columns(pl.lit(0).alias('query_id')), .1, .95)
