"""Hand-enumerable checks of directional conditional randomization tails."""
import polars as pl
import pytest

from vdjmatch.evalue.directional import directional_tails


def test_ties_nohits_and_complete_denominators():
    observed = pl.DataFrame({'query_id': ['a', 'b'], 'group_id': [0, 0], 'score': [2.0, None]})
    reference = pl.DataFrame({'query_id': ['a'] * 3 + ['b'] * 3, 'group_id': [0] * 6,
        'replicate': [0, 1, 2] * 2, 'score': [2.0, 1.0, None, None, 1.0, 2.0]})
    query = reference.with_columns(pl.when(pl.col('query_id') == 'a').then(3.0).otherwise(pl.col('score')).alias('score'))
    result = directional_tails(observed, reference, query)
    assert result['reference_draws'].to_list() == [3, 3]
    assert result['reference_exceedances'].to_list() == [1, 3]
    assert result['p_reference'].to_list() == [0.5, 1.0]
    assert result['p_query'].to_list() == [1.0, 1.0]
    assert result['p_both'].to_list() == [1.0, 1.0]


def test_maximum_needs_no_independence_or_validity_of_other_null():
    # Four equally likely outcomes give a valid uniform marginal reference tail.
    observed = pl.DataFrame({'query_id': list('abcd'), 'group_id': [0] * 4, 'score': [4.0, 3.0, 2.0, 1.0]})
    null = pl.DataFrame({'query_id': [q for q in 'abcd' for _ in range(3)],
        'group_id': [0] * 12, 'replicate': list(range(3)) * 4,
        'score': [1.0, 2.0, 3.0] * 4})
    # Opposite ordering makes the two tails completely dependent.
    other = null.with_columns(pl.Series('score', [4.0] * 3 + [4.0, 4.0, 0.0] + [4.0, 0.0, 0.0] + [0.0] * 3))
    result = directional_tails(observed, null, other)
    assert result['p_reference'].to_list() == [0.25, 0.5, 0.75, 1.0]
    assert result['p_query'].to_list() == [1.0, 0.75, 0.5, 0.25]
    assert result['p_both'].to_list() == [1.0, 0.75, 0.75, 1.0]
    # Even an arbitrary small second tail cannot invalidate the valid first tail.
    arbitrary = null.with_columns(pl.lit(0.0).alias('score'))
    result = directional_tails(observed, null, arbitrary)
    for alpha in [0.25, 0.5, 0.75, 1.0]:
        assert result.filter(pl.col('p_both') <= alpha).height / result.height <= alpha
    swapped = directional_tails(observed, arbitrary, null)
    assert result['p_both'].to_list() == swapped['p_both'].to_list()


def test_invalid_or_incomplete_calibration_fails():
    observed = pl.DataFrame({'query_id': ['a', 'b'], 'group_id': [0, 0], 'score': [1.0, 2.0]})
    valid = observed.with_columns(pl.lit(0).alias('replicate'))
    invalids = [valid.head(1), pl.concat([valid, valid.head(1)]),
        valid.with_columns(pl.lit(1).alias('group_id')),
        valid.with_columns(pl.Series('replicate', [0, 1])),
        valid.with_columns(pl.lit(float('inf')).alias('score')),
        valid.with_columns(pl.lit(float('nan')).alias('score'))]
    for invalid in invalids:
        with pytest.raises(ValueError):
            directional_tails(observed, invalid, valid)
    with pytest.raises(ValueError, match='unique'):
        directional_tails(pl.concat([observed, observed.head(1)]), valid, valid)


def test_exhaustive_reference_tails_equal_hypergeometric_probability():
    from itertools import combinations
    from math import comb

    # Five equally weighted control keys; two match the query. Enumerate every
    # size-two reference, rather than approximating this known null by simulation.
    scores = [float(len(set(reference) & {0, 1})) for reference in combinations(range(5), 2)]
    observed = pl.DataFrame({'query_id': range(len(scores)), 'group_id': [0] * len(scores), 'score': scores})
    rows = [(case, 0, replicate, value) for case in range(len(scores))
            for replicate, value in enumerate(scores[:case] + scores[case + 1:])]
    null = pl.DataFrame(rows, schema=['query_id', 'group_id', 'replicate', 'score'], orient='row')
    result = directional_tails(observed, null, null)
    any_match = 1 - comb(5 - 2, 2) / comb(5, 2)
    assert any_match == 0.7
    assert result.filter(pl.col('score') == 1)['p_reference'].to_list() == [any_match] * 6
    assert result.filter(pl.col('score') == 2)['p_reference'].to_list() == [0.1]
    assert result.filter(pl.col('score') == 0)['p_reference'].to_list() == [1.0] * 3
