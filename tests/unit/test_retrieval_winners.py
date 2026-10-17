"""Packed streaming winners agree exactly with an independent dense oracle."""

import numpy as np
import polars as pl
import pytest

from vdjmatch.runner.retrieval_winners import WinnerAccumulator


def test_unordered_batches_match_dense_oracle_with_every_tie_and_nohit():
    rng = np.random.default_rng(193)
    n_queries, n_targets = 19, 17
    scorers = ['original', 'conditional']
    targets = [f'label-{i}' for i in range(n_targets)]
    accumulator = WinnerAccumulator(n_queries, targets, scorers)
    dense = np.full((2, n_queries, n_targets), -np.inf)
    for s, scorer in enumerate(scorers):
        q = rng.integers(0, n_queries - 1, 400)
        t = rng.integers(0, n_targets, 400)
        values = rng.integers(-3, 4, 400).astype(float)
        values[:3] = [np.nan, np.inf, -np.inf]
        supported = rng.random(400) > .2
        order = rng.permutation(400)
        for indices in np.array_split(order, 11):
            accumulator.consume(scorer, q[indices], t[indices], values[indices], supported[indices])
        for row, target, value, support in zip(q, t, values, supported):
            if support and np.isfinite(value):
                dense[s, row, target] = max(dense[s, row, target], value)
    parts = list(accumulator.frames(max_rows=34))
    assert all(frame.height <= 34 for frame in parts)
    actual = pl.concat(parts)
    assert actual.schema == {'query_id': pl.UInt32, 'epitope': pl.String,
                             'scorer': pl.String, 'score': pl.Float64, 'supported': pl.Boolean}
    expected = []
    for s, scorer in enumerate(scorers):
        for q in range(n_queries):
            maximum = dense[s, q].max()
            if np.isfinite(maximum):
                expected.extend((q, targets[t], scorer, maximum, True)
                                for t in np.flatnonzero(dense[s, q] == maximum))
            else:
                expected.append((q, None, scorer, None, False))
    assert actual.rows() == expected
    assert accumulator.maxima.dtype == np.float64
    assert accumulator.ties.dtype == np.uint8
    assert accumulator.ties.shape == (2, n_queries, 3)


def test_improvement_clears_old_ties_and_repeated_candidates_are_idempotent():
    accumulator = WinnerAccumulator(2, ['A', 'B', 'C'], ['s'])
    accumulator.consume('s', [0, 0], [0, 1], [1., 1.])
    accumulator.consume('s', [0, 0, 0], [2, 2, 1], [2., 2., 1.])
    accumulator.consume('s', [0, 0, 0], [0, 1, 0], [2., 2., 1.])
    accumulator.consume('s', [0, 1], [2, 0], [100., np.inf], [False, True])
    assert pl.concat(accumulator.frames()).rows() == [
        (0, 'A', 's', 2., True), (0, 'B', 's', 2., True),
        (0, 'C', 's', 2., True), (1, None, 's', None, False)]


def test_stable_string_identity_and_negative_finite_scores():
    accumulator = WinnerAccumulator(3, ['A', 'B'], ['s'])
    accumulator.consume('s', [1, 0, 1], [1, 0, 0], [-5., -2., -5.])
    parts = list(accumulator.frames(pl.Series(['second', 'first', 'missing']), max_rows=2))
    assert pl.concat(parts).rows() == [
        ('second', 'A', 's', -2., True), ('first', 'A', 's', -5., True),
        ('first', 'B', 's', -5., True), ('missing', None, 's', None, False)]
    assert all(frame.height <= 2 for frame in parts)


@pytest.mark.parametrize('kwargs', [
    {'n_queries': -1}, {'n_queries': 1.5}, {'n_queries': True}, {'n_queries': 2**32 + 1},
    {'targets': []}, {'targets': ['A', 'A']}, {'targets': ['A', None]},
    {'targets': ['']}, {'scorers': []}, {'scorers': ['s', 's']},
    {'max_bytes': 0}, {'max_bytes': True}, {'max_bytes': 8},
])
def test_invalid_shapes_names_and_budgets(kwargs):
    arguments = {'n_queries': 1, 'targets': ['A'], 'scorers': ['s']}
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        WinnerAccumulator(**arguments)


def test_memory_guard_counts_packed_bits_and_maxima():
    accumulator = WinnerAccumulator(2, [str(i) for i in range(9)], ['s'], max_bytes=20)
    assert accumulator.maxima.nbytes + accumulator.ties.nbytes == 20
    with pytest.raises(ValueError, match='budget'):
        WinnerAccumulator(2, [str(i) for i in range(9)], ['s'], max_bytes=19)
    with pytest.raises(ValueError, match='65536'):
        WinnerAccumulator(0, [str(i) for i in range(65537)], ['s'])


@pytest.mark.parametrize('q,t,values,support', [
    ([2], [0], [1.], None), ([-1], [0], [1.], None),
    ([0], [2], [1.], None), ([0], [-1], [1.], None),
    ([0.], [0], [1.], None), ([True], [0], [1.], None),
    ([0], [0.], [1.], None), ([[0]], [[0]], [[1.]], None),
    ([0], [0, 1], [1.], None), ([0], [0], ['bad'], None),
    ([0], [0], [1.], [1]), ([0], [0], [1.], [True, False]),
])
def test_invalid_candidate_arrays_leave_state_unchanged(q, t, values, support):
    accumulator = WinnerAccumulator(2, ['A', 'B'], ['s'])
    with pytest.raises(ValueError):
        accumulator.consume('s', q, t, values, support)
    assert not accumulator.ties.any()
    assert np.isneginf(accumulator.maxima).all()


def test_invalid_scorer_export_bounds_and_identities():
    accumulator = WinnerAccumulator(2, ['A', 'B'], ['s'])
    with pytest.raises(ValueError, match='scorer'):
        accumulator.consume('absent', [0], [0], [1.])
    for ids in [['same', 'same'], ['only-one'], ['valid', None]]:
        with pytest.raises(ValueError, match='query_ids'):
            list(accumulator.frames(ids))
    for max_rows in [1, 0, 1.5, True]:
        with pytest.raises(ValueError, match='max_rows'):
            list(accumulator.frames(max_rows=max_rows))


def test_empty_queries_and_empty_batch():
    accumulator = WinnerAccumulator(0, ['A'], ['s'])
    accumulator.consume('s', np.array([], dtype=int), np.array([], dtype=int), np.array([]))
    assert list(accumulator.frames()) == []
    assert list(accumulator.frames([])) == []
