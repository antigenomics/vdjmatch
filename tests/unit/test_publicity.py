"""Publicity uses substitutions only, inclusive exact hits, and unique control."""
import pytest

from vdjmatch.aggregate.publicity import junction_publicity


def test_batched_publicity_preserves_query_order_and_counts():
    queries = ['CASSF', 'CASAF', 'CAGGF', 'CASSGF']
    control = ['CASSF', 'CASSF', 'CASAF', 'CASGF', 'CASSGF']
    one = junction_publicity(queries, control, threads=1)
    assert one.equals(junction_publicity(queries, control, threads=2))
    assert one['cdr3'].to_list() == queries
    assert one['control_exact_junctions'].to_list() == [1, 1, 0, 1]
    assert one['control_hamming1_junctions'].to_list() == [3, 3, 1, 1]
    assert one['control_junctions'].to_list() == [4]*4
    assert junction_publicity([], control).is_empty()
    for bad in ([], ['CAXF'], [None]):
        with pytest.raises(ValueError):
            junction_publicity(queries, bad)
