from vdjmatch.aggregate.benchmark import admitted_targets, capped_query_indices


def test_cap_respects_all_target_memberships_and_order_permutation():
    targets = admitted_targets({"A": 31, "B": 52, "C": 30})
    assert targets == {"A", "B"}
    cases = [("one",), ("two",), ("three",), ("four",)]
    labels = [{"A", "B"}, {"A"}, {"B"}, {"C"}]
    first = capped_query_indices(cases, labels, targets, cap=1)
    second = capped_query_indices(cases[::-1], labels[::-1], targets, cap=1)
    assert {cases[i] for i in first} == {cases[::-1][i] for i in second}
    assert 3 not in first
    assert all(sum(target in labels[i] for i in first) <= 1 for target in targets)
