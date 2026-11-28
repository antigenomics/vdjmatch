"""Cohort rules for matched full-catalogue VDJdb benchmarks."""

from collections import Counter
import json
from random import Random


def admitted_targets(source_counts):
    """Admit full pMHC labels with more than 30 valid source records."""
    return {target for target, count in source_counts.items() if count > 30}


def capped_query_indices(case_keys, annotations, targets, *, seed=20261128, cap=300):
    """Select physical queries without dropping any of their eligible labels.

    Call separately for each receptor species, lane and MHC class. Case keys
    must be unique tuples of the observed junctions and V calls. The supplied
    annotations and fixed reference do not change.
    """
    if len(case_keys) != len(annotations) or len(set(case_keys)) != len(case_keys):
        raise ValueError("expected one annotation set per unique physical query key")
    if not isinstance(cap, int) or isinstance(cap, bool) or cap < 1:
        raise ValueError("cap must be a positive integer")

    eligible_labels = [targets.intersection(labels) for labels in annotations]
    order = sorted(
        (i for i, labels in enumerate(eligible_labels) if labels),
        key=lambda i: json.dumps(case_keys[i]),
    )
    Random(seed).shuffle(order)

    occupancy = Counter()
    selected = []
    for i in order:
        labels = eligible_labels[i]
        if all(occupancy[label] < cap for label in labels):
            selected.append(i)
            occupancy.update(labels)
    return sorted(selected)
