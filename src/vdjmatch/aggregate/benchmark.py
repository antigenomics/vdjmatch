"""Cohort rules for matched full-catalogue VDJdb benchmarks."""

from collections import Counter
from hashlib import sha256
import json


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

    def order_key(i):
        key = json.dumps([seed, *case_keys[i]], separators=(",", ":")).encode()
        return sha256(key).digest(), key

    occupancy = Counter()
    selected = []
    for i in sorted(range(len(case_keys)), key=order_key):
        labels = targets.intersection(annotations[i])
        if labels and all(occupancy[label] < cap for label in labels):
            selected.append(i)
            occupancy.update(labels)
    return sorted(selected)
