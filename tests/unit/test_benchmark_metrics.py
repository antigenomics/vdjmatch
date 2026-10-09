"""Distinguish AP from historical trapezoidal PR, including ties."""
import importlib.util
from pathlib import Path
import pytest

def test_average_precision_ties_and_perfect_ranking():
    spec = importlib.util.spec_from_file_location("benchmark_metrics", Path(__file__).parents[2] / "bench/metrics.py")
    metrics = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(metrics)
    tied = [(1, 0), (0, 0), (0, 0), (0, 0)]
    assert metrics.pr_auc_balanced(tied) == 0.75
    assert metrics.average_precision(tied) == 0.25
    assert metrics.average_precision(tied, 0.5) == 0.5
    perfect = [(1, 2), (1, 2), (0, 1)]
    assert metrics.average_precision(perfect) == 1
    assert metrics.average_precision(perfect, 0.5) == 1
    with pytest.raises(ValueError):
        metrics.average_precision(tied, 1)
