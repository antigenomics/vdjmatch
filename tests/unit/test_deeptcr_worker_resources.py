"""Mock the resource boundary without loading DeepTCR or training a model."""
import importlib.util
import json
import multiprocessing
from pathlib import Path
import sys
from types import SimpleNamespace
import pytest


def test_explicit_pool_closes_before_default_training(tmp_path, monkeypatch):
    path = Path(__file__).resolve().parents[2] / "bench/_deeptcr_worker.py"
    spec = importlib.util.spec_from_file_location("resource_worker", path)
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)
    monkeypatch.setattr(worker.sys, "platform", "linux")
    with pytest.raises(ValueError, match="cgroup"):
        worker._threads({})
    assert worker._threads(dict(allocated_cpu_budget=1)) == 1
    for budget in (True, 0, 1.5):
        with pytest.raises(ValueError, match="allocated_cpu_budget"):
            worker._threads(dict(allocated_cpu_budget=budget))
    with pytest.raises(ValueError, match="allocation"):
        worker._threads(dict(threads=2, allocated_cpu_budget=1))
    monkeypatch.setattr(worker.sys, "platform", "darwin")
    monkeypatch.setattr(worker.os, "cpu_count", lambda: 4)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    assert worker._threads(dict(threads=2)) == 2
    for value in (True, 0, -1, 1.5, 3):
        with pytest.raises(ValueError):
            worker._threads(dict(threads=value))
    for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "OMP_THREAD_LIMIT", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
                "TF_NUM_INTRAOP_THREADS", "TF_NUM_INTEROP_THREADS", "NUMBA_NUM_THREADS"):
        monkeypatch.setenv(key, "1")
    events = []
    pool = SimpleNamespace(close=lambda: events.append("close"), join=lambda: events.append("join"))
    def context(kind):
        assert kind == "spawn"
        return SimpleNamespace(Pool=lambda processes: (events.append(("pool", processes)) or pool))
    monkeypatch.setattr(multiprocessing, "get_context", context)
    threading = SimpleNamespace(
        set_intra_op_parallelism_threads=lambda n: events.append(("tf_intra", n)),
        set_inter_op_parallelism_threads=lambda n: events.append(("tf_inter", n)))
    monkeypatch.setitem(sys.modules, "tensorflow", SimpleNamespace(config=SimpleNamespace(threading=threading)))
    class Model:
        def __init__(self, name, tf_verbosity):
            assert tf_verbosity == 3
        def Load_Data(self, **kwargs):
            assert kwargs["p"] is pool
            assert kwargs["beta_sequences"].tolist() == ["placeholder"]
            events.append("load")
        def Train_VAE(self, **kwargs):
            assert events[-2:] == ["close", "join"]
            assert kwargs == dict(graph_seed=0, split_seed=0, suppress_output=False)
            events.append("train")
        def Cluster(self, **kwargs):
            assert kwargs == dict(set="all", clustering_method="phenograph", n_jobs=1)
            self.Cluster_Assignments = [0]
    monkeypatch.setitem(sys.modules, "DeepTCR.DeepTCR", SimpleNamespace(DeepTCR_U=Model))
    payload, output = tmp_path / "input.json", tmp_path / "output.json"
    payload.write_text(json.dumps(dict(beta=["placeholder"])))
    monkeypatch.setattr(sys, "argv", [str(path), str(payload), str(output)])
    monkeypatch.setattr(worker.os, "chdir", lambda path: None)
    worker.main()
    assert ("pool", 1) in events
    assert ("tf_intra", 1) in events and ("tf_inter", 1) in events
    assert json.loads(output.read_text()) == dict(labels=[0], n=1)


def test_caller_validates_and_passes_allocation_before_subprocess(tmp_path, monkeypatch):
    path = Path(__file__).resolve().parents[2] / "bench/deeptcr_cluster.py"
    monkeypatch.setitem(sys.modules, "_cluster_common", SimpleNamespace())
    monkeypatch.setattr(sys, "path", sys.path.copy())
    spec = importlib.util.spec_from_file_location("resource_caller", path)
    caller = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(caller)
    from vdjtools import cores
    monkeypatch.setattr(cores, "available_cores", lambda: 2)
    calls = []
    def run(command, **kwargs):
        received = json.loads(Path(command[-2]).read_text())
        calls.append(received)
        assert received["threads"] == 1
        assert received["allocated_cpu_budget"] == 2
        Path(command[-1]).write_text(json.dumps(dict(labels=[0])))
    monkeypatch.setattr(caller.subprocess, "run", run)
    payload = dict(beta=["placeholder"], allocated_cpu_budget=999)
    assert caller.run_deeptcr(payload) == [0]
    assert payload["allocated_cpu_budget"] == 999
    for threads in (True, 0, -1, 1.5, 3):
        with pytest.raises((ValueError, TypeError)):
            caller.run_deeptcr(dict(payload, threads=threads))
    assert len(calls) == 1
