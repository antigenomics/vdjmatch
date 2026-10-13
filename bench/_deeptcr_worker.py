#!/usr/bin/env python3
"""DeepTCR_U unsupervised clustering worker — runs INSIDE the cmp-deeptcr conda env.

Reads a JSON payload {beta:[...], v_beta:[...], (alpha:[...], v_alpha:[...]), name:...} from argv[1],
trains the DeepTCR variational autoencoder on the CDR3(+V) sequences, clusters the latent features
with phenograph (DeepTCR's default Cluster() backend), and writes the per-sequence cluster
assignment array (load order, ints; phenograph outliers = -1) as JSON to argv[2].

Called by bench/deeptcr_cluster.py (vdjmatch venv) via `conda run -n cmp-deeptcr`.
"""
import json
import multiprocessing
import os
import shutil
import sys
import tempfile


def _threads(payload):
    threads = payload.get("threads", 1)
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be a positive integer")
    allocated = payload.get("allocated_cpu_budget")
    if allocated is None and sys.platform.startswith("linux"):
        raise ValueError("Linux cgroup allocation requires a validated caller CPU budget")
    limits = [os.cpu_count() or 1]
    if allocated is not None:
        if isinstance(allocated, bool) or not isinstance(allocated, int) or allocated < 1:
            raise ValueError("allocated_cpu_budget must be a positive integer")
        limits.append(allocated)
    if hasattr(os, "sched_getaffinity"):
        limits.append(len(os.sched_getaffinity(0)))
    if "SLURM_CPUS_PER_TASK" in os.environ:
        allocated = int(os.environ["SLURM_CPUS_PER_TASK"])
        if allocated < 1:
            raise ValueError("SLURM_CPUS_PER_TASK must be positive")
        limits.append(allocated)
    if threads > min(limits):
        raise ValueError("threads exceeds the available CPU allocation")
    return threads


def main():
    payload = json.load(open(sys.argv[1]))
    threads = _threads(payload)
    # Apply ceilings before importing numerical runtimes in parent or spawn workers.
    for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "OMP_THREAD_LIMIT", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    os.environ["TF_NUM_INTRAOP_THREADS"] = str(threads)
    os.environ["TF_NUM_INTEROP_THREADS"] = "1"
    os.environ["NUMBA_NUM_THREADS"] = str(threads)
    import numpy as np
    import tensorflow as tf
    tf.config.threading.set_intra_op_parallelism_threads(threads)
    tf.config.threading.set_inter_op_parallelism_threads(1)
    outp = sys.argv[2]
    name = payload.get("name", "deeptcr_run")
    beta = payload["beta"]
    v_beta = payload.get("v_beta")
    alpha = payload.get("alpha")
    v_alpha = payload.get("v_alpha")

    # DeepTCR writes models/results under cwd; isolate in a temp dir so reruns are clean (idempotent).
    work = tempfile.mkdtemp(prefix="deeptcr_")
    os.chdir(work)

    from DeepTCR.DeepTCR import DeepTCR_U
    np.random.seed(0)

    dt = DeepTCR_U(name, tf_verbosity=3)
    load_kw = {"beta_sequences": np.array(beta)}
    if v_beta is not None:
        load_kw["v_beta"] = np.array(v_beta)
    if alpha is not None:
        load_kw["alpha_sequences"] = np.array(alpha)
    if v_alpha is not None:
        load_kw["v_alpha"] = np.array(v_alpha)
    # Spawn avoids forking a process after importing TensorFlow. Ordered starmap
    # remains DeepTCR's own embedding implementation through its public pool API.
    pool = multiprocessing.get_context("spawn").Pool(processes=threads)
    try:
        dt.Load_Data(p=pool, **load_kw)
    finally:
        pool.close()
        pool.join()

    # Default-ish VAE training; graph_seed/split_seed for reproducibility. stop_criterion=0.01 is the
    # library default (early-stop on reconstruction-loss plateau), so this is standard unsupervised usage.
    dt.Train_VAE(graph_seed=0, split_seed=0, suppress_output=False)

    # DeepTCR's default clustering backend. phenograph forces graph communities; isolated nodes -> -1.
    dt.Cluster(set="all", clustering_method="phenograph", n_jobs=1)
    idx = np.asarray(dt.Cluster_Assignments).astype(int).tolist()

    json.dump({"labels": idx, "n": len(idx)}, open(outp, "w"))
    shutil.rmtree(work, ignore_errors=True)
    print(f"[deeptcr-worker] {name}: {len(idx)} sequences, "
          f"{len(set(i for i in idx if i >= 0))} phenograph clusters")


if __name__ == "__main__":
    main()
