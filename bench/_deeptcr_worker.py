#!/usr/bin/env python3
"""DeepTCR_U unsupervised clustering worker — runs INSIDE the cmp-deeptcr conda env.

Reads a JSON payload {beta:[...], v_beta:[...], (alpha:[...], v_alpha:[...]), name:...} from argv[1],
trains the DeepTCR variational autoencoder on the CDR3(+V) sequences, clusters the latent features
with phenograph (DeepTCR's default Cluster() backend), and writes the per-sequence cluster
assignment array (load order, ints; phenograph outliers = -1) as JSON to argv[2].

Called by bench/deeptcr_cluster.py (vdjmatch venv) via `conda run -n cmp-deeptcr`.
"""
import json
import os
import shutil
import sys
import tempfile

import numpy as np


def main():
    payload = json.load(open(sys.argv[1]))
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
    dt.Load_Data(**load_kw)

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
