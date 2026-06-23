#!/usr/bin/env python3
"""DeepTCR (unsupervised) -> clustering benchmark predictions on the shared vdjmatch clonotype sets.

Idempotent. Run with the vdjmatch venv:

    ./.venv/bin/python bench/deeptcr_cluster.py

This script (venv) builds the shared sets and ships each to the `cmp-deeptcr` conda env, where
bench/_deeptcr_worker.py trains the DeepTCR_U variational autoencoder on the CDR3-beta (+V_beta)
sequences and clusters its latent features with phenograph (DeepTCR's default Cluster() backend).
Cluster assignments come back aligned to load order (phenograph outliers = -1 -> scored as
singletons), then scored with the canonical metric -> bench/predictions/deeptcr/clustering.tsv.

DeepTCR is the heaviest tool (TensorFlow VAE). It is time-boxed; if a set's VAE does not finish in
TIMEOUT seconds it is recorded as a failure and skipped. TRB is the headline; TRA and paired are run
if time allows (paired uses both alpha+beta chains, which DeepTCR natively supports).
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _cluster_common as C  # noqa: E402

ENV = "cmp-deeptcr"
WORKER = Path(__file__).resolve().parent / "_deeptcr_worker.py"
TIMEOUT = 1500  # seconds per set (hard cap; best-effort)


def run_deeptcr(payload):
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        inp = td / "in.json"; inp.write_text(json.dumps(payload))
        outp = td / "out.json"
        subprocess.run(["conda", "run", "-n", ENV, "python", str(WORKER), str(inp), str(outp)],
                       check=True, timeout=TIMEOUT)
        return json.loads(outp.read_text())["labels"]


def main():
    s = C.sets()
    rows = []
    plan = [
        ("TRB", {"beta": s["TRB"]["cdr3"], "v_beta": [v + "*01" for v in s["TRB"]["v"]],
                 "name": "vdjmatch_TRB"}, s["TRB"]["epi"], "VAE + phenograph (default), beta+V"),
        ("TRA", {"beta": s["TRA"]["cdr3"], "v_beta": [v + "*01" for v in s["TRA"]["v"]],
                 "name": "vdjmatch_TRA"}, s["TRA"]["epi"], "VAE + phenograph (default), alpha+V as beta slot"),
        ("paired", {"beta": s["paired"]["cb"], "v_beta": [v + "*01" for v in s["paired"]["vb"]],
                    "alpha": s["paired"]["ca"], "v_alpha": [v + "*01" for v in s["paired"]["va"]],
                    "name": "vdjmatch_paired"}, s["paired"]["epi"],
         "VAE + phenograph (default), paired alpha+beta"),
    ]
    for name, payload, epi, note in plan:
        try:
            labels = run_deeptcr(payload)
        except subprocess.TimeoutExpired:
            print(f"[deeptcr] {name}: TIMEOUT after {TIMEOUT}s -> skipped")
            rows.append({"set": name, "macro_purity": "NA", "retention": "NA", "n_clusters": "NA",
                         "n_clonotypes": len(epi), "note": f"TIMEOUT >{TIMEOUT}s ({note})"})
            continue
        except subprocess.CalledProcessError as e:
            print(f"[deeptcr] {name}: FAILED ({e}) -> skipped")
            rows.append({"set": name, "macro_purity": "NA", "retention": "NA", "n_clusters": "NA",
                         "n_clonotypes": len(epi), "note": f"FAILED ({note})"})
            continue
        pur, ret, nc = C.score(labels, epi)
        rows.append({"set": name, "macro_purity": pur, "retention": ret, "n_clusters": nc,
                     "n_clonotypes": len(epi), "note": note})
    C.write_tsv("deeptcr", rows)


if __name__ == "__main__":
    main()
