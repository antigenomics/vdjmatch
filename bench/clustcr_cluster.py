#!/usr/bin/env python3
"""clusTCR -> clustering benchmark predictions on the shared vdjmatch clonotype sets.

Idempotent. Run with the vdjmatch venv:

    ./.venv/bin/python bench/clustcr_cluster.py

This script (venv) builds the shared sets via bench/_cluster_common.py, ships each CDR3 list to the
`cmp-clustcr` conda env where clusTCR runs (Clustering().fit on a pandas.Series of CDR3s, default
two-step MCL params), reads back {cdr3 -> cluster}, scores with the canonical metric, and writes
bench/predictions/clustcr/clustering.tsv.

clusTCR clusters a single CDR3 column, so TRA/TRB use that chain directly; paired uses the
concatenated alpha+beta CDR3 string (clusTCR has no native paired mode).
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _cluster_common as C  # noqa: E402

ENV = "cmp-clustcr"
# worker runs inside the conda env: reads a JSON list of strings, returns {string -> cluster_id}
WORKER = r"""
import sys, json
import pandas as pd
from clustcr import Clustering
seqs = json.load(open(sys.argv[1]))
ser = pd.Series(seqs).drop_duplicates()
r = Clustering().fit(ser)
df = r.clusters_df  # columns: junction_aa, cluster
out = {str(a): int(c) for a, c in zip(df['junction_aa'], df['cluster'])}
json.dump(out, open(sys.argv[2], 'w'))
"""


def run_clustcr(seqs):
    """Return {seq -> cluster_id} for the given CDR3 list (deduped inside the worker)."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        wf = td / "worker.py"; wf.write_text(WORKER)
        inp = td / "in.json"; inp.write_text(json.dumps(seqs))
        outp = td / "out.json"
        subprocess.run(["conda", "run", "-n", ENV, "python", str(wf), str(inp), str(outp)],
                       check=True)
        return json.loads(outp.read_text())


def main():
    s = C.sets()
    rows = []
    for name in ("TRB", "TRA", "paired"):
        d = s[name]
        if name == "paired":
            seqs = [a + b for a, b in zip(d["ca"], d["cb"])]   # concatenated alpha+beta
            note = "default two-step; paired=concat(alpha,beta)"
        else:
            seqs = d["cdr3"]
            note = "default two-step"
        groups = run_clustcr(seqs)
        labels = C.labels_from_groups(seqs, groups)
        pur, ret, nc = C.score(labels, d["epi"])
        rows.append({"set": name, "macro_purity": pur, "retention": ret, "n_clusters": nc,
                     "n_clonotypes": len(seqs), "note": note})
    C.write_tsv("clustcr", rows)


if __name__ == "__main__":
    main()
