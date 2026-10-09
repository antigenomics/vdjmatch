"""Standalone TCRGP predictor (runs in cmp-tcrgp; gpflow/TF1). Predicts P(binds epitope) for cdr3b rows.

    conda run -n cmp-tcrgp python bench/_tcrgp_predict.py <input.csv> <model_path> <out.npy>
input.csv must have a 'cdr3b' column. Writes an (n,) array of probabilities (same row order).
"""
import pickle
import sys
from pathlib import Path

import numpy as np

REPO = Path.home() / "vcs/code/TCRGP"
sys.path.insert(0, str(REPO))
import tcrgp                                                         # noqa: E402

inp, model_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
params = pickle.load(open(model_path, "rb"), encoding="bytes")
p = tcrgp.predict(inp, params, cdr3b="cdr3b", alphabet_db_file_path=str(REPO / "data" / "alphabeta_db.tsv"))
np.save(out, np.asarray(p).ravel())
print(f"[tcrgp] {model_path} -> {out}  (n={len(np.asarray(p).ravel())})", file=sys.stderr)
