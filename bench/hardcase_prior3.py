"""Composition prior v3: 1-vs-OTHER-epitope contrast (not vs OLGA).

The decoys in detection are OTHER epitopes' binders, not OLGA randoms. If E's binders and the decoys share
a repertoire-level composition (e.g. LLW vs LLL, both A*02 yellow-fever-context TCRs), an OLGA-null prior
(v2) can't separate them. The discriminative contrast is log P(f | E binders) / P(f | other-epitope
binders), built entirely from the VDJdb reference database (no test-label leakage, epitope-general). This
asks what makes E's binders compositionally distinct from OTHER binders --- the actual detection contrast.

    .venv/bin/python bench/hardcase_prior3.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench                                                        # noqa: E402
import hardcase as H                                                 # noqa: E402
import hardcase_prior2 as P2                                         # noqa: E402
from metrics import roc_auc, pr_auc_balanced                         # noqa: E402
from vdjmatch import db                                              # noqa: E402

CASES = P2.CASES
_OTHER = {}


def other_binders(task, locus, n=15000):
    """pooled unique CDR3 of all human A*02 binders of OTHER epitopes (the detection contrast)."""
    key = (task, locus)
    if key not in _OTHER:
        import benchmark as B
        v26 = db.load(_bench.source(), species="HomoSapiens")
        o = (_bench.valid_cdr3(v26.filter((pl.col("gene") == locus) & pl.col("mhc_a").str.contains(B.A02)
                                          & (pl.col("epitope") != H.FP.E[task])))
             .unique("cdr3", maintain_order=True))
        _OTHER[key] = o.head(n)["cdr3"].to_list()
    return _OTHER[key]


def prior(task, locus):
    ref = H.FP.ref_table(task, locus)["cdr3"].to_list()
    bg = other_binders(task, locus)
    R = np.array([P2.comp(s) for s in ref])
    Bg = np.array([P2.comp(s) for s in bg])
    go = P2.gauss_logodds(R, Bg)
    # k-mer log-odds with the OTHER-epitope pool as background (not the AA marginal)
    kref = H.kmer_logodds(ref, k=3)
    kbg = H.kmer_logodds(bg, k=3)

    def score(s):
        return go(P2.comp(s)) + (kref(s) - kbg(s))
    return score


def ranks(x):
    return np.argsort(np.argsort(x)).astype(float)


if __name__ == "__main__":
    print(f"{'case':10s} {'NED':>14} {'prior(1-v-other)':>18} {'FUSED':>14} {'dROC':>7}")
    deltas = {}
    for t, l in CASES:
        d = H.FP.task_table(t, l)
        y = np.array(d["label"].to_list())
        base = H.FP.baseline_scores(t, l)
        cd = np.array([base.get(s, 0.0) for s in d["cdr3"].to_list()])
        sc = prior(t, l)
        pr = np.array([sc(s) for s in d["cdr3"].to_list()])
        fuse = ranks(cd) + ranks(pr)
        c, p, f = P2.au(y, cd), P2.au(y, pr), P2.au(y, fuse)
        dd = f[0] - c[0]
        deltas[f"{t}-{l}"] = dd
        flag = "  <-- GAIN" if dd > 0.02 else ("  REGRESS" if dd < -0.02 else "")
        print(f"{t+'-'+l:10s} {c[0]:.3f}/{c[1]:.3f} {p[0]:.3f}/{p[1]:.3f}{'':>4} "
              f"{f[0]:.3f}/{f[1]:.3f} {dd:>+7.3f}{flag}")
    print(f"\nmean dROC {np.mean(list(deltas.values())):+.3f}  min {min(deltas.values()):+.3f}  "
          f"LLW dROC {deltas['LLW-TRB']:+.3f}")
