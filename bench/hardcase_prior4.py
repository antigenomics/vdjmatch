"""Composition-contrast prior: compare the BACKGROUND null --- OLGA (v2), other-epitope binders (v3),
and the real post-selection control airr_control (v4, the TCRNET/mirpy choice). Same prior everywhere:
log P(composition|ref) / P(composition|background) over apex hydropathy / aromatic / charge / length +
a k-mer log-odds, rank-fused with NED. The winning background is the one whose NED-fusion gains on the
hard cases (LLW, LLL, GLC-TRA) while regressing NO cell (min dROC >= ~0).

    .venv/bin/python bench/hardcase_prior4.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hardcase as H                                                 # noqa: E402
import hardcase_control as HC                                        # noqa: E402
import hardcase_prior2 as P2                                         # noqa: E402
import hardcase_prior3 as P3                                         # noqa: E402
from metrics import roc_auc, pr_auc_balanced                         # noqa: E402

CASES = P2.CASES


def prior_from_bg(ref, bg):
    """generic composition-contrast prior: log P(comp|ref)/P(comp|bg) + k-mer (ref vs bg) log-odds."""
    go = P2.gauss_logodds(np.array([P2.comp(s) for s in ref]), np.array([P2.comp(s) for s in bg]))
    kref, kbg = H.kmer_logodds(ref, 3), H.kmer_logodds(bg, 3)
    return lambda s: go(P2.comp(s)) + (kref(s) - kbg(s))


def ranks(x):
    return np.argsort(np.argsort(x)).astype(float)


def roc_pr(y, x):
    p = list(zip(y.tolist(), x.tolist()))
    return roc_auc(p), pr_auc_balanced(p)


def main():
    print(f"{'case':10s} {'NED':>11} | {'OLGA':>15} | {'other-bind':>15} | {'airr_control':>15}")
    print(f"{'':10s} {'ROC':>11} | {'fused dROC':>15} | {'fused dROC':>15} | {'fused dROC':>15}")
    agg = {"OLGA": [], "other": [], "control": []}
    for t, l in CASES:
        d = H.FP.task_table(t, l)
        y = np.array(d["label"].to_list())
        qs = d["cdr3"].to_list()
        base = H.FP.baseline_scores(t, l)
        cd = np.array([base.get(s, 0.0) for s in qs])
        ned_roc = roc_pr(y, cd)[0]
        ref = H.FP.ref_table(t, l)["cdr3"].to_list()
        bgs = {"OLGA": P2.olga_cdr3(l), "other": P3.other_binders(t, l), "control": HC.control_cdr3(l)}
        out = {}
        for name, bg in bgs.items():
            sc = prior_from_bg(ref, bg)
            prs = np.array([sc(s) for s in qs])
            froc, fpr = roc_pr(y, ranks(cd) + ranks(prs))
            out[name] = froc - ned_roc
            agg[name].append(froc - ned_roc)
        print(f"{t+'-'+l:10s} {ned_roc:>11.3f} | {out['OLGA']:>+15.3f} | {out['other']:>+15.3f} | "
              f"{out['control']:>+15.3f}")
    print("-" * 74)
    for name in ("OLGA", "other", "control"):
        a = np.array(agg[name])
        print(f"  {name:12s} mean dROC {a.mean():>+.3f}  min {a.min():>+.3f}  "
              f"#regress {(a < -0.02).sum()}  #gain {(a > 0.02).sum()}")


if __name__ == "__main__":
    main()
