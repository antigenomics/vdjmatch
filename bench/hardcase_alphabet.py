"""Integrate the LLW signal as a CHEMICAL-ALPHABET PSSM channel (the user's strong/weak-alphabet idea),
fused with NED --- a cleaner single mechanism than the ad-hoc apex-hydropathy/aromatic/charge features.

Reduce the 20 AA to chemical classes (hydrophobic / aromatic / polar / positive / negative), build a k-mer
log-odds of the query's reduced CDR3 under the epitope reference vs the OTHER-epitope binders, and rank-fuse
with NED. This is one interpretable "physicochemical PSSM" channel alongside the identity-PSSM that NED
already uses --- a dual-PSSM scheme. Held-out (sweep, 50/50) is the leakage-free test; here we sanity-check
the panel + that no cell regresses.

    .venv/bin/python bench/hardcase_alphabet.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hardcase as H                                                 # noqa: E402
import hardcase_prior3 as P3                                         # noqa: E402
from metrics import roc_auc, pr_auc_balanced                        # noqa: E402

# chemical-class alphabet: hydrophobic / aromatic / polar / positive / negative (Kosmrlj strong=H,F)
_CLASS = {}
for cls, aas in [("H", "AILMV"), ("F", "FWY"), ("P", "STNQCGP"), ("+", "KRH"), ("-", "DE")]:
    for a in aas:
        _CLASS[a] = cls


def reduce(s):
    return "".join(_CLASS.get(c, "P") for c in s)


def alpha_prior(task, locus, k=3):
    ref = [reduce(s) for s in H.FP.ref_table(task, locus)["cdr3"].to_list()]
    bg = [reduce(s) for s in P3.other_binders(task, locus)]
    kref, kbg = H.kmer_logodds(ref, k), H.kmer_logodds(bg, k)
    return lambda s: kref(reduce(s)) - kbg(reduce(s))


def ranks(x):
    return np.argsort(np.argsort(x)).astype(float)


def roc_pr(y, x):
    p = list(zip(y.tolist(), x.tolist()))
    return roc_auc(p), pr_auc_balanced(p)


if __name__ == "__main__":
    print(f"{'case':10s} {'NED':>8} {'chem-PSSM alone':>16} {'FUSED ROC/PR':>16} {'dROC':>7}")
    deltas = []
    for t, l in P3.CASES:
        d = H.FP.task_table(t, l)
        y = np.array(d["label"].to_list())
        qs = d["cdr3"].to_list()
        base = H.FP.baseline_scores(t, l)
        cd = np.array([base.get(s, 0.0) for s in qs])
        sc = alpha_prior(t, l)
        pr = np.array([sc(s) for s in qs])
        f = ranks(cd) + ranks(pr)
        nr = roc_pr(y, cd)[0]
        pa = roc_pr(y, pr)
        fr = roc_pr(y, f)
        dd = fr[0] - nr
        deltas.append(dd)
        flag = "  <-- GAIN" if dd > 0.02 else ("  REGRESS" if dd < -0.02 else "")
        print(f"{t+'-'+l:10s} {nr:>8.3f} {pa[0]:.3f}/{pa[1]:.3f}{'':>4} {fr[0]:.3f}/{fr[1]:.3f}{'':>4} "
              f"{dd:>+7.3f}{flag}")
    print(f"\nchem-alphabet PSSM fusion: mean dROC {np.mean(deltas):+.3f}  min {np.min(deltas):+.3f}")
