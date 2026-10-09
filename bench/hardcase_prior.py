"""Test an UNSUPERVISED reference-composition prior as a vdjmatch extension.

LLW-TRB's recoverable signal (hardcase.py) is chemical apex composition, which NED (neighbour density)
misses. A prior that scores a query by how well its CDR3 composition matches the EPITOPE REFERENCE's
composition (apex hydropathy, aromatic fraction, charge, k-mer motif) is reference-based and
epitope-general (the same procedure everywhere, like the germline prior) --- no per-epitope training, no
test-label leakage. We measure: NED alone, the prior alone, and their rank-fusion, on every case, to
check the LLW gain comes without regressing the others (single-parameter-set principle).

    .venv/bin/python bench/hardcase_prior.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hardcase as H                                                 # noqa: E402
from metrics import roc_auc, pr_auc_balanced                         # noqa: E402

CASES = [("NLV", "TRB"), ("LLW", "TRB"), ("LLL", "TRB"), ("YLQ", "TRB"), ("GLC", "TRB"),
         ("GLC", "TRA"), ("YLQ", "TRA")]


def ref_comp_prior(task, locus):
    """Fit Gaussian models of apex-KD / aromatic / charge and a k-mer log-odds from the reference;
    return scorer(cdr3) -> log-likelihood of matching the reference composition (higher = more binder-like)."""
    ref = H.FP.ref_table(task, locus)["cdr3"].to_list()
    kmer = H.kmer_logodds(ref, k=3)

    def feats(s):
        ap = H._apex(s)
        return (np.mean([H.KD[c] for c in ap if c in H.KD] or [0]),
                sum(c in H.AROM for c in s) / max(len(s), 1),
                sum(H.CHG.get(c, 0) for c in s))
    M = np.array([feats(s) for s in ref])
    mu, sd = M.mean(0), M.std(0) + 1e-6

    def score(s):
        f = np.array(feats(s))
        gll = -0.5 * np.sum(((f - mu) / sd) ** 2)                    # Gaussian log-lik of composition match
        return gll + kmer(s)
    return score


def ranks(x):
    """ranks (1..n) of values, ties averaged -> scale-free fusion key."""
    order = np.argsort(np.argsort(x))
    return order.astype(float)


def evaluate(task, locus):
    d = H.FP.task_table(task, locus)
    y = np.array(d["label"].to_list())
    base = H.FP.baseline_scores(task, locus)
    cd = np.array([base.get(s, 0.0) for s in d["cdr3"].to_list()])
    sc = ref_comp_prior(task, locus)
    pr = np.array([sc(s) for s in d["cdr3"].to_list()])
    fuse = ranks(cd) + ranks(pr)                                     # rank-sum fusion (scale-robust)

    def au(x):
        p = list(zip(y.tolist(), x.tolist()))
        return roc_auc(p), pr_auc_balanced(p)
    return {"y": y, "NED": au(cd), "prior": au(pr), "fused": au(fuse)}


if __name__ == "__main__":
    print(f"{'case':10s} {'NED ROC/PR':>16} {'prior ROC/PR':>16} {'FUSED ROC/PR':>16} {'dROC':>7}")
    deltas = []
    for t, l in CASES:
        r = evaluate(t, l)
        c, p, f = r["NED"], r["prior"], r["fused"]
        d = f[0] - c[0]
        deltas.append(d)
        flag = "  <-- gain" if d > 0.02 else ("  REGRESS" if d < -0.02 else "")
        print(f"{t+'-'+l:10s} {c[0]:.3f}/{c[1]:.3f}{'':>4} {p[0]:.3f}/{p[1]:.3f}{'':>4} "
              f"{f[0]:.3f}/{f[1]:.3f}{'':>4} {d:>+7.3f}{flag}")
    print(f"\nmean dROC {np.mean(deltas):+.3f}  min {np.min(deltas):+.3f}  "
          f"(accept only if min >= ~0 and LLW gains)")
