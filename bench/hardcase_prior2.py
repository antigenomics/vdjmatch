"""Composition prior v2: ref-vs-OLGA-background LOG-ODDS (directional), not Gaussian ref-match.

v1 (hardcase_prior.py) failed on LLW because "closeness to the reference mean" has no direction: a decoy
sitting at the ref mean scores as high as a binder. The real signal is monotonic (LLW binders are more
polar/aromatic at the apex than background). The fix: score each composition feature by the log-likelihood
RATIO log P(f|ref)/P(f|OLGA), so a feature shifted away from background gets a signed, directional weight
--- the same null (OLGA) vdjmatch already uses for the E-value. Still unsupervised + epitope-general.

    .venv/bin/python bench/hardcase_prior2.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hardcase as H                                                 # noqa: E402
from compare import TESTDATA                                         # noqa: E402
from metrics import roc_auc, pr_auc_balanced                         # noqa: E402

CASES = [("NLV", "TRB"), ("LLW", "TRB"), ("LLL", "TRB"), ("YLQ", "TRB"), ("GLC", "TRB"),
         ("GLC", "TRA"), ("YLQ", "TRA")]
_OLGA = {}


def olga_cdr3(locus, n=20000):
    if locus not in _OLGA:
        f = "sample4_olga_airr.txt" if locus == "TRB" else "sample5_olga_airr.txt"
        s = (pl.read_csv(TESTDATA / f, separator="\t", infer_schema_length=0)
             .select(cdr3="junction_aa").drop_nulls().unique("cdr3", maintain_order=True))
        _OLGA[locus] = s.head(n)["cdr3"].to_list()
    return _OLGA[locus]


def comp(s):
    ap = H._apex(s)
    return np.array([np.mean([H.KD[c] for c in ap if c in H.KD] or [0]),     # apex hydropathy
                     np.mean([H.KD[c] for c in s if c in H.KD] or [0]),       # whole hydropathy
                     sum(c in H.AROM for c in s) / max(len(s), 1),            # aromatic frac
                     sum(H.CHG.get(c, 0) for c in s),                         # net charge
                     len(s)])                                                 # length


def gauss_logodds(refvals, bgvals):
    """per-feature scorer: log N(f;ref) - log N(f;bg)  (directional ref-vs-background)."""
    mr, sr = refvals.mean(0), refvals.std(0) + 1e-6
    mb, sb = bgvals.mean(0), bgvals.std(0) + 1e-6

    def score(f):
        lr = -0.5 * ((f - mr) / sr) ** 2 - np.log(sr)
        lb = -0.5 * ((f - mb) / sb) ** 2 - np.log(sb)
        return float(np.sum(lr - lb))
    return score


def prior(task, locus):
    ref = H.FP.ref_table(task, locus)["cdr3"].to_list()
    bg = olga_cdr3(locus)
    R = np.array([comp(s) for s in ref])
    B = np.array([comp(s) for s in bg])
    go = gauss_logodds(R, B)
    kmer = H.kmer_logodds(ref, k=3)                                  # already a ref-vs-AA-marginal log-odds

    def score(s):
        return go(comp(s)) + kmer(s)
    return score


def ranks(x):
    return np.argsort(np.argsort(x)).astype(float)


def au(y, x):
    p = list(zip(y.tolist(), x.tolist()))
    return roc_auc(p), pr_auc_balanced(p)


if __name__ == "__main__":
    print(f"{'case':10s} {'NED':>14} {'prior(logodds)':>16} {'FUSED':>14} {'dROC':>7}")
    deltas = {}
    for t, l in CASES:
        d = H.FP.task_table(t, l)
        y = np.array(d["label"].to_list())
        base = H.FP.baseline_scores(t, l)
        cd = np.array([base.get(s, 0.0) for s in d["cdr3"].to_list()])
        sc = prior(t, l)
        pr = np.array([sc(s) for s in d["cdr3"].to_list()])
        fuse = ranks(cd) + ranks(pr)
        c, p, f = au(y, cd), au(y, pr), au(y, fuse)
        dd = f[0] - c[0]
        deltas[f"{t}-{l}"] = dd
        flag = "  <-- GAIN" if dd > 0.02 else ("  REGRESS" if dd < -0.02 else "")
        print(f"{t+'-'+l:10s} {c[0]:.3f}/{c[1]:.3f} {p[0]:.3f}/{p[1]:.3f}{'':>2} "
              f"{f[0]:.3f}/{f[1]:.3f} {dd:>+7.3f}{flag}")
    print(f"\nmean dROC {np.mean(list(deltas.values())):+.3f}  min {min(deltas.values()):+.3f}  "
          f"LLW dROC {deltas['LLW-TRB']:+.3f}")
