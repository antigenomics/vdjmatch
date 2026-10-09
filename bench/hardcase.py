"""Isolated hard-case feature analysis for the poor detection epitopes (GLC-TRA, LLW-TRB, GLC-TRB, ...).

Read-only on the library; writes nothing to committed scoring. For each (epitope, chain) it builds the
binder/decoy feature table, reports per-feature univariate AUC, and a 5-fold cross-validated logistic
"ceiling" AUC over all sequence features --- the key question: *is the binder/decoy signal recoverable
from CDR3 sequence at all*, and which features carry it. Features are built from the epitope REFERENCE
(no test-label leakage); intrinsic features (length, hydropathy, charge) are label-free.

    .venv/bin/python bench/hardcase.py GLC TRA   # one case (verbose)
    .venv/bin/python bench/hardcase.py           # the three headline cases, summary table
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _feat_probe as FP                                              # noqa: E402
from metrics import roc_auc, pr_auc_balanced                         # noqa: E402
from sklearn.linear_model import LogisticRegression                  # noqa: E402
from sklearn.model_selection import StratifiedKFold                  # noqa: E402
from sklearn.preprocessing import StandardScaler                     # noqa: E402

KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2,
      "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}
CHG = {"K": 1, "R": 1, "D": -1, "E": -1}
AROM = set("FWY")


def _apex(s, w=5):
    """central window of a CDR3 (the specificity-determining contact residues)."""
    if len(s) <= w:
        return s
    m = len(s) // 2
    return s[max(0, m - w // 2): m + (w - w // 2)]


def kmer_logodds(refs, k=3):
    """ref position-independent k-mer log-odds vs an AA-marginal background -> scorer(cdr3)->mean log-odds."""
    cnt = Counter()
    aa = Counter()
    for s in refs:
        for i in range(len(s) - k + 1):
            cnt[s[i:i + k]] += 1
        aa.update(s)
    tot = sum(cnt.values()) + 1.0
    naa = sum(aa.values()) or 1
    pa = {a: aa[a] / naa for a in aa}
    V = 20 ** k

    def bg(km):
        p = 1.0
        for c in km:
            p *= pa.get(c, 1e-3)
        return max(p, 1e-12)

    def score(s):
        vals = []
        for i in range(len(s) - k + 1):
            km = s[i:i + k]
            pr = (cnt.get(km, 0) + 1.0) / (tot + V)
            vals.append(np.log(pr / bg(km)))
        return float(np.mean(vals)) if vals else 0.0
    return score


def vj_logodds(ref_col, bg_col):
    """log P(gene|ref) / P(gene|background) -> dict gene->logodds (Laplace)."""
    r, b = Counter(ref_col), Counter(bg_col)
    nr, nb = sum(r.values()) + len(r), sum(b.values()) + len(b)
    keys = set(r) | set(b)
    return {g: np.log(((r[g] + 1) / nr) / ((b[g] + 1) / nb)) for g in keys}


def features(task, locus):
    """Return (X feature matrix, y labels, names, extras dict) for a detection case."""
    d = FP.task_table(task, locus)
    ref = FP.ref_table(task, locus)
    base = FP.baseline_scores(task, locus)
    cdr3 = d["cdr3"].to_list()
    y = np.array(d["label"].to_list())
    reflen = np.array([len(s) for s in ref["cdr3"].to_list()])
    refmode = float(np.median(reflen))
    kmer = kmer_logodds(ref["cdr3"].to_list(), k=3)
    vlo = vj_logodds(ref["v"].to_list(), d["v"].to_list())          # V usage: ref vs this query pool
    jlo = vj_logodds(ref["j"].to_list(), d["j"].to_list())
    rows, names = [], ["len", "len_dev", "kd_mean", "kd_apex", "charge", "arom", "kmer", "vlo", "jlo", "NED"]
    for s, v, j in zip(cdr3, d["v"].to_list(), d["j"].to_list()):
        ap = _apex(s)
        rows.append([
            len(s),
            -abs(len(s) - refmode),                                  # closeness to ref length mode
            np.mean([KD[c] for c in s if c in KD] or [0]),
            np.mean([KD[c] for c in ap if c in KD] or [0]),
            sum(CHG.get(c, 0) for c in s),
            sum(c in AROM for c in s) / max(len(s), 1),
            kmer(s),
            vlo.get(v, 0.0),
            jlo.get(j, 0.0),
            base.get(s, 0.0),
        ])
    X = np.array(rows, float)
    return X, y, names, {"n": len(y), "npos": int(y.sum()), "reflen_mode": refmode, "nref": ref.height}


def cv_ceiling(X, y, seed=0):
    """5-fold CV logistic-regression held-out ROC/PR over all features (the recoverable-signal ceiling)."""
    if y.sum() < 5 or (len(y) - y.sum()) < 5:
        return None
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    oof = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X[tr])
        lr = LogisticRegression(max_iter=2000, class_weight="balanced").fit(sc.transform(X[tr]), y[tr])
        oof[te] = lr.predict_proba(sc.transform(X[te]))[:, 1]
    pairs = list(zip(y.tolist(), oof.tolist()))
    return roc_auc(pairs), pr_auc_balanced(pairs)


def auc1(x, y):
    """univariate ROC-AUC of a single feature (orientation-free: report max(auc, 1-auc) + sign)."""
    a = roc_auc(list(zip(y.tolist(), x.tolist())))
    return (a, "+") if a >= 0.5 else (1 - a, "-")


def diagnose(task, locus, verbose=True):
    X, y, names, ex = features(task, locus)
    base_auc = roc_auc(list(zip(y.tolist(), X[:, names.index("NED")].tolist())))
    base_pr = pr_auc_balanced(list(zip(y.tolist(), X[:, names.index("NED")].tolist())))
    ceil = cv_ceiling(X, y)
    uni = {n: auc1(X[:, i], y) for i, n in enumerate(names)}
    if verbose:
        print(f"\n=== {task} {locus} ===  n={ex['n']} pos={ex['npos']} neg={ex['n']-ex['npos']} "
              f"ref={ex['nref']} reflen_mode={ex['reflen_mode']:.0f}")
        print(f"  NED baseline: ROC={base_auc:.3f}  balPR={base_pr:.3f}")
        if ceil:
            print(f"  CV ceiling (all feats): ROC={ceil[0]:.3f}  balPR={ceil[1]:.3f}   "
                  f"(gain {ceil[0]-base_auc:+.3f} ROC)")
        print("  univariate |AUC| (sign):")
        for n in sorted(uni, key=lambda k: -uni[k][0]):
            print(f"    {n:9s} {uni[n][0]:.3f} {uni[n][1]}")
    return dict(task=task, locus=locus, base_roc=base_auc, base_pr=base_pr,
                ceil_roc=ceil[0] if ceil else None, ceil_pr=ceil[1] if ceil else None, uni=uni, **ex)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        diagnose(sys.argv[1], sys.argv[2])
    else:
        cases = [("GLC", "TRA"), ("LLW", "TRB"), ("GLC", "TRB")]
        rows = [diagnose(t, l) for t, l in cases]
        print("\n=== summary ===")
        print(f"{'case':10s}{'base_ROC':>9}{'ceil_ROC':>9}{'gain':>7}{'top feature':>16}")
        for r in rows:
            top = max(r["uni"], key=lambda k: r["uni"][k][0])
            g = (r["ceil_roc"] - r["base_roc"]) if r["ceil_roc"] else float("nan")
            print(f"{r['task']+'-'+r['locus']:10s}{r['base_roc']:>9.3f}"
                  f"{r['ceil_roc'] or float('nan'):>9.3f}{g:>+7.3f}{top:>16}")
