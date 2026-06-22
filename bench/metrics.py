"""Imbalance-robust retrieval metrics for the VDJdb benchmark.

Epitope classes span ~100x in size (>=30 to the 3000 cap), so micro-pooled precision/PR-AUC/F1/purity
are dominated by the large classes and are not comparable across epitopes. We cannot physically balance
the set (up-sampling a small epitope would fabricate neighbours and distort the distance distribution),
so we re-normalize the *metric* to a fixed reference prevalence ``pi0`` with the standard formula

    precision(pi0) = TPR*pi0 / (TPR*pi0 + FPR*(1-pi0)),     TPR = TP/P,  FPR = FP/N

evaluated along the ranked list. At ``pi0 = 0.5`` this is the balanced precision TPR/(TPR+FPR). ROC-AUC
is already prevalence-invariant (rank-only) and needs no correction; report it as the discriminator.
Pair these per-epitope metrics with macro-averaging across epitopes.
"""
from __future__ import annotations


def _sorted(pairs):
    s = sorted(pairs, key=lambda x: -x[1])
    P = sum(l for l, _ in s)
    return s, P, len(s) - P


def pr_auc(pairs) -> float:
    """Micro PR-AUC at the observed prevalence (average precision)."""
    s, P, _ = _sorted(pairs)
    if P == 0:
        return float("nan")
    tp = fp = 0
    pr = pp = area = 0.0
    pp = 1.0
    for lab, _ in s:
        tp += lab
        fp += 1 - lab
        r, p = tp / P, tp / (tp + fp)
        area += (r - pr) * (p + pp) / 2
        pr, pp = r, p
    return area


def roc_auc(pairs) -> float:
    """ROC-AUC (prevalence-invariant rank discriminator)."""
    s, P, N = _sorted(pairs)
    if P == 0 or N == 0:
        return float("nan")
    tp = fp = 0
    pfpr = ptpr = area = 0.0
    for lab, _ in s:
        tp += lab
        fp += 1 - lab
        tpr, fpr = tp / P, fp / N
        area += (fpr - pfpr) * (tpr + ptpr) / 2
        pfpr, ptpr = fpr, tpr
    return area


def pr_auc_balanced(pairs, pi0: float = 0.5) -> float:
    """PR-AUC with precision re-normalized to reference prevalence ``pi0`` (balanced by default)."""
    s, P, N = _sorted(pairs)
    if P == 0 or N == 0:
        return float("nan")
    tp = fp = 0
    pr = area = 0.0
    pp = 1.0
    for lab, _ in s:
        tp += lab
        fp += 1 - lab
        tpr, fpr = tp / P, fp / N
        denom = tpr * pi0 + fpr * (1 - pi0)
        p0 = (tpr * pi0) / denom if denom else 1.0
        area += (tpr - pr) * (p0 + pp) / 2
        pr, pp = tpr, p0
    return area


def f1_balanced(pairs, pi0: float = 0.5) -> float:
    """Best F1 along the ranked list with precision re-normalized to ``pi0`` (balanced by default)."""
    s, P, N = _sorted(pairs)
    if P == 0 or N == 0:
        return float("nan")
    tp = fp = 0
    best = 0.0
    for lab, _ in s:
        tp += lab
        fp += 1 - lab
        tpr, fpr = tp / P, fp / N
        denom = tpr * pi0 + fpr * (1 - pi0)
        p0 = (tpr * pi0) / denom if denom else 1.0
        if p0 + tpr:
            best = max(best, 2 * p0 * tpr / (p0 + tpr))
    return best


def demo():
    # a ranker that puts most positives first scores well on all; balanced PR >= micro PR when rare
    import random
    rng = random.Random(0)
    pos = [(1, rng.gauss(2, 1)) for _ in range(20)]
    neg = [(0, rng.gauss(0, 1)) for _ in range(2000)]            # 1% prevalence
    pairs = pos + neg
    assert 0.5 < roc_auc(pairs) <= 1.0
    assert pr_auc_balanced(pairs) > pr_auc(pairs)                 # balancing lifts the rare-class PR
    assert 0.0 <= f1_balanced(pairs) <= 1.0
    print(f"micro PR={pr_auc(pairs):.3f} balanced PR={pr_auc_balanced(pairs):.3f} "
          f"ROC={roc_auc(pairs):.3f} balanced F1={f1_balanced(pairs):.3f}")


if __name__ == "__main__":
    demo()
