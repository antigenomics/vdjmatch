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


def _curve(pairs):
    """Cumulative (tp, fp) points after each tied-score GROUP (descending score), plus P, N. Grouping
    equal scores gives them 0.5 credit (Mann-Whitney-consistent) — without it, tied positives that
    happen to precede tied negatives in the input inflate every AUC."""
    s = sorted(pairs, key=lambda x: -x[1])
    P = sum(l for l, _ in s)
    N = len(s) - P
    pts = [(0, 0)]
    tp = fp = i = 0
    while i < len(s):
        sc = s[i][1]
        j = i
        while j < len(s) and s[j][1] == sc:
            tp += s[j][0]
            fp += 1 - s[j][0]
            j += 1
        pts.append((tp, fp))
        i = j
    return pts, P, N


def pr_auc(pairs) -> float:
    """Micro PR-AUC at the observed prevalence (average precision)."""
    pts, P, _ = _curve(pairs)
    if P == 0:
        return float("nan")
    pr, pp, area = 0.0, 1.0, 0.0
    for tp, fp in pts[1:]:
        r, p = tp / P, tp / (tp + fp) if tp + fp else 1.0
        area += (r - pr) * (p + pp) / 2
        pr, pp = r, p
    return area


def roc_auc(pairs) -> float:
    """ROC-AUC (prevalence-invariant rank discriminator), ties grouped -> 0.5 credit."""
    pts, P, N = _curve(pairs)
    if P == 0 or N == 0:
        return float("nan")
    area = 0.0
    for (tp0, fp0), (tp1, fp1) in zip(pts, pts[1:]):
        area += (fp1 - fp0) / N * (tp1 + tp0) / (2 * P)
    return area


def auc01(pairs, max_fpr: float = 0.1) -> float:
    """IMMREP-style early-retrieval metric: McClish-standardised partial ROC-AUC over FPR<=max_fpr,
    normalised to [0.5, 1] (sklearn ``roc_auc_score(max_fpr=...)``). Random/all-equal -> 0.5."""
    from sklearn.metrics import roc_auc_score
    y = [a for a, _ in pairs]
    if len(set(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, [b for _, b in pairs], max_fpr=max_fpr))


def pr_auc_balanced(pairs, pi0: float = 0.5) -> float:
    """PR-AUC with precision re-normalized to reference prevalence ``pi0`` (balanced by default)."""
    pts, P, N = _curve(pairs)
    if P == 0 or N == 0:
        return float("nan")
    pr, pp, area = 0.0, 1.0, 0.0
    for tp, fp in pts[1:]:
        tpr, fpr = tp / P, fp / N
        denom = tpr * pi0 + fpr * (1 - pi0)
        p0 = (tpr * pi0) / denom if denom else 1.0
        area += (tpr - pr) * (p0 + pp) / 2
        pr, pp = tpr, p0
    return area


def f1_balanced(pairs, pi0: float = 0.5) -> float:
    """Best F1 along the ranked list with precision re-normalized to ``pi0`` (balanced by default)."""
    pts, P, N = _curve(pairs)
    if P == 0 or N == 0:
        return float("nan")
    best = 0.0
    for tp, fp in pts[1:]:
        tpr, fpr = tp / P, fp / N
        denom = tpr * pi0 + fpr * (1 - pi0)
        p0 = (tpr * pi0) / denom if denom else 1.0
        if p0 + tpr:
            best = max(best, 2 * p0 * tpr / (p0 + tpr))
    return best


def bootstrap_ci(pairs, metric, B: int = 2000, seed: int = 0, alpha: float = 0.05):
    """Stratified bootstrap CI for a ``(label,score) -> float`` metric. Resamples positives among
    positives and negatives among negatives WITH replacement (P and N held fixed, so prevalence is
    preserved), recomputes the metric B times, and returns ``(estimate, lo, hi)`` at the alpha/2 and
    1-alpha/2 percentiles. Used to colour tables (winner / tie / worst) by CI overlap. NaN-safe."""
    import random
    est = metric(pairs)
    pos = [p for p in pairs if p[0]]
    neg = [p for p in pairs if not p[0]]
    if not pos or not neg:
        return est, float("nan"), float("nan")
    rng = random.Random(seed)
    np_, nn = len(pos), len(neg)
    vals = []
    for _ in range(B):
        rs = [pos[rng.randrange(np_)] for _ in range(np_)] + [neg[rng.randrange(nn)] for _ in range(nn)]
        v = metric(rs)
        if v == v:                                               # drop NaN draws
            vals.append(v)
    if not vals:
        return est, float("nan"), float("nan")
    vals.sort()
    lo = vals[max(0, int((alpha / 2) * len(vals)))]
    hi = vals[min(len(vals) - 1, int((1 - alpha / 2) * len(vals)))]
    return est, lo, hi


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
    est, lo, hi = bootstrap_ci(pairs, roc_auc, B=500)
    assert lo <= est <= hi and 0.0 <= lo <= hi <= 1.0            # CI brackets the estimate
    print(f"micro PR={pr_auc(pairs):.3f} balanced PR={pr_auc_balanced(pairs):.3f} "
          f"ROC={roc_auc(pairs):.3f} [{lo:.3f},{hi:.3f}] balanced F1={f1_balanced(pairs):.3f}")


if __name__ == "__main__":
    demo()
