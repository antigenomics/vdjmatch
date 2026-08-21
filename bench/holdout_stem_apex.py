"""Resolve the V+CAS paradox: is the joint blow-up "just V", or sparse apex k-mers?

V+CAS is essentially V (CAS is the universal germline stem present in ~every TRB CDR3), so log P(V,CAS|ref)/
P(V,CAS|bg) reduces to the V-marginal log-odds and CANNOT carry more than V. We split each query's joint
score by k-mer position into STEM (k-mers touching the first/last 3 residues = germline frame) and APEX
(k-mers fully inside the loop tip), and report each part's ROC + its correlation with the V-only score.
Expectation: stem-joint ~ V (robust ~0.6, corr high); apex-joint = the blow-up (sparse fingerprints, ~0.99
full, collapses cross-study, corr with V low).

    .venv/bin/python bench/holdout_stem_apex.py TRB
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import holdout_eval as HE                                            # noqa: E402
import holdout_features as HF                                        # noqa: E402
from metrics import roc_auc                                          # noqa: E402


def part_score(m, v, s, which, k=3, edge=3):
    """mean joint log-odds over only STEM or APEX k-mers of s."""
    n = len(s)
    vals = []
    for i in range(n - k + 1):
        is_apex = i >= edge and i + k <= n - edge
        if (which == "apex") != is_apex:
            continue
        km = s[i:i + k]
        pr = m.rc.get((v, km), 0) + 1.0
        pb = m.bc.get((v, km), 0) + 1.0
        vals.append(math.log((pr / (m.nr + m.Vsp)) / (pb / (m.nb + m.Vsp))))
    return float(np.mean(vals)) if vals else 0.0


def corr(a, b):
    a, b = np.array(a), np.array(b)
    a, b = a - a.mean(), b - b.mean()
    d = math.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d else float("nan")


def main(locus):
    cache, models = HF.build_models(locus, "airr")
    recs = cache["recs"]
    print(f"\n=== {locus}: joint blow-up decomposed by k-mer position ===")
    print(f"{'epi':5}{'V':>7}{'stem':>7}{'apex':>7}{'full':>7}   {'corr(stem,V)':>13}{'corr(apex,V)':>14}")
    for sh in HE.EPI:
        e = HE.EPI[sh]
        if e not in models or not any(r["true"] == e for r in recs):
            continue
        m = models[e]
        vk = m["Vkmer"]
        y = [int(r["true"] == e) for r in recs]
        sV = [m["V"].get(r["v"], 0.0) for r in recs]
        sStem = [part_score(vk, r["v"], r["cdr3"], "stem") for r in recs]
        sApex = [part_score(vk, r["v"], r["cdr3"], "apex") for r in recs]
        sFull = [vk.score(r["v"], r["cdr3"], loo=r["cdr3"] in m["ref_set"]) for r in recs]
        print(f"{sh:5}{roc_auc(list(zip(y, sV))):>7.3f}{roc_auc(list(zip(y, sStem))):>7.3f}"
              f"{roc_auc(list(zip(y, sApex))):>7.3f}{roc_auc(list(zip(y, sFull))):>7.3f}"
              f"   {corr(sStem, sV):>13.3f}{corr(sApex, sV):>14.3f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TRB")
