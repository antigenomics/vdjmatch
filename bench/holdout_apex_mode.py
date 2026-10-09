"""Test the apex-k-mer 'mode': does adding the leakage-safe V-conditioned apex-loop k-mer channel to the
shipped NED+V+J+len scorer improve ROC / AUC0.1 (esp. NLV, whose signal is the CDR3 apex motif, not V)
WITHOUT regression on the other epitopes? Min-p (Sidak) fusion is non-diluting, so a channel can only help.

    .venv/bin/python bench/holdout_apex_mode.py [TRB|TRA]
"""
from __future__ import annotations

import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
import holdout_eval as HE                                            # noqa: E402
from benchmark import A02, release, vgene                            # noqa: E402
from holdout_controls import _airr                                   # noqa: E402
from holdout_features import _emp_p, _jmap, build_models, ned_sim    # noqa: E402
from metrics import auc01, roc_auc                                   # noqa: E402


def apex_kmers(s, k=3, edge=3):
    return [s[i:i + k] for i in range(edge, len(s) - k + 1 - edge)] if len(s) > 2 * edge + k else []


class ApexVK:
    """V-conditioned apex-loop k-mer log-odds (corrected tuple order; leakage-safe, survives cross-study)."""
    def __init__(self, ref_vc, bg_vc, k=3):
        self.k, self.Vsp = k, 20 ** k
        self.rc, self.bc = Counter(), Counter()
        for s, v in ref_vc:
            for km in apex_kmers(s, k):
                self.rc[(v, km)] += 1
        for s, v in bg_vc:
            for km in apex_kmers(s, k):
                self.bc[(v, km)] += 1
        self.nr, self.nb = sum(self.rc.values()) + 1.0, sum(self.bc.values()) + 1.0

    def score(self, v, s, loo=False):
        d = 1 if loo else 0
        vals = []
        for km in apex_kmers(s, self.k):
            pr = max(self.rc.get((v, km), 0) - d, 0) + 1.0
            vals.append(math.log((pr / (self.nr - d + self.Vsp)) /
                                 ((self.bc.get((v, km), 0) + 1.0) / (self.nb + self.Vsp))))
        return float(np.mean(vals)) if vals else 0.0


def comb_minp(colmap, null):
    nulls = {c: sorted(colmap[c][null]) for c in colmap}
    K = len(colmap)
    return [-math.log10(max(1.0 - (1.0 - min(_emp_p(colmap[c][i], nulls[c]) for c in colmap)) ** K, 1e-12))
            for i in range(len(colmap[next(iter(colmap))]))]


def main(locus, ref="full"):
    cache, models = build_models(locus, "airr", ref)
    recs = cache["recs"]; jm = _jmap(locus)
    ned = [ned_sim(r) for r in recs]
    rdf = _bench.valid_cdr3(release("vdjdb2026").filter(pl.col("mhc_a").str.contains(A02) & (pl.col("gene") == locus)))
    bgd = _airr("human", locus, 60000); bg_vc = list(zip(bgd["cdr3"].to_list(), bgd["v"].to_list()))
    apex = {}
    for sh, e in HE.EPI.items():
        if e not in models:
            continue
        sub = rdf.filter(pl.col("epitope") == e).group_by("cdr3").agg(pl.col("v").first())
        if sub.height < 20:
            continue
        rvc = list(zip(sub["cdr3"].to_list(), [vgene(x) for x in sub["v"]]))
        avk = ApexVK(rvc, bg_vc)
        # unsupervised reliability: AUC(apex scores of E ref-binders vs airr_control) — no test labels
        rs = [avk.score(v, s, loo=True) for s, v in rvc[:800]]
        bs = [avk.score(v, s) for s, v in bg_vc[:800]]
        reliab = roc_auc([(1, x) for x in rs] + [(0, x) for x in bs])
        apex[e] = (avk, set(sub["cdr3"].to_list()), reliab)
    null = [i for i, r in enumerate(recs) if r["true"] is None]
    present = [sh for sh in HE.EPI if HE.EPI[sh] in models and any(r["true"] == HE.EPI[sh] for r in recs)]
    TAU = 0.60                                                       # include apex only where reliable
    print(f"\n=== {locus}/{ref}: base (NED+V+J+len) vs +apex-kmer (gated at reliab>={TAU}) ===")
    print(f"{'epi':5}{'reliab':>7}{'use':>5}{'ROC_base':>9}{'ROC_apex':>9}{'dROC':>7}{'A01_base':>9}{'A01_apex':>9}{'dA01':>7}")
    db, da, db1, da1 = [], [], [], []
    for sh in present:
        e = HE.EPI[sh]; m = models[e]
        cols = {"NED": np.array([ned[i].get(e, 0.0) for i in range(len(recs))]),
                "V": np.array([m["V"].get(r["v"], 0.0) for r in recs]),
                "J": np.array([m["J"].get(jm.get(r["cdr3"], ""), 0.0) for r in recs]),
                "len": np.array([m["len"].get(len(r["cdr3"]), 0.0) for r in recs])}
        avk, aset, reliab = apex[e]
        use = reliab >= TAU
        acol = np.array([avk.score(r["v"], r["cdr3"], loo=r["cdr3"] in aset) for r in recs])
        y = [int(r["true"] == e) for r in recs]
        base = comb_minp(cols, null)
        apx = comb_minp({**cols, "apex": acol}, null) if use else base
        rb, ra = roc_auc(list(zip(y, base))), roc_auc(list(zip(y, apx)))
        ab, aa = auc01(list(zip(y, base))), auc01(list(zip(y, apx)))
        db.append(rb); da.append(ra); db1.append(ab); da1.append(aa)
        print(f"{sh:5}{reliab:>7.2f}{'Y' if use else '.':>5}{rb:>9.3f}{ra:>9.3f}{ra - rb:>+7.3f}"
              f"{ab:>9.3f}{aa:>9.3f}{aa - ab:>+7.3f}")
    print(f"{'MEAN':5}{np.mean(db):>9.3f}{np.mean(da):>9.3f}{np.mean(da) - np.mean(db):>+7.3f}"
          f"{np.nanmean(db1):>9.3f}{np.nanmean(da1):>9.3f}{np.nanmean(da1) - np.nanmean(db1):>+7.3f}")
    print(f"  min dROC {min(a - b for a, b in zip(da, db)):+.3f}  (>=0 => no regression)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TRB")
