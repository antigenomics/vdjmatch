"""Reference-vs-background FEATURE channels for the held-out benchmark (diverse-epitope recovery).

The diverse epitopes (ELA/LLW/LLL) have no neighbour-density signal: every query has equally-many far
(edit 4-5) reference neighbours, binder or not. Their signal lives in repertoire-level FEATURES learned
unsupervised from the reference DB vs a real post-selection background (airr_control):
  V    : gene-level V-gene usage log-odds (allele stripped; TRBV7-9 != TRBV7-2)
  kmer : CDR3 3-mer log-odds
  Vkmer: JOINT (V-gene, 3-mer) log-odds  -- TRBV7-9+RSG, not P(V)*P(kmer)
  comp : Gaussian log-odds over [apex hydropathy, whole hydropathy, aromatic frac, net charge,
         MJ strong-residue frac, length]  (hydrophobicity / strong-weak modes)
Each is log P(feature|E binders)/P(feature|background); the query's own CDR3 is held out of the model.
Per epitope E the reference size is constant across queries, so these size-aware LLRs are valid ranking
scores; we report per-epitope ROC for each channel + the neighbour-density NED + a rank fusion.

    .venv/bin/python bench/holdout_features.py [TRB|TRA] [airr|other]
"""
from __future__ import annotations

import bisect
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
import holdout_eval as HE                                            # noqa: E402
from benchmark import A02, PSSM_SCALE, SOFTV_BETA, release, shortlist, vgene   # noqa: E402
from hardcase import AROM, CHG, KD, _apex, kmer_logodds, vj_logodds  # noqa: E402
from holdout_controls import _airr                                   # noqa: E402
from metrics import auc01, roc_auc                                   # noqa: E402
from vdjmatch.match import vgene as _vg                              # noqa: E402

# Miyazawa-Jernigan "strong" interacting residues (large hydrophobic/aromatic core); the strong/weak mode.
MJ_STRONG = set("FMILYWVC")


def _comp(s):
    ap = _apex(s)
    return np.array([np.mean([KD[c] for c in ap if c in KD] or [0]),       # apex hydropathy
                     np.mean([KD[c] for c in s if c in KD] or [0]),        # whole hydropathy
                     sum(c in AROM for c in s) / max(len(s), 1),           # aromatic frac
                     sum(CHG.get(c, 0) for c in s),                        # net charge
                     sum(c in MJ_STRONG for c in s) / max(len(s), 1),      # MJ strong-residue frac
                     float(len(s))])


def _gauss_logodds(R, B):
    mr, sr = R.mean(0), R.std(0) + 1e-6
    mb, sb = B.mean(0), B.std(0) + 1e-6

    def score(f):
        return float(np.sum((-0.5 * ((f - mr) / sr) ** 2 - np.log(sr)) -
                            (-0.5 * ((f - mb) / sb) ** 2 - np.log(sb))))
    return score


def _kmers(s, k=3):
    return [s[i:i + k] for i in range(len(s) - k + 1)]


class VKmer:
    """JOINT (V-gene, k-mer) log-odds with Laplace smoothing; supports leave-one-out by count subtraction."""

    def __init__(self, ref_vc, bg_vc, k=3):
        self.k, self.Vsp = k, 20 ** k
        self.rc, self.bc = Counter(), Counter()
        for s, v in ref_vc:                                        # ref_vc is (cdr3, vgene)
            for km in _kmers(s, k):
                self.rc[(v, km)] += 1
        for s, v in bg_vc:
            for km in _kmers(s, k):
                self.bc[(v, km)] += 1
        self.nr, self.nb = sum(self.rc.values()) + 1.0, sum(self.bc.values()) + 1.0

    def score(self, v, s, loo=False):
        d = 1 if loo else 0                                         # drop q's own contribution if in ref
        vals = []
        for km in _kmers(s, self.k):
            pr = max(self.rc.get((v, km), 0) - d, 0) + 1.0
            vals.append(math.log((pr / (self.nr - d + self.Vsp)) /
                                 ((self.bc.get((v, km), 0) + 1.0) / (self.nb + self.Vsp))))
        return float(np.mean(vals)) if vals else 0.0


class Kmer:
    """CDR3 k-mer log-odds vs an AA-marginal background; LOO by count subtraction."""

    def __init__(self, refs, k=3):
        self.k, self.V = k, 20 ** k
        self.cnt, aa = Counter(), Counter()
        for s in refs:
            for km in _kmers(s, k):
                self.cnt[km] += 1
            aa.update(s)
        self.tot = sum(self.cnt.values()) + 1.0
        naa = sum(aa.values()) or 1
        self.pa = {a: aa[a] / naa for a in aa}

    def _bg(self, km):
        p = 1.0
        for c in km:
            p *= self.pa.get(c, 1e-3)
        return max(p, 1e-12)

    def score(self, s, loo=False):
        d = 1 if loo else 0
        vals = []
        for km in _kmers(s, self.k):
            pr = (max(self.cnt.get(km, 0) - d, 0) + 1.0) / (self.tot - d + self.V)
            vals.append(math.log(pr / self._bg(km)))
        return float(np.mean(vals)) if vals else 0.0


def _len_logodds(ref_len, bg_len):
    from collections import Counter as _C
    rr, bb = _C(ref_len), _C(bg_len)
    nr, nb = sum(rr.values()) + len(rr), sum(bb.values()) + len(bb)
    return {L: math.log(((rr[L] + 1) / nr) / ((bb[L] + 1) / nb)) for L in set(rr) | set(bb)}


def build_models(locus, bg="airr", ref="full"):
    """Per-epitope reference feature models (V/J/length/kmer/comp) vs a background; LOO per query."""
    cache = HE.build(locus, ref)
    rdf = release("vdjdb2026").filter(pl.col("mhc_a").str.contains(A02))
    if ref == "shortlist":
        rdf = shortlist(rdf, min_refs=2)
    r = (_bench.valid_cdr3(rdf.filter(pl.col("gene") == locus)).group_by("cdr3")
         .agg(pl.col("epitope").first(), pl.col("v").first(), pl.col("j").first()))
    by = defaultdict(list)
    for c, e, v, j in zip(r["cdr3"], r["epitope"], r["v"], r["j"]):
        by[e].append((c, vgene(v), vgene(j)))
    bgd = _airr("human", locus, 60000)
    bg_vc = list(zip(bgd["cdr3"].to_list(), bgd["v"].to_list()))
    bg_j, bg_len = bgd["j"].to_list(), [len(c) for c in bgd["cdr3"]]
    models = {}
    for sh, e in HE.EPI.items():
        ref_vcj = by.get(e, [])
        if len(ref_vcj) < 20:
            continue
        ref_cdr3 = [c for c, v, j in ref_vcj]
        ref_vc = [(c, v) for c, v, j in ref_vcj]
        ref_set = set(ref_cdr3)
        contrast = bg_vc if bg == "airr" else [(c, v) for ee, lst in by.items() if ee != e for c, v, j in lst]
        bcd = [c for c, v in contrast]
        models[e] = {
            "ref_set": ref_set,
            "V": vj_logodds([v for c, v, j in ref_vcj], [v for c, v in contrast]),
            "J": vj_logodds([j for c, v, j in ref_vcj], bg_j),
            "len": _len_logodds([len(c) for c in ref_cdr3], bg_len),
            "kmer": Kmer(ref_cdr3, k=3),
            "Vkmer": VKmer(ref_vc, contrast),
            "comp": _gauss_logodds(np.array([_comp(s) for s in ref_cdr3]),
                                   np.array([_comp(s) for s in bcd])),
        }
    return cache, models


def ned_sim(rec):
    """size-invariant neighbour-density per epitope (the baseline argmax score)."""
    sim = defaultdict(float)
    cs = rec["ctrl"]
    for ps, ed, he, hv in rec["hits"]:
        w = 1.0 if _vg.gene_family(rec["v"]) == _vg.gene_family(hv) else SOFTV_BETA * _vg.vsim(rec["v"], hv)
        if w > 0:
            sim[he] += w * math.exp(-ps / PSSM_SCALE) / (bisect.bisect_right(cs, ed) + 1)
    return sim


def main(locus, bg="airr"):
    cache, models = build_models(locus, bg)
    recs = cache["recs"]
    chans = ["NED", "V", "kmer", "Vkmer", "comp"]
    present = [sh for sh in HE.EPI if HE.EPI[sh] in models and any(r["true"] == HE.EPI[sh] for r in recs)]
    # precompute per-query scores per channel per epitope
    ned = [ned_sim(r) for r in recs]
    print(f"\n=== {locus}: per-epitope ROC by FEATURE channel (bg={bg}, n={len(recs)}) ===")
    print(f"{'epitope':8}" + "".join(f"{c:>8}" for c in chans) + f"{'FUSE':>8}{'ORACLE':>8}")
    means = defaultdict(list)
    for sh in present:
        e = HE.EPI[sh]
        m = models[e]
        y = [int(r["true"] == e) for r in recs]
        sc = {
            "NED": [ned[i].get(e, 0.0) for i in range(len(recs))],
            "V": [m["V"].get(r["v"], 0.0) for r in recs],
            "kmer": [m["kmer"].score(r["cdr3"], loo=r["cdr3"] in m["ref_set"]) for r in recs],
            "Vkmer": [m["Vkmer"].score(r["v"], r["cdr3"], loo=r["cdr3"] in m["ref_set"]) for r in recs],
            "comp": [m["comp"](_comp(r["cdr3"])) for r in recs],
        }
        line = f"{sh:8}"
        rocs = {}
        for c in chans:
            rocs[c] = roc_auc(list(zip(y, sc[c])))
            means[c].append(rocs[c])
            line += f"{rocs[c]:>8.3f}"
        rk = {c: np.argsort(np.argsort(sc[c])) for c in chans}
        fuse = sum(rk[c] for c in chans)
        rf = roc_auc(list(zip(y, fuse.tolist())))
        means["FUSE"].append(rf)
        orc = max(rocs.values())                                   # best single channel = achievable ceiling
        means["ORACLE"].append(orc)
        line += f"{rf:>8.3f}{orc:>8.3f}"
        print(line)
    print(f"{'MEAN':8}" + "".join(f"{np.mean(means[c]):>8.3f}" for c in chans)
          + f"{np.mean(means['FUSE']):>8.3f}{np.mean(means['ORACLE']):>8.3f}")


def _emp_p(score, null_sorted):
    """empirical upper-tail p: fraction of the airr_control null at or above `score` (Laplace-smoothed)."""
    import bisect as _b
    ge = len(null_sorted) - _b.bisect_left(null_sorted, score)
    return (1 + ge) / (1 + len(null_sorted))


def _jmap(locus):
    """cdr3 -> J gene from the held-out TSVs (the query cache stores no J)."""
    m = {}
    for name in HE.DATASETS[locus]:
        f = HE.HD / f"{name}.tsv"
        if f.exists():
            d = pl.read_csv(f, separator="\t", infer_schema_length=0)
            for c, j in zip(d["cdr3"], d["j"]):
                m.setdefault(c, vgene(j) if j else "")
    return m


def robust_eval(locus, ref="full", alpha=1e-3):
    """Leakage-robust scorer: NED + V + J + length log-odds, combined by calibrated min-p (Sidak),
    non-diluting; first-hit Poisson gate. Per-epitope null = the airr_control negatives. ROC + confusion."""
    import numpy as np
    cache, models = build_models(locus, "airr", ref)
    recs = cache["recs"]
    jm = _jmap(locus)
    ned = [ned_sim(r) for r in recs]
    present = [sh for sh in HE.EPI if HE.EPI[sh] in models and any(r["true"] == HE.EPI[sh] for r in recs)]
    null = [i for i, r in enumerate(recs) if r["true"] is None]
    comb = {}
    for sh in present:
        e = HE.EPI[sh]
        m = models[e]
        cols = {
            "NED": np.array([ned[i].get(e, 0.0) for i in range(len(recs))]),
            "V": np.array([m["V"].get(r["v"], 0.0) for r in recs]),
            "J": np.array([m["J"].get(jm.get(r["cdr3"], ""), 0.0) for r in recs]),
            "len": np.array([m["len"].get(len(r["cdr3"]), 0.0) for r in recs]),
        }
        nulls = {c: sorted(cols[c][null]) for c in cols}
        K = len(cols)
        cp = []
        for i in range(len(recs)):
            pmin = min(_emp_p(cols[c][i], nulls[c]) for c in cols)
            cp.append(-math.log10(max(1.0 - (1.0 - pmin) ** K, 1e-12)))    # Sidak over K channels
        comb[e] = cp
    # gate: continuous first-hit control-calibrated Poisson E-value (V-agnostic, not floored), at alpha
    from vdjmatch.evalue import first_hit                            # noqa: E402
    Ntot = sum(cache["n_epi"].values())
    gated = []
    for r in recs:
        t1 = [(c, e) for c, e in r["first"] if c <= 1]
        gated.append(first_hit.pvalue(t1, r["ctrl"], Ntot, cache["M"])["p_enrichment"] < alpha)
    assigns = []
    for i in range(len(recs)):
        best = max(present, key=lambda sh: comb[HE.EPI[sh]][i])
        assigns.append(HE.EPI[best] if gated[i] else None)
    print(f"\n=== {locus} / {ref}: leakage-robust (NED+V+J+len rank, first-hit gate) ===")
    print(f"{'epi':5}{'n+':>5}{'TP':>5}{'FN':>5}{'FP':>4}{'prec':>7}{'rec':>7}{'ROC':>7}{'AUC0.1':>8}{'dBASE':>7}")
    BASE = {"NLV": 0.690, "LLW": 0.511, "LLL": 0.622, "ELA": 0.522, "YLQ": 0.954, "GLC": 0.860,
            "NLV_TRA": 0.674, "LLL_TRA": 0.710, "GLC_TRA": 0.876, "YLQ_TRA": 0.889}
    rocs, a01s = [], []
    for sh in present:
        e = HE.EPI[sh]
        y = [int(r["true"] == e) for r in recs]
        roc = roc_auc(list(zip(y, comb[e])))
        a01 = auc01(list(zip(y, comb[e])))
        a01s.append(a01)
        tp = sum(1 for i, r in enumerate(recs) if r["true"] == e and assigns[i] == e)
        fn = sum(1 for i, r in enumerate(recs) if r["true"] == e and assigns[i] != e)
        fp = sum(1 for i, r in enumerate(recs) if r["true"] != e and assigns[i] == e)
        npos = tp + fn
        prec = tp / (tp + fp) if tp + fp else float("nan")
        b = BASE.get(sh if locus == "TRB" else sh + "_TRA", float("nan"))
        rocs.append(roc)
        print(f"{sh:5}{npos:>5}{tp:>5}{fn:>5}{fp:>4}{prec:>7.3f}{tp/npos if npos else 0:>7.3f}"
              f"{roc:>7.3f}{a01:>8.3f}{roc-b:>+7.3f}")
    print(f"  macro ROC {np.mean(rocs):.3f}  macro AUC0.1 {np.nanmean(a01s):.3f}  "
          f"(baseline ROC {np.mean([BASE.get(sh if locus=='TRB' else sh+'_TRA',0) for sh in present]):.3f})")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "robust":
        robust_eval(sys.argv[2] if len(sys.argv) > 2 else "TRB",
                    sys.argv[3] if len(sys.argv) > 3 else "full")
    else:
        main(sys.argv[1] if len(sys.argv) > 1 else "TRB", sys.argv[2] if len(sys.argv) > 2 else "airr")
