"""What exactly drives the joint V+k-mer signal: PCR-error trails, V-protocol bias, or genuine motif?

Mechanistic decomposition (no raw clonotypes printed; k-mers + V genes are aggregate features):
  trail   : split held-out positives by whether they have an edit<=1 same-epitope reference neighbour
            (a near-duplicate / PCR-error trail). Report Vkmer ROC for the WITH-trail vs NO-trail subset.
            If NO-trail positives still score high -> the signal is a genuine shared motif, not a trail.
  Vvskmer : cross-study ROC for V-only, k-mer-only (V-agnostic), and joint -- isolates V-protocol bias
            (V usage can be protocol-driven) from CDR3 motif (no known batch mechanism but PCR trails).
  motifs  : top (V, 3-mer) log-odds contributors aggregated over each epitope's held-out positives.

    .venv/bin/python bench/holdout_kmer_probe.py TRB
"""
from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
import holdout_audit as HA                                           # noqa: E402
import holdout_eval as HE                                            # noqa: E402
import holdout_features as HF                                        # noqa: E402
from benchmark import A02, release, vgene                            # noqa: E402
from holdout_controls import _airr                                   # noqa: E402
from metrics import roc_auc                                          # noqa: E402
from vdjmatch.match import vgene as _vg                              # noqa: E402


def ntrail(rec, e):
    """# same-epitope reference neighbours at edit<=1 (near-duplicate / PCR-error-trail count)."""
    return sum(1 for ps, ed, he, hv in rec["hits"] if he == e and ed <= 1)


def main(locus):
    cache, models = HF.build_models(locus, "airr")
    recs = cache["recs"]

    # --- trail split: does the signal survive when positives have NO edit<=1 same-E neighbour? ---
    print(f"\n=== {locus}: Vkmer ROC split by near-duplicate (edit<=1 same-E) presence ===")
    print(f"{'epitope':8}{'n_pos':>6}{'%w/trail':>9}{'ROC_trail':>10}{'n':>5}{'ROC_notrail':>12}{'n':>5}")
    for sh in HE.EPI:
        e = HE.EPI[sh]
        if e not in models or not any(r["true"] == e for r in recs):
            continue
        m = models[e]
        def vk(r):
            return m["Vkmer"].score(r["v"], r["cdr3"], loo=r["cdr3"] in m["ref_set"])
        neg = [(0, vk(r)) for r in recs if r["true"] != e]
        pos = [(r, ntrail(r, e), vk(r)) for r in recs if r["true"] == e]
        wt = [(1, s) for r, nt, s in pos if nt >= 1]
        nt0 = [(1, s) for r, nt, s in pos if nt == 0]
        frac = 100 * len(wt) / max(len(pos), 1)
        roc_wt = roc_auc(wt + neg) if wt else float("nan")
        roc_nt0 = roc_auc(nt0 + neg) if nt0 else float("nan")
        print(f"{sh:8}{len(pos):>6}{frac:>8.0f}%{roc_wt:>10.3f}{len(wt):>5}{roc_nt0:>12.3f}{len(nt0):>5}")

    # --- V-only vs kmer-only vs joint, cross-study (exclude held-out source study) ---
    print(f"\n=== {locus}: cross-study ROC by channel (source study removed from model) ===")
    print(f"{'epitope':8}{'V_full':>8}{'V_xs':>8}{'kmer_full':>10}{'kmer_xs':>9}{'joint_full':>11}{'joint_xs':>9}")
    rdf = _bench.valid_cdr3(release("vdjdb2026").filter(pl.col("mhc_a").str.contains(A02)
                                                        & (pl.col("gene") == locus)))
    bgd = _airr("human", locus, 60000)
    bg_v = bgd["v"].to_list()
    bg_vc = list(zip(bgd["cdr3"].to_list(), bg_v))
    from hardcase import kmer_logodds, vj_logodds                    # noqa: E402
    for sh, src in HA.SRC.items():
        e = HE.EPI[sh]
        sub = rdf.filter(pl.col("epitope") == e)
        if sub.height < 20 or not any(r["true"] == e for r in recs):
            continue
        sid = src.split("/")[-1]
        kept = sub.filter(~pl.col("reference_id").str.contains(sid, literal=True))
        fv = sub.group_by("cdr3").agg(pl.col("v").first())
        kv = kept.group_by("cdr3").agg(pl.col("v").first())
        ok = kv.height >= 20
        fset, kset = set(fv["cdr3"].to_list()), set(kv["cdr3"].to_list())
        fvg, kvg = [vgene(x) for x in fv["v"]], [vgene(x) for x in kv["v"]]
        # models
        Vf = vj_logodds(fvg, bg_v); Vx = vj_logodds(kvg, bg_v) if ok else None
        Kf = HF.Kmer(fv["cdr3"].to_list()); Kx = HF.Kmer(kv["cdr3"].to_list()) if ok else None
        Jf = HF.VKmer(list(zip(fv["cdr3"], fvg)), bg_vc)
        Jx = HF.VKmer(list(zip(kv["cdr3"], kvg)), bg_vc) if ok else None
        y = [int(r["true"] == e) for r in recs]
        def roc_for(getter):
            return roc_auc(list(zip(y, [getter(r) for r in recs])))
        vF = roc_for(lambda r: Vf.get(r["v"], 0.0))
        vX = roc_for(lambda r: Vx.get(r["v"], 0.0)) if ok else float("nan")
        kF = roc_for(lambda r: Kf.score(r["cdr3"], loo=r["cdr3"] in fset))
        kX = roc_for(lambda r: Kx.score(r["cdr3"], loo=r["cdr3"] in kset)) if ok else float("nan")
        jF = roc_for(lambda r: Jf.score(r["v"], r["cdr3"], loo=r["cdr3"] in fset))
        jX = roc_for(lambda r: Jx.score(r["v"], r["cdr3"], loo=r["cdr3"] in kset)) if ok else float("nan")
        print(f"{sh:8}{vF:>8.3f}{vX:>8.3f}{kF:>10.3f}{kX:>9.3f}{jF:>11.3f}{jX:>9.3f}")

    # --- APEX-restricted: exclude the germline C...F frame (first 3 / last 3 residues) ---
    # If the joint blow-up is germline-stem fingerprinting, apex-only joint should NOT inflate and
    # should reflect genuine loop-tip motif (high for convergent, chance for diverse), honestly.
    def apex_kmers(s, k=3, edge=3):
        return [(i, s[i:i + k]) for i in range(edge, len(s) - k + 1 - edge)] if len(s) > 2 * edge + k else []

    class ApexVKmer:
        def __init__(self, ref_vc, bg_vc, k=3):
            self.k = k; self.Vsp = 20 ** k
            self.rc, self.bc = Counter(), Counter()
            for v, s in ref_vc:
                for _, km in apex_kmers(s, k):
                    self.rc[(v, km)] += 1
            for v, s in bg_vc:
                for _, km in apex_kmers(s, k):
                    self.bc[(v, km)] += 1
            self.nr, self.nb = sum(self.rc.values()) + 1.0, sum(self.bc.values()) + 1.0

        def score(self, v, s, loo=False):
            d = 1 if loo else 0
            vals = []
            for _, km in apex_kmers(s, self.k):
                pr = max(self.rc.get((v, km), 0) - d, 0) + 1.0
                vals.append(np.log((pr / (self.nr - d + self.Vsp)) /
                                   ((self.bc.get((v, km), 0) + 1.0) / (self.nb + self.Vsp))))
            return float(np.mean(vals)) if vals else 0.0

    print(f"\n=== {locus}: APEX-only joint V+kmer ROC (germline frame excluded), full vs cross-study ===")
    print(f"{'epitope':8}{'apexJ_full':>11}{'apexJ_xs':>10}  (compare joint_full 0.99 / joint_xs above)")
    for sh, src in HA.SRC.items():
        e = HE.EPI[sh]
        sub = rdf.filter(pl.col("epitope") == e)
        if sub.height < 20 or not any(r["true"] == e for r in recs):
            continue
        sid = src.split("/")[-1]
        kept = sub.filter(~pl.col("reference_id").str.contains(sid, literal=True))
        fv = sub.group_by("cdr3").agg(pl.col("v").first())
        kv = kept.group_by("cdr3").agg(pl.col("v").first())
        ok = kv.height >= 20
        fset, kset = set(fv["cdr3"].to_list()), set(kv["cdr3"].to_list())
        Af = ApexVKmer(list(zip(fv["cdr3"], [vgene(x) for x in fv["v"]])), bg_vc)
        Ax = ApexVKmer(list(zip(kv["cdr3"], [vgene(x) for x in kv["v"]])), bg_vc) if ok else None
        y = [int(r["true"] == e) for r in recs]
        aF = roc_auc(list(zip(y, [Af.score(r["v"], r["cdr3"], loo=r["cdr3"] in fset) for r in recs])))
        aX = (roc_auc(list(zip(y, [Ax.score(r["v"], r["cdr3"], loo=r["cdr3"] in kset) for r in recs])))
              if ok else float("nan"))
        print(f"{sh:8}{aF:>11.3f}{aX:>10.3f}")

    # --- top (V,3-mer) motifs by aggregate log-odds over held-out positives ---
    print(f"\n=== {locus}: top (V, 3-mer) contributors among held-out positives ===")
    for sh in HE.EPI:
        e = HE.EPI[sh]
        if e not in models or not any(r["true"] == e for r in recs):
            continue
        m = models[e]
        agg = Counter()
        for r in (r for r in recs if r["true"] == e):
            for i in range(len(r["cdr3"]) - 2):
                km = r["cdr3"][i:i + 3]
                pr = max(m["Vkmer"].rc.get((r["v"], km), 0) - (1 if r["cdr3"] in m["ref_set"] else 0), 0) + 1.0
                lo = np.log((pr / (m["Vkmer"].nr + m["Vkmer"].Vsp)) /
                            ((m["Vkmer"].bc.get((r["v"], km), 0) + 1.0) / (m["Vkmer"].nb + m["Vkmer"].Vsp)))
                agg[(r["v"], km)] += lo
        top = ", ".join(f"{v}+{km}({s:.0f})" for (v, km), s in agg.most_common(6))
        print(f"  {sh}: {top}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TRB")
