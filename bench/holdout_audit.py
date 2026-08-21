"""Adversarial leakage audit of the joint V+k-mer channel's near-perfect ROC.

Two decisive tests:
  split   : per-epitope Vkmer ROC computed separately against (a) the airr_control negatives only
            [E-binder vs random repertoire -- the easy contrast] and (b) other-epitope binders only
            [E vs other selected TCRs -- the hard, argmax-relevant contrast].
  xstudy  : rebuild the Vkmer model EXCLUDING the held-out source study's sequences; re-score the
            held-out positives. If ROC holds -> genuine epitope motif; if it collapses -> study-batch
            memorisation. Only meaningful for epitopes with independent studies (NLV/YLQ/GLC).

    .venv/bin/python bench/holdout_audit.py TRB
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
import holdout_eval as HE                                            # noqa: E402
import holdout_features as HF                                        # noqa: E402
from benchmark import A02, release, vgene                            # noqa: E402
from holdout_controls import _airr                                   # noqa: E402
from metrics import roc_auc                                          # noqa: E402

# held-out source studies (reference_id) per epitope, to exclude for the cross-study test
SRC = {"NLV": "https://github.com/antigenomics/vdjdb-db/issues/252",
       "LLW": "https://github.com/antigenomics/vdjdb-db/issues/193",
       "LLL": "https://github.com/antigenomics/vdjdb-db/issues/193",
       "ELA": "https://github.com/antigenomics/vdjdb-db/issues/193"}


def split(locus):
    cache, models = HF.build_models(locus, "airr")
    recs = cache["recs"]
    print(f"\n=== {locus}: Vkmer ROC by negative type ===")
    print(f"{'epitope':8}{'n_pos':>7}{'vs_airr':>9}{'n_air':>7}{'vs_other':>10}{'n_oth':>7}")
    for sh in HE.EPI:
        e = HE.EPI[sh]
        if e not in models or not any(r["true"] == e for r in recs):
            continue
        m = models[e]
        pos = [(r, m["Vkmer"].score(r["v"], r["cdr3"], loo=r["cdr3"] in m["ref_set"]))
               for r in recs if r["true"] == e]
        air = [(r, m["Vkmer"].score(r["v"], r["cdr3"], loo=r["cdr3"] in m["ref_set"]))
               for r in recs if r["true"] is None]
        oth = [(r, m["Vkmer"].score(r["v"], r["cdr3"], loo=r["cdr3"] in m["ref_set"]))
               for r in recs if r["true"] is not None and r["true"] != e]
        roc_air = roc_auc([(1, s) for _, s in pos] + [(0, s) for _, s in air])
        roc_oth = roc_auc([(1, s) for _, s in pos] + [(0, s) for _, s in oth])
        print(f"{sh:8}{len(pos):>7}{roc_air:>9.3f}{len(air):>7}{roc_oth:>10.3f}{len(oth):>7}")


def xstudy(locus):
    """Rebuild Vkmer excluding each epitope's held-out source study; re-score held-out positives vs all neg."""
    cache = HE.build(locus, "full")
    recs = cache["recs"]
    rdf = release("vdjdb2026").filter(pl.col("mhc_a").str.contains(A02) & (pl.col("gene") == locus))
    rdf = _bench.valid_cdr3(rdf)
    bgd = _airr("human", locus, 60000)
    bg_vc = list(zip(bgd["cdr3"].to_list(), bgd["v"].to_list()))
    print(f"\n=== {locus}: Vkmer ROC with held-out study EXCLUDED from the model ===")
    print(f"{'epitope':8}{'study_refs':>11}{'kept_refs':>10}{'roc_full':>9}{'roc_xstudy':>11}")
    for sh, src in SRC.items():
        e = HE.EPI[sh]
        sub = rdf.filter(pl.col("epitope") == e)
        if sub.height < 20 or not any(r["true"] == e for r in recs):
            continue
        in_study = sub.filter(pl.col("reference_id").str.contains(src.split("/")[-1], literal=True))
        kept = sub.filter(~pl.col("reference_id").str.contains(src.split("/")[-1], literal=True))
        # dedup by cdr3 for each model
        kv = kept.group_by("cdr3").agg(pl.col("v").first())
        fv = sub.group_by("cdr3").agg(pl.col("v").first())
        ref_set = set(fv["cdr3"].to_list())
        m_full = HF.VKmer(list(zip(fv["cdr3"], [vgene(x) for x in fv["v"]])), bg_vc)
        m_xs = HF.VKmer(list(zip(kv["cdr3"], [vgene(x) for x in kv["v"]])), bg_vc) if kv.height >= 20 else None
        y = [int(r["true"] == e) for r in recs]
        sf = [m_full.score(r["v"], r["cdr3"], loo=r["cdr3"] in ref_set) for r in recs]
        rf = roc_auc(list(zip(y, sf)))
        if m_xs is None:
            rx = float("nan")
        else:
            xset = set(kv["cdr3"].to_list())
            sx = [m_xs.score(r["v"], r["cdr3"], loo=r["cdr3"] in xset) for r in recs]
            rx = roc_auc(list(zip(y, sx)))
        print(f"{sh:8}{in_study.select('cdr3').unique().height:>11}{kv.height:>10}{rf:>9.3f}{rx:>11.3f}")


if __name__ == "__main__":
    loc = sys.argv[1] if len(sys.argv) > 1 else "TRB"
    split(loc)
    xstudy(loc)
