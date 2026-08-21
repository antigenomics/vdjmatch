"""Sweep VDJdb2026 A*02 TRB epitopes: which does NED detect poorly (1-vs-other ROC/PR < 0.6), and why.

Per epitope: 50/50 seed-42 split (reference half builds the index; positive-query half is scored against
it), negatives = seed-42 sample of OTHER epitopes' binders. NED 1-vs-other ROC/balPR. For the flagged
(<0.6) epitopes, the same feature diagnosis as hardcase.py (univariate AUCs + 5-fold CV logistic ceiling)
classifies the failure: unrecoverable (ceiling<0.60) / recoverable-composition (ceiling-NED>+0.08 and a
composition feature tops) / tie-artifact (ROC ok, balPR collapses) / other. READ-ONLY library; new files only.

    .venv/bin/python bench/hardcase_sweep.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench                                                        # noqa: E402
import benchmark as B                                               # noqa: E402
import hardcase as H                                                 # noqa: E402
import hardcase_prior2 as P2                                         # noqa: E402
import hardcase_control as HC                                        # noqa: E402
from metrics import roc_auc, pr_auc_balanced                         # noqa: E402
from vdjmatch import db                                              # noqa: E402
from vdjmatch.evalue import background                               # noqa: E402
from vdjmatch.evalue.first_hit import scope                          # noqa: E402

LOCUS = "TRB"
MIN_BIND = 60
MAX_POS, NEG_MULT = 200, 3
OUT = Path(__file__).resolve().parent / "_hardcase"
OUT.mkdir(exist_ok=True)


def binders(v26):
    a = _bench.valid_cdr3(v26.filter((pl.col("gene") == LOCUS) & pl.col("mhc_a").str.contains(B.A02)))
    return a.with_columns(v=pl.col("v").map_elements(B.vgene, return_dtype=pl.Utf8))


def split(df, frac=0.5, seed=42):
    d = df.unique("cdr3", maintain_order=True).sort("cdr3").sample(fraction=1.0, shuffle=True, seed=seed)
    k = int(d.height * frac)
    return d.head(k), d.tail(d.height - k)


def feature_diag(ref_df, qpool, y, NED):
    """univariate AUCs + CV ceiling for a flagged epitope (qpool: cdr3,v rows aligned with y/NED)."""
    refc = ref_df["cdr3"].to_list()
    kmer = H.kmer_logodds(refc, 3)
    vlo = H.vj_logodds(ref_df["v"].to_list(), qpool["v"].to_list())
    reflen = float(np.median([len(s) for s in refc]))
    rows = []
    for s, v in zip(qpool["cdr3"].to_list(), qpool["v"].to_list()):
        ap = H._apex(s)
        rows.append([len(s), -abs(len(s) - reflen),
                     np.mean([H.KD[c] for c in ap if c in H.KD] or [0]),
                     sum(c in H.AROM for c in s) / max(len(s), 1),
                     sum(H.CHG.get(c, 0) for c in s), kmer(s), vlo.get(v, 0.0)])
    names = ["len", "len_dev", "kd_apex", "arom", "charge", "kmer", "vlo", "NED"]
    X = np.column_stack([np.array(rows, float), NED])
    uni = {n: H.auc1(X[:, i], y)[0] for i, n in enumerate(names)}
    ceil = H.cv_ceiling(X, y)
    return uni, (ceil[0] if ceil else None)


def classify(roc, pr, ceil, uni):
    if ceil is not None and ceil < 0.60:
        return "unrecoverable"
    comp = max(("kd_apex", "arom", "kmer"), key=lambda k: uni.get(k, 0))
    top = max(uni, key=lambda k: uni[k])
    if ceil and roc and ceil - roc > 0.08 and top in ("kd_apex", "arom", "kmer"):
        return "recoverable-composition"
    if roc and roc >= 0.6 and pr < 0.6:
        return "tie-artifact"
    return "other"


def main():
    import argparse
    bgmode = argparse.ArgumentParser()
    bgmode.add_argument("--bg", default="control", choices=["olga", "other", "control"],
                        help="prior background null (default: airr_control post-selection repertoire)")
    BGM = bgmode.parse_args().bg
    v26 = db.load(_bench.source(), species="HomoSapiens")
    b = binders(v26)
    sizes = (b.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
             .filter(pl.col("n") >= MIN_BIND).sort("n", descending=True))
    eps = sizes["epitope"].to_list()
    ctrl = background(LOCUS)
    # global prior background (olga/control are the same for every epitope -> precompute; other is per-epitope)
    GBG = HC.control_cdr3(LOCUS, 60000) if BGM == "control" else (P2.olga_cdr3(LOCUS) if BGM == "olga" else None)
    GBG_comp = np.array([P2.comp(s) for s in GBG]) if GBG else None
    GBG_kmer = H.kmer_logodds(GBG, 3) if GBG else None
    print(f"{len(eps)} epitopes >={MIN_BIND} binders; bg={BGM}; scoring 1-vs-other ...", file=sys.stderr)
    rows = []
    for E in eps:
        pos = b.filter(pl.col("epitope") == E)
        ref_df, posq = split(pos)
        if posq.height > MAX_POS:
            posq = posq.sample(MAX_POS, seed=42)
        oth = b.filter(pl.col("epitope") != E).unique("cdr3", maintain_order=True).sort("cdr3")
        nneg = min(posq.height * NEG_MULT, oth.height)
        neg = oth.sample(nneg, seed=42)
        qpool = pl.concat([posq.select("cdr3", "v"), neg.select("cdr3", "v")])
        y = np.array([1] * posq.height + [0] * neg.height)
        tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v = B.ref_index(ref_df, LOCUS)
        qs, qv = qpool["cdr3"].to_list(), qpool["v"].to_list()
        sc, _ = B.vdjmatch_classify(tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v, ctrl, qs, qv, [E],
                                    1e-3, True, params=scope(5, 2, 2))
        cd = np.array([sc[q][E][0] for q in qs])
        pairs = list(zip(y.tolist(), cd.tolist()))
        roc, pr = roc_auc(pairs), pr_auc_balanced(pairs)
        # composition-contrast prior fused with NED (rank-sum); background per --bg
        if GBG is not None:
            bg_comp, kbg = GBG_comp, GBG_kmer
        else:
            bgc = oth.head(8000)["cdr3"].to_list()
            bg_comp, kbg = np.array([P2.comp(s) for s in bgc]), H.kmer_logodds(bgc, 3)
        go = P2.gauss_logodds(np.array([P2.comp(s) for s in ref_df["cdr3"].to_list()]), bg_comp)
        kref = H.kmer_logodds(ref_df["cdr3"].to_list(), 3)
        prs = np.array([go(P2.comp(s)) + (kref(s) - kbg(s)) for s in qs])
        fuse = np.argsort(np.argsort(cd)).astype(float) + np.argsort(np.argsort(prs)).astype(float)
        fp = list(zip(y.tolist(), fuse.tolist()))
        froc, fpr = roc_auc(fp), pr_auc_balanced(fp)
        cls, ceil, uni = "ok", None, {}
        if roc < 0.6 or pr < 0.6:
            uni, ceil = feature_diag(ref_df, qpool, y, cd)
            cls = classify(roc, pr, ceil, uni)
        mhc = v26.filter(pl.col("epitope") == E)["mhc_a"][0]
        rows.append(dict(epitope=E, mhc=mhc, n_pos=int(posq.height), n_neg=int(neg.height),
                         NED_roc=round(roc, 3), NED_pr=round(pr, 3),
                         fused_roc=round(froc, 3), fused_pr=round(fpr, 3), dROC=round(froc - roc, 3),
                         ceiling_roc=round(ceil, 3) if ceil else None,
                         kd_apex_auc=round(uni.get("kd_apex", 0), 3), arom_auc=round(uni.get("arom", 0), 3),
                         kmer_auc=round(uni.get("kmer", 0), 3), vlo_auc=round(uni.get("vlo", 0), 3),
                         failure_class=cls))
        print(f"  {E:12s} ROC={roc:.3f}->{froc:.3f} ({froc-roc:+.3f}) PR={pr:.3f} {cls}", file=sys.stderr)
    df = pl.DataFrame(rows).sort("NED_roc")
    df.write_csv(OUT / "sweep.tsv", separator="\t")
    print(f"\nwrote {OUT/'sweep.tsv'} ({df.height} epitopes)")
    d = df["dROC"].to_numpy()
    print(f"\nPRIOR FUSION across all {df.height} epitopes: mean dROC {d.mean():+.3f}  "
          f"min {d.min():+.3f}  #regress(<-0.02) {(d < -0.02).sum()}  #gain(>+0.02) {(d > 0.02).sum()}")
    poor = df.filter((pl.col("NED_roc") < 0.6) | (pl.col("NED_pr") < 0.6))
    print(f"\n<0.6 cases: {poor.height}")
    print(poor.group_by("failure_class").len().sort("len", descending=True))
    print(poor.select("epitope", "NED_roc", "fused_roc", "dROC", "ceiling_roc", "failure_class").head(30))


if __name__ == "__main__":
    main()
