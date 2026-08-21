"""Argmax multi-epitope detection (authoritative scheme, EXCLUSION_POLICY.md / latest decision).

Each held-out test query q (sample1 NLV, sample2 LLW/LLL, sample6 GLC/YLQ) is scored against the FULL
vdjdb2026 A*02 reference with its own exact match removed (exclude_exact LOO), and assigned its argmax
epitope across all candidate epitopes (if its first-hit enrichment is significant). Per epitope E the
confusion is: q true-E & assign==E -> TP; true-E & assign!=E -> FN; not-E & assign==E -> FP; else TN.
ROC/PR per E come from the per-E NED score (threshold-free); precision/recall/F1 from the argmax call.

    .venv/bin/python bench/detection_argmax.py [TRB|TRA]
"""
from __future__ import annotations

import bisect
import math
import sys
from collections import defaultdict
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
from benchmark import (A02, EPI, PSSM_SCALE, SOFTV_BETA, _pssm_targets, ref_index,   # noqa: E402
                       release, vgene)
from compare import TESTDATA                                         # noqa: E402
from metrics import pr_auc_balanced, roc_auc                         # noqa: E402
from vdjmatch.evalue import background, first_hit                    # noqa: E402
from vdjmatch.match import vgene as _vg                              # noqa: E402

ALPHA = 1e-3
EPIS = ["NLV", "LLW", "LLL", "YLQ", "GLC"]


def query_pool(locus):
    """Held-out test TCRs for a locus -> list of (cdr3, v, true_epitope_aa or None). LLW/LLL are TRB-only."""
    rows = []
    s1 = _bench.valid_cdr3(pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
                           .filter(pl.col("gene") == locus).select("cdr3", lab="type", v="v.segm")).unique("cdr3")
    rows += [(c, vgene(v), EPI["NLV"] if lab == "cmv" else None)
             for c, lab, v in zip(s1["cdr3"], s1["lab"], s1["v"])]
    if locus == "TRB":
        s2 = _bench.valid_cdr3(pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t")
                               .rename({"antigen.epitope": "lab"}).select("cdr3", "lab", v="v.segm")).unique("cdr3")
        rows += [(c, vgene(v), lab if lab in (EPI["LLW"], EPI["LLL"]) else None)
                 for c, lab, v in zip(s2["cdr3"], s2["lab"], s2["v"])]
    cc, vc = ("cdr3_beta_aa", "TRBV") if locus == "TRB" else ("cdr3_alpha_aa", "TRAV")
    t6 = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv").with_columns(pos=pl.col("padj") < 1e-5)
    t6 = (_bench.valid_cdr3(t6.select(cdr3=cc, v=vc, epitope="epitope_aa", pos="pos"))
          .filter(pl.col("epitope").is_in([EPI["GLC"], EPI["YLQ"]])).unique("cdr3"))   # GLC/YLQ test set only
    rows += [(c, vgene(v), ep if p else None)        # padj<1e-5 -> positive of its epitope, else negative
             for c, v, ep, p in zip(t6["cdr3"], t6["v"], t6["epitope"], t6["pos"])]
    seen = {}
    for c, v, e in rows:                                            # dedup by cdr3, keep first
        seen.setdefault(c, (v, e))
    return [(c, v, e) for c, (v, e) in seen.items()]


def annotate(locus, shortlist=False):
    """Return per-query (true, assign, ned_by_epi) for the test pool vs the vdjdb2026 A*02 reference.
    shortlist=True restricts candidate epitopes to clonotype-pMHC pairs seen in >=2 references (drops the
    singleton epitopes that otherwise steal the NED argmax via their tiny reference size)."""
    from benchmark import shortlist as _short
    ref = release("vdjdb2026").filter(pl.col("mhc_a").str.contains(A02))
    if shortlist:
        ref = _short(ref, min_refs=2)
    tgt, ref_epi, ref_v, n_epi, _, _ = ref_index(ref, locus)
    ctrl = background(locus)
    M, Ntot = len(ctrl), len(ref_epi)
    pool = query_pool(locus)
    qs, qv, true = [r[0] for r in pool], [r[1] for r in pool], [r[2] for r in pool]
    params = first_hit.scope(5, 2, 2)
    th, cc = first_hit.scan(tgt, ref_epi, ctrl, qs, target_v=ref_v, params=params, exclude_exact=True)
    pt = _pssm_targets(tgt, ref_epi, ref_v, qs, 5)
    recs = []
    for (q, vqg, t, c), tr in zip(zip(qs, qv, th, cc), true):
        cs = sorted(c)
        t1 = [(cost, e) for cost, e, *_ in t if cost <= 1]
        sig = first_hit.pvalue(t1, c, Ntot, M)["p_enrichment"] < ALPHA
        # ARGMAX score = reference-size-invariant similarity: control-calibrated closeness density with the
        # +1 control pseudocount, but NO N_e factor -> equal similarity to two epitopes scores equally, and
        # a singleton epitope (one neighbour) cannot out-score a large one with many close neighbours.
        # Spurious hits are filtered by `sig` (the enrichment significance), not by epitope size.
        dens = defaultdict(float)
        for ps, ed, he, hv in pt[q]:
            w = 1.0 if _vg.gene_family(vqg) == _vg.gene_family(hv) else SOFTV_BETA * _vg.vsim(vqg, hv)
            if w <= 0:
                continue
            nc = bisect.bisect_right(cs, ed) + 1
            dens[he] += w * math.exp(-ps / PSSM_SCALE) / nc
        a = max(dens, key=dens.get) if dens else None
        recs.append((tr, a if (sig and a) else None, dict(dens)))
    return recs


def report(locus, shortlist=False):
    recs = annotate(locus, shortlist=shortlist)
    n = len(recs)
    tag = "≥2-ref shortlist" if shortlist else "full A*02"
    print(f"\n=== argmax detection, {locus} (n={n} queries; {tag} reference) ===")
    print(f"{'epitope':8} {'n_pos':>6} {'TP':>4} {'FN':>4} {'FP':>4} {'TN':>5} "
          f"{'prec':>6} {'recall':>7} {'F1':>6} {'ROC':>6} {'balPR':>6}")
    for sh in EPIS:
        e = EPI[sh]
        pairs, tp = [], 0
        fn = fp = tn = npos = fn_none = fn_other = 0
        for true, assign, dens in recs:
            y = int(true == e)
            npos += y
            sc = dens.get(e, 0.0)
            pairs.append((y, sc))
            called = int(assign == e)
            tp += y & called
            if y and not called:
                fn += 1
                fn_none += assign is None
                fn_other += assign is not None
            fp += (not y) & called
            tn += (not y) & (not called)
        if npos == 0:
            continue
        prec = tp / (tp + fp) if tp + fp else float("nan")
        rec = tp / (tp + fn) if tp + fn else float("nan")
        f1 = 2 * prec * rec / (prec + rec) if prec + rec and prec == prec and rec == rec else float("nan")
        print(f"{sh:8} {npos:>6} {tp:>4} {fn:>4} {fp:>4} {tn:>5} {prec:>6.3f} {rec:>7.3f} {f1:>6.3f} "
              f"{roc_auc(pairs):>6.3f} {pr_auc_balanced(pairs):>6.3f}   FN: {fn_none} none / {fn_other} other-epi")


if __name__ == "__main__":
    loci = [sys.argv[1]] if len(sys.argv) > 1 and sys.argv[1] in ("TRA", "TRB") else ["TRB", "TRA"]
    for loc in loci:
        report(loc, shortlist=False)
        report(loc, shortlist=True)
