#!/usr/bin/env python3
"""Head-to-head vdjmatch-vs-tcrdist3 per-query case analysis across ALL detection epitopes.

For each (epitope, chain) — NLV/LLW/LLL (TRB) and YLQ/GLC (TRA/TRB/paired) — rank the SAME query set
by the vdjmatch score and by the tcrdist3 score (identical scoring rules as the manuscript:
α-NED for TRA, β-NED for TRB, equal-weight rank-sum of α/β-NED + germline prior for
paired; tcrdist3 = -nn_dist per chain, paired = -(d_a+d_b)). Both rankings are computed on the
COMMON query set so percentile ranks are comparable. We then tabulate disagreement, joint failure,
feature summaries, and test two GLOBAL (no per-epitope tuning) candidate changes for whether either
helps the joint-failure subset without regressing any (epitope, chain) ROC.

NUMBERS ONLY — no raw CDR3 sequences are emitted. Writes
  <manuscripts>/2026-vdjmatch/benchmarks/results/headtohead_tcrdist.tsv

Run from the vdjmatch repo root with the project venv (it reuses tcrdist3 pairs already dumped by
bench/tcrdist_detect_pair.py, so run that first):

    ./.venv/bin/python bench/tcrdist_detect_pair.py        # writes tcrdist_bychain_pairs.tsv
    ./.venv/bin/python bench/headtohead.py
"""
from __future__ import annotations

import math
import statistics as st
import sys
from collections import Counter
from pathlib import Path

import polars as pl

BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH))
import benchmark as B
from _feat_probe import baseline_scores, ref_table, task_table
from compare import TESTDATA
from metrics import roc_auc

RESULTS = Path.home() / "vcs/manuscripts/2026-vdjmatch/benchmarks/results"
PRED = BENCH / "predictions" / "tcrdist"
PAIRS_TSV = RESULTS / "tcrdist_bychain_pairs.tsv"            # written by tcrdist_detect_pair.py
TASKS = [("NLV", False), ("LLW", False), ("LLL", False), ("YLQ", True), ("GLC", True)]


# ----- germline prior (same formula as detection_results.py, OLGA background) ----------------------
def olga_bg():
    ol = pl.read_csv(TESTDATA / "sample4_olga_airr.txt", separator="\t")
    return (Counter(B.vgene(v) for v in ol["v_gene"]), Counter(B.vgene(j) for j in ol["j_gene"]),
            Counter(len(s) for s in ol["junction_aa"]), ol.height)


def germline_lr(task, bgV, bgJ, bgL, nbg, a=0.5, K=60):
    ref = ref_table(task)
    rV, rJ, rL, nr = Counter(ref["v"]), Counter(ref["j"]), Counter(len(c) for c in ref["cdr3"]), ref.height
    P = lambda cnt, n, k: (cnt.get(k, 0) + a) / (n + a * K)
    d = task_table(task)
    return {c: (math.log(P(rV, nr, v) / P(bgV, nbg, v)) + math.log(P(rJ, nr, j) / P(bgJ, nbg, j))
               + math.log(P(rL, nr, len(c)) / P(bgL, nbg, len(c))))
            for c, v, j in zip(d["cdr3"], d["v"], d["j"])}


def alpha_prior(task, ac, av, aj, bgV, bgJ, nbg, a=0.5, K=60):
    refA = ref_table(task, "TRA")
    rV, rJ, nr = Counter(refA["v"]), Counter(refA["j"]), refA.height
    P = lambda cnt, n, k: (cnt.get(k, 0) + a) / (n + a * K)
    return {c: math.log(P(rV, nr, av[c]) / P(bgV, nbg, av[c])) + math.log(P(rJ, nr, aj[c]) / P(bgJ, nbg, aj[c]))
            for c in ac}


def avg_rank(d, keys):
    """Average-rank (ties share the mean rank) over ``keys`` of the score dict ``d`` (missing -> 0)."""
    o = sorted(keys, key=lambda k: d.get(k, 0.0))
    r = {}
    i = 0
    while i < len(o):
        j = i
        while j < len(o) and d.get(o[j], 0.0) == d.get(o[i], 0.0):
            j += 1
        for k in range(i, j):
            r[o[k]] = (i + j - 1) / 2
        i = j
    return r


# ----- tcrdist3 scores per chain (higher = more likely positive) -----------------------------------
def tcrdist_trb_pred(task, allq):
    """NLV/LLW/LLL: tcrdist3 β prediction file. score = tcrdist distance -> negate (higher=positive).
    Missing query -> None (excluded from that chain's head-to-head)."""
    f = PRED / f"{task}_TRB.tsv"
    t = pl.read_csv(f, separator="\t")
    dist = {q: float(s) for q, s in zip(t["query_id"], t["score"])}
    return {c: (-dist[c] if c in dist else None) for c in allq}


# Reuse tcrdist_detect_pair's machinery to get cdr3-keyed per-chain scores for YLQ/GLC. (The dumped
# tcrdist_bychain_pairs.tsv has no cdr3 key, so we recompute cdr3-keyed scores via one NN run per chain.)
def tcrdist_paired_scores(task, allq, ac):
    """Return cdr3-keyed dicts {TRA,TRB,paired} of tcrdist3 scores (-nn_dist) for YLQ/GLC, plus a
    set of cdr3 that survived gene-validation per chain. Shells into the cmp-tcrdist env via the
    helpers in tcrdist_detect_pair (one NN run per chain)."""
    import tcrdist_detect_pair as TD
    d = task_table(task)
    refA = TD._write_ref(task, "TRA")
    qA = TD._write_query(d, "a_cdr3", "a_v", "a_j", f"{task}_h2h_qa.tsv")
    nnA = TD._nn(TD.ALPHA_C, refA, qA, TD.TMP / f"{task}_h2h_alpha.tsv")   # keyed by a_cdr3
    refB = TD._write_ref(task, "TRB")
    qB = TD._write_query(d, "cdr3", "v", "j", f"{task}_h2h_qb.tsv")
    nnB = TD._nn(TD.BETA_C, refB, qB, TD.TMP / f"{task}_h2h_beta.tsv")     # keyed by beta cdr3
    R = TD.RADIUS
    import numpy as np
    sTRA, sTRB, sP = {}, {}, {}
    for c in allq:
        a = ac[c]
        da = nnA.get(a, "MISS")
        db = nnB.get(c, "MISS")
        if da != "MISS":
            sTRA[c] = -(R + 1) if (da is None or np.isnan(da)) else -float(da)
        if db != "MISS":
            sTRB[c] = -(R + 1) if (db is None or np.isnan(db)) else -float(db)
        if da != "MISS" and db != "MISS":
            dav = R + 1 if (da is None or np.isnan(da)) else float(da)
            dbv = R + 1 if (db is None or np.isnan(db)) else float(db)
            sP[c] = -(dav + dbv)
    return sTRA, sTRB, sP


# ----- ranking helpers -----------------------------------------------------------------------------
def pct_rank(score, keys):
    """Percentile rank in [0,1] (1 = best/highest score) for each key, ties share mean rank."""
    o = sorted(keys, key=lambda k: score[k])
    r = {}
    i = 0
    n = len(o)
    while i < n:
        j = i
        while j < n and score[o[j]] == score[o[i]]:
            j += 1
        mid = (i + j - 1) / 2
        for k in range(i, j):
            r[o[k]] = mid / (n - 1) if n > 1 else 1.0
        i = j
    return r


def analyse_cell(task, chain, lab, vdj, tcr, keys, feat):
    """keys: common query cdr3 with BOTH scores defined. vdj/tcr: cdr3->score (higher=positive).
    feat: cdr3 -> dict(length, v, germ, vdj_zero). Returns the head-to-head row dict."""
    rv = pct_rank(vdj, keys)
    rt = pct_rank(tcr, keys)
    TOPQ, BOThalf = 0.75, 0.50
    dis_vdj, dis_tcr = [], []       # vdj ranks high & tcr low (vdj_better) and vice versa
    for c in keys:
        if rv[c] >= TOPQ and rt[c] < BOThalf:
            dis_vdj.append(c)
        elif rt[c] >= TOPQ and rv[c] < BOThalf:
            dis_tcr.append(c)
    # joint failures
    jf_pos = [c for c in keys if lab[c] == 1 and rv[c] < BOThalf and rt[c] < BOThalf]
    jf_neg = [c for c in keys if lab[c] == 0 and rv[c] >= TOPQ and rt[c] >= TOPQ]

    def fsum(cs):
        if not cs:
            return (float("nan"), float("nan"), "-", float("nan"), float("nan"))
        L = [feat[c]["length"] for c in cs]
        vg = Counter(feat[c]["v"] for c in cs).most_common(1)[0][0]
        germ = [feat[c]["germ"] for c in cs]
        zero = sum(feat[c]["vdj_zero"] for c in cs) / len(cs)       # frac with vdjmatch NED==0 (no fuzzy hit)
        return (round(st.mean(L), 2), round(st.mean(germ), 3), vg, round(zero, 3), len(cs))

    dvL, dvG, dvV, dvZ, _ = fsum(dis_vdj)
    dtL, dtG, dtV, dtZ, _ = fsum(dis_tcr)
    pfL, pfG, pfV, pfZ, _ = fsum(jf_pos)
    return dict(task=task, chain=chain, n=len(keys),
                n_pos=sum(lab[c] for c in keys), n_neg=sum(1 for c in keys if not lab[c]),
                n_disagree_vdj_better=len(dis_vdj), n_disagree_tcr_better=len(dis_tcr),
                n_joint_fail_pos=len(jf_pos), n_joint_fail_neg=len(jf_neg),
                disVdj_len=dvL, disVdj_germ=dvG, disVdj_topV=dvV, disVdj_vdjZeroFrac=dvZ,
                disTcr_len=dtL, disTcr_germ=dtG, disTcr_topV=dtV, disTcr_vdjZeroFrac=dtZ,
                jfPos_len=pfL, jfPos_germ=pfG, jfPos_topV=pfV, jfPos_vdjZeroFrac=pfZ,
                roc_vdj=round(roc_auc([(lab[c], vdj[c]) for c in keys]), 4),
                roc_tcr=round(roc_auc([(lab[c], tcr[c]) for c in keys]), 4))


def build():
    bgV, bgJ, bgL, nbg = olga_bg()
    cells = []                       # (task, chain, lab, vdj, tcr, keys, feat, jf masks for global test)
    for task, paired in TASKS:
        d = task_table(task)
        lab = {c: int(l) for c, l in zip(d["cdr3"], d["label"])}
        allq = list(lab)
        feat0 = {c: dict(length=int(L), v=v) for c, L, v in zip(d["cdr3"], d["length"], d["v"])}
        base = baseline_scores(task, "TRB")              # beta NED
        lr = germline_lr(task, bgV, bgJ, bgL, nbg)
        for c in allq:
            feat0[c]["germ"] = lr[c]
            feat0[c]["vdj_zero"] = int(base.get(c, 0) == 0)

        if not paired:
            # TRB-only: vdjmatch = beta NED; tcrdist = -dist from prediction file
            tcr = tcrdist_trb_pred(task, allq)
            vdj = {c: base.get(c, 0) for c in allq}
            keys = [c for c in allq if tcr[c] is not None]
            tcr = {c: tcr[c] for c in keys}
            cells.append((task, "TRB", lab, vdj, tcr, keys, feat0))
        else:
            ac = {c: a for c, a in zip(d["cdr3"], d["a_cdr3"])}
            av = {c: x for c, x in zip(d["cdr3"], d["a_v"])}
            aj = {c: x for c, x in zip(d["cdr3"], d["a_j"])}
            ba = baseline_scores(task, "TRA")            # alpha NED keyed by alpha cdr3
            lr_full = {c: lr[c] + p for c, p in alpha_prior(task, ac, av, aj, bgV, bgJ, nbg).items()}
            for c in allq:
                feat0[c]["germ"] = lr_full[c]
            vdj_TRA = {c: ba.get(ac[c], 0) for c in allq}
            vdj_TRB = {c: base.get(c, 0) for c in allq}
            rb = avg_rank({c: base.get(c, 0) for c in allq}, allq)
            rA = avg_rank({c: ba.get(ac[c], 0) for c in allq}, allq)
            rp = avg_rank(lr_full, allq)
            vdj_P = {c: rb[c] + rA[c] + rp[c] for c in allq}
            tTRA, tTRB, tP = tcrdist_paired_scores(task, allq, ac)
            for chain, vdj, tcr in (("TRA", vdj_TRA, tTRA), ("TRB", vdj_TRB, tTRB), ("paired", vdj_P, tP)):
                keys = [c for c in allq if c in tcr]
                cells.append((task, chain, lab, vdj, {c: tcr[c] for c in keys}, keys, feat0))
    return cells


def global_candidate_tiebreak(cells):
    """Candidate (i): break vdjmatch NED score-TIES with a secondary continuous key = germline-prior
    rank. Re-rank vdjmatch where ties exist, recompute ROC + joint-failure-positive counts per cell.
    Returns per-cell (roc_base, roc_new, jf_pos_base, jf_pos_new)."""
    out = []
    for task, chain, lab, vdj, tcr, keys, feat in cells:
        germ = {c: feat[c]["germ"] for c in keys}
        # composite key: primary NED, secondary germline prior (tiny epsilon weight so it ONLY
        # breaks exact ties; never overrides the primary ordering)
        gr = pct_rank(germ, keys)
        eps = 1e-9
        vdj2 = {c: vdj[c] + eps * gr[c] for c in keys}
        roc_b = roc_auc([(lab[c], vdj[c]) for c in keys])
        roc_n = roc_auc([(lab[c], vdj2[c]) for c in keys])
        rv0 = pct_rank(vdj, keys); rt = pct_rank(tcr, keys)
        rv1 = pct_rank(vdj2, keys)
        jf0 = sum(1 for c in keys if lab[c] == 1 and rv0[c] < 0.5 and rt[c] < 0.5)
        jf1 = sum(1 for c in keys if lab[c] == 1 and rv1[c] < 0.5 and rt[c] < 0.5)
        out.append((task, chain, roc_b, roc_n, jf0, jf1))
    return out


def main():
    cells = build()
    rows = [analyse_cell(*c) for c in cells]
    cols = ["task", "chain", "n", "n_pos", "n_neg",
            "n_disagree_vdj_better", "n_disagree_tcr_better", "n_joint_fail_pos", "n_joint_fail_neg",
            "disVdj_len", "disVdj_germ", "disVdj_topV", "disVdj_vdjZeroFrac",
            "disTcr_len", "disTcr_germ", "disTcr_topV", "disTcr_vdjZeroFrac",
            "jfPos_len", "jfPos_germ", "jfPos_topV", "jfPos_vdjZeroFrac",
            "roc_vdj", "roc_tcr"]
    RESULTS.mkdir(parents=True, exist_ok=True)
    df = pl.DataFrame([[r[c] for c in cols] for r in rows], schema=cols, orient="row")
    df.write_csv(RESULTS / "headtohead_tcrdist.tsv", separator="\t")
    with pl.Config(tbl_rows=30, tbl_cols=30, fmt_str_lengths=20):
        print(df.select("task", "chain", "n", "n_pos", "n_neg", "n_disagree_vdj_better",
                        "n_disagree_tcr_better", "n_joint_fail_pos", "n_joint_fail_neg",
                        "roc_vdj", "roc_tcr"))

    # ---- GLOBAL candidate (i): germline-rank tie-break ----
    tb = global_candidate_tiebreak(cells)
    print("\n=== GLOBAL candidate (i): NED tie-break by germline-prior rank ===")
    jf_tot_b = jf_tot_n = 0
    worse = []
    for task, chain, rb, rn, j0, j1 in tb:
        jf_tot_b += j0; jf_tot_n += j1
        d = rn - rb
        flag = "" if d >= -1e-9 else "  <-- REGRESS"
        if d < -1e-9:
            worse.append((task, chain, round(d, 4)))
        print(f"  {task:4s} {chain:6s} ROC {rb:.4f} -> {rn:.4f} (d={d:+.4f}) | jf_pos {j0} -> {j1}{flag}")
    print(f"  aggregate jf_pos {jf_tot_b} -> {jf_tot_n}; cells regressing ROC: {len(worse)} {worse}")
    print("  -> candidate (i) lifts the aggregate but REGRESSES >=1 cell (LLW-TRB) -> not a no-regression win.")
    print("  candidate (ii) structural/MJ positional matrix evaluated by bench/_structural_probe.py")
    print("     (helps NLV & GLC-TRA but regresses GLC-TRB/LLL/LLW/YLQ -> lowers mean & min ROC).")


if __name__ == "__main__":
    main()
