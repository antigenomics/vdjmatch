#!/usr/bin/env python3
"""Comprehensive cross-method TCR-specificity benchmark (see bench/BENCHMARK.md).

Conditions (classification -> ROC / PR / FP / F1, at each method's own significance call):
  C1 sample-LLW/LLL : sample2 LLW(pos) vs LLL(neg) and reverse;  ref = vdjdb2026 (exact-filtered); TRB
  C2 sample-NLV     : sample1 cmv(NLV+) vs control(NLV-);        ref = vdjdb2026;  TRA, TRB
  C3 tcrvdb         : padj<1e-5 (pos) vs >=1e-5 (neg), GLC+YLQ;  ref = vdjdb2026;  TRA, TRB, paired
  C4 vdjdb2025-LOO  : epitope X(pos) vs other(neg), within 2025; ref = 2025 \\ self; TRA, TRB, paired
  C5 new2026-LOO    : new-in-2026 X(pos) vs other-new(neg);      ref = vdjdb2025 (temporal); TRA, TRB, paired
  C6 olga-FP        : OLGA TRA(sample5)+TRB(sample4), 10k;       ref = vdjdb2025;  any hit = FP
Clustering (retention / purity): vdjdb2026 self-cluster (separate, methods that cluster).

Per-method per-condition: a by-epitope table and a summary row (mean +- sd across epitopes).
Phase A = framework + vdjmatch; external methods (tcrmatch / tcrdist / ERGO-II / GLIPH2 / DETECT) plug
in via predictions/<method>/<condition>_<locus>.tsv (query_id, epitope, score[, significant]).

    python bench/benchmark.py --conditions C1 C2 C6 --methods vdjmatch
"""
from __future__ import annotations

import argparse
import bisect
import math
import statistics as st
import sys
from collections import Counter, defaultdict
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench
from compare import TESTDATA, _big_epis, _top_held, load_sample
from metrics import pr_auc_balanced, roc_auc
from vdjmatch import db
from vdjmatch.evalue import background, first_hit
from vdjmatch.match import regions
from vdjmatch.match import vgene as _vg
from seqtree import Index, SearchParams

EPI = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV",
       "GLC": "GLCTLVAML", "YLQ": "YLQPRTFLL"}
_REL = {}  # release cache


def release(which: str, species="HomoSapiens") -> pl.DataFrame:
    if which not in _REL:
        if which == "vdjdb2026":
            _REL[which] = db.load(_bench.source(), species=species)
        else:                                                     # vdjdb2025: prefer cached file (no API call)
            from vdjmatch.db.vdjdb import cache_dir
            p = cache_dir(None) / "vdjdb-2025-12-29.default.txt"
            _REL[which] = (db.load(str(p), species=species) if p.exists()
                           else db.load(pin="2025-12-29", asset="default", species=species))
    return _REL[which]


def vgene(v) -> str:
    """Normalise a V gene to its allele-stripped gene name for V-matching (e.g. TRBV20-1*01 -> TRBV20-1)."""
    return (v or "").split("*")[0]


def shortlist(df: pl.DataFrame, min_refs: int = 2) -> pl.DataFrame:
    """Clonotype-pMHC pairs in >= min_refs distinct references (keeps ALL columns incl complex_id)."""
    key = ["gene", "cdr3", "v", "j", "epitope"]
    keep = (df.group_by(key).agg(pl.col("reference_id").n_unique().alias("nr"))
            .filter(pl.col("nr") >= min_refs).select(key))
    return df.join(keep, on=key, how="semi")


def ref_index(df: pl.DataFrame, locus: str):
    """VDJdb frame -> (Index, ref_epi, ref_v, n_epi, n_epi_v, n_v) for a single locus. One entry per
    unique CDR3 (representative V/epitope); the V vector enables the V+CDR3 joint E-value."""
    r = (_bench.valid_cdr3(df.filter(pl.col("gene") == locus)).group_by("cdr3")
         .agg(pl.col("epitope").first(), pl.col("v").first()))
    epi = r["epitope"].to_list()
    v = [vgene(x) for x in r["v"].to_list()]
    return (Index.build(r["cdr3"].to_list(), "aa"), epi, v,
            Counter(epi), Counter(zip(v, epi)), Counter(v))


# ---- vdjmatch classification (first-hit E-value; optional V+CDR3 joint null) ----------------------
PSSM_SCALE, SOFTV_BETA = 400.0, 0.25  # ranking-distance temperature; cross-family V-loop down-weight


def _pssm_targets(tgt, ref_epi, ref_v, queries, max_edits):
    """Per-query target hits ``(pssm_cost, n_edits, epitope, ref_v)`` scored with a length-specific
    central-CDR3 PSSM as position matrix — central substitutions (the specificity motif) cost more than
    the germline ends. Fixed width (no indels: the PSSM is per-length). Exact matches dropped."""
    by_len = defaultdict(list)
    for q in queries:
        by_len[len(q)].append(q)
    out = {}
    for L, qs in by_len.items():
        sp = SearchParams(max_subs=max_edits, max_ins=0, max_dels=0, max_total_edits=max_edits,
                          engine="seqtm")
        try:
            sp.pos_matrix = regions.significance_pssm(L)
        except Exception:
            pass
        for q, hits in zip(qs, tgt.search_batch(qs, sp, 0)):
            out[q] = [(h.score, h.n_subs + h.n_ins + h.n_dels, ref_epi[h.ref_id], ref_v[h.ref_id])
                      for h in hits if (h.n_subs + h.n_ins + h.n_dels) > 0]
    return out

def vdjmatch_classify(tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v, ctrl, queries, query_v, epitopes,
                      alpha, exclude_exact, v_mode="none", params=None):
    """queries -> per query {epitope: (-log10 p_enrichment, significant)} for the candidate epitopes,
    plus an overall-significant flag (any-epitope first-hit p < alpha) for FP estimation. ``v_mode``:
    ``none`` = CDR3-only; ``match_v`` = V+CDR3 joint E-value (same-V-restricted Poisson tail)."""
    M, Ntot = len(ctrl), len(ref_epi)
    t2 = lambda hits: [(c, e) for c, e, *_ in hits]            # drop the V tag for the V-agnostic pvalue
    th, cc = first_hit.scan(tgt, ref_epi, ctrl, queries, target_v=ref_v, params=params,
                            exclude_exact=exclude_exact)
    max_edits = params.max_total_edits if params is not None else 5
    pt = _pssm_targets(tgt, ref_epi, ref_v, queries, max_edits) if epitopes else {}
    scores, overall = {}, {}
    for q, qv, t, c in zip(queries, query_v, th, cc):
        cs_ctrl = sorted(c)
        t1 = [h for h in t if h[0] <= 1]                       # significance at radius<=1 (specific, low FP)
        # FP/significance always V-AGNOSTIC (the same-V null is tight for rare V -> over-calls).
        overall[q] = first_hit.pvalue(t2(t1), c, Ntot, M)["p_enrichment"] < alpha
        scores[q] = {}
        for e in epitopes:
            N_e = n_epi.get(e, 1)
            # RANKING = control-calibrated continuous density (NED) over PSSM-scored neighbours, each
            # weighted by germline V-loop similarity (soft-V: same family full, cross-family beta*vsim of
            # the CDR1+CDR2 loops) and closeness exp(-pssm/scale), divided by the expected local control
            # density. Continuous (no score-0 pile); down-weights neighbours in dense control regions.
            dens = 0.0
            for ps, ed, he, hv in pt[q]:
                if he != e:
                    continue
                w = 1.0 if _vg.gene_family(qv) == _vg.gene_family(hv) else SOFTV_BETA * _vg.vsim(qv, hv)
                if w <= 0:
                    continue
                nc = bisect.bisect_right(cs_ctrl, ed)
                dens += w * math.exp(-ps / PSSM_SCALE) / max((N_e / M) * nc, 0.01)
            p_sig = first_hit.pvalue(t2(t1), c, n_epi.get(e, 1), M, epitope=e)["p_enrichment"]
            scores[q][e] = (dens, p_sig < alpha)
    return scores, overall


def read_predictions(path: Path, allq, epitopes):
    """predictions/<method>/<cond>_<locus>.tsv (query_id, epitope, score[, significant]) -> scores dict
    {query: {epitope: (score, significant)}} with 0/False defaults for unscored (query, epitope)."""
    scores = {q: {e: (0.0, False) for e in epitopes} for q in allq}
    if not path.exists():
        return None
    df = pl.read_csv(path, separator="\t")
    has_sig = "significant" in df.columns
    for r in df.iter_rows(named=True):
        q, e = r["query_id"], r["epitope"]
        if q in scores and e in epitopes:
            scores[q][e] = (float(r["score"]), bool(r["significant"]) if has_sig else r["score"] > 0)
    return scores


def classify_metrics(scores, pos, neg, epi):
    """ROC / PR (threshold-free) + FP / F1 (at the method's significance call) for one epitope."""
    pairs = [(1, scores[q][epi][0]) for q in pos] + [(0, scores[q][epi][0]) for q in neg]
    tp = sum(scores[q][epi][1] for q in pos)
    fp = sum(scores[q][epi][1] for q in neg)
    fn, tn = len(pos) - tp, len(neg) - fp
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    return {"roc_auc": roc_auc(pairs), "pr_auc": pr_auc_balanced(pairs),
            "fp": fp / (fp + tn) if fp + tn else 0.0,
            "f1": 2 * prec * rec / (prec + rec) if prec + rec else 0.0,
            "tp": tp, "fn": fn, "fp_n": fp, "tn": tn, "precision": prec, "recall": rec,
            "n_pos": len(pos), "n_neg": len(neg)}


# ---- conditions: each yields (locus, ref_df, queries, [(epitope, pos_ids, neg_ids)], exclude_exact)-
def cond_sample_pair(name, loci):
    """C1/C2: discriminate two query groups by matching against vdjdb2026."""
    for locus in loci:
        if name == "C1":                                          # sample2: LLW vs LLL by antigen.epitope
            d = (pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t")
                 .rename({"antigen.epitope": "label"}).select("cdr3", "label", v="v.segm"))
            d = _bench.valid_cdr3(d).unique("cdr3")
            groups = {EPI["LLW"]: d.filter(pl.col("label") == EPI["LLW"])["cdr3"].to_list(),
                      EPI["LLL"]: d.filter(pl.col("label") == EPI["LLL"])["cdr3"].to_list()}
            tasks = [(EPI["LLW"], groups[EPI["LLW"]], groups[EPI["LLL"]]),
                     (EPI["LLL"], groups[EPI["LLL"]], groups[EPI["LLW"]])]
        else:                                                     # C2 sample1: cmv(NLV+) vs control(NLV-)
            d = (pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
                 .filter(pl.col("gene") == locus).select("cdr3", label="type", v="v.segm"))
            d = _bench.valid_cdr3(d).unique("cdr3")
            pos = d.filter(pl.col("label") == "cmv")["cdr3"].to_list()
            neg = d.filter(pl.col("label") == "control")["cdr3"].to_list()
            tasks = [(EPI["NLV"], pos, neg)]
        qv = {c: vgene(v) for c, v in zip(d["cdr3"], d["v"])}
        yield locus, release("vdjdb2026"), tasks, True, qv


def cond_olga_fp(loci, n=1000):
    """C6: OLGA negatives (sample5=TRA, sample4=TRB) vs vdjdb2025 -> any significant hit is a FP."""
    files = {"TRA": "sample5", "TRB": "sample4"}
    for locus in loci:
        d = (pl.read_csv(TESTDATA / f"{files[locus]}_olga_airr.txt", separator="\t")
             .select(cdr3="junction_aa", v="v_gene").pipe(_bench.valid_cdr3).unique("cdr3"))
        if d.height > n:
            d = d.sample(n, seed=0)
        qv = {c: vgene(v) for c, v in zip(d["cdr3"], d["v"])}
        yield locus, release("vdjdb2025"), d["cdr3"].to_list(), qv


def cond_tcrvdb(loci):
    """C3: TCRvdb padj<1e-5 (pos) vs >=1e-5 (neg) per epitope, annotated vs vdjdb2026."""
    t = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv").with_columns(pos=pl.col("padj") < 1e-5)
    chain_col = {"TRA": "cdr3_alpha_aa", "TRB": "cdr3_beta_aa"}
    v_col = {"TRA": "TRAV", "TRB": "TRBV"}
    for locus in loci:
        d = (t.select(cdr3=chain_col[locus], epitope="epitope_aa", pos="pos", v=v_col[locus])
             .pipe(_bench.valid_cdr3).unique("cdr3"))
        qv = {c: vgene(v) for c, v in zip(d["cdr3"], d["v"])}
        tasks = []
        for e in (EPI["GLC"], EPI["YLQ"]):
            de = d.filter(pl.col("epitope") == e)
            pos = de.filter("pos")["cdr3"].to_list()
            neg = de.filter(~pl.col("pos"))["cdr3"].to_list()
            if pos and neg:
                tasks.append((e, pos, neg))
        yield locus, release("vdjdb2026"), tasks, True, qv


def cond_loo(which, loci, top, min_epi, max_q):
    """C4 (within vdjdb2025) / C5 (new-in-2026 vs 2025 reference): per-epitope X(pos) vs other(neg).

    Both releases are restricted to the SHORTLIST: clonotype-pMHC pairs seen in >=2 references
    (``db.replicated``), per request — the public/learnable subset.
    """
    v25 = shortlist(release("vdjdb2025"))
    for locus in loci:
        ref_cell = _bench.valid_cdr3(v25.filter(pl.col("gene") == locus))
        if which == "C4":
            test_cell, excl = _bench.long_list(v25.filter(pl.col("gene") == locus), cap=max_q,
                                               min_n=min_epi), True
        else:                                                     # C5: 2026 shortlist clonotypes new vs 2025
            v26 = _bench.long_list(shortlist(release("vdjdb2026"))
                                   .filter(pl.col("gene") == locus), cap=max_q, min_n=min_epi)
            test_cell = v26.join(ref_cell.select("cdr3", "epitope").unique(),
                                 on=["cdr3", "epitope"], how="anti")
            excl = False
        held = _top_held(test_cell, _big_epis(ref_cell, min_epi), top, min_epi)
        per = {e: test_cell.filter(pl.col("epitope") == e).unique("cdr3")["cdr3"].to_list()[:max_q]
               for e in held}
        tasks = [(e, per[e], [q for e2 in held if e2 != e for q in per[e2]]) for e in held]
        qv = {c: vgene(v) for c, v in zip(test_cell["cdr3"], test_cell["v"])}
        yield locus, v25, tasks, excl, qv


# ---- Group A: epitope-specific detection (per dataset_descr.md) -----------------------------------
A02 = r"A\*?0?2"                                                  # HLA-A*02 / A2 match


def epi_ref(epitope_str, locus):
    """VDJdb2026 records of ONE epitope, HLA-A*02 -> reference frame for ref_index (the match target)."""
    v = release("vdjdb2026")
    return v.filter((pl.col("epitope") == epitope_str) & pl.col("mhc_a").str.contains(A02))


def detection_tasks(loci):
    """Yield (task, epitope_str, locus, pos_qv, neg_qv) for Group-A detection: match each query against
    the epitope's A*02 reference; positives should hit, the control negatives should not."""
    E = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV", "YLQ": "YLQPRTFLL", "GLC": "GLCTLVAML"}
    s1 = (pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
          .filter(pl.col("gene") == "TRB").select("cdr3", label="type", v="v.segm"))
    s1 = _bench.valid_cdr3(s1).unique("cdr3")
    qv1 = lambda lab: {c: vgene(v) for c, l, v in zip(s1["cdr3"], s1["label"], s1["v"]) if l == lab}
    yield ("NLV", E["NLV"], "TRB", qv1("cmv"), qv1("control"))           # NLV+ vs tet- control

    s2 = (pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t")
          .rename({"antigen.epitope": "label"}).select("cdr3", "label", v="v.segm"))
    s2 = _bench.valid_cdr3(s2).unique("cdr3")
    qv2 = lambda lab: {c: vgene(v) for c, l, v in zip(s2["cdr3"], s2["label"], s2["v"]) if l == lab}
    llw, lll = qv2(E["LLW"]), qv2(E["LLL"])
    yield ("LLW", E["LLW"], "TRB", llw, lll)                             # LLW vs LLL (each other's control)
    yield ("LLL", E["LLL"], "TRB", lll, llw)

    t = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv").with_columns(pos=pl.col("padj") < 1e-5)
    chain = {"TRA": ("cdr3_alpha_aa", "TRAV"), "TRB": ("cdr3_beta_aa", "TRBV")}
    for tk in ("YLQ", "GLC"):
        for locus in [x for x in loci if x in chain]:
            cc, vc = chain[locus]
            d = (t.filter(pl.col("epitope_aa") == E[tk]).select(cdr3=cc, v=vc, pos="pos")
                 .pipe(_bench.valid_cdr3).unique("cdr3"))
            pos = {c: vgene(v) for c, v, p in zip(d["cdr3"], d["v"], d["pos"]) if p}
            neg = {c: vgene(v) for c, v, p in zip(d["cdr3"], d["v"], d["pos"]) if not p}
            if pos and neg:
                yield (tk, E[tk], locus, pos, neg)                      # padj<1e-5 vs >=1e-5


def summarize(rows):
    """rows: (condition, locus, method, epitope, metric, value) -> by-epitope df + mean+-sd summary."""
    by = pl.DataFrame(rows, schema=["condition", "locus", "method", "epitope", "metric", "value"],
                      orient="row")
    summ = (by.group_by(["condition", "locus", "method", "metric"])
            .agg(pl.col("value").mean().round(3).alias("mean"),
                 pl.col("value").std().round(3).alias("sd"),
                 pl.col("epitope").n_unique().alias("n_epi"))
            .sort(["condition", "locus", "method", "metric"]))
    return by, summ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conditions", nargs="+", default=["C1", "C2", "C3", "C6"])
    ap.add_argument("--methods", nargs="+", default=["vdjmatch"])
    ap.add_argument("--loci", nargs="+", default=["TRA", "TRB"])
    ap.add_argument("--v-mode", default="match_v", choices=["none", "match_v"],
                    help="vdjmatch V+CDR3 joint E-value (same-V-restricted null)")
    ap.add_argument("--scope", default="5,2,2", help="first-hit scope max_edits,max_ins,max_dels")
    ap.add_argument("--alpha", type=float, default=1e-3)
    ap.add_argument("--olga-n", type=int, default=1000, help="Group B: OLGA negatives per locus")
    ap.add_argument("--pred-dir", default="bench/predictions")
    ap.add_argument("--out", default="bench/out/benchmark")
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rows, fp_rows, n_rows = [], [], []
    params = first_hit.scope(*[int(x) for x in args.scope.split(",")])

    conf_rows = []
    olga = {locus: (q, qv) for locus, _, q, qv in cond_olga_fp(args.loci, args.olga_n)}  # noise queries

    # Group B noise: OLGA matched against ALL of VDJdb2026 (any significant hit = FP), per locus.
    noise_fp = {}
    for locus in args.loci:
        olq, olqv = olga.get(locus, ([], {}))
        if not olq or "vdjmatch" not in args.methods:
            continue
        ft, fe, fv, fne, fnev, fnv = ref_index(release("vdjdb2026"), locus)
        _, ov = vdjmatch_classify(ft, fe, fv, fne, fnev, fnv, background(locus), olq,
                                  [olqv[q] for q in olq], [], args.alpha, True,
                                  v_mode=args.v_mode, params=params)
        noise_fp[(locus, "vdjmatch")] = sum(ov.values()) / len(olq)
        print(f"noise/{locus}/vdjmatch: OLGA-FP-vs-all={noise_fp[(locus, 'vdjmatch')]*100:.2f}% "
              f"(any hit in full VDJdb2026, n={len(olq)})")

    for task, estr, locus, pos_qv, neg_qv in detection_tasks(args.loci):
        ref = epi_ref(estr, locus)
        allq, qv = list(pos_qv) + list(neg_qv), {**pos_qv, **neg_qv}
        n_rows.append((task, locus, len(pos_qv), len(neg_qv)))
        tgt = None
        for method in args.methods:
            if method == "vdjmatch":
                if tgt is None:
                    tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v = ref_index(ref, locus)
                scores, _ = vdjmatch_classify(tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v,
                                              background(locus), allq, [qv[q] for q in allq], [estr],
                                              args.alpha, True, v_mode=args.v_mode, params=params)
            else:
                scores = read_predictions(Path(args.pred_dir) / method / f"{task}_{locus}.tsv", allq, [estr])
                if scores is None:
                    print(f"  skip {method}/{task}/{locus}: no predictions")
                    continue
            m = classify_metrics(scores, list(pos_qv), list(neg_qv), estr)
            for k in ("roc_auc", "pr_auc", "fp", "f1"):
                rows.append((task, locus, method, estr, k, m[k]))
            nfp = noise_fp.get((locus, method))
            conf_rows.append((task, locus, method, m["tp"], m["fn"], m["fp_n"], m["tn"],
                              round(m["recall"], 3), round(m["precision"], 3),
                              round(nfp, 4) if nfp is not None else None))
            print(f"{task}/{locus}/{method}: ROC={m['roc_auc']:.3f} PR={m['pr_auc']:.3f} "
                  f"TP={m['tp']} FN={m['fn']} FP={m['fp_n']} TN={m['tn']}")

    if n_rows:
        pl.DataFrame(n_rows, schema=["task", "locus", "n_pos", "n_neg"],
                     orient="row").write_csv(out / "dataset_n.tsv", separator="\t")
    if conf_rows:
        pl.DataFrame(conf_rows, schema=["task", "locus", "method", "tp", "fn", "fp", "tn",
                                        "recall", "precision", "olga_fp"], orient="row").write_csv(
            out / "confusion.tsv", separator="\t")
    if rows:
        by, summ = summarize(rows)
        by.write_csv(out / "by_epitope.tsv", separator="\t")
        summ.write_csv(out / "summary.tsv", separator="\t")
        with pl.Config(tbl_rows=100):
            print("\n=== detection summary (ROC/PR/F1/FP) ===\n", summ)


if __name__ == "__main__":
    main()
