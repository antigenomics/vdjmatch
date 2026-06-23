#!/usr/bin/env python3
"""Case-by-case vdjmatch vs tcrdist on the C1/C2 discrimination tasks.

For every query, classify both methods as right/wrong and bucket: both_right, both_wrong,
vdjmatch_only, tcrdist_only. For the cases tcrdist gets right but vdjmatch misses (and vice versa),
dump the TCR (cdr3, V) plus vdjmatch's first-hit diagnostics — nearest same-(true-)epitope target hit
radius, how many same-V same-epitope targets and control sequences sit at that radius, and the E/p — so
we can see *why* vdjmatch missed (no same-V neighbour / control too dense / competing epitope nearer).

    python bench/failure_analysis.py --condition C1 --locus TRB --v-mode match_v
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark as B
from vdjmatch.evalue import background, first_hit


def tcrdist_scores(cond, locus, epitopes):
    p = Path("bench/predictions/tcrdist") / f"{cond}_{locus}.tsv"
    sc = {}
    if p.exists():
        for r in pl.read_csv(p, separator="\t").iter_rows(named=True):
            sc.setdefault(r["query_id"], {})[r["epitope"]] = float(r["score"])
    return sc


def diag(t, c, qv, e, n_epi_v, n_v, M, match_v):
    """vdjmatch first-hit diagnostics for query toward true epitope e."""
    N = n_epi_v.get((qv, e), 1) if match_v else 1
    d = first_hit.pvalue_v(t, c, qv, N, M, epitope=e, match_v=match_v)
    # nearest target of ANY epitope (what vdjmatch's first hit actually is) and its epitope
    near = t[0] if t else None
    return d, (near[1] if near else None), (near[0] if near else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--condition", default="C1")
    ap.add_argument("--locus", default="TRB")
    ap.add_argument("--v-mode", default="match_v", choices=["none", "match_v"])
    ap.add_argument("--scope", default="5,2,2")
    ap.add_argument("--show", type=int, default=25, help="rows to print per bucket")
    args = ap.parse_args()
    params = first_hit.scope(*[int(x) for x in args.scope.split(",")])
    match_v = args.v_mode == "match_v"

    cells = (B.cond_sample_pair(args.condition, [args.locus]) if args.condition in ("C1", "C2")
             else B.cond_tcrvdb([args.locus]))
    locus, ref_df, tasks, excl, qv = next(cells)
    tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v = B.ref_index(ref_df, locus)
    ctrl = background(locus); M = len(ctrl)
    epitopes = sorted({e for e, _, _ in tasks})
    allq = sorted({q for _, p, n in tasks for q in (*p, *n)})
    # truth: C1 -> epitope (which pos-group); C2 -> "pos"/"neg" (cmv vs control) on the single epitope
    single = len(tasks) == 1
    if single:
        e0, pos0, neg0 = tasks[0]
        truth = {q: "pos" for q in pos0} | {q: "neg" for q in neg0}
    else:
        truth = {q: e for e, pos, _ in tasks for q in pos}

    vscores, _ = B.vdjmatch_classify(tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v, ctrl, allq,
                                     [qv[q] for q in allq], epitopes, 1e-3, excl, v_mode=args.v_mode,
                                     params=params)
    th, cc = first_hit.scan(tgt, ref_epi, ctrl, allq, target_v=ref_v, params=params, exclude_exact=excl)
    thb, ccb = dict(zip(allq, th)), dict(zip(allq, cc))
    td = tcrdist_scores(args.condition, locus, epitopes)

    rows = []
    for q in allq:
        tr = truth[q]
        if single:                                                # C2: significance call on the single epitope
            v_call = "pos" if vscores[q][e0][1] else "neg"
            td_s = td.get(q, {}).get(e0, 0.0)
            td_call = "pos" if td_s > 0 else "neg"
        else:                                                     # C1: argmax epitope
            v_call = max(epitopes, key=lambda e: vscores[q][e][0])
            td_call = max(epitopes, key=lambda e: td.get(q, {}).get(e, 0.0)) if q in td else "?"
        vr, tdr = v_call == tr, td_call == tr
        bucket = ("both_right" if vr and tdr else "both_wrong" if not vr and not tdr
                  else "vdjmatch_only" if vr else "tcrdist_only")
        rows.append((q, qv[q], tr, v_call, td_call, bucket))

    df = pl.DataFrame(rows, schema=["cdr3", "v", "truth", "vdj", "tcrdist", "bucket"], orient="row")
    print(f"\n{args.condition}/{locus}  v_mode={args.v_mode}  (n={df.height})")
    print(df.group_by("bucket").len().sort("len", descending=True))

    # the cases tcrdist gets right and vdjmatch misses -> dump with diagnostics
    for b in ("tcrdist_only", "vdjmatch_only"):
        sub = df.filter(pl.col("bucket") == b)
        print(f"\n=== {b}: {sub.height} cases (showing {min(args.show, sub.height)}) ===")
        for r in sub.head(args.show).iter_rows(named=True):
            q = r["cdr3"]
            true_e = e0 if single else r["truth"]
            d, near_e, near_r = diag(thb[q], ccb[q], r["v"], true_e, n_epi_v, n_v, M, match_v)
            print(f"  {q:20} V={r['v']:10} truth={r['truth']:9} vdj={r['vdj']:9} td={r['tcrdist']:9} "
                  f"| nearest-true r={d['radius']} n_t={d['n_target']} n_ctrl={d['n_control']} "
                  f"E={d['E']:.2g} p={d['p_enrichment']:.1e} | actual-1st-hit={near_e}@{near_r}")


if __name__ == "__main__":
    main()
