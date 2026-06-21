#!/usr/bin/env python3
"""Leave-one-out-by-epitope evaluation of VDJAM (non-circular).

VDJAM is derived from same-antigen substitution pairs, so testing it on the same epitopes would be
circular. For each held-out epitope E*, we re-derive the NDN substitution matrix from **all other**
epitopes, then run VDJdb-vs-VDJdb retrieval on E*: query each E* CDR3 against the whole chain (exact
self excluded), score hits with each matrix, and measure how well same-E* hits rank above
other-epitope hits (micro PR-AUC over (query, hit) pairs). We compare VDJAM vs BLOSUM62 vs unit cost.

    python bench/loo_vdjam.py --chain TRB --top 8 --scope 2,0,0,2
"""
from __future__ import annotations

import argparse
import math

import polars as pl
from seqtree import Index, SearchParams, SubstitutionMatrix, amino_acids

from vdjmatch import db
from gen_vdjam import hamming1_events, position_background, estimate_matrix, AA  # noqa: E402


def subst_from_scores(scores: dict[tuple[str, str], float], scale: int = 100) -> SubstitutionMatrix:
    """Build a seqtree SubstitutionMatrix from learned off-diagonal log-odds (dominant diagonal)."""
    order = amino_acids(); n = len(order)
    off = [v for (a, b), v in scores.items() if a != b]  # diagonal log-odds are meaningless
    hi = round(scale * (max(off) + 1.0)); lo = round(scale * (min(off) - 1.0))
    grid = [[lo for _ in range(n)] for _ in range(n)]
    for i in range(n):
        grid[i][i] = hi                                   # dominant, uniform self-similarity
    idx = {a: i for i, a in enumerate(order)}
    for (a, b), s in scores.items():
        if a != b:
            grid[idx[a]][idx[b]] = grid[idx[b]][idx[a]] = round(scale * s)
    return SubstitutionMatrix.from_similarity(grid)


def pr_auc(labels_scores: list[tuple[int, float]]) -> float:
    """Micro PR-AUC; score higher = more confident. labels in {0,1}."""
    s = sorted(labels_scores, key=lambda x: -x[1])
    P = sum(l for l, _ in s)
    if P == 0:
        return float("nan")
    tp = fp = 0
    prev_r, prev_p, area = 0.0, 1.0, 0.0
    for lab, _ in s:
        if lab:
            tp += 1
        else:
            fp += 1
        r, p = tp / P, tp / (tp + fp)
        area += (r - prev_r) * (p + prev_p) / 2
        prev_r, prev_p = r, p
    return area


def retrieval_prauc(index: Index, ep_of_ref: list[str], queries: list[str], true_epi: str,
                    matrix, scope: str) -> float:
    s, _, _, t = [int(x) for x in (scope.split(",") + ["0", "0", "0"])][:4] or (2, 0, 0, 2)
    # bound the search purely by substitution count; max_penalty huge so the (scaled) matrix only
    # *ranks* hits and never prunes them — all three matrices then see the same candidate set.
    params = SearchParams(max_subs=s, max_total_edits=max(s, t), engine="seqtm",
                          matrix=matrix or "", max_penalty=10**9)
    res = index.search_batch(queries, params, 0)
    pairs = []
    for qi, hl in enumerate(res):
        for h in hl:
            if h.score == 0 and ep_of_ref[h.ref_id] == true_epi:
                continue  # drop exact self / identical same-epitope clone
            pairs.append((1 if ep_of_ref[h.ref_id] == true_epi else 0, -float(h.score)))
    return pr_auc(pairs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default="/Users/mikesh/vcs/manuscripts/2026-vdjmatch/test_data/sample3_vdjdb.txt")
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--chain", default="TRB")
    ap.add_argument("--min-epi", type=int, default=50)
    ap.add_argument("--top", type=int, default=8, help="held-out epitopes (largest)")
    ap.add_argument("--scope", default="2,0,0,2")
    args = ap.parse_args()

    vdj = db.load(args.pmhc, asset="full", species=args.species).filter(pl.col("gene") == args.chain)
    uc = vdj.select("cdr3", "epitope").unique()
    all_cdr3 = uc["cdr3"].unique().to_list()
    bg = position_background(all_cdr3)

    # one retrieval index over the whole chain; ref_id -> epitope (first label)
    refs = uc.group_by("cdr3").agg(pl.col("epitope").first())
    ref_list = refs["cdr3"].to_list()
    ep_of_ref = refs["epitope"].to_list()
    index = Index.build(ref_list, "aa")
    blosum = SubstitutionMatrix.blosum62()
    unit = SubstitutionMatrix.unit(len(amino_acids()))

    sizes = (uc.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
               .filter(pl.col("n") >= args.min_epi).sort("n", descending=True))
    held = sizes["epitope"].to_list()[:args.top]
    print(f"chain={args.chain}; held-out epitopes={len(held)}; scope={args.scope}")
    print(f"{'epitope':14}{'n':>5}{'PR_VDJAM':>10}{'PR_BLOSUM':>11}{'PR_unit':>9}")
    dV, dU = [], []
    for epi in held:
        ev_train = []
        for e in sizes["epitope"].to_list():
            if e == epi:
                continue
            cds = uc.filter(pl.col("epitope") == e)["cdr3"].unique().to_list()
            ev_train.extend([x for x in hamming1_events(cds, args.min_epi)])
        # NDN-region events only (central positions), per the gen_vdjam finding
        ev_ndn = [x for x in ev_train if 4 <= x[2] < x[3] - 6]
        vdjam = subst_from_scores(estimate_matrix(ev_ndn, bg))
        q = uc.filter(pl.col("epitope") == epi)["cdr3"].unique().to_list()
        a = retrieval_prauc(index, ep_of_ref, q, epi, vdjam, args.scope)
        b = retrieval_prauc(index, ep_of_ref, q, epi, blosum, args.scope)
        u = retrieval_prauc(index, ep_of_ref, q, epi, unit, args.scope)
        dV.append(a - b); dU.append(a - u)
        print(f"{epi:14}{len(q):>5}{a:>10.3f}{b:>11.3f}{u:>9.3f}")
    import statistics as st
    print(f"\nmean dPR-AUC  VDJAM-BLOSUM = {st.mean(dV):+.3f}   VDJAM-unit = {st.mean(dU):+.3f}")
    print(f"epitopes where VDJAM>BLOSUM: {sum(1 for d in dV if d>0)}/{len(dV)}")


if __name__ == "__main__":
    main()
