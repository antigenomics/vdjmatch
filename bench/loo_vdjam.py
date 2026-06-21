#!/usr/bin/env python3
"""Leave-one-out-by-epitope evaluation of VDJAM and region-aware (NDN-weighted) scoring.

VDJAM is derived from same-antigen substitution pairs, so testing on the same epitopes is circular.
For each held-out epitope E* we re-derive the NDN substitution matrix from **all other** epitopes,
then run VDJdb-vs-VDJdb retrieval on E*: a single unit-cost search (scope 2, no indels) fixes the
candidate set, and each (query, hit) pair is scored four ways and ranked. Micro PR-AUC for
"hit shares E*" measures each scoring's specificity resolution:

  unit          : Hamming (substitution count)
  BLOSUM62      : sum of per-substitution dissimilarity
  VDJAM (flat)  : sum of NDN-VDJAM dissimilarity, all positions equal
  VDJAM region  : NDN-VDJAM dissimilarity weighted by (1 - germline retention) per position
                  -> germline V/J flanks discounted, NDN core emphasized

    python bench/loo_vdjam.py --chain TRB --top 8 --max-queries 600
"""
from __future__ import annotations

import argparse
import statistics as st

import polars as pl
from Bio.Align import substitution_matrices
from seqtree import Index, SearchParams

from vdjmatch import db
from vdjmatch.match import regions
from gen_vdjam import hamming1_events, position_background, estimate_matrix, AA  # noqa: E402


def dissim_from_scores(scores: dict[tuple[str, str], float]) -> dict[tuple[str, str], float]:
    """similarity log-odds -> non-negative dissimilarity (most-similar pair = 0)."""
    off = {k: v for k, v in scores.items() if k[0] != k[1]}
    hi = max(off.values())
    return {k: hi - v for k, v in off.items()}


def blosum_dissim() -> dict[tuple[str, str], float]:
    b = substitution_matrices.load("BLOSUM62")
    vals = {(a, c): float(b[a, c]) for a in AA for c in AA if a != c}
    hi = max(vals.values())
    return {k: hi - v for k, v in vals.items()}


def pr_auc(pairs: list[tuple[int, float]]) -> float:
    s = sorted(pairs, key=lambda x: -x[1])
    P = sum(l for l, _ in s)
    if P == 0:
        return float("nan")
    tp = fp = 0
    pr, pp, area = 0.0, 1.0, 0.0
    for lab, _ in s:
        tp += lab; fp += 1 - lab
        r, p = tp / P, tp / (tp + fp)
        area += (r - pr) * (p + pp) / 2
        pr, pp = r, p
    return area


def score_pairs(qseqs, qv, qj, cand, ref_cdr3, ref_epi, true_epi, chain, ret, dis, weighted):
    """(label, -score) pairs over all candidate hits; score = sum of (weighted) dissimilarity."""
    out = []
    for i, hl in enumerate(cand):
        w = regions.position_weights(len(qseqs[i]), qv[i], qj[i], chain, ret) if weighted else None
        for ri in hl:
            r = ref_cdr3[ri]
            if len(r) != len(qseqs[i]):
                continue  # subs-only candidates are same length
            if r == qseqs[i]:
                continue  # exact self / identical clone -- not an informative retrieval
            sc = 0.0
            for p, (a, b) in enumerate(zip(qseqs[i], r)):
                if a != b:
                    sc += (w[p] if w else 1.0) * dis.get((a, b), 0.0)
            out.append((1 if ref_epi[ri] == true_epi else 0, -sc))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default="/Users/mikesh/vcs/manuscripts/2026-vdjmatch/test_data/sample3_vdjdb.txt")
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--chain", default="TRB")
    ap.add_argument("--min-epi", type=int, default=50)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--max-queries", type=int, default=600)
    ap.add_argument("--subs", type=int, default=2, help="candidate search substitution budget")
    args = ap.parse_args()

    vdj = db.load(args.pmhc, asset="full", species=args.species).filter(pl.col("gene") == args.chain)
    uc = vdj.select("cdr3", "v", "j", "epitope").unique()
    bg = position_background(uc["cdr3"].unique().to_list())
    ret = regions.load_retention()
    blo = blosum_dissim()
    unit = {(a, c): 1.0 for a in AA for c in AA if a != c}

    refs = uc.group_by("cdr3").agg(pl.col("epitope").first())
    ref_cdr3 = refs["cdr3"].to_list()
    ref_epi = refs["epitope"].to_list()
    index = Index.build(ref_cdr3, "aa")
    cand_params = SearchParams(max_subs=args.subs, max_total_edits=args.subs, engine="seqtm")

    sizes = (uc.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
               .filter(pl.col("n") >= args.min_epi).sort("n", descending=True))
    held = sizes["epitope"].to_list()[:args.top]
    print(f"chain={args.chain}; held-out epitopes={len(held)}")
    print(f"{'epitope':13}{'n':>5}{'unit':>8}{'BLOSUM':>8}{'VDJAM':>8}{'VDJAM_reg':>10}")
    acc = {k: [] for k in ("unit", "blosum", "vdjam", "region")}
    for epi in held:
        ev = []
        for e in sizes["epitope"].to_list():
            if e == epi:
                continue
            cds = uc.filter(pl.col("epitope") == e)["cdr3"].unique().to_list()
            ev += [x for x in hamming1_events(cds, args.min_epi) if 4 <= x[2] < x[3] - 6]  # NDN
        vdis = dissim_from_scores(estimate_matrix(ev, bg))

        q = uc.filter(pl.col("epitope") == epi).unique("cdr3").head(args.max_queries)
        qs, qv, qj = q["cdr3"].to_list(), q["v"].to_list(), q["j"].to_list()
        cand = [[h.ref_id for h in hl] for hl in index.search_batch(qs, cand_params, 0)]
        row = {
            "unit": pr_auc(score_pairs(qs, qv, qj, cand, ref_cdr3, ref_epi, epi, args.chain, ret, unit, False)),
            "blosum": pr_auc(score_pairs(qs, qv, qj, cand, ref_cdr3, ref_epi, epi, args.chain, ret, blo, False)),
            "vdjam": pr_auc(score_pairs(qs, qv, qj, cand, ref_cdr3, ref_epi, epi, args.chain, ret, vdis, False)),
            "region": pr_auc(score_pairs(qs, qv, qj, cand, ref_cdr3, ref_epi, epi, args.chain, ret, vdis, True)),
        }
        for k in acc:
            acc[k].append(row[k])
        print(f"{epi:13}{len(qs):>5}{row['unit']:>8.3f}{row['blosum']:>8.3f}"
              f"{row['vdjam']:>8.3f}{row['region']:>10.3f}")
    print(f"\nmean PR-AUC:  unit {st.mean(acc['unit']):.3f}  BLOSUM {st.mean(acc['blosum']):.3f}  "
          f"VDJAM {st.mean(acc['vdjam']):.3f}  VDJAM-region {st.mean(acc['region']):.3f}")
    print(f"region vs flat VDJAM: {st.mean(acc['region']) - st.mean(acc['vdjam']):+.3f}  | "
          f"region vs BLOSUM: {st.mean(acc['region']) - st.mean(acc['blosum']):+.3f}  | "
          f"region>BLOSUM in {sum(1 for a,b in zip(acc['region'],acc['blosum']) if a>b)}/{len(held)}")


if __name__ == "__main__":
    main()
