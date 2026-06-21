#!/usr/bin/env python3
"""V-gene-stratified scoring: does down-weighting CDR3 borders pay off more for CROSS-V comparisons?

The V-gene lever (the hypothesis): same-V same-length CDR3 pairs share the germline V prefix, so all
their mismatches already sit in the NDN core -- position weighting is then a near no-op (this is why
region weighting barely moved overall retrieval). CROSS-V pairs additionally differ at the
germline-fixed V-flank, which is gene-identity noise rather than specificity signal. So down-weighting
borders (the central-significance PSSM) or germline flanks (the retention weight) should improve
specificity resolution *specifically* among cross-V pairs.

We re-use the leave-one-out candidate machinery of loo_vdjam.py, then split every scored (query, hit)
pair into same-V vs cross-V (by V gene family) and pool across held-out epitopes into one micro PR-AUC
per stratum, under flat BLOSUM62, possig-weighted BLOSUM62, and germline-retention-weighted BLOSUM62.

    python bench/vgene_strat.py --chain TRB --subs 4
"""
from __future__ import annotations

import argparse
import os

import polars as pl
from seqtree import Index, SearchParams

from vdjmatch import db
from vdjmatch.match import regions, vgene
from loo_vdjam import named_dissim, pr_auc  # noqa: E402
from gen_vdjam import AA  # noqa: E402  (kept for parity / future matrix arms)

VSIM_BINS = [0.0, 0.3, 0.5, 0.7, 0.9, 1.0]  # germline CDR1+CDR2 similarity buckets for cross-V pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default=os.environ.get("VDJDB_SAMPLE", "test_data/sample3_vdjdb.txt"),
                    help="VDJdb export TSV (default $VDJDB_SAMPLE or test_data/sample3_vdjdb.txt)")
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--chain", default="TRB")
    ap.add_argument("--min-epi", type=int, default=50)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--max-queries", type=int, default=600)
    ap.add_argument("--subs", type=int, default=4)
    args = ap.parse_args()

    vdj = db.load(args.pmhc, asset="full", species=args.species).filter(pl.col("gene") == args.chain)
    uc = vdj.select("cdr3", "v", "j", "epitope").unique()
    ret = regions.load_retention()
    blo = named_dissim("BLOSUM62")

    refs = uc.group_by("cdr3").agg(pl.col("v").first(), pl.col("j").first(), pl.col("epitope").first())
    ref_cdr3 = refs["cdr3"].to_list()
    ref_v = refs["v"].to_list()
    ref_epi = refs["epitope"].to_list()
    index = Index.build(ref_cdr3, "aa")
    cand_params = SearchParams(max_subs=args.subs, max_total_edits=args.subs, engine="seqtm")

    sizes = (uc.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
               .filter(pl.col("n") >= args.min_epi).sort("n", descending=True))
    held = sizes["epitope"].to_list()[:args.top]
    print(f"species={args.species} chain={args.chain}; held-out epitopes={len(held)}; subs={args.subs}")

    # pooled (label, -score) pairs per V-stratum per scoring; plus cross-V co-specificity by V-similarity
    pools = {s: {m: [] for m in ("flat", "possig", "region")} for s in ("same", "cross")}
    vbin_pos = [0] * (len(VSIM_BINS) - 1)
    vbin_tot = [0] * (len(VSIM_BINS) - 1)
    for epi in held:
        q = uc.filter(pl.col("epitope") == epi).unique("cdr3").head(args.max_queries)
        qs, qv, qj = q["cdr3"].to_list(), q["v"].to_list(), q["j"].to_list()
        cand = index.search_batch(qs, cand_params, 0)
        for i, hl in enumerate(cand):
            L = len(qs[i])
            wsig = regions.significance_weights(L)
            wreg = regions.position_weights(L, qv[i], qj[i], args.chain, ret)
            qfam = regions.gene_family(qv[i])
            for h in hl:
                ri = h.ref_id
                r = ref_cdr3[ri]
                if len(r) != L or r == qs[i]:
                    continue
                lab = 1 if ref_epi[ri] == epi else 0
                same = regions.gene_family(ref_v[ri]) == qfam
                strat = "same" if same else "cross"
                flat = possig = region = 0.0
                for p, (a, b) in enumerate(zip(qs[i], r)):
                    if a != b:
                        d = blo.get((a, b), 0.0)
                        flat += d
                        possig += wsig[p] * d
                        region += wreg[p] * d
                pools[strat]["flat"].append((lab, -flat))
                pools[strat]["possig"].append((lab, -possig))
                pools[strat]["region"].append((lab, -region))
                if not same:  # does germline V-similarity recover the same-V co-specificity prior?
                    s = vgene.vsim(qv[i], ref_v[ri])
                    b = min(len(VSIM_BINS) - 2, max(0, sum(s >= e for e in VSIM_BINS[1:]) ))
                    vbin_tot[b] += 1
                    vbin_pos[b] += lab

    print(f"{'stratum':8}{'pairs':>9}{'pos%':>7}{'flat':>9}{'possig':>9}{'region':>9}"
          f"{'d(possig)':>11}{'d(region)':>11}")
    for s in ("same", "cross"):
        n = len(pools[s]["flat"])
        pos = sum(l for l, _ in pools[s]["flat"]) / n if n else 0.0
        a = {m: pr_auc(pools[s][m]) for m in ("flat", "possig", "region")}
        print(f"{s:8}{n:>9}{pos*100:>6.1f}%{a['flat']:>9.3f}{a['possig']:>9.3f}{a['region']:>9.3f}"
              f"{a['possig']-a['flat']:>+11.3f}{a['region']-a['flat']:>+11.3f}")
    print("\ncross-V co-specificity by germline CDR1+CDR2 similarity (does soft-V grouping recover the "
          "same-V prior?):")
    print(f"{'vsim bin':12}{'pairs':>9}{'co-specific%':>14}")
    for b in range(len(vbin_tot)):
        lo, hi = VSIM_BINS[b], VSIM_BINS[b + 1]
        pct = 100 * vbin_pos[b] / vbin_tot[b] if vbin_tot[b] else 0.0
        print(f"[{lo:.1f},{hi:.1f})   {vbin_tot[b]:>9}{pct:>13.1f}%")


if __name__ == "__main__":
    main()
