#!/usr/bin/env python3
"""Which germline FR/CDR region (if any) carries the V-gene co-specificity prior?

The V match is a strong prior, but coarse CDR1+CDR2 similarity does not recover it (vgene_scan.py).
Two hypotheses remain: (i) the decisive germline signal is in the framework (FR), not CDR1/CDR2; or
(ii) it IS CDR1/CDR2 but -- like CDR3 -- needs near-exact identity, which a coarse similarity washed
out. We decide this by exact-match decomposition: among cross-V CDR3-neighbour pairs (human TRB MHC-I),
stratify by WHICH germline regions are identical between the two V genes, and ask whether any region or
combination, when identical, lifts cross-V co-specificity toward the same-V level.

  - per region R in {FR1,CDR1,FR2,CDR2,FR3}: P(co-specific | R identical) and lift vs cross-V baseline
  - combinations: CDR1&CDR2 identical; any-CDR identical; all-FR identical
  - and by CDR1+CDR2 edit distance (0,1-2,3-5,6+): does "0-1 mismatch" (CDR3-like tolerance) suffice?

    python bench/vregion_decompose.py --chain TRB --class MHCI --subs 2
"""
from __future__ import annotations

import argparse
import os

import polars as pl
from seqtree import Index, SearchParams

from vdjmatch import db
from vdjmatch.match import regions, vgene

REG = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3"]


def lev(a: str, b: str) -> int:
    if a == b:
        return 0
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default=os.environ.get("VDJDB_SAMPLE", "test_data/sample3_vdjdb.txt"),
                    help="VDJdb export TSV (default $VDJDB_SAMPLE or test_data/sample3_vdjdb.txt)")
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--chain", default="TRB")
    ap.add_argument("--mhc-class", default="MHCI")
    ap.add_argument("--min-epi", type=int, default=50)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--max-queries", type=int, default=600)
    ap.add_argument("--subs", type=int, default=2)
    args = ap.parse_args()

    vreg = vgene.load_v_regions(args.chain)
    vdj = (db.load(args.pmhc, asset="full", species=args.species)
             .filter((pl.col("gene") == args.chain) & (pl.col("mhc_class") == args.mhc_class)))
    uc = (vdj.select("cdr3", "v", "epitope").unique()
             .filter(pl.col("cdr3").str.contains("^[ACDEFGHIKLMNPQRSTVWY]+$")))
    refs = uc.group_by("cdr3").agg(pl.col("v").first(), pl.col("epitope").first())
    ref_cdr3, ref_v, ref_epi = refs["cdr3"].to_list(), refs["v"].to_list(), refs["epitope"].to_list()
    index = Index.build(ref_cdr3, "aa")
    params = SearchParams(max_subs=args.subs, max_total_edits=args.subs, engine="seqtm")
    sizes = (uc.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
               .filter(pl.col("n") >= args.min_epi).sort("n", descending=True))
    held = sizes["epitope"].to_list()[:args.top]

    same = [0, 0]
    cross = [0, 0]
    by_reg = {r: [0, 0] for r in REG}                  # region identical across the two V genes
    combos = {"CDR1&CDR2": [0, 0], "anyCDR": [0, 0], "allFR": [0, 0]}
    dist = {"0": [0, 0], "1-2": [0, 0], "3-5": [0, 0], "6+": [0, 0]}  # CDR1+CDR2 edit distance
    print(f"species={args.species} chain={args.chain} class={args.mhc_class} subs={args.subs} "
          f"epitopes={len(held)}")
    for epi in held:
        q = uc.filter(pl.col("epitope") == epi).unique("cdr3").head(args.max_queries)
        qs, qv = q["cdr3"].to_list(), q["v"].to_list()
        cand = index.search_batch(qs, params, 0)
        for i, hl in enumerate(cand):
            qfam = regions.gene_family(qv[i])
            ra = vreg.get(qfam)
            for h in hl:
                ri = h.ref_id
                r = ref_cdr3[ri]
                if len(r) != len(qs[i]) or r == qs[i]:
                    continue
                lab = 1 if ref_epi[ri] == epi else 0
                rfam = regions.gene_family(ref_v[ri])
                if rfam == qfam:
                    same[0] += lab; same[1] += 1
                    continue
                cross[0] += lab; cross[1] += 1
                rb = vreg.get(rfam)
                if ra is None or rb is None:
                    continue
                ident = {x: ra[x] == rb[x] for x in REG}
                for x in REG:
                    if ident[x]:
                        by_reg[x][0] += lab; by_reg[x][1] += 1
                if ident["cdr1"] and ident["cdr2"]:
                    combos["CDR1&CDR2"][0] += lab; combos["CDR1&CDR2"][1] += 1
                if ident["cdr1"] or ident["cdr2"]:
                    combos["anyCDR"][0] += lab; combos["anyCDR"][1] += 1
                if ident["fwr1"] and ident["fwr2"] and ident["fwr3"]:
                    combos["allFR"][0] += lab; combos["allFR"][1] += 1
                d = lev(ra["cdr1"], rb["cdr1"]) + lev(ra["cdr2"], rb["cdr2"])
                k = "0" if d == 0 else "1-2" if d <= 2 else "3-5" if d <= 5 else "6+"
                dist[k][0] += lab; dist[k][1] += 1

    def pct(pt):
        return 100 * pt[0] / pt[1] if pt[1] else float("nan")
    sp, cr = pct(same), pct(cross)
    print(f"\nsame-V: {sp:.1f}%  (n={same[1]})   cross-V baseline: {cr:.1f}%  (n={cross[1]})\n")
    print(f"{'cross-V stratum':28}{'co-spec%':>10}{'n':>9}{'lift vs cross':>15}")
    for x in REG:
        print(f"  {x+' identical':26}{pct(by_reg[x]):>10.1f}{by_reg[x][1]:>9}"
              f"{pct(by_reg[x])/cr:>14.2f}x")
    for k, pt in combos.items():
        print(f"  {k+' identical':26}{pct(pt):>10.1f}{pt[1]:>9}{pct(pt)/cr:>14.2f}x")
    print(f"\ncross-V co-specificity by CDR1+CDR2 edit distance (0 = identical loops):")
    for k in ("0", "1-2", "3-5", "6+"):
        print(f"  edit {k:5}{pct(dist[k]):>10.1f}%   (n={dist[k][1]})")
    print(f"\ndecisive? a region/combination is decisive only if its co-spec approaches same-V "
          f"({sp:.0f}%). cross-V floor is {cr:.0f}%.")


if __name__ == "__main__":
    main()
