#!/usr/bin/env python3
"""V pseudosequence test (b): which germline CDR1/CDR2 *positions*, shared across different V genes,
predict CDR3 co-specificity?

The region-level test (vregion_decompose.py) shows cross-V co-specificity rises as germline CDR1+CDR2
approach identity. Here we go to single-residue resolution: among cross-V CDR3-neighbour pairs whose
V genes both carry the modal-length contacting loops (human TRB: CDR1 len 5, CDR2 len 6), we ask, for
each of the 11 loop positions, P(co-specific | the two V genes share that residue) and its lift over
the cross-V baseline. Positions with high lift form a candidate V "pseudosequence" (the MHC-pseudoseq
analogue); a flat profile would mean position identity is irrelevant. Full IMGT-numbered coverage
(FR + the DE-loop) is the arda milestone; this is the length-constant CDR1/CDR2 first cut.

    python bench/vpseudo.py --chain TRB --mhc-class MHCI --subs 2
"""
from __future__ import annotations

import argparse

import polars as pl
from seqtree import Index, SearchParams

import _bench
from vdjmatch import db
from vdjmatch.match import regions, vgene

L1, L2 = 5, 6  # modal human-TRB CDR1, CDR2 lengths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default=None)
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--chain", default="TRB")
    ap.add_argument("--mhc-class", default="MHCI")
    ap.add_argument("--min-epi", type=int, default=30)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--max-queries", type=int, default=600)
    ap.add_argument("--subs", type=int, default=2)
    args = ap.parse_args()

    vreg = vgene.load_v_regions(args.chain)
    loops = {v: d["cdr1"] + d["cdr2"] for v, d in vreg.items()
             if len(d["cdr1"]) == L1 and len(d["cdr2"]) == L2}
    vdj = (db.load(_bench.source(args.pmhc), species=args.species)
             .filter((pl.col("gene") == args.chain) & (pl.col("mhc_class") == args.mhc_class)))
    uc = _bench.long_list(vdj, cap=3000, min_n=args.min_epi)
    refs = uc.group_by("cdr3").agg(pl.col("v").first(), pl.col("epitope").first())
    ref_cdr3, ref_v, ref_epi = refs["cdr3"].to_list(), refs["v"].to_list(), refs["epitope"].to_list()
    index = Index.build(ref_cdr3, "aa")
    params = SearchParams(max_subs=args.subs, max_total_edits=args.subs, engine="seqtm")
    sizes = (uc.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
               .filter(pl.col("n") >= args.min_epi).sort(["n", "epitope"], descending=[True, False]))
    held = sizes["epitope"].to_list()[:args.top]

    base = [0, 0]
    match = [[0, 0] for _ in range(L1 + L2)]   # per position: [co-spec, total] when residue shared
    for epi in held:
        q = uc.filter(pl.col("epitope") == epi).unique("cdr3").sort("cdr3").head(args.max_queries)
        qs, qv = q["cdr3"].to_list(), q["v"].to_list()
        cand = index.search_batch(qs, params, 0)
        for i, hl in enumerate(cand):
            qfam = regions.gene_family(qv[i])
            ql = loops.get(qfam)
            if ql is None:
                continue
            for h in hl:
                r = ref_cdr3[h.ref_id]
                if len(r) != len(qs[i]) or r == qs[i]:
                    continue                                  # same length, never exact self
                rfam = regions.gene_family(ref_v[h.ref_id])
                if rfam == qfam:
                    continue                                  # cross-V only
                rl = loops.get(rfam)
                if rl is None:
                    continue
                lab = 1 if ref_epi[h.ref_id] == epi else 0
                base[0] += lab; base[1] += 1
                for p in range(L1 + L2):
                    if ql[p] == rl[p]:
                        match[p][0] += lab; match[p][1] += 1

    b = base[0] / base[1] if base[1] else float("nan")
    print(f"species={args.species} {args.chain} {args.mhc_class} subs={args.subs}; "
          f"cross-V pairs (modal loops) = {base[1]}; baseline co-spec = {b*100:.1f}%\n")
    print(f"{'loop pos':12}{'residue-match co-spec':>22}{'n_match':>9}{'lift':>8}")
    for p in range(L1 + L2):
        lab = "CDR1." + str(p + 1) if p < L1 else "CDR2." + str(p - L1 + 1)
        pm = match[p][0] / match[p][1] if match[p][1] else float("nan")
        print(f"{lab:12}{pm*100:>21.1f}%{match[p][1]:>9}{(pm/b if b else 0):>7.2f}x")
    print("\nlift = P(co-spec | residue shared at this position) / cross-V baseline. A V pseudosequence "
          "would show a few positions with lift >> 1; a flat ~1 profile means position identity is "
          "uninformative beyond whole-gene identity.")


if __name__ == "__main__":
    main()
