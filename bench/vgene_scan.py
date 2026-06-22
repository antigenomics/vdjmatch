#!/usr/bin/env python3
"""Comprehensive V-gene / CDR1-2 scan across species x chain x MHC class.

Two questions, swept over every {HomoSapiens,MusMusculus} x {TRA,TRB} x {MHCI,MHCII} cell of VDJdb:

  1. Is a V match a co-specificity prior? (same-V vs cross-V neighbours sharing the held-out epitope)
  2. Does germline CDR1+CDR2 similarity *interpolate* that prior? -- i.e. is there a CDR1-2 effect
     that would justify soft V-clustering (fuzzy V matching)?

For each cell we run leave-one-out-style VDJdb-vs-VDJdb retrieval (unit-cost, same-length subs-only
candidate pairs), pool scored neighbour pairs, and report same-V/cross-V co-specificity, their ratio,
and the point-biserial correlation r(vsim, co-specific) among cross-V pairs. A CDR1-2 effect would make
cross-V co-specificity climb toward the same-V level as V-similarity -> 1 (r > 0, high-sim bin elevated).
A flat gradient (r ~ 0) in every cell refutes focusing on CDR1/CDR2 for fuzzy V matching.

    python bench/vgene_scan.py --subs 2 --top 6 --min-epi 30
"""
from __future__ import annotations

import argparse
import os
from math import sqrt

import polars as pl
from seqtree import Index, SearchParams

import _bench
from vdjmatch import db
from vdjmatch.match import regions, vgene

VSIM_BINS = [0.0, 0.3, 0.5, 0.7, 0.9, 1.0]


def scan_cell(uc: pl.DataFrame, chain: str, top: int, min_epi: int, max_q: int, subs: int):
    """One species/chain/class cell. Returns a summary dict or None if too little data."""
    uc = uc.filter(pl.col("cdr3").str.contains("^[ACDEFGHIKLMNPQRSTVWY]+$"))  # drop non-std residues
    refs = uc.group_by("cdr3").agg(pl.col("v").first(), pl.col("epitope").first())
    ref_cdr3 = refs["cdr3"].to_list()
    ref_v = refs["v"].to_list()
    ref_epi = refs["epitope"].to_list()
    sizes = (uc.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
               .filter(pl.col("n") >= min_epi).sort(["n", "epitope"], descending=[True, False]))
    held = sizes["epitope"].to_list()[:top]
    if len(held) < 2:
        return None
    index = Index.build(ref_cdr3, "aa")
    params = SearchParams(max_subs=subs, max_total_edits=subs, engine="seqtm")

    same = [0, 0]   # [pos, tot]
    cross = [0, 0]
    vbin = [[0, 0] for _ in VSIM_BINS[:-1]]
    # point-biserial accumulators over cross-V pairs: corr(vsim, label)
    n = sx = sx2 = sy = sxy = 0.0
    per_epi_same = []  # same-V co-specificity per held epitope (is the V effect universal?)
    nq = 0
    for epi in held:
        q = uc.filter(pl.col("epitope") == epi).unique("cdr3").head(max_q)
        qs, qv = q["cdr3"].to_list(), q["v"].to_list()
        nq += len(qs)
        cand = index.search_batch(qs, params, 0)
        e_same = [0, 0]
        for i, hl in enumerate(cand):
            L = len(qs[i])
            qfam = regions.gene_family(qv[i])
            for h in hl:
                ri = h.ref_id
                r = ref_cdr3[ri]
                if len(r) != L or r == qs[i]:
                    continue
                lab = 1 if ref_epi[ri] == epi else 0
                if regions.gene_family(ref_v[ri]) == qfam:
                    same[0] += lab; same[1] += 1
                    e_same[0] += lab; e_same[1] += 1
                else:
                    cross[0] += lab; cross[1] += 1
                    s = vgene.vsim(qv[i], ref_v[ri])
                    b = min(len(VSIM_BINS) - 2, sum(s >= e for e in VSIM_BINS[1:]))
                    vbin[b][0] += lab; vbin[b][1] += 1
                    n += 1; sx += s; sx2 += s * s; sy += lab; sxy += s * lab
        if e_same[1]:
            per_epi_same.append(e_same[0] / e_same[1])
    if cross[1] == 0 or same[1] == 0:
        return None
    # point-biserial r
    num = n * sxy - sx * sy
    den = sqrt(max(1e-9, (n * sx2 - sx * sx) * (n * sy - sy * sy)))
    r = num / den if den else 0.0
    return dict(epi=len(held), nq=nq, same=same, cross=cross, vbin=vbin, r=r,
                per_epi=per_epi_same)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default=None,
                    help="VDJdb export TSV (default: $VDJDB_SAMPLE or the HF-pinned release)")
    ap.add_argument("--top", type=int, default=6)
    ap.add_argument("--min-epi", type=int, default=30)
    ap.add_argument("--max-queries", type=int, default=400)
    ap.add_argument("--subs", type=int, default=2)
    ap.add_argument("--out", default="appendix/figures/vgene_scan.dat")
    args = ap.parse_args()

    org = {"HomoSapiens": "human", "MusMusculus": "mouse"}
    short = {"HomoSapiens": "h", "MusMusculus": "m"}
    rows = []
    bars = []  # (short_label, same%, cross%) per cell for the summary figure
    print(f"subs={args.subs} top={args.top} min-epi={args.min_epi}\n")
    print(f"{'cell':22}{'epi':>4}{'sameV n':>9}{'same%':>7}{'crossV n':>10}{'cross%':>8}"
          f"{'ratio':>7}{'r(vsim,co)':>11}{'lo(n)':>13}{'hi(n)':>13}")
    src = _bench.source(args.pmhc)
    for sp in ("HomoSapiens", "MusMusculus"):
        vdj = db.load(src, species=sp)
        for chain in ("TRA", "TRB"):
            for cls in ("MHCI", "MHCII"):
                uc = _bench.long_list(vdj.filter((pl.col("gene") == chain)
                                                  & (pl.col("mhc_class") == cls)),
                                      cap=3000, min_n=args.min_epi)  # composition-controlled
                res = scan_cell(uc, chain, args.top, args.min_epi, args.max_queries, args.subs)
                cell = f"{org[sp]} {chain} {cls}"
                if res is None:
                    print(f"{cell:22}  -- too little data --")
                    continue
                sp_pct = 100 * res["same"][0] / res["same"][1]
                cr_pct = 100 * res["cross"][0] / res["cross"][1]
                ratio = sp_pct / cr_pct if cr_pct else float("nan")
                lo = res["vbin"][0]                                  # vsim [0,0.3)
                hi = [res["vbin"][3][0] + res["vbin"][4][0],         # vsim >=0.7 (pooled, robust)
                      res["vbin"][3][1] + res["vbin"][4][1]]
                lo_pct = 100 * lo[0] / lo[1] if lo[1] else 0.0
                hi_pct = 100 * hi[0] / hi[1] if hi[1] else 0.0
                print(f"{cell:22}{res['epi']:>4}{res['same'][1]:>9}{sp_pct:>6.1f}%"
                      f"{res['cross'][1]:>10}{cr_pct:>7.1f}%{ratio:>7.1f}{res['r']:>11.3f}"
                      f"{lo_pct:>8.1f}%({lo[1]}){hi_pct:>7.1f}%({hi[1]})")
                bars.append((f"{short[sp]}{chain[-1]}-{'I' if cls=='MHCI' else 'II'}",
                             sp_pct, cr_pct))
                for b, (p, t) in enumerate(res["vbin"]):
                    rows.append(f"{org[sp]}\t{chain}\t{cls}\t{VSIM_BINS[b]:.1f}\t{VSIM_BINS[b+1]:.1f}"
                                f"\t{t}\t{(p/t if t else 0):.4f}\t{sp_pct/100:.4f}\t{cr_pct/100:.4f}")
    from pathlib import Path
    import subprocess
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("org\tchain\tclass\tvsim_lo\tvsim_hi\tn\tco_spec\tsameV\tcrossV\n"
                   + "\n".join(rows) + "\n")
    print(f"\nlo% = cross-V co-specificity at CDR1+CDR2 similarity [0,0.3); hi% at >=0.7 (pooled).")
    print("CDR1-2 effect would need hi% -> same%, r > 0. Wrote per-cell bins to", args.out)

    # summary figure: same-V vs cross-V co-specificity per cell (the V prior is universal; CDR1-2 not)
    bdat = out.parent / "vgene_scan_bars.dat"
    bdat.write_text("cell same cross\n" + "".join(f"{c} {s:.1f} {x:.1f}\n" for c, s, x in bars))
    gp = out.parent / "_vg.gp"
    gp.write_text(f"""
set terminal svg size 680,360 font 'Helvetica,12' background rgb 'white'
set output 'vgene_scan.svg'
set style data histogram; set style histogram cluster gap 1
set style fill solid 0.85 border -1; set boxwidth 0.9
set ylabel 'co-specificity of CDR3 neighbours (%)'; set yrange [0:100]
set grid ytics lc rgb '#e5e7eb'; set xtics rotate by -25; set key top right
set title 'A V match is a universal co-specificity prior across species/chain/class'
plot '{bdat.name}' using 2:xtic(1) title 'same V' lc rgb '#2563eb', \\
     '' using 3 title 'cross V' lc rgb '#94a3b8'
""")
    subprocess.run(["gnuplot", gp.name], cwd=out.parent, check=True)
    gp.unlink()
    print("wrote", out.parent / "vgene_scan.svg")


if __name__ == "__main__":
    main()
