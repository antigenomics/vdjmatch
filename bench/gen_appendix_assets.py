#!/usr/bin/env python3
"""Generate figures + example alignments for the scoring appendix.

Emits into appendix/:
  figures/region_corr.svg   — per-region (V/NDN/J) correlation of the learned matrix with BLOSUM62
  figures/loo_prauc.svg     — leave-one-out retrieval PR-AUC (VDJAM / BLOSUM62 / unit) per epitope
  aln/<epi>_msa.fasta       — same-length CDR3 blocks per epitope for texshade (multiple alignment)
  aln/pairwise.tex          — seqtree pairwise alignments of different-length same-epitope CDR3s
"""
from __future__ import annotations

import os
import subprocess
from collections import Counter
from pathlib import Path

import polars as pl
from seqtree import Index, SearchParams

import _bench
from vdjmatch import db
from vdjmatch.match import cigar, load_vdjam, search_params

OUT = Path("appendix")
# results from bench/gen_vdjam.py (human) and bench/loo_vdjam.py (TRB, scope 2)
REGION_CORR = {"V": (0.08, 0.07), "NDN": (0.42, 0.35), "J": (0.22, 0.11)}  # (TRB, TRA) r(BLOSUM62)
# (epitope, unit, BLOSUM62, VDJAM-flat, VDJAM-region) retrieval PR-AUC, leave-one-out by epitope
LOO = [("GILGFVFTL", 0.579, 0.657, 0.578, 0.585), ("NLVPMVATV", 0.182, 0.268, 0.180, 0.192),
       ("GLCTLVAML", 0.312, 0.354, 0.231, 0.231), ("YLQPRTFLL", 0.679, 0.712, 0.712, 0.724),
       ("RAKFKQLL", 0.145, 0.167, 0.095, 0.095), ("KLGGALQAK", 0.285, 0.287, 0.287, 0.289)]


def _gnuplot(script: str, cwd: Path):
    (cwd / "_tmp.gp").write_text(script)
    subprocess.run(["gnuplot", "_tmp.gp"], cwd=cwd, check=True)
    (cwd / "_tmp.gp").unlink()


def region_corr_fig(figdir: Path):
    rows = "\n".join(f"{i} {reg} {trb} {tra}" for i, (reg, (trb, tra)) in enumerate(REGION_CORR.items()))
    (figdir / "region_corr.dat").write_text("i region TRB TRA\n" + rows + "\n")
    _gnuplot(f"""
set terminal svg size 560,360 font 'Helvetica,13' background rgb 'white'
set output 'region_corr.svg'
set style data histogram
set style histogram cluster gap 1
set style fill solid 0.85 border -1
set boxwidth 0.9
set ylabel 'Pearson r vs BLOSUM62'
set yrange [0:0.5]
set grid ytics lc rgb '#e5e7eb'
set xtics nomirror
set key top left
set title 'Substitution-matrix agreement with BLOSUM62 by CDR3 region'
plot 'region_corr.dat' using 3:xtic(2) title 'TRB' lc rgb '#2563eb', \
     '' using 4 title 'TRA' lc rgb '#f59e0b'
""", figdir)


def retention_fig(figdir: Path):
    """Germline-retention profiles (OLGA-derived, via mirpy) for a few V/J genes: probability each
    CDR3 offset from the anchor is germline-encoded — the basis for the region weight 1 - retention."""
    tsv = "src/vdjmatch/resources/trimming/human_vj_retention.tsv"
    want = [("V", "TRBV19"), ("V", "TRBV9"), ("J", "TRBJ2-7"), ("J", "TRBJ1-1")]
    prof = {k: {} for k in want}
    with open(tsv) as fh:
        next(fh)
        for line in fh:
            chain, seg, gene, off, p, _aa = line.rstrip("\n").split("\t")
            from vdjmatch.match.regions import gene_family
            key = (seg, gene_family(gene))
            if chain == "TRB" and key in prof:
                prof[key][int(off)] = float(p)
    cols = [f"{seg}:{g}" for seg, g in want]
    rows = ["off " + " ".join(c.replace(" ", "") for c in cols)]
    for off in range(12):
        rows.append(str(off) + " " + " ".join(f"{prof[k].get(off, 0.0):.4f}" for k in want))
    (figdir / "retention.dat").write_text("\n".join(rows) + "\n")
    plots = ", ".join(
        f"'retention.dat' using 1:{i+2} with linespoints lw 2 pt 7 ps 0.5 title '{cols[i]}'"
        for i in range(len(want)))
    _gnuplot(f"""
set terminal svg size 560,360 font 'Helvetica,12' background rgb 'white'
set output 'retention.svg'
set xlabel 'CDR3 offset from anchor (aa)'
set ylabel 'P(position is germline-encoded)'
set yrange [0:1.02]
set grid lc rgb '#e5e7eb'
set key bottom left
set title 'OLGA-derived germline-retention (V from N-anchor, J from C-anchor)'
plot {plots}
""", figdir)


def loo_fig(figdir: Path):
    rows = "\n".join(f"{i} {e} {u} {b} {v} {r}" for i, (e, u, b, v, r) in enumerate(LOO))
    (figdir / "loo_prauc.dat").write_text("i epi unit BLOSUM VDJAM region\n" + rows + "\n")
    _gnuplot("""
set terminal svg size 760,400 font 'Helvetica,12' background rgb 'white'
set output 'loo_prauc.svg'
set style data histogram
set style histogram cluster gap 1
set style fill solid 0.85 border -1
set boxwidth 0.9
set ylabel 'retrieval PR-AUC (held-out epitope)'
set yrange [0:0.8]
set grid ytics lc rgb '#e5e7eb'
set xtics nomirror rotate by -30
set key top right
set title 'Leave-one-out-by-epitope retrieval (TRB, scope 2)'
plot 'loo_prauc.dat' using 3:xtic(2) title 'unit' lc rgb '#d1d5db', \\
     '' using 4 title 'BLOSUM62' lc rgb '#9ca3af', \\
     '' using 5 title 'VDJAM (flat)' lc rgb '#93c5fd', \\
     '' using 6 title 'VDJAM (region)' lc rgb '#2563eb'
""", figdir)


def _msa_svg(block: list[str], out: Path, title: str):
    """Render a same-length CDR3 block as a conservation-shaded letter grid (self-contained SVG):
    columns shaded by per-column consensus fraction (conserved flanks dark, variable NDN light)."""
    rows, L = len(block), len(block[0])
    cons = []  # per-column (consensus residue, conservation fraction)
    for p in range(L):
        col = Counter(s[p] for s in block)
        aa, c = col.most_common(1)[0]
        cons.append((aa, c / rows))
    cw, ch, x0, y0 = 22, 17, 8, 50      # y0 = grid top; title + conservation track sit above it
    track_base = y0 - 6                  # bottom of conservation bars (clear gap below the title)
    W, H = x0 + L * cw + 8, y0 + (rows + 1) * ch + 26
    def fill(f):  # white -> blue by conservation
        r = int(255 - f * (255 - 37)); g = int(255 - f * (255 - 99)); b = int(255 - f * (255 - 235))
        return f"rgb({r},{g},{b})"
    e = [f"<svg xmlns='http://www.w3.org/2000/svg' width='{W}' height='{H}' font-family='monospace'>",
         f"<rect width='{W}' height='{H}' fill='white'/>",
         f"<text x='{x0}' y='16' font-size='13' font-family='Helvetica' font-weight='bold'>{title}</text>"]
    for p in range(L):  # conservation track (bars grow up from track_base, below the title)
        bh = cons[p][1] * 14
        e.append(f"<rect x='{x0+p*cw+3}' y='{track_base-bh:.1f}' width='{cw-6}' height='{bh:.1f}' fill='#94a3b8'/>")
    for i, s in enumerate(block):
        for p, aa in enumerate(s):
            x, y = x0 + p * cw, y0 + i * ch
            tc = "white" if cons[p][1] > 0.6 else "#334155"
            e.append(f"<rect x='{x}' y='{y}' width='{cw}' height='{ch}' fill='{fill(cons[p][1])}' stroke='white'/>")
            e.append(f"<text x='{x+cw/2}' y='{y+ch-4}' font-size='12' text-anchor='middle' fill='{tc}'>{aa}</text>")
    cy = y0 + rows * ch + 2  # consensus row
    for p, (aa, _) in enumerate(cons):
        e.append(f"<text x='{x0+p*cw+cw/2}' y='{cy+ch-4}' font-size='12' text-anchor='middle' font-weight='bold'>{aa}</text>")
    e.append(f"<text x='{x0}' y='{cy+ch+14}' font-size='10' fill='#64748b'>consensus; column shade = conservation (dark = conserved flank, light = variable NDN core)</text>")
    e.append("</svg>")
    out.write_text("\n".join(e))


def msa_blocks(vdj: pl.DataFrame, alndir: Path, epis: list[str], chain="TRB", n=14):
    """Same-length CDR3 block per epitope -> a conservation-shaded SVG (+ FASTA for reference)."""
    sub = _bench.valid_cdr3(vdj.filter(pl.col("gene") == chain))
    for epi in epis:
        cds = sub.filter(pl.col("epitope") == epi)["cdr3"].unique().to_list()
        L = Counter(len(c) for c in cds).most_common(1)[0][0]  # modal length
        block = [c for c in cds if len(c) == L][:n]
        (alndir / f"{epi}.fasta").write_text("".join(f">{epi}_{i}\n{c}\n" for i, c in enumerate(block)))
        _msa_svg(block, alndir / f"{epi}_msa.svg", f"{epi}  (TRB CDR3, length {L}, n={len(block)})")


def pairwise_alns(vdj: pl.DataFrame, alndir: Path, epi="GILGFVFTL", chain="TRB"):
    """seqtree pairwise alignments of same-epitope CDR3s, incl. different lengths -> a verbatim block.
    Uses unit (edit-distance) cost for an interpretable display: substitutions stay substitutions,
    indels appear only for genuine length differences."""
    cds = (_bench.valid_cdr3(vdj).filter((pl.col("gene") == chain) & (pl.col("epitope") == epi))
           ["cdr3"].unique().to_list())
    ref = next(c for c in cds if len(c) == 13)
    idx = Index.build(cds, "aa")
    p = search_params("2,2,2,4", engine="seqtm", gap_open=1, gap_extend=1)  # unit cost
    hits = [h for h in idx.search(ref, p) if h.score > 0]
    # pick partners spanning lengths: a couple same-length (subs) + a couple length-changing (indels)
    chosen, lens = [], set()
    for h in sorted(hits, key=lambda h: (abs(len(cds[h.ref_id]) - 13), h.score)):
        L = len(cds[h.ref_id])
        if L not in lens or len([1 for x in chosen if len(cds[x.ref_id]) == L]) < 1:
            chosen.append(h); lens.add(L)
        if len(chosen) >= 5:
            break
    lines = [f"query (ref): {ref}   epitope {epi}   (unit / edit-distance alignment)", ""]
    for h in chosen:
        aln = idx.align(h.ref_id, ref, p)
        lines += [f"vs {cds[h.ref_id]}  (len {len(cds[h.ref_id])}; subs={h.n_subs} ins={h.n_ins} "
                  f"dels={h.n_dels}; CIGAR {cigar.to_cigar(aln.ops)})",
                  f"  {aln.aligned_query}", f"  {cigar.match_line(aln.ops)}", f"  {aln.aligned_ref}", ""]
    (alndir / "pairwise.txt").write_text("\n".join(lines) + "\n")


def main():
    figdir = OUT / "figures"; alndir = OUT / "aln"
    figdir.mkdir(parents=True, exist_ok=True); alndir.mkdir(parents=True, exist_ok=True)
    region_corr_fig(figdir)
    retention_fig(figdir)
    loo_fig(figdir)
    vdj = db.load(_bench.source(), species="HomoSapiens")
    msa_blocks(vdj, alndir, ["GILGFVFTL", "NLVPMVATV"])
    pairwise_alns(vdj, alndir)
    print("wrote appendix figures + alignments")


if __name__ == "__main__":
    main()
