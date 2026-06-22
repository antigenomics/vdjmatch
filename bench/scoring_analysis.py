#!/usr/bin/env python3
"""CDR3 scoring analyses for the appendix: signal:noise vs edit distance, positional substitution
significance (centre vs V/J borders), and the TRA-vs-TRB central composition (D-gene glycine).

Emits SVGs into appendix/figures/: purity_vs_distance.svg, position_significance.svg, tra_trb_gly.svg.

    python bench/scoring_analysis.py
"""
from __future__ import annotations

import os
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import polars as pl
from seqtree import Index, SearchParams

import _bench
from vdjmatch import db

OUT = Path("appendix/figures")
AA = "ACDEFGHIKLMNPQRSTVWY"


def _gnuplot(script: str):
    (OUT / "_t.gp").write_text(script)
    subprocess.run(["gnuplot", "_t.gp"], cwd=OUT, check=True)
    (OUT / "_t.gp").unlink()


def cdr3_epitopes(df: pl.DataFrame) -> dict[str, frozenset]:
    """{cdr3 -> frozenset(epitopes)} from a clonotype frame (already composition-controlled)."""
    d = defaultdict(set)
    for c, e in zip(df["cdr3"], df["epitope"]):
        d[c].add(e)
    return {c: frozenset(s) for c, s in d.items()}


def purity_vs_distance(ll, chain="TRB", n_queries=500, dmax=5, seed=0, org="human"):
    """Macro purity (per-epitope average of P[neighbour shares epitope]) and lift (purity / chance
    pi(E) -- enrichment over the admixed random-control rate, the E-value view), by Hamming distance,
    on the composition-controlled long list. Macro + lift are robust to the ~100x epitope-size range
    that makes pooled (micro) purity read ~constant on the 10x-dominated raw release."""
    import statistics as _st
    ce = cdr3_epitopes(ll)
    cdrs = list(ce)
    N = len(cdrs)
    epi_n = Counter(e for s in ce.values() for e in s)
    idx = Index.build(cdrs, "aa")
    rng = __import__("random").Random(seed)
    q = rng.sample(cdrs, min(n_queries, N))
    res = idx.search_batch(q, SearchParams(max_subs=dmax, max_total_edits=dmax, engine="seqtm"), 0)
    macro, lift = defaultdict(list), defaultdict(list)
    for qi, hl in zip(q, res):
        qe = ce[qi]
        s, n = Counter(), Counter()
        for h in hl:
            k = h.n_subs
            if k == 0 or cdrs[h.ref_id] == qi:
                continue                                          # never count exact self-hits
            n[k] += 1
            if ce[cdrs[h.ref_id]] & qe:
                s[k] += 1
        pi = max(epi_n[e] for e in qe) / N                        # chance a random clonotype shares qe
        for k in n:
            p = s[k] / n[k]
            macro[k].append(p); lift[k].append(p / pi if pi else float("nan"))
    pur = {k: _st.mean(macro[k]) for k in sorted(macro)}
    lif = {k: _st.mean(lift[k]) for k in sorted(lift)}
    print(f"[{org} {chain}] macro purity:", {k: round(v, 3) for k, v in pur.items()},
          "| lift:", {k: round(v, 1) for k, v in lif.items()})
    rows = "\n".join(f"{k} {pur[k]:.4f} {lif[k]:.2f}" for k in sorted(pur))
    (OUT / "purity.dat").write_text("d purity lift\n" + rows + "\n")
    _gnuplot(f"""
set terminal svg size 540,360 font 'Helvetica,12' background rgb 'white'
set output 'purity_vs_distance.svg'
set xlabel 'CDR3 Hamming distance (substitutions)'; set ylabel 'macro purity (per-epitope mean)'
set xrange [0.5:{dmax+0.5}]; set yrange [0:*]; set grid lc rgb '#e5e7eb'; set xtics 1
set key top right; set title 'Signal:noise vs edit distance ({org} {chain})'
plot 'purity.dat' using 1:2 with linespoints lw 2.5 pt 7 lc rgb '#2563eb' notitle
""")
    return pur


def _eb_shrink(binned):
    """Beta-Binomial empirical-Bayes posterior-mean estimator for per-bin rates. Fits a global
    Beta(a,b) prior by moments (mean + between-bin variance) and returns f(same,tot)=(same+a)/(tot+a+b),
    so sparse offset bins are shrunk toward the global rate. Returns (global_mean, f)."""
    pop = [(s, t) for s, t in binned if t > 0]
    tot = sum(t for _, t in pop)
    m = sum(s for s, _ in pop) / tot
    var = sum(t * (s / t - m) ** 2 for s, t in pop) / tot
    k = max(5.0, min(1000.0, m * (1 - m) / var - 1)) if var > 1e-9 else 200.0
    a, b = m * k, (1 - m) * k
    return m, (lambda s, t: (s + a) / (t + a + b))


def position_significance(ll, chain="TRB", n_queries=4000, maxd=8, seed=0, org="human"):
    """P(neighbours share epitope | a single substitution at offset d from the V/J anchors), end-
    anchored (not relative position; the germline anchors sit at fixed absolute offsets) and
    Beta-Binomial-smoothed. Bundles the two end profiles as the positional-significance PSSM factor and
    reports the BLOSUM-severity x position interaction. On the composition-controlled long list."""
    from Bio.Align import substitution_matrices
    blo = substitution_matrices.load("BLOSUM62")
    ce = cdr3_epitopes(ll)
    cdrs = list(ce)
    idx = Index.build(cdrs, "aa")
    rng = __import__("random").Random(seed)
    q = rng.sample(cdrs, min(n_queries, len(cdrs)))
    res = idx.search_batch(q, SearchParams(max_subs=1, max_total_edits=1, engine="seqtm"), 0)
    sv = [[0, 0] for _ in range(maxd + 1)]   # [same,tot] by V-offset (>=maxd pooled = NDN core)
    sj = [[0, 0] for _ in range(maxd + 1)]   # by J-offset
    grid = {}                                # (offset-zone, BLOSUM-severity) -> [same,tot]
    for qi, hl in zip(q, res):
        qe = ce[qi]; L = len(qi)
        for h in hl:
            if h.n_subs != 1 or len(cdrs[h.ref_id]) != L or L < 2:
                continue
            r = cdrs[h.ref_id]
            diff = [p for p in range(L) if qi[p] != r[p]]
            if len(diff) != 1:
                continue
            p = diff[0]; lab = 1 if ce[r] & qe else 0
            dv, dj = min(p, maxd), min(L - 1 - p, maxd)
            sv[dv][0] += lab; sv[dv][1] += 1
            sj[dj][0] += lab; sj[dj][1] += 1
            d = min(dv, dj)
            zone = "anchor" if d <= 1 else "mid" if d <= 3 else "core"
            s = blo[qi[p], r[p]]
            sev = "conservative" if s >= 1 else "neutral" if s == 0 else "radical"
            g = grid.setdefault((zone, sev), [0, 0]); g[0] += lab; g[1] += 1
    _, fv = _eb_shrink(sv); _, fj = _eb_shrink(sj)
    pv = [fv(*sv[d]) for d in range(maxd + 1)]
    pj = [fj(*sj[d]) for d in range(maxd + 1)]
    print(f"[{org} {chain}] P(same|V-offset 0..{maxd}):", [round(x, 2) for x in pv])
    print(f"[{org} {chain}] P(same|J-offset 0..{maxd}):", [round(x, 2) for x in pj])
    print("BLOSUM-severity x position  P(same) (n):")
    for zone in ("anchor", "mid", "core"):
        print(f"  {zone:7}", {sev: (round(grid[zone, sev][0] / grid[zone, sev][1], 2),
                                     grid[zone, sev][1]) for sev in ("conservative", "neutral", "radical")
                              if (zone, sev) in grid})
    res = Path("src/vdjmatch/resources/trimming/position_significance.tsv")
    res.write_text("side\toffset\tn\tp_same\n" + "".join(
        f"V\t{d}\t{sv[d][1]}\t{pv[d]:.4f}\nJ\t{d}\t{sj[d][1]}\t{pj[d]:.4f}\n" for d in range(maxd + 1)))
    (OUT / "possig.dat").write_text("offset psame_v psame_j\n"
                                    + "\n".join(f"{d} {pv[d]:.4f} {pj[d]:.4f}" for d in range(maxd + 1)) + "\n")
    _gnuplot(f"""
set terminal svg size 560,360 font 'Helvetica,12' background rgb 'white'
set output 'position_significance.svg'
set xlabel 'offset from anchor (residues)'; set ylabel 'P(neighbour shares epitope | sub here)'
set xrange [-0.3:{maxd}]; set yrange [0:*]; set grid lc rgb '#e5e7eb'; set key bottom right
set title 'Single-substitution significance vs anchor offset ({org} {chain})'
plot 'possig.dat' using 1:2 with linespoints lw 2.5 pt 7 lc rgb '#7c3aed' title 'from V anchor', \\
     '' using 1:3 with linespoints lw 2.5 pt 5 lc rgb '#059669' title 'from J anchor'
""")
    return pv, pj


def tra_trb_gly(vdj, bins=10, org="human"):
    """Glycine fraction by relative CDR3 position, TRA vs TRB (TRB D-gene -> central poly-G).
    On composition-controlled long lists per chain."""
    def gly_profile(chain):
        cdrs = _bench.long_list(vdj.filter(pl.col("gene") == chain), cap=3000, min_n=30)["cdr3"].to_list()
        g = [0] * bins; n = [0] * bins
        for s in cdrs:
            if len(s) < 2:
                continue
            for i, a in enumerate(s):
                b = min(bins - 1, int(i / (len(s) - 1) * bins))
                n[b] += 1; g[b] += (a == "G")
        return [g[b] / n[b] if n[b] else 0.0 for b in range(bins)]
    pa, pb = gly_profile("TRA"), gly_profile("TRB")
    print("[TRA] central Gly frac:", round(pa[bins // 2], 3), "| [TRB]:", round(pb[bins // 2], 3))
    rows = "\n".join(f"{(b+0.5)/bins:.3f} {pa[b]:.4f} {pb[b]:.4f}" for b in range(bins))
    (OUT / "gly.dat").write_text("relpos TRA TRB\n" + rows + "\n")
    _gnuplot("""
set terminal svg size 560,360 font 'Helvetica,12' background rgb 'white'
set output 'tra_trb_gly.svg'
set xlabel 'relative CDR3 position'; set ylabel 'glycine (G) fraction'
set xrange [0:1]; set yrange [0:*]; set grid lc rgb '#e5e7eb'; set key top right
set title 'Central glycine enrichment ({org}): TRB (D gene) vs TRA'
plot 'gly.dat' using 1:3 with linespoints lw 2.5 pt 7 lc rgb '#2563eb' title 'TRB', \
     '' using 1:2 with linespoints lw 2.5 pt 5 lc rgb '#f59e0b' title 'TRA'
""")
    return pa, pb


def matrices_fig():
    """Summary bar chart (2026-06-11-ZENODO, composition-controlled long list, balanced PR-AUC, dist
    <=2): the data-derived VDJAM trails; the genetic-code null VDJAMr ties BLOSUM62; only BLOSUM62
    reweighted by the central-position factor (possig) clearly leads. Means from bench/loo_vdjam.py."""
    data = [("VDJAM", 0.524), ("VDJAM-region", 0.529), ("unit", 0.543), ("BLOSUM62", 0.564),
            ("structural", 0.564), ("PAM250", 0.572), ("VDJAMr", 0.581), ("BLOSUM+possig", 0.598)]
    win = len(data) - 1  # the position-weighted bar, drawn in a contrasting colour
    rows = "\n".join(f"{i} {n} {v}" for i, (n, v) in enumerate(data))
    (OUT / "matrices.dat").write_text("i name prauc\n" + rows + "\n")
    tics = ", ".join(f"'{n}' {i}" for i, (n, _) in enumerate(data))
    _gnuplot(f"""
set terminal svg size 620,360 font 'Helvetica,12' background rgb 'white'
set output 'matrices.svg'
set style fill solid 0.85 border -1; set boxwidth 0.6
set ylabel 'mean balanced PR-AUC (leave-one-out)'; set yrange [0.5:0.62]
set grid ytics lc rgb '#e5e7eb'; set xtics ({tics}) rotate by -25; unset key
set title 'Only central-position weighting clearly beats BLOSUM62 (human TRB, dist <= 2)'
plot 'matrices.dat' using 1:3 with boxes lc rgb '#94a3b8' notitle, \\
     'matrices.dat' using 1:($1=={win}?$3:1/0) with boxes lc rgb '#059669' notitle
""")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Example case: human TRB/TRA. The same pipeline applies to human/mouse and TRA/TRB -- point
    # VDJDB_SAMPLE at another VDJdb export and set VDJDB_SPECIES to regenerate for a different case.
    species = os.environ.get("VDJDB_SPECIES", "HomoSapiens")
    org = {"HomoSapiens": "human", "MusMusculus": "mouse"}.get(species, species)
    vdj = db.load(_bench.source(), species=species)
    ll = _bench.long_list(vdj.filter(pl.col("gene") == "TRB"), cap=3000, min_n=30)
    purity_vs_distance(ll, "TRB", org=org)
    position_significance(ll, "TRB", org=org)
    tra_trb_gly(vdj, org=org)
    matrices_fig()
    print("wrote purity_vs_distance.svg, position_significance.svg, tra_trb_gly.svg, matrices.svg")


if __name__ == "__main__":
    main()
