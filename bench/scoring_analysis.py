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


def position_significance(ll, chain="TRB", n_queries=2000, bins=8, seed=0, org="human"):
    """P(neighbours share epitope | the single substitution falls at relative CDR3 position).
    Centre (NDN/contact) vs V/J borders (germline-fixed). On the composition-controlled long list."""
    ce = cdr3_epitopes(ll)
    cdrs = list(ce)
    idx = Index.build(cdrs, "aa")
    rng = __import__("random").Random(seed)
    q = rng.sample(cdrs, min(n_queries, len(cdrs)))
    res = idx.search_batch(q, SearchParams(max_subs=1, max_total_edits=1, engine="seqtm"), 0)
    same = [0] * bins; tot = [0] * bins
    for qi, hl in zip(q, res):
        qe = ce[qi]
        for h in hl:
            if h.n_subs != 1:
                continue
            r = cdrs[h.ref_id]
            if len(r) != len(qi):
                continue
            diff = [p for p in range(len(qi)) if qi[p] != r[p]]
            if len(diff) != 1 or len(qi) < 2:
                continue
            b = min(bins - 1, int(diff[0] / (len(qi) - 1) * bins))
            tot[b] += 1
            if ce[r] & qe:
                same[b] += 1
    prob = [same[b] / tot[b] if tot[b] else 0.0 for b in range(bins)]
    print(f"[{chain}] P(same epitope | sub at rel-pos bin):", [round(p, 3) for p in prob],
          "| n:", [tot[b] for b in range(bins)])
    # bundle the well-populated profile as a positional informativeness factor for scoring:
    # weight(relpos) = 1 - P(same | sub here), normalised to mean 1 (centre > borders).
    good = [(b, prob[b]) for b in range(bins) if tot[b] >= 20]
    raw = [1.0 - p for _, p in good]
    mean = sum(raw) / len(raw)
    res = Path("src/vdjmatch/resources/trimming/position_significance.tsv")
    res.write_text("relpos\tp_same\tweight\n" + "".join(
        f"{(b+0.5)/bins:.4f}\t{p:.4f}\t{(1.0-p)/mean:.4f}\n" for b, p in good))
    # only plot well-populated bins (the conserved flanks rarely vary -> tiny, noisy counts)
    rows = "\n".join(f"{(b+0.5)/bins:.3f} {prob[b]:.4f} {tot[b]}" for b in range(bins) if tot[b] >= 20)
    (OUT / "possig.dat").write_text("relpos psame n\n" + rows + "\n")
    _gnuplot(f"""
set terminal svg size 560,360 font 'Helvetica,12' background rgb 'white'
set output 'position_significance.svg'
set xlabel 'relative CDR3 position (0 = V/Cys anchor, 1 = J anchor)'
set ylabel 'P(neighbour shares epitope | sub here)'
set xrange [0:1]; set yrange [0:*]; set grid lc rgb '#e5e7eb'
set title 'Single-substitution significance by CDR3 position ({org} {chain})'; unset key
set label 'V border' at 0.06,graph 0.08 tc rgb '#9ca3af'; set label 'NDN core' at 0.42,graph 0.08 tc rgb '#9ca3af'
set label 'J border' at 0.80,graph 0.08 tc rgb '#9ca3af'
plot 'possig.dat' using 1:2 with linespoints lw 2.5 pt 7 lc rgb '#7c3aed' notitle
""")
    return prob


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
    data = [("VDJAM", 0.519), ("VDJAM-region", 0.525), ("unit", 0.545), ("structural", 0.572),
            ("BLOSUM62", 0.572), ("PAM250", 0.575), ("VDJAMr", 0.585), ("BLOSUM+possig", 0.623)]
    win = len(data) - 1  # the position-weighted bar, drawn in a contrasting colour
    rows = "\n".join(f"{i} {n} {v}" for i, (n, v) in enumerate(data))
    (OUT / "matrices.dat").write_text("i name prauc\n" + rows + "\n")
    tics = ", ".join(f"'{n}' {i}" for i, (n, _) in enumerate(data))
    _gnuplot(f"""
set terminal svg size 620,360 font 'Helvetica,12' background rgb 'white'
set output 'matrices.svg'
set style fill solid 0.85 border -1; set boxwidth 0.6
set ylabel 'mean balanced PR-AUC (leave-one-out)'; set yrange [0.45:0.65]
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
