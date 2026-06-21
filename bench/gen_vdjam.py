#!/usr/bin/env python3
"""Derive a data-driven TCR substitution matrix (VDJAM) from VDJdb.

A BLOSUM/PAM analogue for CDR3s, which cannot be multiple-aligned and have over-conserved
flanks (every TRB CDR3 starts ``CASS`` ...). The observation unit is a **same-antigen
Hamming-1 pair**: two CDR3s annotated to the same epitope, equal length, differing at exactly
one position. The substitution log-odds use a **position-specific background** so over-conserved
columns contribute ~no signal (substitutions there are rare *and* their background expectation is
already large), which removes the germline-conservation bias analytically.

Estimator (per chain c, optionally per region):
    f(a,b)  = count of observed {a,b} substitution events
    e(a,b)  = sum over event slots (L,p) of  n_events(L,p) * 2 * bg(a|L,p) * bg(b|L,p)
    VDJAM(a,b) = log2( (f+alpha) / (e+alpha) )                         (a != b)
where bg(.|L,p) is the residue frequency at position p among all length-L CDR3s of chain c.

Regions (V-germline / NDN / J-germline) are split by position; v1 uses a conserved-flank
heuristic (V = first ``--vlen``, J = last ``--jlen``, NDN = middle). The germline+trimming model
for V/J ("auto-computed from V/J trimming probabilities") is the planned refinement.

    python bench/gen_vdjam.py --pmhc /path/to/vdjdb_full.txt
"""
from __future__ import annotations

import argparse
import math
import os
from collections import Counter, defaultdict

import polars as pl
from seqtree import Index, SearchParams

from vdjmatch import db

AA = "ACDEFGHIKLMNPQRSTVWY"


def hamming1_events(cdr3s: list[str], min_epi: int) -> list[tuple[str, str, int, int]]:
    """All same-set Hamming-1 substitution events (a, b, pos, L) among the CDR3s of one epitope."""
    uc = list(dict.fromkeys(cdr3s))
    if len(uc) < 2:
        return []
    idx = Index.build(uc, "aa")
    p = SearchParams(max_subs=1, max_total_edits=1, engine="seqtm")
    res = idx.search_batch(uc, p, 0)
    out = []
    for i, hl in enumerate(res):
        x = uc[i]
        for h in hl:
            j = h.ref_id
            if j <= i or h.n_subs != 1 or h.n_ins or h.n_dels:
                continue
            y = uc[j]
            if len(x) != len(y):
                continue
            diff = [k for k in range(len(x)) if x[k] != y[k]]
            if len(diff) == 1 and x[diff[0]] in AA and y[diff[0]] in AA:
                out.append((x[diff[0]], y[diff[0]], diff[0], len(x)))
    return out


def position_background(cdr3s: list[str]) -> dict[tuple[int, int], dict[str, float]]:
    """bg[(L, pos)][aa] = residue frequency at position pos among length-L CDR3s."""
    counts: dict[tuple[int, int], Counter] = defaultdict(Counter)
    for s in cdr3s:
        L = len(s)
        for p, a in enumerate(s):
            if a in AA:
                counts[(L, p)][a] += 1
    return {k: {a: c / sum(v.values()) for a, c in v.items()} for k, v in counts.items()}


def estimate_matrix(events, bg, alpha=1.0) -> dict[tuple[str, str], float]:
    """Position-debiased substitution log-odds from events (a,b,pos,L)."""
    f = Counter()
    slots = Counter()  # (L,pos) -> n events there
    for a, b, p, L in events:
        key = tuple(sorted((a, b)))
        f[key] += 1
        slots[(L, p)] += 1
    e = defaultdict(float)
    for (L, p), n in slots.items():
        freq = bg.get((L, p), {})
        items = list(freq.items())
        for ai in range(len(items)):
            a, fa = items[ai]
            for bj in range(ai, len(items)):
                b, fb = items[bj]
                key = tuple(sorted((a, b)))
                e[key] += n * (fa * fb if a == b else 2 * fa * fb)
    score = {}
    for a in AA:
        for b in AA:
            if a < b or a == b:
                key = (a, b)
                score[key] = math.log2((f.get(key, 0) + alpha) / (e.get(key, 0.0) + alpha))
    return score


def corr(m1: dict, m2: dict) -> float:
    """Pearson correlation of off-diagonal substitution scores shared by both matrices."""
    xs, ys = [], []
    for a in AA:
        for b in AA:
            if a < b:
                k = (a, b)
                if k in m1 and k in m2:
                    xs.append(m1[k]); ys.append(m2[k])
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return cov / (sx * sy) if sx and sy else float("nan")


def blosum62_offdiag() -> dict[tuple[str, str], float]:
    from Bio.Align import substitution_matrices
    b = substitution_matrices.load("BLOSUM62")
    return {tuple(sorted((a, bb))): float(b[a, bb]) for a in AA for bb in AA if a < bb}


def legacy_vdjam_offdiag(path) -> dict[tuple[str, str], float]:
    out = {}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            a, bb, s = ln.split()
            if a < bb and a in AA and bb in AA:
                out[(a, bb)] = float(s)
    return out


def region_of(pos: int, L: int, vlen: int, jlen: int) -> str:
    if pos < vlen:
        return "V"
    if pos >= L - jlen:
        return "J"
    return "NDN"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default=os.environ.get("VDJDB_SAMPLE", "test_data/sample3_vdjdb.txt"),
                    help="VDJdb export TSV (default $VDJDB_SAMPLE or test_data/sample3_vdjdb.txt)")
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--min-epi", type=int, default=30)
    ap.add_argument("--vlen", type=int, default=4)
    ap.add_argument("--jlen", type=int, default=6)
    ap.add_argument("--legacy", default="src/vdjmatch/resources/vdjam.txt")
    args = ap.parse_args()

    vdj = db.load(args.pmhc, asset="full", species=args.species)
    blosum = blosum62_offdiag()
    legacy = legacy_vdjam_offdiag(args.legacy)

    print(f"{'chain':5} {'region':5} {'events':>8} {'r_BLOSUM':>9} {'r_legacy':>9}")
    for chain in ("TRA", "TRB"):
        sub = vdj.filter(pl.col("gene") == chain)
        # epitopes with >= min_epi unique CDR3s
        big = (sub.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
                  .filter(pl.col("n") >= args.min_epi)["epitope"].to_list())
        all_cdr3 = sub["cdr3"].unique().to_list()
        bg = position_background(all_cdr3)
        events = []
        for epi in big:
            cdr3s = sub.filter(pl.col("epitope") == epi)["cdr3"].unique().to_list()
            events.extend(hamming1_events(cdr3s, args.min_epi))
        by_region = {"all": events}
        for reg in ("V", "NDN", "J"):
            by_region[reg] = [ev for ev in events if region_of(ev[2], ev[3], args.vlen, args.jlen) == reg]
        for reg, evs in by_region.items():
            if not evs:
                print(f"{chain:5} {reg:5} {0:>8}")
                continue
            m = estimate_matrix(evs, bg)
            print(f"{chain:5} {reg:5} {len(evs):>8} {corr(m, blosum):>9.3f} {corr(m, legacy):>9.3f}")


if __name__ == "__main__":
    main()
