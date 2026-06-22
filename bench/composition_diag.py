#!/usr/bin/env python3
"""Benchmark composition diagnostics for the dense 2026 VDJdb release.

The raw release is dominated by a few 10x mega-studies (SLLMWITQV ~30k clonotypes) and contains
phage-display / length-degenerate epitopes, which makes naive purity (P[neighbour shares epitope])
read ~0.9 regardless of distance. This script validates the composition controls:

  1. detect phage-display (length-degenerate) epitopes;
  2. build a composition-controlled long list (cap mega-epitopes to a random N, quarantine phage);
  3. purity vs CDR3 edit distance three ways:
       - micro  : pooled P(neighbour shares epitope)            (mega-epitope dominated)
       - macro  : per-epitope purity, averaged over epitopes    (size-balanced)
       - lift   : purity / chance pi(E)  == enrichment over an admixed random control (the E-value view)
  4. rarefaction: how deep must you sample a mega-epitope to cover its specificity space?

    python bench/composition_diag.py --chain TRB --cap 3000
"""
from __future__ import annotations

import argparse
import random
import statistics as st
from collections import Counter, defaultdict

import polars as pl
from seqtree import Index, SearchParams

import _bench
from vdjmatch import db


def purity_three_ways(ll: pl.DataFrame, chain: str, nq: int, dmax: int, seed: int):
    ce: dict[str, set] = defaultdict(set)
    for c, e in zip(ll["cdr3"], ll["epitope"]):
        ce[c].add(e)
    cdrs = list(ce)
    N = len(cdrs)
    epi_n = Counter(e for s in ce.values() for e in s)           # clonotypes per epitope (long list)
    idx = Index.build(cdrs, "aa")
    rng = random.Random(seed)
    q = rng.sample(cdrs, min(nq, N))
    res = idx.search_batch(q, SearchParams(max_subs=dmax, max_total_edits=dmax, engine="seqtm"), 0)
    micro_s, micro_n = Counter(), Counter()
    macro, lift = defaultdict(list), defaultdict(list)
    for qi, hl in zip(q, res):
        qe = ce[qi]
        per_k_s, per_k_n = Counter(), Counter()
        for h in hl:
            k = h.n_subs
            if k == 0 or cdrs[h.ref_id] == qi:
                continue                                          # never count exact self-hits
            per_k_n[k] += 1
            if ce[cdrs[h.ref_id]] & qe:
                per_k_s[k] += 1
        pi = max(epi_n[e] for e in qe) / N                       # chance a random clonotype shares qe
        for k in per_k_n:
            micro_s[k] += per_k_s[k]; micro_n[k] += per_k_n[k]
            p = per_k_s[k] / per_k_n[k]
            macro[k].append(p)
            lift[k].append(p / pi if pi else float("nan"))
    return {k: (round(micro_s[k] / micro_n[k], 3), round(st.mean(macro[k]), 3),
                round(st.mean(lift[k]), 1)) for k in sorted(micro_n)}


def rarefaction(full: list[str], depths, seed=0):
    """Coverage = fraction of held-out clonotypes within edit distance 1 of a depth-d random sample
    (how much of the specificity space a sample of size d already 'reaches')."""
    rng = random.Random(seed)
    pool = full[:]
    rng.shuffle(pool)
    out = []
    holdout = pool[:4000]                                        # fixed held-out probe set
    rest = pool[4000:]
    p1 = SearchParams(max_subs=1, max_total_edits=1, engine="seqtm")
    for d in depths:
        if d > len(rest):
            break
        idx = Index.build(rest[:d], "aa")
        hit = sum(1 for hl in idx.search_batch(holdout, p1, 0) if any(h.n_subs <= 1 for h in hl))
        out.append((d, round(hit / len(holdout), 3)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRB")
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--cap", type=int, default=3000)
    ap.add_argument("--nq", type=int, default=500)
    ap.add_argument("--dmax", type=int, default=5)
    args = ap.parse_args()

    vdj = db.load(_bench.source(), species=args.species).filter(pl.col("gene") == args.chain)
    uc = _bench.valid_cdr3(vdj).select("cdr3", "v", "j", "epitope", "mhc_class").unique()
    print(f"{args.species} {args.chain}: {uc.height} unique clonotypes, "
          f"{uc['epitope'].n_unique()} epitopes")

    spec = _bench.spectratype(vdj).filter(pl.col("n") >= 30)
    anom = _bench.spectratype_anomalies(vdj)
    print(f"\nspectratype (n>=30 epitopes): global eff_len median "
          f"{spec['eff_len'].median():.1f}; anomalous (eff_len<2 or modal_frac>=0.9): {len(anom)}")
    if anom:
        print("  e.g.", spec.filter(pl.col("epitope").is_in(list(anom)))
              .sort("eff_len").head(6).to_dicts())

    mplx = _bench.multiplex_studies(vdj, max_epitopes=100)
    print(f"\nmultiplex studies (reference reporting >100 epitopes): {mplx.height}")
    print("  ", mplx.head(5).to_dicts())

    sizes = uc.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
    mega = sizes.filter(pl.col("n") > args.cap).sort("n", descending=True)
    print(f"\nmega-epitopes (>{args.cap} clonotypes): {mega.height}")
    print("  ", mega.head(6).to_dicts())

    ll = _bench.long_list(vdj, cap=args.cap, min_n=30)
    print(f"\nlong list (cap {args.cap}, min_n 30, spectratype-anomalous quarantined): {ll.height} "
          f"clonotypes, {ll['epitope'].n_unique()} epitopes")

    print(f"\npurity vs CDR3 edit distance  (dist: micro, macro, lift):")
    for k, (mi, ma, li) in purity_three_ways(ll, args.chain, args.nq, args.dmax, seed=0).items():
        print(f"  d={k}:  micro {mi:.3f}   macro {ma:.3f}   lift {li:.1f}x")

    # rarefaction on the single largest epitope (sampling depth -> coverage)
    top = sizes.sort("n", descending=True).row(0)
    full = uc.filter(pl.col("epitope") == top[0])["cdr3"].unique().to_list()
    print(f"\nrarefaction on {top[0]} (n={len(full)}): depth -> coverage (held-out within dist 1)")
    for d, cov in rarefaction(full, [100, 300, 1000, 3000, 10000, 20000]):
        print(f"  depth {d:>6}: {cov*100:5.1f}%")
    print("\nlift = purity / chance pi(E): enrichment over an admixed random control (the random tail "
          "of the epitope). lift>>1 = real motif signal; lift~1 = same-epitope-by-size-chance.")


if __name__ == "__main__":
    main()
