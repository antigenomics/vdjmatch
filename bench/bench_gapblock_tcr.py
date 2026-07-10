#!/usr/bin/env python
"""Which gap layout should win: the best score, or the best score plus a prior?

``gapblock_score`` already enumerates every one of the ``L+1`` block positions and keeps the
best. The open question is what "best" means. Four rungs, differing *only* in how the winning
layout is chosen:

    fixed centre   one layout, the centre. The score never gets a vote.
    candidates     five layouts -- after the first 3 or 4 residues, before the last 3 or 4, or
                   the centre -- and the score picks among them. (The originally proposed rule.)
    central prior  every layout, scored, with lam * |midpoint - centre| added before the argmin.
    flat           every layout, scored, and the lowest score wins outright.

They are compared at a **matched false-positive rate**, not a matched score budget. A freer rung
finds lower scores for the same pair, so a fixed budget would hand it more hits and call that a
win. ``seqtree.threshold_for_evalue`` gives each rung, and each query, the cutoff at which its
ball admits ``E*`` expected chance neighbours from the same 250k control. Only then are
precision and recall comparable.

Positives are same-epitope pairs on the VDJdb >=2-reference shortlist; the query's own sequence is
excluded. Everything is reported twice: over all pairs, and over the pairs where the two sequences
differ in length -- the only pairs where the choice of layout can possibly matter.

Usage:
    python bench/bench_gapblock_tcr.py [--queries 2000] [--e-targets 1.0 0.1 0.01]

2026-07-10
"""
from __future__ import annotations

import argparse
import collections
import gzip
import random
import statistics as stt
import sys
import time
from importlib import resources

import polars as pl

import seqtree
from seqtree.evalue import thetas_from_scores
from seqtree.gapblock import GapBlockIndex, central_prior

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from _bench import source, valid_cdr3  # noqa: E402

from vdjmatch import db  # noqa: E402

HUGE = 10 ** 6
D_MAX = 3


def candidate_prior(lam: int, offsets=(3, 4)):
    """Restrict the block to a handful of plausible starts; let the score pick among them.

    Candidates, in the longer sequence's coordinates: ``o`` residues in from the Cys anchor,
    ``o`` residues in from the Phe anchor, and the centre. With ``lam`` large this is a hard
    restriction; the score still chooses. Non-negative and zero at ``d == 0``, so the ball holds.
    """
    def prior(i: int, d: int, m: int) -> int:
        if d == 0:
            return 0
        shorter = m - d
        cands = {shorter // 2}
        for o in offsets:
            cands.add(min(o, shorter))
            cands.add(max(0, shorter - o))
        return lam * min(abs(i - c) for c in cands)
    return prior


def rungs(lam: int):
    return [
        ("fixed centre", central_prior(HUGE)),
        ("candidates (3,4,mid)", candidate_prior(HUGE)),
        (f"central prior lam={lam}", central_prior(lam)),
        ("flat (score alone)", None),
    ]


def shortlist(src, species, gene, min_refs, min_size, cap, seed):
    vdj = valid_cdr3(db.load(src, species=species)).filter(pl.col("gene") == gene)
    keep = (vdj.group_by("epitope")
            .agg(pl.col("reference_id").n_unique().alias("n"))
            .filter(pl.col("n") >= min_refs)["epitope"])
    sub = vdj.filter(pl.col("epitope").is_in(keep))
    by_ep = collections.defaultdict(set)
    for ep, c in zip(sub["epitope"].to_list(), sub["cdr3"].to_list()):
        by_ep[ep].add(c)
    groups = {}
    for ep in sorted(by_ep):
        seqs = sorted(by_ep[ep])
        if len(seqs) < min_size:
            continue
        groups[ep] = random.Random(f"{seed}:{ep}").sample(seqs, cap) if len(seqs) > cap else seqs
    return groups


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pmhc", default=None)
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--gene", default="TRB")
    ap.add_argument("--min-refs", type=int, default=2)
    ap.add_argument("--min-size", type=int, default=20)
    ap.add_argument("--cap", type=int, default=600)
    ap.add_argument("--queries", type=int, default=2000)
    ap.add_argument("--theta-max", type=int, default=50)
    ap.add_argument("--gap-open", type=int, default=None,
                    help="default 2*scale = 28. At 28 a single indel costs more than a "
                         "typical calibrated cutoff, so gapped hits are nearly unreachable; "
                         "pass 14 to let them compete.")
    ap.add_argument("--e-targets", type=float, nargs="*", default=[1.0, 0.1, 0.01])
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    t_start = time.perf_counter()
    with resources.files("seqtree").joinpath("data/control_human_trb_aa.txt.gz").open("rb") as fh:
        pool = [ln.strip() for ln in gzip.open(fh, "rt") if ln.strip()]
    M = len(pool)

    groups = shortlist(source(args.pmhc), args.species, args.gene,
                       args.min_refs, args.min_size, args.cap, args.seed)
    # a CDR3 may serve several epitopes; a hit is a true positive if the epitope sets intersect
    eps_of = collections.defaultdict(set)
    for ep, seqs in groups.items():
        for s in seqs:
            eps_of[s].add(ep)
    refs = sorted(eps_of)
    rid_of = {s: i for i, s in enumerate(refs)}
    N = len(refs)
    print(f"{len(groups)} epitopes, {N:,} unique CDR3 (control M = {M:,})")

    mat = seqtree.SubstitutionMatrix.blosum62()
    scale = mat.scale()
    go = args.gap_open if args.gap_open is not None else 2 * scale
    lam = int(1.5 * scale)

    t0 = time.perf_counter()
    ctrl = GapBlockIndex(pool, "aa", d_max=D_MAX)
    tgt = GapBlockIndex(refs, "aa", d_max=D_MAX)
    print(f"indices built in {time.perf_counter()-t0:.0f}s  "
          f"(gap_open={go}, lam={lam}, theta_max={args.theta_max}, d_max={D_MAX})")

    rng = random.Random(args.seed + 7)
    queries = rng.sample(refs, min(args.queries, N))

    # co-specific partner counts, for the recall denominators
    denom_all, denom_gap = {}, {}
    for q in queries:
        mates = {r for ep in eps_of[q] for r in groups[ep] if r != q}
        denom_all[q] = len(mates)
        denom_gap[q] = sum(1 for r in mates if 1 <= abs(len(r) - len(q)) <= D_MAX)

    results: dict[tuple[str, float, str], list] = collections.defaultdict(list)
    for name, prior in rungs(lam):
        t0 = time.perf_counter()
        ctrl_scores, tgt_hits = [], []
        for q in queries:
            ctrl_scores.append([h[1] for h in ctrl.search(q, args.theta_max, mat, gap_open=go,
                                                          gap_prior=prior)])
            tgt_hits.append([(h[0], h[1]) for h in tgt.search(q, args.theta_max, mat, gap_open=go,
                                                              gap_prior=prior)])
        for e in args.e_targets:
            thetas = thetas_from_scores(ctrl_scores, N, M, e, args.theta_max, exclude_exact=True)
            for stratum in ("all", "gapped"):
                tp = fp = rec_num = rec_den = 0
                usable = 0
                for q, th, hits in zip(queries, thetas, tgt_hits):
                    if th < 0:
                        continue
                    usable += 1
                    den = denom_all[q] if stratum == "all" else denom_gap[q]
                    rec_den += den
                    for rid, s in hits:
                        r = refs[rid]
                        if r == q or s > th:
                            continue
                        d = abs(len(r) - len(q))
                        if stratum == "gapped" and d == 0:
                            continue
                        good = bool(eps_of[q] & eps_of[r])
                        tp += good
                        fp += not good
                        rec_num += good
                prec = tp / (tp + fp) if tp + fp else float("nan")
                rec = rec_num / rec_den if rec_den else float("nan")
                f1 = 2 * prec * rec / (prec + rec) if prec + rec else float("nan")
                results[(name, e, stratum)] = [prec, rec, f1, tp, fp, usable,
                                               stt.median([t for t in thetas if t >= 0] or [-1])]
        print(f"  {name:<24} {time.perf_counter()-t0:5.0f}s", flush=True)

    for stratum in ("all", "gapped"):
        title = ("ALL same-epitope pairs" if stratum == "all"
                 else f"ONLY pairs differing in length (1 <= |dlen| <= {D_MAX})")
        print(f"\n=== {title} ===")
        for e in args.e_targets:
            print(f"\n  E* = {e}   (queries with a usable cutoff / median theta)")
            print(f"  {'rung':<24}{'precision':>11}{'recall':>9}{'F1':>8}{'TP':>8}{'FP':>8}"
                  f"{'usable':>9}{'theta':>7}")
            best_p = max((results[(n, e, stratum)][0] for n, _ in rungs(lam)), default=0)
            best_f = max((results[(n, e, stratum)][2] for n, _ in rungs(lam)), default=0)
            for name, _ in rungs(lam):
                p, r, f1, tp, fp, u, th = results[(name, e, stratum)]
                mp = "*" if p == best_p else " "
                mf = "*" if f1 == best_f else " "
                print(f"  {name:<24}{p:>10.3f}{mp}{r:>9.3f}{f1:>7.3f}{mf}{tp:>8,}{fp:>8,}"
                      f"{u:>9,}{th:>7.0f}")
        print("\n  * marks the best precision and the best F1 in each block.")

    print(f"\ntotal wall {time.perf_counter()-t_start:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
