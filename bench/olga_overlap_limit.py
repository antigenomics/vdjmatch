#!/usr/bin/env python3
"""OLGA-overlap theoretical limit: P(a Pgen-drawn TCR coincides with the epitope reference).

The OLGA 'spurious-hit' rate is an estimator of 1 - P(no hit) = the probability that a TCR drawn
from OLGA's generative Pgen distribution falls within the first-hit search scope of *some* reference
sequence. This is a property of (reference, scope, Pgen) -- a coincidental-collision floor set by how
much Pgen mass sits in the reference's edit-neighbourhood, not a method artifact -- and is estimated
directly by generating OLGA TCRs (sample5).

For each OLGA query we count in-scope reference hits (cost>=1) within edit radius k:
    P_overlap(k) = fraction with >=1 hit within k edits      = 1 - P(no hit within k)
    Lambda(k)    = mean #hits within k edits                 = sum_r Pgen-mass of r's k-edit ball
    Poisson cross-check: P_overlap(k) ~= 1 - exp(-Lambda(k)) (independent rare collisions)
The number of reference neighbours of a Pgen draw is the Poisson rate Lambda, so the empirical
neighbour count IS the Pgen integral. Contrast with the significant-call rate (p_enrichment < alpha),
which is what vdjmatch actually flags after control calibration.

    python bench/olga_overlap_limit.py --olga-n 20000
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench
from compare import load_sample
from vdjmatch import db
from vdjmatch.evalue import background, first_hit
from seqtree import Index


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--olga-n", type=int, default=5000, help="OLGA negatives to scan (0 = all; full ~250k is slow)")
    ap.add_argument("--alpha", type=float, default=1e-3, help="significance cutoff for contrast")
    ap.add_argument("--species", default="HomoSapiens")
    args = ap.parse_args()

    vdj = db.load(_bench.source(), species=args.species).filter(pl.col("gene") == "TRB")
    ref = _bench.valid_cdr3(vdj).group_by("cdr3").agg(pl.col("epitope").first())
    ref_cdr3, ref_epi = ref["cdr3"].to_list(), ref["epitope"].to_list()
    tgt = Index.build(ref_cdr3, "aa")
    ctrl = background("TRB")
    N, M = len(tgt), len(ctrl)

    q5 = load_sample("sample5")
    if args.olga_n and q5.height > args.olga_n:
        q5 = q5.sample(args.olga_n, seed=0)
    qlist = q5["cdr3"].to_list()
    print(f"reference N={N} unique TRB CDR3; control M={M}; OLGA queries={len(qlist)}")

    th, cc = first_hit.scan(tgt, ref_epi, ctrl, qlist, exclude_exact=True, progress=True)

    print(f"\nOLGA-overlap theoretical limit  (n={len(qlist)} Pgen draws, scope <=5 edits / <=2 ins / <=2 del)")
    print(f"{'radius k':>8} {'P_overlap=1-P(no hit)':>22} {'Lambda (mean #hits)':>20} {'1-exp(-Lambda)':>15}")
    for k in (1, 2, 3, 4, 5):
        nhits = [sum(1 for c, _ in t if c <= k) for t in th]
        n = len(nhits)
        p_overlap = sum(1 for x in nhits if x > 0) / n
        lam = sum(nhits) / n
        print(f"{k:>8} {p_overlap:>22.4f} {lam:>20.4f} {1 - math.exp(-lam):>15.4f}")

    sig = sum(1 for t, c in zip(th, cc)
              if first_hit.pvalue(t, c, N, M)["p_enrichment"] < args.alpha) / len(qlist)
    print(f"\nfor contrast, vdjmatch significant-call rate (p_enrichment < {args.alpha}): {sig:.4f}")
    print("  -> raw overlap is the Pgen-set coincidental ceiling; the E-value test flags only the "
          "control-enriched subset.")


if __name__ == "__main__":
    main()
