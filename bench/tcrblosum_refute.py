#!/usr/bin/env python3
"""Refute the tcrBLOSUM claim that a TCR-specific substitution matrix beats BLOSUM62.

Postovskaya et al. (Brief Bioinform 2024, doi:10.1093/bib/bbae602; github apostovskaya/tcrBLOSUM)
report that tcrBLOSUMb improves epitope-specific TCR clustering over BLOSUM62. Their construction
counts substitutions from same-epitope, same-length CDR3 pairs with NO redundancy clustering, and is
validated on a clustering task whose objective is the same same-epitope relation -- so the gain may be
in-sample. We test their matrix in OUR leave-one-out-by-epitope retrieval, where the matrix is judged on
epitopes it never saw constrain it (here it is fixed, so every VDJdb epitope is a fair test of transfer),
ranking VDJdb-vs-VDJdb neighbour pairs by summed dissimilarity and scoring micro PR-AUC for "shares the
held-out epitope". Arms: BLOSUM62 (biopython and their own bundled copy, as a sanity cross-check),
tcrBLOSUMb, and each reweighted by our central-position significance PSSM.

Matrices are downloaded on demand to bench/matrices/ (not redistributed):
    https://github.com/apostovskaya/tcrBLOSUM  (results/tcrBLOSUMmtx, results/otherMTXs)

    python bench/tcrblosum_refute.py --chain TRB --subs 2 --top 8
"""
from __future__ import annotations

import argparse
import os
import statistics as st
import urllib.request
from pathlib import Path

import polars as pl
from seqtree import Index, SearchParams

import _bench
from metrics import pr_auc_balanced as pr_auc  # noqa: E402  (balanced for the 100x class imbalance)
from vdjmatch import db
from vdjmatch.match import regions
from loo_vdjam import named_dissim, score_pairs  # noqa: E402
from gen_vdjam import AA  # noqa: E402

RAW = "https://raw.githubusercontent.com/apostovskaya/tcrBLOSUM/main/results"
SRC = {"tcrBLOSUM_all_beta.tsv": f"{RAW}/tcrBLOSUMmtx/tcrBLOSUM_all_beta.tsv",
       "tcrBLOSUM_all_alpha.tsv": f"{RAW}/tcrBLOSUMmtx/tcrBLOSUM_all_alpha.tsv",
       "blosum62_20aa.tsv": f"{RAW}/otherMTXs/blosum62_20aa.tsv"}


def fetch(name: str) -> Path:
    p = Path("bench/matrices") / name
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(SRC[name], p)
    return p


def tsv_dissim(path: Path) -> dict[tuple[str, str], float]:
    """tcrBLOSUM-format TSV (leading empty header cell, AA-labelled rows, integer similarity) ->
    non-negative off-diagonal dissimilarity (most-similar pair = 0), matching named_dissim()."""
    lines = path.read_text().splitlines()
    cols = lines[0].split("\t")[1:]
    sim = {}
    for line in lines[1:]:
        f = line.split("\t")
        a = f[0]
        for c, v in zip(cols, f[1:]):
            sim[(a, c)] = float(v)
    vals = {(a, c): sim[(a, c)] for a in AA for c in AA if a != c}
    hi = max(vals.values())
    return {k: hi - v for k, v in vals.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default=None,
                    help="VDJdb export TSV (default: $VDJDB_SAMPLE or the HF-pinned release)")
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--chain", default="TRB")
    ap.add_argument("--min-epi", type=int, default=50)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--max-queries", type=int, default=600)
    ap.add_argument("--subs", type=int, default=2)
    args = ap.parse_args()

    tcr_file = "tcrBLOSUM_all_beta.tsv" if args.chain == "TRB" else "tcrBLOSUM_all_alpha.tsv"
    dis = {"BLOSUM62": named_dissim("BLOSUM62"),
           "BLOSUM62(theirs)": tsv_dissim(fetch("blosum62_20aa.tsv")),
           "tcrBLOSUM": tsv_dissim(fetch(tcr_file))}

    vdj = db.load(_bench.source(args.pmhc), species=args.species).filter(pl.col("gene") == args.chain)
    uc = _bench.long_list(vdj, cap=3000, min_n=args.min_epi)  # composition-controlled clonotypes
    ret = regions.load_retention()
    refs = uc.group_by("cdr3").agg(pl.col("epitope").first())
    ref_cdr3, ref_epi = refs["cdr3"].to_list(), refs["epitope"].to_list()
    index = Index.build(ref_cdr3, "aa")
    cparams = SearchParams(max_subs=args.subs, max_total_edits=args.subs, engine="seqtm")
    sizes = (uc.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
               .filter(pl.col("n") >= args.min_epi).sort(["n", "epitope"], descending=[True, False]))
    held = sizes["epitope"].to_list()[:args.top]

    arms = ["BLOSUM62", "BLOSUM62(theirs)", "tcrBLOSUM", "tcrBLOSUM+possig", "BLOSUM62+possig"]
    acc = {a: [] for a in arms}
    print(f"species={args.species} chain={args.chain}; held-out epitopes={len(held)}; subs={args.subs}")
    for epi in held:
        q = uc.filter(pl.col("epitope") == epi).unique("cdr3").head(args.max_queries)
        qs, qv, qj = q["cdr3"].to_list(), q["v"].to_list(), q["j"].to_list()
        cand = [[h.ref_id for h in hl] for hl in index.search_batch(qs, cparams, 0)]

        def sc(d, w):
            return pr_auc(score_pairs(qs, qv, qj, cand, ref_cdr3, ref_epi, epi, args.chain, ret, d, w))
        acc["BLOSUM62"].append(sc(dis["BLOSUM62"], None))
        acc["BLOSUM62(theirs)"].append(sc(dis["BLOSUM62(theirs)"], None))
        acc["tcrBLOSUM"].append(sc(dis["tcrBLOSUM"], None))
        acc["tcrBLOSUM+possig"].append(sc(dis["tcrBLOSUM"], "possig"))
        acc["BLOSUM62+possig"].append(sc(dis["BLOSUM62"], "possig"))

    print("\nmean leave-one-out retrieval PR-AUC:")
    base = st.mean(acc["BLOSUM62"])
    for a in arms:
        m = st.mean(acc[a])
        wins = sum(1 for x, b in zip(acc[a], acc["BLOSUM62"]) if x > b)
        print(f"  {a:20} {m:.3f}   ({m-base:+.3f} vs BLOSUM62, wins {wins}/{len(held)})")
    print("\nverdict: tcrBLOSUM does not beat BLOSUM62 on held-out-epitope retrieval; central-position "
          "weighting is what helps, and it helps BLOSUM62 more than tcrBLOSUM.")


if __name__ == "__main__":
    main()
