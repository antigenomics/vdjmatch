#!/usr/bin/env python3
"""Annotation accuracy on the high-confidence VDJdb shortlist (>=2 distinct references).

VDJdb labels vary in confidence. The cleanest gold standard is a clonotype-epitope association reported
independently in two or more references (``db.replicated``, n_refs>=2) -- analogous to the mhcmatch
shortlist. We measure vdjmatch's epitope-annotation accuracy on this shortlist by leave-one-out
nearest-neighbour prediction against the rest of the latest VDJdb release: for each shortlist clonotype
we search at the signal:noise-optimal scope (subs 1), exclude its own exact CDR3, and predict the
epitope by distance-weighted neighbour vote. Two rules: CDR3-only, and V-restricted (neighbours must
share the V family -- the strong V prior). Reported against a single-reference control of the same size.

    python bench/shortlist_accuracy.py --species HomoSapiens --subs 1
"""
from __future__ import annotations

import argparse
import os
import random
from collections import defaultdict

import polars as pl
from seqtree import Index, SearchParams

from vdjmatch import db
from vdjmatch.match import regions

VALID = "^[ACDEFGHIKLMNPQRSTVWY]+$"


def predict(cand_ids, ref_cdr3, ref_epi_sets, ref_vfam, qcdr3, qvfam, v_restrict):
    """Distance-weighted neighbour vote -> (predicted_epitope|None, true_reachable_set)."""
    votes = defaultdict(float)
    reach = set()
    for h in cand_ids:
        rc = ref_cdr3[h.ref_id]
        if rc == qcdr3:
            continue                      # exclude the clonotype's own exact CDR3 (self / replicate)
        if v_restrict and qvfam not in ref_vfam[h.ref_id]:
            continue
        w = 1.0 / h.n_subs                # closer neighbours weigh more
        for e in ref_epi_sets[h.ref_id]:
            votes[e] += w
            reach.add(e)
    if not votes:
        return None, reach
    return max(votes, key=votes.get), reach


def accuracy(queries, index, params, ref_cdr3, ref_epi_sets, ref_vfam):
    qc = [c for c, _, _ in queries]
    res = index.search_batch(qc, params, 0)
    n = top1 = top1v = rec = 0
    for (qcdr3, qvfam, true_epi), hl in zip(queries, res):
        n += 1
        p, reach = predict(hl, ref_cdr3, ref_epi_sets, ref_vfam, qcdr3, qvfam, False)
        pv, _ = predict(hl, ref_cdr3, ref_epi_sets, ref_vfam, qcdr3, qvfam, True)
        top1 += (p == true_epi)
        top1v += (pv == true_epi)
        rec += (true_epi in reach)
    return n, top1 / n, top1v / n, rec / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default=os.environ.get("VDJDB_SAMPLE", "test_data/sample3_vdjdb.txt"),
                    help="VDJdb export TSV (default $VDJDB_SAMPLE or test_data/sample3_vdjdb.txt)")
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--subs", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    vdj = db.load(args.pmhc, asset="full", species=args.species)
    short = db.replicated(vdj, min_refs=2)
    rng = random.Random(args.seed)
    print(f"species={args.species}  subs={args.subs}  (latest-release benchmark)\n")
    print(f"{'chain':6}{'set':16}{'n':>7}{'top1 CDR3':>11}{'top1 V-restr':>14}{'recall':>9}")
    for chain in ("TRA", "TRB"):
        cv = vdj.filter((pl.col("gene") == chain) & pl.col("cdr3").str.contains(VALID))
        # reference DB: unique cdr3 -> epitope set, vfam set
        ec = cv.group_by("cdr3").agg(pl.col("epitope").unique(), pl.col("v").unique())
        ref_cdr3 = ec["cdr3"].to_list()
        ref_epi_sets = [set(e) for e in ec["epitope"].to_list()]
        ref_vfam = [{regions.gene_family(x) for x in vs} for vs in ec["v"].to_list()]
        index = Index.build(ref_cdr3, "aa")
        params = SearchParams(max_subs=args.subs, max_total_edits=args.subs, engine="seqtm")

        sl = short.filter((pl.col("gene") == chain) & pl.col("cdr3").str.contains(VALID))
        sl_q = [(c, regions.gene_family(v), e)
                for c, v, e in zip(sl["cdr3"], sl["v"], sl["epitope"])]
        # single-reference control of the same epitopes, matched in size
        sl_epi = set(sl["epitope"].to_list())
        ctrl_all = (cv.filter(pl.col("epitope").is_in(list(sl_epi)))
                      .select("cdr3", "v", "epitope").unique())
        ctrl_pairs = [(c, regions.gene_family(v), e)
                      for c, v, e in zip(ctrl_all["cdr3"], ctrl_all["v"], ctrl_all["epitope"])]
        sl_keys = {(c, e) for c, _, e in sl_q}
        ctrl_pairs = [t for t in ctrl_pairs if (t[0], t[2]) not in sl_keys]
        rng.shuffle(ctrl_pairs)
        ctrl_q = ctrl_pairs[:len(sl_q)]

        for label, qs in (("shortlist (>=2 ref)", sl_q), ("control (1 ref)", ctrl_q)):
            if not qs:
                continue
            n, t1, t1v, r = accuracy(qs, index, params, ref_cdr3, ref_epi_sets, ref_vfam)
            print(f"{chain:6}{label:16}{n:>7}{t1*100:>10.1f}%{t1v*100:>13.1f}%{r*100:>8.1f}%")
    print("\ntop1 = nearest-neighbour-vote epitope == true; V-restr requires same V family; "
          "recall = true epitope reachable within scope (exact self excluded).")


if __name__ == "__main__":
    main()
