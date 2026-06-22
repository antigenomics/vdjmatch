#!/usr/bin/env python3
"""Is the OLGA spurious rate control-mismatch, or genuine coincidental overlap? (chain-consistent)

CHAIN-CONSISTENT design: queries = sample4 (TRB OLGA), reference = VDJdb-beta (TRB), all controls TRB.
(An earlier version used sample5 as queries — sample5 is ALPHA chain, so its 0% with a "matched" control
was a same-CHAIN artifact, not a real result.) We compare the TRB-OLGA significant-call rate under three
TRB controls at matched size: the real selected repertoire (`human_trb_aa`), the SAME OLGA run
(sample4 disjoint half), and an INDEPENDENT freshly-generated OLGA set (different seed) — the latter
resolves the fixed-seed caveat. If the independent control matches the real control (and the same-run
control is the outlier), the spurious rate is genuine coincidental overlap.

    python bench/control_pgen_test.py --n-query 5000 --n-control 50000
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench
from compare import load_sample
from vdjmatch import db
from vdjmatch.evalue import background, first_hit
from seqtree import Index


def load_fresh_olga(path: str) -> pl.DataFrame:
    """A raw olga-generate_sequences TSV (headerless: nt, cdr3_aa, V, J) -> valid unique CDR3."""
    d = pl.read_csv(path, separator="\t", has_header=False,
                    new_columns=["nt", "cdr3", "v", "j"]).select("cdr3")
    return _bench.valid_cdr3(d).unique("cdr3").sample(fraction=1.0, shuffle=True, seed=0)


def cost_lists(idx, queries, params, desc, prog=True):
    return first_hit._cost_lists(idx, queries, params, 0, True, 10000, desc, prog)


def sig_rate(th, cc_raw, N, M, alpha):
    cc = [[c for c, _ in t] for t in cc_raw]
    hits = sum(1 for t, c in zip(th, cc) if first_hit.pvalue(t, c, N, M)["p_enrichment"] < alpha)
    return hits / len(th)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-query", type=int, default=5000)
    ap.add_argument("--n-control", type=int, default=50000, help="OLGA & real control size (matched); full ~250k is slow")
    ap.add_argument("--alpha", type=float, default=1e-3)
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--fresh-olga", default="bench/out/fresh_trb_olga.tsv",
                    help="independent TRB OLGA set (olga-generate_sequences --humanTRB, different seed)")
    args = ap.parse_args()

    # CHAIN-CONSISTENT: queries = sample4 (TRB OLGA); controls all TRB.
    olga = load_sample("sample4").sample(fraction=1.0, shuffle=True, seed=0)
    nq, nc = args.n_query, args.n_control
    if olga.height < nq + nc:
        nc = olga.height - nq
        print(f"sample4 only {olga.height} unique; using n_control={nc}")
    qlist = olga[:nq]["cdr3"].to_list()
    samerun_cdr3 = olga[nq:nq + nc]["cdr3"].to_list()           # sample4 (same OLGA run as queries)
    indep_cdr3 = load_fresh_olga(args.fresh_olga)[:nc]["cdr3"].to_list()  # fresh OLGA, independent seed

    vdj = db.load(_bench.source(), species=args.species).filter(pl.col("gene") == "TRB")
    ref = _bench.valid_cdr3(vdj).group_by("cdr3").agg(pl.col("epitope").first())
    ref_epi = ref["epitope"].to_list()
    tgt = Index.build(ref["cdr3"].to_list(), "aa")
    real_ctrl = background("TRB", size=nc)
    samerun_ctrl = Index.build(samerun_cdr3, "aa")
    indep_ctrl = Index.build(indep_cdr3, "aa")
    N, M_real, M_same, M_indep = len(tgt), len(real_ctrl), len(samerun_ctrl), len(indep_ctrl)
    print(f"target N={N}; controls M real={M_real} same-run={M_same} fresh={M_indep}; queries={len(qlist)} (TRB)")

    params = first_hit.scope()
    th_raw = cost_lists(tgt, qlist, params, "target")
    th = [[(c, ref_epi[r]) for c, r in t] for t in th_raw]
    cc_real = cost_lists(real_ctrl, qlist, params, "real-control")
    cc_same = cost_lists(samerun_ctrl, qlist, params, "samerun-control")
    cc_indep = cost_lists(indep_ctrl, qlist, params, "fresh-control")

    r_real = sig_rate(th, cc_real, N, M_real, args.alpha)
    r_same = sig_rate(th, cc_same, N, M_same, args.alpha)
    r_indep = sig_rate(th, cc_indep, N, M_indep, args.alpha)
    print(f"\nTRB OLGA (sample4) significant-call rate (p_enrichment < {args.alpha}) — all controls TRB:")
    print(f"  real control (human_trb_aa, selected):       {r_real*100:6.3f}%")
    print(f"  OLGA control, same run (sample4 half):       {r_same*100:6.3f}%")
    print(f"  OLGA control, INDEPENDENT run (fresh olga):  {r_indep*100:6.3f}%   <- fixed-seed caveat check")
    print(f"  alpha (nominal):                             {args.alpha*100:6.3f}%")


if __name__ == "__main__":
    main()
