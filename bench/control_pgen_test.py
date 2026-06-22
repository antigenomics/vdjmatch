#!/usr/bin/env python3
"""Is the OLGA spurious rate control-mismatch, or genuine coincidental overlap?

Hypothesis (from bench/olga_overlap_limit.py): the OLGA significant-call rate (~9.3%) >> alpha because
the bundled control `human_trb_aa` is a *real, thymically-selected* repertoire, whereas the OLGA query
draws are raw Pgen (no selection). The enrichment test then sees OLGA queries hit the (Pgen-biased,
public-TCR-rich) epitope target more than the selected control predicts.

Direct test: split sample5 (OLGA) into disjoint query / control halves. The OLGA control is i.i.d.
with the OLGA queries, so a perfectly matched null -> enrichment should vanish -> significant rate
-> alpha. We compare the significant-call rate under (A) the real bundled control vs (B) the OLGA
control, at matched control size.

    python bench/control_pgen_test.py --n-query 20000 --n-control 200000
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench
from compare import load_sample
from compare import TESTDATA
from vdjmatch import db
from vdjmatch.evalue import background, first_hit
from seqtree import Index


def load_olga(name: str) -> pl.DataFrame:
    """An OLGA AIRR file -> valid unique CDR3 (junction_aa)."""
    d = pl.read_csv(TESTDATA / f"{name}_olga_airr.txt", separator="\t").select(cdr3="junction_aa")
    return _bench.valid_cdr3(d).unique("cdr3")


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
    args = ap.parse_args()

    olga = load_sample("sample5").sample(fraction=1.0, shuffle=True, seed=0)
    nq, nc = args.n_query, args.n_control
    if olga.height < nq + nc:
        nc = olga.height - nq
        print(f"sample5 only {olga.height} unique; using n_control={nc}")
    qlist = olga[:nq]["cdr3"].to_list()
    olga_ctrl_cdr3 = olga[nq:nq + nc]["cdr3"].to_list()

    indep_cdr3 = load_olga("sample4").sample(fraction=1.0, shuffle=True, seed=0)[:nc]["cdr3"].to_list()

    vdj = db.load(_bench.source(), species=args.species).filter(pl.col("gene") == "TRB")
    ref = _bench.valid_cdr3(vdj).group_by("cdr3").agg(pl.col("epitope").first())
    ref_epi = ref["epitope"].to_list()
    tgt = Index.build(ref["cdr3"].to_list(), "aa")
    real_ctrl = background("TRB", size=nc)
    olga_ctrl = Index.build(olga_ctrl_cdr3, "aa")               # sample5 (same OLGA run as queries)
    indep_ctrl = Index.build(indep_cdr3, "aa")                  # sample4 (independent OLGA run)
    N, M_real, M_olga, M_indep = len(tgt), len(real_ctrl), len(olga_ctrl), len(indep_ctrl)
    print(f"target N={N}; controls M real={M_real} sample5={M_olga} sample4={M_indep}; queries={len(qlist)}")

    params = first_hit.scope()
    th_raw = cost_lists(tgt, qlist, params, "target")
    th = [[(c, ref_epi[r]) for c, r in t] for t in th_raw]
    cc_real = cost_lists(real_ctrl, qlist, params, "real-control")
    cc_olga = cost_lists(olga_ctrl, qlist, params, "olga-control-sample5")
    cc_indep = cost_lists(indep_ctrl, qlist, params, "olga-control-sample4")

    r_real = sig_rate(th, cc_real, N, M_real, args.alpha)
    r_olga = sig_rate(th, cc_olga, N, M_olga, args.alpha)
    r_indep = sig_rate(th, cc_indep, N, M_indep, args.alpha)
    print(f"\nOLGA (sample5) significant-call rate (p_enrichment < {args.alpha}):")
    print(f"  real control (human_trb_aa, selected):      {r_real*100:6.3f}%")
    print(f"  OLGA control, same run (sample5 half):      {r_olga*100:6.3f}%")
    print(f"  OLGA control, INDEPENDENT run (sample4):     {r_indep*100:6.3f}%   <- fixed-seed caveat check")
    print(f"  alpha (nominal):                            {args.alpha*100:6.3f}%")


if __name__ == "__main__":
    main()
