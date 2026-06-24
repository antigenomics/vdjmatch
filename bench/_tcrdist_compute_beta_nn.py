#!/usr/bin/env python3
"""tcrdist3 BETA-chain 1-NN distance compute (runs in the `cmp-tcrdist` conda env).

Beta-chain NN-distance analogue of bench/_tcrdist_compute.py (which emits the samples-prediction
contract). Reads a reference TSV (cdr3, v, j) and a query TSV (cdr3, v, j) for the BETA chain, builds
TCRrep reps (chains=["beta"]), computes sparse rectangular TCRdist (query x reference) within
``--radius`` from ``trR.rw_beta``, and writes per-query nearest-neighbour distance (query_cdr3,
nn_dist). Exact-CDR3 self-matches are dropped (mirrors vdjmatch exclude_exact) unless --keep-exact.
Queries with no reference within radius get nn_dist = NaN. Invalid V/J genes are dropped (logged).
Driven by tcrdist_detect_pair.py. Kept separate from _tcrdist_compute.py so alpha/beta/paired share an
identical query universe and NaN-handling.
"""
import argparse
import os
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.sparse as sp
from tcrdist.repertoire import TCRrep

CPUS = max(1, (os.cpu_count() or 2) - 1)


def _rep(df, organism, valid_v, valid_j):
    """cols cdr3,v,j -> (TCRrep, n_dropped) with invalid V/J genes dropped (beta)."""
    d = df.rename(columns={"cdr3": "cdr3_b_aa", "v": "v_b_gene", "j": "j_b_gene"}).copy()
    d["v_b_gene"] = d["v_b_gene"].where(d["v_b_gene"].isin(valid_v))
    d["j_b_gene"] = d["j_b_gene"].where(d["j_b_gene"].isin(valid_j))
    n0 = len(d)
    d = d.dropna(subset=["v_b_gene", "j_b_gene", "cdr3_b_aa"])
    d["count"] = 1
    tr = TCRrep(cell_df=d, organism=organism, chains=["beta"], compute_distances=False)
    tr.cpus = CPUS
    return tr, n0 - len(d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--queries", required=True)
    ap.add_argument("--radius", type=int, default=90)
    ap.add_argument("--keep-exact", action="store_true",
                    help="keep exact-CDR3 self-hits (default: drop, mirroring vdjmatch exclude_exact)")
    ap.add_argument("--out", required=True, help="TSV: query_cdr3, nn_dist (NaN if none within radius)")
    ap.add_argument("--organism", default="human")
    args = ap.parse_args()

    ref = pd.read_csv(args.ref, sep="\t")
    qry = pd.read_csv(args.queries, sep="\t")
    from tcrdist.repertoire_db import RefGeneSet
    genes = set(RefGeneSet("alphabeta_gammadelta_db.tsv").all_genes[args.organism].keys())
    vv, jj = {g for g in genes if g.startswith("TRBV")}, {g for g in genes if g.startswith("TRBJ")}

    trR, dR = _rep(ref, args.organism, vv, jj)
    trQ, dQ = _rep(qry, args.organism, vv, jj)
    ref_cdr3 = trR.clone_df["cdr3_b_aa"].to_numpy()
    q_cdr3 = trQ.clone_df["cdr3_b_aa"].to_numpy()
    print(f"[beta] ref clones={len(ref_cdr3)} (dropped {dR}); query clones={len(q_cdr3)} "
          f"(dropped {dQ}); radius={args.radius}")

    trR.compute_sparse_rect_distances(df=trQ.clone_df, df2=trR.clone_df,
                                      radius=args.radius, chunk_size=100)
    M = trR.rw_beta
    M = M.tocsr() if sp.issparse(M) else sp.csr_matrix(M)

    rows = []
    for i in range(M.shape[0]):
        s, e = M.indptr[i], M.indptr[i + 1]
        cols, vals = M.indices[s:e], M.data[s:e]
        keep = vals != 0
        cols, dist = cols[keep], np.where(vals[keep] == -1, 0, vals[keep])
        qc = q_cdr3[i]
        if not args.keep_exact and len(cols):
            notself = ref_cdr3[cols] != qc
            cols, dist = cols[notself], dist[notself]
        nn = float(np.min(dist)) if len(dist) else float("nan")
        rows.append((qc, nn))

    pd.DataFrame(rows, columns=["query_cdr3", "nn_dist"]).to_csv(args.out, sep="\t", index=False)
    print(f"[beta] wrote {args.out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
