#!/usr/bin/env python3
"""tcrdist3 ALPHA-chain 1-NN distance compute (runs in the `cmp-tcrdist` conda env).

Alpha-chain analogue of bench/_tcrdist_compute.py: reads a reference TSV (cdr3, v, j) and a query TSV
(cdr3, v, j) for the ALPHA chain, builds TCRrep reps (chains=["alpha"], renaming to cdr3_a_aa /
v_a_gene / j_a_gene), computes sparse rectangular TCRdist (query x reference) within ``--radius`` from
``trR.rw_alpha``, and writes per-query nearest-neighbour distance (query_cdr3, nn_dist). Exact-CDR3
self-matches are dropped (mirrors vdjmatch exclude_exact) unless --keep-exact. Queries with no reference
within radius get nn_dist = NaN. Invalid V/J genes are dropped (logged). Driven by tcrdist_detect_pair.py.
"""
import argparse
import os
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.sparse as sp
from tcrdist.repertoire import TCRrep

CPUS = max(1, (os.cpu_count() or 2) - 1)                          # parallel parasail; leave one core


def _rep(df, organism, valid_v, valid_j):
    """cols cdr3,v,j -> (TCRrep, kept clone_df, n_dropped) with invalid V/J genes dropped (alpha)."""
    d = df.rename(columns={"cdr3": "cdr3_a_aa", "v": "v_a_gene", "j": "j_a_gene"}).copy()
    d["v_a_gene"] = d["v_a_gene"].where(d["v_a_gene"].isin(valid_v))
    d["j_a_gene"] = d["j_a_gene"].where(d["j_a_gene"].isin(valid_j))
    n0 = len(d)
    d = d.dropna(subset=["v_a_gene", "j_a_gene", "cdr3_a_aa"])
    d["count"] = 1
    tr = TCRrep(cell_df=d, organism=organism, chains=["alpha"], compute_distances=False)
    tr.cpus = CPUS
    return tr, d, n0 - len(d)


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
    vv, jj = {g for g in genes if g.startswith("TRAV")}, {g for g in genes if g.startswith("TRAJ")}

    trR, _, dR = _rep(ref, args.organism, vv, jj)
    trQ, _, dQ = _rep(qry, args.organism, vv, jj)
    ref_cdr3 = trR.clone_df["cdr3_a_aa"].to_numpy()
    q_cdr3 = trQ.clone_df["cdr3_a_aa"].to_numpy()
    print(f"[alpha] ref clones={len(ref_cdr3)} (dropped {dR}); query clones={len(q_cdr3)} "
          f"(dropped {dQ}); radius={args.radius}")

    trR.compute_sparse_rect_distances(df=trQ.clone_df, df2=trR.clone_df,
                                      radius=args.radius, chunk_size=100)
    M = trR.rw_alpha
    M = M.tocsr() if sp.issparse(M) else sp.csr_matrix(M)

    rows = []
    for i in range(M.shape[0]):
        s, e = M.indptr[i], M.indptr[i + 1]
        cols, vals = M.indices[s:e], M.data[s:e]
        # decode sparse: 0 = beyond radius (absent); -1 = true distance 0; v>0 = distance v
        keep = vals != 0
        cols, dist = cols[keep], np.where(vals[keep] == -1, 0, vals[keep])
        qc = q_cdr3[i]
        if not args.keep_exact and len(cols):                    # drop exact-CDR3 self-hits (leakage)
            notself = ref_cdr3[cols] != qc
            cols, dist = cols[notself], dist[notself]
        nn = float(np.min(dist)) if len(dist) else float("nan")
        rows.append((qc, nn))

    pd.DataFrame(rows, columns=["query_cdr3", "nn_dist"]).to_csv(args.out, sep="\t", index=False)
    print(f"[alpha] wrote {args.out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
