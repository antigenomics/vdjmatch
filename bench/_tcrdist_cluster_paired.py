#!/usr/bin/env python3
"""tcrdist3 NATIVE paired within-set distance for single-linkage clustering (runs in `cmp-tcrdist`).

Reads a paired clonotype TSV (cdr3_a, v_a, j_a, cdr3_b, v_b, j_b), builds a PAIRED TCRrep
(chains=['alpha','beta']), computes the square sparse self-distance, and emits an edge list (i, j, dist)
for off-diagonal pairs whose JOINT tcrdist (rw_alpha + rw_beta) <= ``--radius``. tcrdist3 stores the two
chains separately and the paired distance is their sum; ``compute_sparse_rect_distances`` thresholds each
chain at ``--radius`` (per-chain), so to capture every pair with SUM <= R we run the compute at radius R
on EACH chain (a pair with sum <= R has each chain <= R) and then keep those whose summed distance <= R.
Indices i/j are 0-based input-row indices (dedup-collapsed clones expanded back to all their rows).
Numbers-only output. Driven by tcrdist_cluster_perepi.py.
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
KEY = ["cdr3_a_aa", "v_a_gene", "j_a_gene", "cdr3_b_aa", "v_b_gene", "j_b_gene"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True,
                    help="TSV: cdr3_a, v_a, j_a, cdr3_b, v_b, j_b")
    ap.add_argument("--radius", type=int, required=True, help="JOINT (alpha+beta) tcrdist cutoff")
    ap.add_argument("--out", required=True, help="edge list TSV: i, j, dist (i<j, joint dist<=radius)")
    ap.add_argument("--organism", default="human")
    args = ap.parse_args()

    df = pd.read_csv(args.inp, sep="\t").reset_index(drop=True)
    df = df.rename(columns={"cdr3_a": "cdr3_a_aa", "v_a": "v_a_gene", "j_a": "j_a_gene",
                            "cdr3_b": "cdr3_b_aa", "v_b": "v_b_gene", "j_b": "j_b_gene"})
    df["_row"] = np.arange(len(df))

    from tcrdist.repertoire_db import RefGeneSet
    genes = set(RefGeneSet("alphabeta_gammadelta_db.tsv").all_genes[args.organism].keys())
    for col, pref in (("v_a_gene", "TRAV"), ("j_a_gene", "TRAJ"),
                      ("v_b_gene", "TRBV"), ("j_b_gene", "TRBJ")):
        df[col] = df[col].where(df[col].isin({g for g in genes if g.startswith(pref)}))
    n0 = len(df)
    df = df.dropna(subset=KEY)
    df["count"] = 1
    print(f"[paired] input rows={n0}; kept={len(df)} (dropped {n0 - len(df)} unknown-gene); "
          f"joint radius={args.radius}")

    rows_by_key = df.groupby(KEY)["_row"].apply(list).to_dict()
    tr = TCRrep(cell_df=df.drop(columns=["_row"]), organism=args.organism, chains=["alpha", "beta"],
                compute_distances=False)
    tr.cpus = CPUS
    clone = tr.clone_df.reset_index(drop=True)
    clone_to_rows = [rows_by_key[tuple(r)] for r in clone[KEY].itertuples(index=False, name=None)]

    # per-chain threshold = radius (necessary condition for sum<=radius); we then filter on the SUM.
    tr.compute_sparse_rect_distances(df=clone, df2=clone, radius=args.radius, chunk_size=100)
    Ma = tr.rw_alpha.tocsr() if sp.issparse(tr.rw_alpha) else sp.csr_matrix(tr.rw_alpha)
    Mb = tr.rw_beta.tocsr() if sp.issparse(tr.rw_beta) else sp.csr_matrix(tr.rw_beta)

    def decode_row(M, i):
        s, e = M.indptr[i], M.indptr[i + 1]
        cols, vals = M.indices[s:e], M.data[s:e]
        keep = vals != 0
        cols = cols[keep]
        dist = np.where(vals[keep] == -1, 0, vals[keep])
        return dict(zip(cols.tolist(), dist.tolist()))

    edges = []
    nclone = Ma.shape[0]
    for i in range(nclone):
        ra = decode_row(Ma, i)
        rb = decode_row(Mb, i)
        for cj in (ra.keys() & rb.keys()):                     # within radius on BOTH chains
            if cj <= i:
                continue
            joint = ra[cj] + rb[cj]
            if joint > args.radius:
                continue
            for x in clone_to_rows[i]:
                for y in clone_to_rows[cj]:
                    a, b = (x, y) if x < y else (y, x)
                    if a != b:
                        edges.append((a, b, int(joint)))

    pd.DataFrame(edges, columns=["i", "j", "dist"]).to_csv(args.out, sep="\t", index=False)
    print(f"[paired] wrote {args.out} ({len(edges)} edges)")


if __name__ == "__main__":
    main()
