#!/usr/bin/env python3
"""tcrdist3 WITHIN-SET all-pairs sparse distance for single-linkage clustering (runs in `cmp-tcrdist`).

Reads one clonotype TSV (cdr3, v, j) for a single chain, builds a TCRrep, computes the SQUARE sparse
self-distance (set x set) within ``--radius`` via compute_sparse_rect_distances(df=clone_df,
df2=clone_df), and writes an edge list TSV (i, j, dist) for every off-diagonal pair with tcrdist <=
radius. Indices i/j are 0-based ROW indices into the *input* TSV (NOT tcrdist's internal clone_df order),
so the caller can align labels 1:1 to the shared clonotype list. Rows whose V or J gene is unknown to
tcrdist3 are dropped from the distance computation (they simply produce no edges -> singletons); their
input-row indices never appear in the edge list. Numbers-only output. Driven by tcrdist_cluster_perepi.py.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="TSV: cdr3, v, j (one chain)")
    ap.add_argument("--chain", choices=["alpha", "beta"], required=True)
    ap.add_argument("--radius", type=int, required=True)
    ap.add_argument("--out", required=True, help="edge list TSV: i, j, dist (i<j, dist<=radius)")
    ap.add_argument("--organism", default="human")
    args = ap.parse_args()

    pre = "a" if args.chain == "alpha" else "b"
    df = pd.read_csv(args.inp, sep="\t").reset_index(drop=True)
    df["_row"] = np.arange(len(df))                              # original input-row index (the alignment key)

    from tcrdist.repertoire_db import RefGeneSet
    genes = set(RefGeneSet("alphabeta_gammadelta_db.tsv").all_genes[args.organism].keys())
    vp, jp = ("TRAV", "TRAJ") if args.chain == "alpha" else ("TRBV", "TRBJ")
    vv = {g for g in genes if g.startswith(vp)}
    jj = {g for g in genes if g.startswith(jp)}

    d = df.rename(columns={"cdr3": f"cdr3_{pre}_aa", "v": f"v_{pre}_gene", "j": f"j_{pre}_gene"}).copy()
    d[f"v_{pre}_gene"] = d[f"v_{pre}_gene"].where(d[f"v_{pre}_gene"].isin(vv))
    d[f"j_{pre}_gene"] = d[f"j_{pre}_gene"].where(d[f"j_{pre}_gene"].isin(jj))
    n0 = len(d)
    d = d.dropna(subset=[f"v_{pre}_gene", f"j_{pre}_gene", f"cdr3_{pre}_aa"])
    d["count"] = 1
    print(f"[{args.chain}] input rows={n0}; kept={len(d)} (dropped {n0 - len(d)} unknown-gene); "
          f"radius={args.radius}")

    # TCRrep collapses identical clones into clone_df ONLY when the cell_df has no per-row distinguishing
    # column — so we must NOT hand it _row (that would force one clone per input row, i.e. the full
    # 9214x9214 distance for the paired set instead of ~3050x3050 on unique alpha clones). Drop _row for
    # the rep; recover the clone->input-rows mapping ourselves by grouping the kept rows on (cdr3,v,j).
    keycols = [f"cdr3_{pre}_aa", f"v_{pre}_gene", f"j_{pre}_gene"]
    rows_by_key = d.groupby(keycols)["_row"].apply(list).to_dict()
    tr = TCRrep(cell_df=d.drop(columns=["_row"]), organism=args.organism, chains=[args.chain],
                compute_distances=False)
    tr.cpus = CPUS
    clone = tr.clone_df.reset_index(drop=True)
    clone_to_rows = [rows_by_key[tuple(r)] for r in clone[keycols].itertuples(index=False, name=None)]

    tr.compute_sparse_rect_distances(df=clone, df2=clone, radius=args.radius, chunk_size=100)
    M = tr.rw_alpha if args.chain == "alpha" else tr.rw_beta
    M = M.tocsr() if sp.issparse(M) else sp.csr_matrix(M)

    edges = []
    for i in range(M.shape[0]):
        s, e = M.indptr[i], M.indptr[i + 1]
        cols, vals = M.indices[s:e], M.data[s:e]
        keep = vals != 0                                        # 0 = beyond radius (absent in csr)
        cols, vals = cols[keep], vals[keep]
        dist = np.where(vals == -1, 0, vals)                    # -1 sentinel = true distance 0
        for cj, dd in zip(cols.tolist(), dist.tolist()):
            if cj <= i:                                         # off-diagonal, emit each clone pair once
                continue
            if dd > args.radius:
                continue
            # expand clone-pair to all input-row pairs it represents (dedup collapse, usually 1x1)
            for ri in clone_to_rows[i]:
                for rj in clone_to_rows[cj]:
                    a, b = (ri, rj) if ri < rj else (rj, ri)
                    if a != b:
                        edges.append((a, b, int(dd)))

    pd.DataFrame(edges, columns=["i", "j", "dist"]).to_csv(args.out, sep="\t", index=False)
    print(f"[{args.chain}] wrote {args.out} ({len(edges)} edges)")


if __name__ == "__main__":
    main()
