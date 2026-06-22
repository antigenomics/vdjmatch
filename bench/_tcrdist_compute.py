#!/usr/bin/env python3
"""tcrdist3 label-transfer compute step (runs in the `cmp-tcrdist` conda env).

Reads a reference TSV (cdr3, v, j, epitope) and a query TSV (cdr3, v, j), builds TCRrep reps, computes
sparse rectangular TCRdist (query x reference) within ``--radius``, and writes predictions in the
compare.py samples contract (query_id, epitope, score, significant). Two arms: 1-NN (score = radius -
nearest same-epitope dist) and k-NN (score = #epitope votes among the k nearest). Invalid V/J genes are
dropped (logged). Driven by bench/tcrdist_samples.py — not run directly.
"""
import argparse
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scipy.sparse as sp
from tcrdist.repertoire import TCRrep


def _rep(df, organism, valid_v, valid_j):
    """cols cdr3,v,j[,epitope] -> (TCRrep, kept clone_df) with invalid V/J genes dropped."""
    d = df.rename(columns={"cdr3": "cdr3_b_aa", "v": "v_b_gene", "j": "j_b_gene"}).copy()
    d["v_b_gene"] = d["v_b_gene"].where(d["v_b_gene"].isin(valid_v))
    d["j_b_gene"] = d["j_b_gene"].where(d["j_b_gene"].isin(valid_j))
    n0 = len(d)
    d = d.dropna(subset=["v_b_gene", "j_b_gene", "cdr3_b_aa"])
    d["count"] = 1
    tr = TCRrep(cell_df=d, organism=organism, chains=["beta"], compute_distances=False)
    tr.cpus = 1
    return tr, n0 - len(d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--queries", required=True)
    ap.add_argument("--targets", required=True, help="comma-separated target epitopes")
    ap.add_argument("--radius", type=int, default=90)
    ap.add_argument("--knn", type=int, default=5)
    ap.add_argument("--sig-radius", type=int, default=24, help="nearest-dist cutoff for 'significant'")
    ap.add_argument("--keep-exact", action="store_true",
                    help="keep exact-CDR3 self-hits (default: drop, mirroring vdjmatch exclude_exact)")
    ap.add_argument("--out-1nn", required=True)
    ap.add_argument("--out-knn", required=True)
    ap.add_argument("--organism", default="human")
    args = ap.parse_args()
    targets = set(args.targets.split(","))

    ref = pd.read_csv(args.ref, sep="\t")
    qry = pd.read_csv(args.queries, sep="\t")
    # valid genes = those tcrdist3 knows (intersect with its all_genes for this organism)
    from tcrdist.repertoire_db import RefGeneSet
    genes = set(RefGeneSet("alphabeta_gammadelta_db.tsv").all_genes[args.organism].keys())
    vv, jj = {g for g in genes if g.startswith("TRBV")}, {g for g in genes if g.startswith("TRBJ")}

    trR, dR = _rep(ref, args.organism, vv, jj)
    trQ, dQ = _rep(qry, args.organism, vv, jj)
    ref_epi = trR.clone_df["epitope"].to_numpy()
    ref_cdr3 = trR.clone_df["cdr3_b_aa"].to_numpy()
    q_cdr3 = trQ.clone_df["cdr3_b_aa"].to_numpy()
    print(f"ref clones={len(ref_cdr3)} (dropped {dR}); query clones={len(q_cdr3)} (dropped {dQ}); "
          f"radius={args.radius}")

    trR.compute_sparse_rect_distances(df=trQ.clone_df, df2=trR.clone_df,
                                      radius=args.radius, chunk_size=100)
    M = trR.rw_beta
    M = M.tocsr() if sp.issparse(M) else sp.csr_matrix(M)

    rows_1nn, rows_knn = [], []
    for i in range(M.shape[0]):
        s, e = M.indptr[i], M.indptr[i + 1]
        cols, vals = M.indices[s:e], M.data[s:e]
        # decode: 0 = beyond radius (absent in csr anyway); -1 = true distance 0; v>0 = distance v
        keep = vals != 0
        cols, dist = cols[keep], np.where(vals[keep] == -1, 0, vals[keep])
        qc = q_cdr3[i]
        if not args.keep_exact and len(cols):                  # drop exact-CDR3 self-hits (leakage)
            notself = ref_cdr3[cols] != qc
            cols, dist = cols[notself], dist[notself]
        if len(cols) == 0:
            continue
        order = np.argsort(dist)
        cols, dist = cols[order], dist[order]
        epis = ref_epi[cols]
        sig = int(dist[0] <= args.sig_radius)
        # 1-NN per target: score = radius - nearest same-epitope distance
        for e_t in targets:
            m = epis == e_t
            if m.any():
                rows_1nn.append((qc, e_t, float(args.radius - dist[m][0]), sig))
        # k-NN: vote among k nearest (any epitope); emit votes for target epitopes
        kk = epis[:args.knn]
        for e_t in targets:
            v = int((kk == e_t).sum())
            if v:
                rows_knn.append((qc, e_t, float(v), sig))
        if sig and not any(r[0] == qc for r in rows_1nn[-len(targets):] if r[3]):
            rows_1nn.append((qc, "_OTHER_", 0.0, 1))            # carry significance even off-target
            rows_knn.append((qc, "_OTHER_", 0.0, 1))

    for path, rows in ((args.out_1nn, rows_1nn), (args.out_knn, rows_knn)):
        pd.DataFrame(rows, columns=["query_id", "epitope", "score", "significant"]).to_csv(
            path, sep="\t", index=False)
        print(f"wrote {path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
