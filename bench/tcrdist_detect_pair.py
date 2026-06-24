#!/usr/bin/env python3
"""tcrdist3 ALPHA / BETA / PAIRED detection for the paired tasks (GLC, YLQ).

Single-chain 1-NN tcrdist of each sample6 query to the epitope's A*02 VDJdb2026 reference, per chain:
  alpha query = (a_cdr3, a_v, a_j) vs ref_table(task,'TRA');  beta query = (cdr3, v, j) vs (...,'TRB').
Score sign: SMALLER tcrdist = more likely positive, so score = -nn_dist. Paired = -(dist_a + dist_b).
Exact-CDR3 self-matches dropped (mirrors vdjmatch exclude_exact). Genes unknown to tcrdist3 -> clone
dropped (counted). ROC via vdjmatch metrics.roc_auc (ties grouped). Temp ref/query TSVs live under /tmp
ONLY (never a git repo). Run from repo root with the project venv:

    ./.venv/bin/python bench/tcrdist_detect_pair.py
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _feat_probe import ref_table, task_table
from metrics import pr_auc_balanced, roc_auc

ENV = "cmp-tcrdist"
BENCH = Path(__file__).resolve().parent
ALPHA_C = BENCH / "_tcrdist_compute_alpha.py"
BETA_C = BENCH / "_tcrdist_compute_beta_nn.py"
RADIUS = 90
TMP = Path(tempfile.mkdtemp(prefix="tcrdist_pair_", dir="/tmp"))


def norm_gene(col: pl.Expr) -> pl.Expr:
    """IMGT-ify: strip, default allele *01 when missing (tcrdist3 wants e.g. TRAV12-1*01)."""
    g = col.cast(pl.Utf8).str.strip_chars()
    return pl.when(g.str.contains(r"\*")).then(g).otherwise(g + "*01")


def _write_ref(task: str, locus: str) -> Path:
    r = (ref_table(task, locus).select("cdr3", v=norm_gene(pl.col("v")), j=norm_gene(pl.col("j")))
         .sort(["cdr3", "v", "j"]).unique("cdr3", keep="first", maintain_order=True))
    p = TMP / f"{task}_{locus}_ref.tsv"
    r.write_csv(p, separator="\t")
    return p


def _write_query(d: pl.DataFrame, cdr3, v, j, name) -> Path:
    q = (d.select(cdr3=pl.col(cdr3), v=norm_gene(pl.col(v)), j=norm_gene(pl.col(j)))
         .sort(["cdr3", "v", "j"]).unique("cdr3", keep="first", maintain_order=True))
    p = TMP / name
    q.write_csv(p, separator="\t")
    return p


def _nn(compute: Path, ref: Path, qry: Path, out: Path) -> dict:
    """Run a chain compute in the conda env; return {query_cdr3 -> nn_dist (NaN if no hit)}."""
    subprocess.run(["conda", "run", "-n", ENV, "python", str(compute),
                    "--ref", str(ref), "--queries", str(qry), "--radius", str(RADIUS),
                    "--out", str(out)], check=True)
    t = pl.read_csv(out, separator="\t")
    return {c: dd for c, dd in zip(t["query_cdr3"], t["nn_dist"])}


def run_task(task: str) -> list[dict]:
    d = task_table(task)                                          # one row per beta-cdr3 query
    labels = {c: int(l) for c, l in zip(d["cdr3"], d["label"])}

    # ALPHA: query (a_cdr3,a_v,a_j) vs TRA ref
    refA = _write_ref(task, "TRA")
    qA = _write_query(d, "a_cdr3", "a_v", "a_j", f"{task}_qa.tsv")
    nnA = _nn(ALPHA_C, refA, qA, TMP / f"{task}_alpha_nn.tsv")    # keyed by a_cdr3
    # BETA: query (cdr3,v,j) vs TRB ref
    refB = _write_ref(task, "TRB")
    qB = _write_query(d, "cdr3", "v", "j", f"{task}_qb.tsv")
    nnB = _nn(BETA_C, refB, qB, TMP / f"{task}_beta_nn.tsv")      # keyed by beta cdr3

    # Map alpha NN (keyed by a_cdr3) back to per-query rows (keyed by beta cdr3).
    a_for = {c: a for c, a in zip(d["cdr3"], d["a_cdr3"])}

    # Per-query distances. A query is in a chain's ROC iff its clone survived gene-validation
    # (i.e. appears in that chain's NN output). NaN nn_dist = no neighbour within radius (kept,
    # ranked worst). Missing key = clone dropped (unknown gene) -> excluded from that chain.
    pairs_a, pairs_b, pairs_p = [], [], []
    drop_a = drop_b = 0
    for c, lab in labels.items():
        ac = a_for[c]
        da = nnA.get(ac, None)                                    # alpha NN dist for this query
        db = nnB.get(c, None)
        if da is None:
            drop_a += 1
        else:
            sa = -(RADIUS + 1) if (da is None or np.isnan(da)) else -float(da)
            pairs_a.append((lab, sa))
        if db is None:
            drop_b += 1
        else:
            sb = -(RADIUS + 1) if (db is None or np.isnan(db)) else -float(db)
            pairs_b.append((lab, sb))
        # paired: both chains' clones must survive gene-validation
        if da is not None and db is not None:
            dav = RADIUS + 1 if np.isnan(da) else float(da)
            dbv = RADIUS + 1 if np.isnan(db) else float(db)
            pairs_p.append((lab, -(dav + dbv)))

    out = []
    fd = Path.home() / "vcs/manuscripts/2026-vdjmatch/figures/data"            # Fig 5 PR curves
    for chain, pairs, dropped in (("alpha", pairs_a, drop_a),
                                  ("beta", pairs_b, drop_b),
                                  ("paired", pairs_p, None)):
        npos = sum(l for l, _ in pairs)
        nneg = len(pairs) - npos
        out.append(dict(task=task, chain=chain,
                        roc_auc=round(roc_auc(pairs), 4),
                        pr_auc=round(pr_auc_balanced(pairs), 4),
                        n_pos=npos, n_neg=nneg, dropped=dropped))
        if fd.is_dir():                                                        # dump balanced PR curve
            P = npos or 1; N = nneg or 1; xs, ys = [0.0], [1.0]
            for thr in sorted({s for _, s in pairs}, reverse=True):
                tp = sum(1 for l, s in pairs if l and s >= thr); fp = sum(1 for l, s in pairs if not l and s >= thr)
                tpr = tp / P; fpr = fp / N
                xs.append(tpr); ys.append(tpr / (tpr + fpr) if tpr + fpr else 1.0)
            cn = {"alpha": "TRA", "beta": "TRB", "paired": "paired"}[chain]
            (fd / f"pr5_{task}_{cn}_tcrdist.dat").write_text(
                "\n".join(f"{x:.4f} {y:.4f}" for x, y in zip(xs, ys)) + "\n")
    print(f"== {task} ==  alpha_dropped={drop_a} beta_dropped={drop_b} "
          f"(paired n={len(pairs_p)})")
    for r in out:
        print(f"  {r['chain']:6s} ROC={r['roc_auc']:.4f} PR={r['pr_auc']:.4f} "
              f"pos={r['n_pos']} neg={r['n_neg']} dropped={r['dropped']}")
    # raw per-query (label, score) pairs for downstream bootstrap CIs (chain names TRA/TRB/paired)
    raw = {"TRA": pairs_a, "TRB": pairs_b, "paired": pairs_p}
    return out, raw


RESULTS = Path.home() / "vcs/manuscripts/2026-vdjmatch/benchmarks/results"


def main():
    allr = []
    rows = ["task\tchain\tlabel\tscore"]
    for task in ("YLQ", "GLC"):
        out, raw = run_task(task)
        allr += out
        for chain in ("TRA", "TRB", "paired"):
            for lab, sc in raw[chain]:
                rows.append(f"{task}\t{chain}\t{int(lab)}\t{sc:.6g}")
    if RESULTS.is_dir():
        p = RESULTS / "tcrdist_bychain_pairs.tsv"
        p.write_text("\n".join(rows) + "\n")
        print(f"\nwrote {p}  ({len(rows) - 1} pairs)")
    print("\nTMP dir:", TMP)
    import json
    print(json.dumps(allr, indent=2))


if __name__ == "__main__":
    main()
