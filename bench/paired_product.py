#!/usr/bin/env python3
"""Paired-by-product detection predictions for methods that score each chain independently.

For a paired clonotype (alpha CDR3, beta CDR3), combine the two single-chain scores into one paired
score: paired_score = score_TRA * score_TRB. This is the natural independence combination for methods
that emit a per-chain confidence (imw-DETECT, ERGO-II) but have no native paired model. The paired
clonotype set comes from _feat_probe.task_table(task, 'TRB'), which carries both cdr3 (beta) and a_cdr3
(alpha) plus the validated label; this is how alpha<->beta are aligned.

Emits predictions/{imw-detect,ergo}/{YLQ,GLC}_paired.tsv, keyed by beta CDR3 (query_id), with
significant = 1 iff BOTH single-chain calls were significant. A chain with no single-chain score (e.g.
DETECT predicted a non-target epitope) contributes score 0 -> paired score 0, matching the detection
reader's no-call convention. YLQ/GLC only (sample6 is the sole paired dataset).

    python bench/paired_product.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _feat_probe import task_table
from metrics import roc_auc

EPI = {"YLQ": "YLQPRTFLL", "GLC": "GLCTLVAML"}
PRED = Path("bench/predictions")
METHODS = ("imw-detect", "ergo")


def load_chain(method: str, task: str, locus: str) -> dict[str, tuple[float, int]]:
    """Read predictions/<method>/<task>_<locus>.tsv -> {query_id: (score, significant)}."""
    p = PRED / method / f"{task}_{locus}.tsv"
    if not p.exists():
        return {}
    t = pl.read_csv(p, separator="\t")
    return {q: (float(s), int(sig)) for q, s, sig in zip(t["query_id"], t["score"], t["significant"])}


def main():
    for method in METHODS:
        for task, pep in EPI.items():
            a_sc = load_chain(method, task, "TRA")                 # keyed by alpha CDR3
            b_sc = load_chain(method, task, "TRB")                 # keyed by beta CDR3
            if not a_sc or not b_sc:
                print(f"{method} {task}/paired: SKIP (missing "
                      f"{'TRA' if not a_sc else 'TRB'} single-chain file)")
                continue
            d = task_table(task, "TRB")                            # cdr3=beta, a_cdr3=alpha, label
            rows, pairs = [], []
            for beta, alpha, lab in zip(d["cdr3"], d["a_cdr3"], d["label"]):
                sa, siga = a_sc.get(alpha, (0.0, 0))
                sb, sigb = b_sc.get(beta, (0.0, 0))
                ps = sa * sb
                rows.append((beta, pep, ps, int(bool(siga and sigb))))
                pairs.append((int(lab), ps))
            pl.DataFrame(rows, schema=["query_id", "epitope", "score", "significant"],
                         orient="row").write_csv(PRED / method / f"{task}_paired.tsv", separator="\t")
            auc = roc_auc(pairs) if pairs else float("nan")
            npos = sum(p for p, _ in pairs)
            print(f"{method} {task}/paired: {len(rows)} scored | "
                  f"ROC-AUC {auc:.3f} (n_pos={npos} n={len(pairs)})")


if __name__ == "__main__":
    main()
