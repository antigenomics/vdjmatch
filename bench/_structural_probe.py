#!/usr/bin/env python3
"""GLOBAL candidate (ii) probe: swap the NED PSSM base matrix to the seqtree ``structural``
(Miyazawa-Jernigan interaction-strength) matrix for ALL epitopes, and measure the per-(epitope,chain)
detection ROC and the joint-failure-positive count vs the committed BLOSUM62 base. No per-epitope
tuning — the same global base swap is applied everywhere.

The swap is done by monkeypatching ``regions.significance_pssm`` to default ``base='structural'`` so
``benchmark._pssm_targets`` (and thus ``baseline_scores``) build the position matrix on the MJ matrix.

NUMBERS ONLY. Run from repo root with the project venv:
    ./.venv/bin/python bench/_structural_probe.py
"""
from __future__ import annotations

import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH))
from vdjmatch.match import regions
from metrics import roc_auc

TASKS = [("NLV", False), ("LLW", False), ("LLL", False), ("YLQ", True), ("GLC", True)]
_orig_pssm = regions.significance_pssm


def _patched(length, base="structural", scale=100):       # global base swap
    return _orig_pssm(length, base=base, scale=scale)


def per_cell_roc(base_name):
    """Return {(task,chain): (roc, pairs)} with NED built on ``base_name`` PSSM base.
    Reloads _feat_probe each call so its B.* caches re-score with the patched/unpatched matrix."""
    # Force a fresh import so cached PSSMs / scores are recomputed under the active patch.
    for m in ("_feat_probe", "benchmark"):
        sys.modules.pop(m, None)
    import benchmark as B  # noqa
    import _feat_probe as fp  # noqa
    out = {}
    for task, paired in TASKS:
        d = fp.task_table(task)
        lab = {c: int(l) for c, l in zip(d["cdr3"], d["label"])}
        allq = list(lab)
        base = fp.baseline_scores(task, "TRB")
        if not paired:
            sc = {c: base.get(c, 0) for c in allq}
            out[(task, "TRB")] = ([(lab[c], sc[c]) for c in allq])
        else:
            ac = {c: a for c, a in zip(d["cdr3"], d["a_cdr3"])}
            ba = fp.baseline_scores(task, "TRA")
            sTRA = {c: ba.get(ac[c], 0) for c in allq}
            sTRB = {c: base.get(c, 0) for c in allq}
            out[(task, "TRA")] = [(lab[c], sTRA[c]) for c in allq]
            out[(task, "TRB")] = [(lab[c], sTRB[c]) for c in allq]
    return out


def main():
    # 1) committed BLOSUM62 NED (no patch)
    regions.significance_pssm = _orig_pssm
    blo = per_cell_roc("blosum62")
    # 2) structural (MJ) NED (patched)
    regions.significance_pssm = _patched
    struc = per_cell_roc("structural")
    regions.significance_pssm = _orig_pssm

    print("=== GLOBAL candidate (ii): NED PSSM base BLOSUM62 -> structural(MJ), per chain ===")
    rocs_b, rocs_s, worse = [], [], []
    for k in sorted(blo):
        rb = roc_auc(blo[k]); rs = roc_auc(struc[k])
        rocs_b.append(rb); rocs_s.append(rs)
        d = rs - rb
        flag = "  <-- REGRESS" if d < -1e-9 else ""
        if d < -1e-9:
            worse.append((k[0], k[1], round(d, 4)))
        print(f"  {k[0]:4s} {k[1]:6s} ROC {rb:.4f} -> {rs:.4f} (d={d:+.4f}){flag}")
    mb = sum(rocs_b) / len(rocs_b); ms = sum(rocs_s) / len(rocs_s)
    print(f"  mean ROC {mb:.4f} -> {ms:.4f} (d={ms-mb:+.4f}); min ROC {min(rocs_b):.4f} -> {min(rocs_s):.4f}")
    print(f"  cells regressing: {len(worse)} {worse}")


if __name__ == "__main__":
    main()
