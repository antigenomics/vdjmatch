#!/usr/bin/env python3
"""TCRMatch (beta-only) -> TASK-NAMED detection prediction files {NLV,LLW,LLL,YLQ,GLC}_TRB.tsv.

bench/tcrmatch_bench.py emits condition-coded files (C1/C2/C3); the detection reader expects per-task
files keyed by query CDR3beta. Here we score each task's TRB query set (from _feat_probe.task_table)
against that epitope's VDJdb2026-A*02 TRB reference (from _feat_probe.ref_table) with the TCRMatch
BLOSUM62 kernel, keeping the max kernel score per query (exact self-hits dropped, mirroring
exclude_exact). significant = score >= 0.97 (TCRMatch's own near-identity call). TCRMatch is beta-only,
so only TRB files are produced. Data read from the manuscript test_data at runtime (never copied).

    python bench/tcrmatch_detect.py
"""
from __future__ import annotations

import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _feat_probe import ref_table, task_table
from metrics import roc_auc
from tcrmatch_samples import MAXLEN, run_tcrmatch, trim

EPI = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV",
       "YLQ": "YLQPRTFLL", "GLC": "GLCTLVAML"}
SIG = 0.97                                                          # TCRMatch near-identity call threshold


def write_ref(task: str, path: Path) -> int:
    """Epitope's A*02 VDJdb2026 TRB reference -> 6-col TCRMatch DB (col0 trimmed key, col3 epitope)."""
    ref = ref_table(task, "TRB")
    kept = 0
    with open(path, "w") as f:
        f.write("trimmed_seq\toriginal_seq\treceptor_group\tepitopes\tsource_organisms\tsource_antigens\n")
        for cdr3 in ref["cdr3"].to_list():
            t = trim(cdr3)
            if len(t) > MAXLEN:
                continue
            f.write(f"{t}\t{cdr3}\t{task}\t{EPI[task]}\t\t\n")
            kept += 1
    return kept


def main():
    pdir = Path("bench/predictions/tcrmatch"); pdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        for task, pep in EPI.items():
            d = task_table(task, "TRB")
            qs = [q for q in d["cdr3"].to_list() if len(trim(q)) <= MAXLEN]
            t2o = defaultdict(list)
            for q in qs:
                t2o[trim(q)].append(q)
            ref = Path(td) / f"{task}_ref.tsv"
            n_ref = write_ref(task, ref)
            qfile = Path(td) / f"{task}_q.txt"
            qfile.write_text("\n".join(t2o) + "\n")
            out = run_tcrmatch(qfile, ref, 16, 0.84)               # low threshold -> score gradient
            best: dict[str, float] = defaultdict(float)
            for line in out.splitlines():
                f = line.split("\t")
                if len(f) < 5:
                    continue
                try:
                    score = float(f[2])
                except ValueError:
                    continue                                       # header / malformed
                ti, tm = f[0], f[1]
                if ti == tm:
                    continue                                       # drop exact self
                for orig in t2o.get(ti, []):
                    if score > best[orig]:
                        best[orig] = score
            rows = [(q, pep, best.get(q, 0.0), int(best.get(q, 0.0) >= SIG)) for q in d["cdr3"].to_list()]
            pl.DataFrame(rows, schema=["query_id", "epitope", "score", "significant"],
                         orient="row").write_csv(pdir / f"{task}_TRB.tsv", separator="\t")
            sc = {q: s for q, _, s, _ in rows}
            pairs = [(int(lab), sc[c]) for c, lab in zip(d["cdr3"], d["label"]) if c in sc]
            auc = roc_auc(pairs) if pairs else float("nan")
            npos = sum(p for p, _ in pairs)
            print(f"tcrmatch {task}/TRB: {len(rows)} scored (ref={n_ref}) | "
                  f"ROC-AUC {auc:.3f} (n_pos={npos} n={len(pairs)})")


if __name__ == "__main__":
    main()
