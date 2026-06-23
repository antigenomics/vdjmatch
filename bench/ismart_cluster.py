#!/usr/bin/env python3
"""iSMART -> clustering benchmark predictions on the shared vdjmatch clonotype sets.

Idempotent. Run with the vdjmatch venv:

    ./.venv/bin/python bench/ismart_cluster.py

iSMART clusters CDR3-beta (+ V gene) by k-mer-guided pairwise alignment. This script (venv) builds
the shared sets, writes a CDR3<TAB>V input, runs iSMARTf3.py in the `cmp-ismart` conda env with
DEFAULT params (threshold 7.5, KmerSize 6, V-gene on), parses the `*_clustered_v3.txt` output
(rows: <orig cols><TAB><Group>; Group is the LAST column), scores, and writes
bench/predictions/ismart/clustering.tsv.

iSMART is beta-only (human TRBV scoring), so only TRB is run. TRA/paired are not supported.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _cluster_common as C  # noqa: E402

ENV = "cmp-ismart"
ISMART_DIR = Path(__file__).resolve().parent / "external" / "iSMART"
WORK = Path(__file__).resolve().parent / "out" / "_ismart"


def run_ismart(cdr3, v):
    WORK.mkdir(parents=True, exist_ok=True)
    inp = WORK / "trb_in.txt"
    # iSMART input (column order): CDR3, V gene, Frequency, Other. Header optional (not starting with C-..F).
    # iSMART's VScore table is keyed on IMGT allele names (e.g. 'TRBV7-9*01'); an absent name makes
    # falign() return 0 and drops the cluster edge. The shared set is subgroup-level (B.vgene); append
    # '*01' so it matches Imgt_Human_TRBV.fasta (2785/2786) and V-gene scoring works as designed.
    with inp.open("w") as fh:
        fh.write("CDR3\tVgene\n")
        for c, vg in zip(cdr3, v):
            fh.write(f"{c}\t{vg}*01\n")
    # default params: -t 7.5, -K 6, V-gene on. iSMART reads Imgt_Human_TRBV.fasta / VgeneScores.txt
    # by relative path, so run from its own dir.
    subprocess.run(["conda", "run", "-n", ENV, "python", "iSMARTf3.py",
                    "-f", str(inp.resolve()), "-o", str(WORK.resolve())],
                   check=True, cwd=ISMART_DIR)
    out = WORK / "trb_in_clustered_v3.txt"
    if not out.exists():
        cands = list(WORK.glob("*_clustered_v3.txt"))
        if not cands:
            raise FileNotFoundError(f"iSMART output not found in {WORK}")
        out = cands[0]
    groups = {}
    for line in out.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        # header line (col0 not a CDR3) -> skip
        if not (parts[0].startswith("C") and len(parts) >= 2):
            continue
        # last column is the Group id; col0 is CDR3
        try:
            cid = int(parts[-1])
        except ValueError:
            continue
        groups[parts[0]] = cid
    return groups


def main():
    s = C.sets()
    d = s["TRB"]
    groups = run_ismart(d["cdr3"], d["v"])
    labels = C.labels_from_groups(d["cdr3"], groups)
    pur, ret, nc = C.score(labels, d["epi"])
    rows = [{"set": "TRB", "macro_purity": pur, "retention": ret, "n_clusters": nc,
             "n_clonotypes": len(d["cdr3"]), "note": "default thr=7.5, K=6, V-gene on"}]
    C.write_tsv("ismart", rows)


if __name__ == "__main__":
    main()
