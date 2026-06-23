#!/usr/bin/env python3
"""GIANA -> clustering benchmark predictions on the shared vdjmatch clonotype sets.

Idempotent. Run with the vdjmatch venv:

    ./.venv/bin/python bench/giana_cluster.py

GIANA (Geometric Isometry based ANtigen-specific tcr Alignment) clusters CDR3-beta with V-gene.
This script (venv) builds the shared sets, writes a CDR3<TAB>V input file, runs GIANA4.py in the
`cmp-giana` conda env with DEFAULT params (isometric thr=7, SW exact=on, V-gene=on), parses the
`*--RotationEncodingBL62.txt` output (lines: CDR3<TAB>cluster_id<TAB>...), scores, and writes
bench/predictions/giana/clustering.tsv.

GIANA is beta-only (uses human TRBV scoring), so only TRB is run. TRA/paired are not supported.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _cluster_common as C  # noqa: E402

ENV = "cmp-giana"
GIANA_DIR = Path(__file__).resolve().parent / "external" / "GIANA"
WORK = Path(__file__).resolve().parent / "out" / "_giana"


def full_v(vfam):
    """GIANA's VScore table is keyed on IMGT allele names from Imgt_Human_TRBV.fasta (e.g.
    'TRBV7-9*01'); a name absent from the table makes falign() return 0 and drops the cluster edge.
    The shared set gives subgroup-level names (B.vgene, e.g. 'TRBV7-9'); appending '*01' matches the
    fasta for 2785/2786 clonotypes, so the V-gene scoring works as designed (standard GIANA usage)."""
    return vfam + "*01"


def run_giana(cdr3, v):
    WORK.mkdir(parents=True, exist_ok=True)
    inp = WORK / "trb_in.txt"
    # GIANA input: col0 CDR3 (starts with C), col1 V gene. Header line optional (must not start with C).
    with inp.open("w") as fh:
        fh.write("CDR3\tV\n")
        for c, vg in zip(cdr3, v):
            fh.write(f"{c}\t{full_v(vg)}\n")
    # default params: thr=7 (-t), exact SW on, V-gene on. Output dir = WORK.
    # GIANA must run from its own dir (reads Imgt_Human_TRBV.fasta / VgeneScores.txt by relative path).
    subprocess.run(["conda", "run", "-n", ENV, "python", "GIANA4.py",
                    "-f", str(inp.resolve()), "-o", str(WORK.resolve()), "-t", "7"],
                   check=True, cwd=GIANA_DIR)
    # output filename: <inputstem>--RotationEncodingBL62.txt in -o dir
    out = WORK / "trb_in--RotationEncodingBL62.txt"
    if not out.exists():
        cands = list(WORK.glob("*RotationEncodingBL62.txt"))
        if not cands:
            raise FileNotFoundError(f"GIANA output not found in {WORK}")
        out = cands[0]
    groups = {}
    for line in out.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        c, cid = parts[0], parts[1]
        groups[c] = int(cid)
    return groups


def main():
    s = C.sets()
    d = s["TRB"]
    groups = run_giana(d["cdr3"], d["v"])
    labels = C.labels_from_groups(d["cdr3"], groups)
    pur, ret, nc = C.score(labels, d["epi"])
    rows = [{"set": "TRB", "macro_purity": pur, "retention": ret, "n_clusters": nc,
             "n_clonotypes": len(d["cdr3"]), "note": "default thr=7, SW exact, V-gene on"}]
    C.write_tsv("giana", rows)


if __name__ == "__main__":
    main()
