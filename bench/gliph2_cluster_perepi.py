#!/usr/bin/env python3
"""GLIPH2 (immGLIPH-equivalent) -> per-epitope clustering on the shared TRB clonotype set.

Runs under the vdjmatch venv. GLIPH2 = the Huang 2020 ``irtools`` binary (the same engine the manuscript
calls "immGLIPH") with the bundled human ``CD48_v2.0`` background (ref_CD48_v2.0.fa / ref_V / ref_L).
GLIPH2 is beta-only, so only TRB is produced.

Pipeline (all local, macOS Mach-O binary via Rosetta):
  1. build the EXACT shared TRB set (CR.single_clonotypes(sl, "TRB"));
  2. write a GLIPH2 metarepertoire (CDR3b<TAB>TRBV<TAB>TRBJ<TAB>CDR3a<TAB>subject:condition<TAB>count)
     into a per-run work dir alongside the binary + reference + parameter file;
  3. run ``irtools.osx -c parameters`` (algorithm=GLIPH2, CD48 reference, defaults from the GLIPH2 readme /
     clusTCR template);
  4. parse the ``*_cluster.txt`` convergence groups -> {cdr3 -> cluster_id}; a CDR3 in several motif groups
     takes its first (largest-listed) group; CDR3s in no >=2 group stay singletons;
  5. build labels aligned 1:1 to the shared set and report CR.per_epitope_scores + aggregate (anchor
     0.971 / 0.444).

    ./.venv/bin/python bench/gliph2_cluster_perepi.py
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parent
sys.path.insert(0, "/Users/mikesh/vcs/manuscripts/2026-vdjmatch/benchmarks/scripts")
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(BENCH.parent / "src"))
import _cluster_common as C  # noqa: E402
import cluster_results as CR  # noqa: E402

LIB = BENCH / "external" / "gliph2" / "lib"          # irtools.osx + CD48 reference live here
WORK = BENCH / "out" / "_gliph2"
BIN = "irtools.osx"

PARAMS = """# vdjmatch GLIPH2 benchmark — Huang 2020 CD48_v2.0 background, default thresholds
out_prefix=trb_out
cdr3_file=trb_in.txt
hla_file=hla.txt
refer_file=ref_CD48_v2.0.fa
v_usage_freq_file=ref_V_CD48_v2.0.txt
cdr3_length_freq_file=ref_L_CD48_v2.0.txt
local_min_pvalue=0.001
p_depth = 1000
global_convergence_cutoff = 1
simulation_depth=1000
kmer_min_depth=3
local_min_OVE=10
algorithm=GLIPH2
all_aa_interchangeable=1
"""


def _setup_work():
    """Fresh work dir with the binary, the three CD48 reference files, and the parameter file."""
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    for f in (BIN, "ref_CD48_v2.0.fa", "ref_V_CD48_v2.0.txt", "ref_L_CD48_v2.0.txt"):
        shutil.copy(LIB / f, WORK / f)
    (WORK / BIN).chmod(0o755)
    (WORK / "parameters").write_text(PARAMS)
    # this irtools build requires an hla_file; single subject "01" with a generic A*02 line. HLA only
    # feeds the hla_score, not cluster membership, so it does not change the convergence groups.
    (WORK / "hla.txt").write_text("01\tA*02:01\tA*02:01\n")


def _write_input(cdr3, v, j):
    """GLIPH2 metarepertoire: CDR3b, TRBV, TRBJ, CDR3a, subject:condition, count (tab-delimited)."""
    with (WORK / "trb_in.txt").open("w") as fh:
        for c, vg, jg in zip(cdr3, v, j):
            fh.write(f"{c}\t{vg}\t{jg}\tNA\t01:bench\t1\n")


def _parse_clusters(path):
    """Parse the GLIPH2 *_cluster.txt convergence-group file -> {cdr3 -> cluster_id}.

    GLIPH2 reports OVERLAPPING motif groups (one per line; the same CDR3 recurs across nested motifs,
    e.g. ``lSIRS`` ⊂ ``lSIR``). To turn this into a partition for purity/retention we assign each CDR3
    to its MOST SPECIFIC (smallest) group — the standard GLIPH2 specificity-group reading. The CDR3b
    members on a line are the whitespace-separated tokens that look like a CDR3 (start 'C', end 'F'/'W',
    length >= 8, all upper-case AA). Only groups with >= 2 members are clusters.
    """
    AA = set("ACDEFGHIKLMNPQRSTVWY")

    def is_cdr3(tok):
        return (len(tok) >= 8 and tok[0] == "C" and tok[-1] in "FW"
                and all(ch in AA for ch in tok))

    grps = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        members = [t for t in line.split() if is_cdr3(t)]
        members = list(dict.fromkeys(members))                  # unique within the group
        if len(members) >= 2:
            grps.append(members)
    groups = {}                                                 # smallest group wins (most specific motif)
    for cid, k in enumerate(sorted(range(len(grps)), key=lambda i: len(grps[i]))):
        for m in grps[k]:
            groups.setdefault(m, cid)
    return groups


def run_gliph2(cdr3, v, j):
    _setup_work()
    _write_input(cdr3, v, j)
    # This irtools build SIGSEGVs on the FINAL hla-scoring/CSV step AFTER the cluster file is fully
    # written (trb_out_cluster.txt is complete; only the empty trb_out_cluster.csv is the casualty).
    # So we do not check the return code — we validate the cluster output instead.
    subprocess.run([f"./{BIN}", "-c", "parameters"], cwd=WORK)
    out = WORK / "trb_out_cluster.txt"
    if not out.exists() or out.stat().st_size == 0:
        cands = [p for p in sorted(WORK.glob("*cluster.txt")) if p.stat().st_size > 0]
        if not cands:
            raise FileNotFoundError(f"GLIPH2 cluster output missing/empty in {WORK}: "
                                    f"{[p.name for p in WORK.glob('trb_out*')]}")
        out = cands[0]
    return _parse_clusters(out)


def produce():
    s = C.sets()
    d = s["TRB"]
    cdr3, epi = d["cdr3"], d["epi"]
    # GLIPH2 wants subgroup TRBV names (e.g. TRBV7-9); the shared v is already B.vgene subgroup-level.
    groups = run_gliph2(cdr3, d["v"], _jvec(cdr3))
    labels = C.labels_from_groups(cdr3, groups)
    return {"TRB": (labels, epi)}


def _jvec(cdr3):
    """Recover the TRBJ gene per shared-set CDR3 (first j per cdr3), aligned to `cdr3`."""
    import polars as pl
    import _bench  # noqa: PLC0415
    import benchmark as B  # noqa: PLC0415
    sl = B.shortlist(B.release("vdjdb2026"))
    m = (_bench.valid_cdr3(sl.filter(pl.col("gene") == "TRB")).group_by("cdr3").agg(pl.col("j").first()))
    jmap = dict(zip(m["cdr3"], m["j"]))
    return [str(jmap.get(c, "NA")).split("*")[0] for c in cdr3]   # GLIPH2 wants gene-level TRBJ


if __name__ == "__main__":
    res = produce()
    for chain, (labels, epi) in res.items():
        agg = CR.score_labels(labels, epi)
        print(f"GLIPH2 {chain}: aggregate purity={agg[0]} retention={agg[1]} clusters={agg[2]} "
              f"(n={len(labels)})")
