#!/usr/bin/env python3
"""Produce predictions/<method>/samples.tsv for TCRMatch, consumed by bench/compare.py (samples mode).

Two reference arms (EXTERNAL_TOOLS.md Q1 = "use both"):
  tcrmatch        -> bundled IEDB reference (external/tcrmatch-src/data/IEDB_data.tsv)
  tcrmatch-vdjdb  -> leakage-removed VDJdb-beta reference, rebuilt as the 6-col TCRMatch DB

Queries are built identically to compare.run_samples (same sample1/2/5 + --olga-n) so query_id
(original CDR3) aligns. Exact self-matches (match_sequence == trimmed query) are dropped in both
arms, mirroring the vdjmatch arm's exclude_exact=True. Run from repo root with the project venv:

    python bench/tcrmatch_samples.py --olga-n 5000
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))   # bench/ -> compare, _bench
import _bench
from compare import SAMPLE_EPI, load_sample
from vdjmatch import db

SRC = Path(__file__).resolve().parent / "external" / "tcrmatch-src"
BIN = SRC / "tcrmatch"
IEDB_REF = SRC / "data" / "IEDB_data.tsv"
MAXLEN = 30   # TCRMatch's k3_sum has a fixed [31][31][31] stack array; trimmed len >=32 segfaults.


def trim(s: str) -> str:
    """Replicate TCRMatch's flank trimming (src/tcrmatch.cpp::trim)."""
    return s[1:-1] if len(s) > 1 and s[0] == "C" and s[-1] in "FW" else s


def build_queries(olga_n: int) -> list[str]:
    q1, q2, q5 = load_sample("sample1"), load_sample("sample2"), load_sample("sample4")  # TRB OLGA
    if olga_n and q5.height > olga_n:
        q5 = q5.sample(olga_n, seed=0)
    return pl.concat([q1, q2, q5]).unique("cdr3")["cdr3"].to_list()


def write_vdjdb_ref(species: str, path: Path):
    """VDJdb-beta unique-CDR3 -> 6-col TCRMatch DB (col0=trimmed key, col3=epitope)."""
    vdj = db.load(_bench.source(), species=species).filter(pl.col("gene") == "TRB")
    ref = _bench.valid_cdr3(vdj).group_by("cdr3").agg(pl.col("epitope").first())
    kept = dropped = 0
    with open(path, "w") as f:
        f.write("trimmed_seq\toriginal_seq\treceptor_group\tepitopes\tsource_organisms\tsource_antigens\n")
        for cdr3, epi in zip(ref["cdr3"], ref["epitope"]):
            t = trim(cdr3)
            if len(t) > MAXLEN:
                dropped += 1
                continue
            f.write(f"{t}\t{cdr3}\tvdjdb\t{epi}\t\t\n")
            kept += 1
    print(f"  vdjdb ref: {kept} entries (dropped {dropped} with trimmed len > {MAXLEN})")


def run_tcrmatch(qfile: Path, ref: Path, threads: int, threshold: float) -> str:
    cmd = [str(BIN), "-i", str(qfile), "-t", str(threads), "-d", str(ref), "-s", str(threshold)]
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def parse(stdout: str, trimmed2orig: dict[str, list[str]]):
    """-> (scores[(qid, epitope)] = max score, significant qids). Exact self-matches dropped."""
    scores: dict[tuple[str, str], float] = {}
    sig: set[str] = set()
    for line in stdout.splitlines():
        if not line or line.startswith("trimmed_input_sequence"):
            continue
        p = line.split("\t")
        ti, match, score = p[0], p[1], float(p[2])
        if match == ti:                                    # exact self -> mirror exclude_exact
            continue
        epis = p[4] if len(p) > 4 else ""
        for qid in trimmed2orig.get(ti, ()):
            sig.add(qid)
            for e in set(epis.split(",")):
                if e:
                    k = (qid, e)
                    scores[k] = max(scores.get(k, 0.0), score)
    return scores, sig


def emit(scores, sig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    scored_q = set()
    with open(path, "w") as f:
        f.write("query_id\tepitope\tscore\tsignificant\n")
        for (qid, epi), sc in scores.items():
            f.write(f"{qid}\t{epi}\t{sc:.6f}\t1\n")
            scored_q.add(qid)
        for qid in sig - scored_q:                         # matched but no epitope label -> still significant
            f.write(f"{qid}\t_OTHER_\t1.0\t1\n")
    return len(scored_q | sig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--olga-n", type=int, default=0, help="subsample OLGA negatives (match compare.py)")
    ap.add_argument("--threshold", type=float, default=0.97,
                    help="tcrmatch -s reporting threshold; 0.97 = binary call (precision), "
                         "lower = nearest-neighbour score gradient for ROC/PR")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--pred-dir", default="bench/predictions")
    ap.add_argument("--arms", nargs="+", default=["iedb", "vdjdb"], choices=["iedb", "vdjdb"])
    args = ap.parse_args()
    if not BIN.exists():
        sys.exit(f"tcrmatch binary missing: {BIN}\nbuild it (see bench/external/tcrmatch-src).")

    qlist = build_queries(args.olga_n)
    trimmed2orig: dict[str, list[str]] = {}
    dropped = 0
    for q in qlist:
        t = trim(q)
        if len(t) > MAXLEN:                                # would segfault tcrmatch -> no-call (score 0)
            dropped += 1
            continue
        trimmed2orig.setdefault(t, []).append(q)
    feed = [origs[0] for origs in trimmed2orig.values()]   # one original per trimmed group (tcrmatch re-trims)
    print(f"queries: {len(qlist)} unique CDR3 -> {len(feed)} unique trimmed fed "
          f"(dropped {dropped} with trimmed len > {MAXLEN})")

    pred = Path(args.pred_dir)
    with tempfile.TemporaryDirectory() as td:
        qfile = Path(td) / "queries.txt"
        qfile.write_text("\n".join(feed) + "\n")
        arms = {"iedb": ("tcrmatch", IEDB_REF),
                "vdjdb": ("tcrmatch-vdjdb", Path(td) / "vdjdb_ref.tsv")}
        for arm in args.arms:
            method, ref = arms[arm]
            if arm == "vdjdb":
                print(f"building VDJdb-beta reference ({args.species})...")
                write_vdjdb_ref(args.species, ref)
            print(f"[{method}] running tcrmatch (-s {args.threshold}, {args.threads} threads) vs {ref.name}...")
            scores, sig = parse(run_tcrmatch(qfile, ref, args.threads, args.threshold), trimmed2orig)
            n = emit(scores, sig, pred / method / "samples.tsv")
            print(f"[{method}] {len(sig)} queries matched; wrote {pred / method / 'samples.tsv'} "
                  f"({n} query ids, {len(scores)} (query,epitope) scores)")


if __name__ == "__main__":
    main()
