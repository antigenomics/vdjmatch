#!/usr/bin/env python3
"""TCRMatch (beta-only) -> benchmark prediction files for the sample conditions C1/C2/C3 (TRB).

Each sample CDR3-beta query is matched (TCRMatch BLOSUM62 kernel) against the VDJdb2026-beta reference;
score(query, E) = max kernel score over matched references labelled with epitope E (exact self-hits
dropped, mirroring exclude_exact). Read from manuscript test_data at runtime (never copied).

    python bench/tcrmatch_bench.py
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench
from compare import TESTDATA
from tcrmatch_samples import BIN, MAXLEN, run_tcrmatch, trim, write_vdjdb_ref

EPI = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV", "GLC": "GLCTLVAML", "YLQ": "YLQPRTFLL"}
PLAN = {  # condition -> (query loader, candidate epitopes); all TRB
    "C1": ("sample2", [EPI["LLW"], EPI["LLL"]]),
    "C2": ("sample1", [EPI["NLV"]]),
    "C3": ("sample6", [EPI["GLC"], EPI["YLQ"]]),
}


def queries(cond: str) -> list[str]:
    if cond == "C1":
        d = pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t").select("cdr3")
    elif cond == "C2":
        d = (pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
             .filter(pl.col("gene") == "TRB").select("cdr3"))
    else:
        d = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv").select(cdr3="cdr3_beta_aa")
    return _bench.valid_cdr3(d).unique("cdr3")["cdr3"].to_list()


def main():
    pdir = Path("bench/predictions/tcrmatch"); pdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        ref = Path(td) / "ref.tsv"
        write_vdjdb_ref("HomoSapiens", ref)                       # 6-col vdjdb-beta DB w/ epitope col
        for cond, (_, cand) in PLAN.items():
            qs = [q for q in queries(cond) if len(trim(q)) <= MAXLEN]
            t2o = defaultdict(list)
            for q in qs:
                t2o[trim(q)].append(q)
            qfile = Path(td) / f"{cond}_q.txt"
            qfile.write_text("\n".join(t2o) + "\n")
            out = run_tcrmatch(qfile, ref, 16, 0.84)              # low threshold -> score gradient
            best: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
            for line in out.splitlines():
                f = line.split("\t")
                if len(f) < 5:
                    continue
                try:
                    score = float(f[2])
                except ValueError:
                    continue                                      # header / malformed line
                ti, tm, epi = f[0], f[1], f[4]
                if ti == tm or epi not in cand:
                    continue                                      # drop exact self; keep candidate epis
                for orig in t2o.get(ti, []):
                    if score > best[orig][epi]:
                        best[orig][epi] = score
            rows = [(q, e, s, int(s >= 0.97)) for q, es in best.items() for e, s in es.items()]
            pl.DataFrame(rows, schema=["query_id", "epitope", "score", "significant"],
                         orient="row").write_csv(pdir / f"{cond}_TRB.tsv", separator="\t")
            print(f"tcrmatch {cond}/TRB: {len(qs)} queries -> {len(rows)} (query,epitope) scores")


if __name__ == "__main__":
    main()
