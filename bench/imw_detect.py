#!/usr/bin/env python3
"""imw-DETECT (ImmuneWatch DETECT) predictions -> detection-benchmark prediction files.

Per dataset_descr.md: DETECT supplies a top-1 epitope + Score per query TCR. For detection task E, a
query is a "hit" iff DETECT predicted E (score is its confidence); predicting any other epitope = no hit
(score 0). Tasks: NLV (sample1), LLW/LLL (sample2), YLQ/GLC (sample6, TRA+TRB flattened). Read from the
manuscript test_data at runtime (never copied).

    python bench/imw_detect.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import TESTDATA

ODS = TESTDATA / "immunedetect_results"
THRESH = 0.5                                                      # DETECT call threshold for the confusion matrix
# task -> (ods sample, target epitope, loci)
PLAN = {
    "NLV": ("sample1", "NLVPMVATV", ["TRB"]),
    "LLW": ("sample2", "LLWNGPMAV", ["TRB"]),
    "LLL": ("sample2", "LLLGIGILV", ["TRB"]),
    "YLQ": ("sample6", "YLQPRTFLL", ["TRA", "TRB"]),
    "GLC": ("sample6", "GLCTLVAML", ["TRA", "TRB"]),
}


def main():
    pdir = Path("bench/predictions/imw-detect"); pdir.mkdir(parents=True, exist_ok=True)
    for task, (sample, target, loci) in PLAN.items():
        d = pd.read_excel(ODS / f"predictions_{sample}.tsv.ods", engine="odf")
        d = d.rename(columns={"junction_aa": "cdr3", "Epitope": "epitope", "Score": "score"})
        d["chain"] = d["v_call"].astype(str).str[:3].map({"TRA": "TRA", "TRB": "TRB"})
        for locus in loci:
            sub = d[d["chain"] == locus]
            rows = []
            for cdr3, epi, score in zip(sub["cdr3"], sub["epitope"], sub["score"]):
                if str(epi) == target:                            # DETECT predicted the target epitope
                    rows.append((cdr3, target, float(score), int(float(score) >= THRESH)))
            pl.DataFrame(rows, schema=["query_id", "epitope", "score", "significant"],
                         orient="row").write_csv(pdir / f"{task}_{locus}.tsv", separator="\t")
            print(f"imw-detect {task}/{locus}: {len(rows)} predicted-{target} of {sub.shape[0]} TCRs")


if __name__ == "__main__":
    main()
