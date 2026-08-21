"""Pgen-matched mock negatives for the LLW/LLL detection tasks (follows mir/comparative/vdjbet.py).

For each positive set we build a negative set of CONTROL clonotypes with (i) exactly the same number of
clonotypes and (ii) the same generation-probability distribution: bin every CDR3 by round(log2 Pgen))
(OLGA), and for each Pgen bin sample as many control sequences (from an OLGA-generated pool) as there are
positives in that bin. This makes the negatives "as generatable as the binders", so detection measures
antigen-specificity, not Pgen. Deterministic (seed 42). LLW/LLL are TRB-only (sample2 has no TRA).

    .venv/bin/python bench/mock_negatives.py        # writes bench/_mock/mock_<epi>_<locus>.tsv
"""
from __future__ import annotations

import math
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
from benchmark import EPI, vgene                                     # noqa: E402
from compare import TESTDATA                                         # noqa: E402

OUT = Path(__file__).resolve().parent / "_mock"
OUT.mkdir(exist_ok=True)
POOL_N = 40000                                                       # OLGA control sequences to Pgen-bin
_OLGA_SUB = {"TRB": "human_T_beta", "TRA": "human_T_alpha"}
_OLGA_FILE = {"TRB": "sample4_olga_airr.txt", "TRA": "sample5_olga_airr.txt"}
_MODEL: dict = {}


def pgen_model(locus):
    if locus not in _MODEL:
        import olga.load_model as lm
        import olga.generation_probability as gp
        d = Path(os.path.dirname(__import__("olga").__file__)) / "default_models" / _OLGA_SUB[locus]
        gm, gd = lm.GenerativeModelVDJ(), lm.GenomicDataVDJ()
        gd.load_igor_genomic_data(str(d / "model_params.txt"), str(d / "V_gene_CDR3_anchors.csv"),
                                  str(d / "J_gene_CDR3_anchors.csv"))
        gm.load_and_process_igor_model(str(d / "model_marginals.txt"))
        _MODEL[locus] = gp.GenerationProbabilityVDJ(gm, gd)
    return _MODEL[locus]


def log2_bin(cdr3, model):
    try:
        p = model.compute_aa_CDR3_pgen(cdr3)
    except Exception:
        return None
    return round(math.log2(p)) if p and p > 0 else None


def control_bins(locus):
    """OLGA control pool grouped by round(log2 Pgen) bin: {bin: list of (cdr3, v, j)}. Cached."""
    cf = OUT / f"control_bins_{locus}_{POOL_N}.parquet"
    if not cf.exists():
        m = pgen_model(locus)
        d = (pl.read_csv(TESTDATA / _OLGA_FILE[locus], separator="\t", infer_schema_length=0)
             .select(cdr3="junction_aa", v="v_gene", j="j_gene").pipe(_bench.valid_cdr3)
             .unique("cdr3", maintain_order=True).head(POOL_N))
        rows = []
        for c, v, j in zip(d["cdr3"], d["v"], d["j"]):
            b = log2_bin(c, m)
            if b is not None:
                rows.append((b, c, vgene(v), vgene(j)))
        pl.DataFrame(rows, schema=["bin", "cdr3", "v", "j"], orient="row").write_parquet(cf)
        print(f"[control {locus}] binned {len(rows)}/{d.height} OLGA seqs", file=sys.stderr)
    p = pl.read_parquet(cf)
    bins = defaultdict(list)
    for b, c, v, j in p.iter_rows():
        bins[b].append((c, v, j))
    return bins


def positives(epi, locus):
    """sample2 LLW/LLL positive CDR3s (TRB only)."""
    d = (pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t", infer_schema_length=0)
         .rename({"antigen.epitope": "lab"}).select("cdr3", "lab", v="v.segm", j="j.segm"))
    return _bench.valid_cdr3(d).unique("cdr3").filter(pl.col("lab") == EPI[epi])


def make_mock(epi, locus, seed=42):
    m = pgen_model(locus)
    pos = positives(epi, locus)
    hist = defaultdict(int)
    skipped = 0
    for c in pos["cdr3"]:
        b = log2_bin(c, m)
        if b is None:
            skipped += 1
        else:
            hist[b] += 1
    pool = control_bins(locus)
    avail = sorted(pool)
    rng = np.random.default_rng(seed)
    mock = []
    for b, count in sorted(hist.items()):
        src = b if pool.get(b) else min(avail, key=lambda x: abs(x - b))   # nearest non-empty bin
        recs = pool[src]
        idx = rng.choice(len(recs), count, replace=len(recs) < count)
        mock += [recs[int(i)] for i in idx]
    df = pl.DataFrame(mock, schema=["cdr3", "v", "j"], orient="row")
    df.write_csv(OUT / f"mock_{epi}_{locus}.tsv", separator="\t")
    print(f"{epi}-{locus}: {pos.height} positives ({skipped} no-Pgen) -> {df.height} Pgen-matched mock "
          f"negatives across {len(hist)} bins")
    return df


if __name__ == "__main__":
    for epi in ("LLW", "LLL"):
        make_mock(epi, "TRB")
