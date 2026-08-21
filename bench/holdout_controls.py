"""Pgen-matched airr_control negatives for the held-out LLW/LLL/ELA datasets.

A negative set of REAL post-selection control TCRs (isalgo/airr_control) matched to the positives in
(i) exactly the same number of clonotypes and (ii) the same generation-probability distribution
(round(log2 Pgen)), Pgen computed with OLGA. Build a precomputed log2-Pgen-binned airr_control pool per
(species, locus) once (parallel), then sample per epitope. Deterministic (seed 42).

    .venv/bin/python bench/holdout_controls.py        # builds pools + writes negatives into hold_out_data
"""
from __future__ import annotations

import math
import multiprocessing as mp
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                       # noqa: E402
from benchmark import vgene                                        # noqa: E402


def _airr(species, locus, n):
    """isalgo/airr_control productive (cdr3, v, j) for a species/locus, deterministically sampled to n."""
    from huggingface_hub import hf_hub_download
    f = hf_hub_download("isalgo/airr_control", repo_type="dataset",
                        filename=f"{species}.{locus.lower()}.aa.vdjtools.tsv.gz")
    d = (pl.read_csv(f, separator="\t", columns=["cdr3aa", "v", "j"], infer_schema_length=0)
         .rename({"cdr3aa": "cdr3"}).pipe(_bench.valid_cdr3)
         .filter(~pl.col("cdr3").str.contains(r"[*_~]")).unique("cdr3", maintain_order=True)
         .with_columns(v=pl.col("v").map_elements(vgene, return_dtype=pl.Utf8),
                       j=pl.col("j").map_elements(vgene, return_dtype=pl.Utf8)))
    return d.sort("cdr3").sample(n, seed=42) if d.height > n else d

HD = Path.home() / "vcs/manuscripts/2026-vdjmatch/hold_out_data"
POOLDIR = Path.home() / "vcs/manuscripts/2026-vdjmatch/extreme-optimization/pgen_pools"
POOLDIR.mkdir(parents=True, exist_ok=True)
POOL_N = 80000
_SUB = {("human", "TRB"): "human_T_beta", ("human", "TRA"): "human_T_alpha",
        ("mouse", "TRB"): "mouse_T_beta", ("mouse", "TRA"): "mouse_T_alpha"}
_WMODEL = {}


def _winit(sub):
    import olga.generation_probability as gp
    import olga.load_model as lm
    d = Path(os.path.dirname(__import__("olga").__file__)) / "default_models" / sub
    vj = "alpha" in sub                                            # alpha chain = V-J recombination (no D)
    gm = lm.GenerativeModelVJ() if vj else lm.GenerativeModelVDJ()
    gd = lm.GenomicDataVJ() if vj else lm.GenomicDataVDJ()
    gd.load_igor_genomic_data(str(d / "model_params.txt"), str(d / "V_gene_CDR3_anchors.csv"),
                              str(d / "J_gene_CDR3_anchors.csv"))
    gm.load_and_process_igor_model(str(d / "model_marginals.txt"))
    _WMODEL["m"] = (gp.GenerationProbabilityVJ if vj else gp.GenerationProbabilityVDJ)(gm, gd)


def _wbin(cdr3):
    try:
        p = _WMODEL["m"].compute_aa_CDR3_pgen(cdr3)
    except Exception:
        return None
    return round(math.log2(p)) if p and p > 0 else None


def pgen_pool(species, locus, n=POOL_N):
    """log2-Pgen-binned airr_control pool: parquet (bin, cdr3, v, j). Built once per (species, locus)."""
    pf = POOLDIR / f"pool_{species}_{locus}_{n}.parquet"
    if not pf.exists():
        d = _airr(species, locus, n)            # airr_control productive (cdr3, v, j)
        cd = d["cdr3"].to_list()
        with mp.Pool(max(1, (os.cpu_count() or 2) - 1), initializer=_winit, initargs=(_SUB[(species, locus)],)) as pool:
            bins = pool.map(_wbin, cd, chunksize=200)
        rows = [(b, c, v, j) for b, c, v, j in zip(bins, cd, d["v"].to_list(), d["j"].to_list()) if b is not None]
        pl.DataFrame(rows, schema=["bin", "cdr3", "v", "j"], orient="row").write_parquet(pf)
        print(f"[pool {species} {locus}] {len(rows)}/{len(cd)} binned", file=sys.stderr)
    return pl.read_parquet(pf)


def matched_negatives(pos_cdr3, species, locus, seed=42):
    """Sample airr_control negatives matching the positives' log2-Pgen histogram (count-exact)."""
    _winit(_SUB[(species, locus)])                                  # single-process model for the (few) positives
    hist = defaultdict(int)
    for c in pos_cdr3:
        b = _wbin(c)
        if b is not None:
            hist[b] += 1
    pool = pgen_pool(species, locus)
    bybin = defaultdict(list)
    for b, c, v, j in pool.iter_rows():
        bybin[b].append((c, v, j))
    avail = sorted(bybin)
    rng = np.random.default_rng(seed)
    out = []
    for b, cnt in sorted(hist.items()):
        src = b if bybin.get(b) else min(avail, key=lambda x: abs(x - b))
        recs = bybin[src]
        for i in rng.choice(len(recs), cnt, replace=len(recs) < cnt):
            out.append(recs[int(i)])
    return out


def add_negatives(name, species="human"):
    locus = name.split("_")[1]
    d = pl.read_csv(HD / f"{name}.tsv", separator="\t", infer_schema_length=0)
    pos = d.filter(pl.col("label") == "1")
    neg = matched_negatives(pos["cdr3"].to_list(), species, locus)
    negdf = pl.DataFrame(neg, schema=["cdr3", "v", "j"], orient="row").with_columns(
        label=pl.lit("0"), source=pl.lit(f"airr_control_pgen_{locus}"))
    out = pl.concat([pos.with_columns(pl.col("label").cast(pl.Utf8)), negdf.select(pos.columns)])
    out.write_csv(HD / f"{name}.tsv", separator="\t")
    print(f"  {name}: {pos.height} pos + {negdf.height} Pgen-matched airr_control neg")


if __name__ == "__main__":
    for nm in ("LLW_TRB", "LLL_TRB", "ELA_TRB", "LLL_TRA"):
        add_negatives(nm)
