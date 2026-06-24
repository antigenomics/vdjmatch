#!/usr/bin/env python3
"""FAIR cross-tool noise false-positive benchmark (replaces the broken default-cutoff Fig 3a).

Each tool's score is on its own scale (imw probability 0-0.78, vdjmatch NED / E-value, tcrdist3
distance, TCRMatch similarity, ERGO/NetTCR P(bind)), so a fixed score cutoff is meaningless: imw's
"0% noise FP" at its default 0.5 is an artifact (at 0.5 it also detects only ~3.5% of REAL binders).
The fair comparison sets EACH tool's threshold to where it recalls 50% of true binders, then measures
the noise FP at that same operating point.

For each (tool, locus):
  1. confidence score per receptor = best-over-panel confidence (higher = more confident):
       vdjmatch  -> max NED score over the 5 epitopes (committed detection score, exact-excluded)
       tcrdist3  -> -min nearest-neighbour distance over the union of the 5 epitope references
       TCRMatch  -> max BLOSUM62 kernel score over the 5 epitope references (TRB only, beta kernel)
       imw-DETECT-> the .ods top-1 Score per receptor (precomputed, no model call)
       ERGO-II   -> max P(bind) over the 5 peptides (TRA via the beta-encoder, per ergo_detect_alpha)
       NetTCR-2.0-> max P(bind) over the 5 peptides (TRB only, beta model)
  2. score (a) the TRUE-BINDER set and (b) the 10k NOISE set, matched BY CDR3 SEQUENCE.
  3. recall-r threshold s* = quantile(true_binder_scores, 1-r); noise FP = fraction noise >= s*.
     Reported at r = 50% (the headline), 25%, 10%.

Panel = 5 HLA-A*02:01 epitopes: NLVPMVATV LLWNGPMAV LLLGIGILV YLQPRTFLL GLCTLVAML.

True binders (per locus): sample3_vdjdb.txt, species=HomoSapiens, gene=TR{A,B}, antigen.epitope in
panel, unique cdr3. NOTE: sample3 is the full VDJdb export, NOT the vdjdb2026 >=2-ref shortlist used
elsewhere -- a known limitation; future work is to recalibrate on the vdjdb2026 shortlist.

Noise (per locus): 10,000 unique CDR3 sampled at seed=42 from sample4 (TRB) / sample5 (TRA) OLGA AIRR.
``.unique('cdr3', maintain_order=True)`` is required for REPRODUCIBILITY: plain ``.unique`` does not
fix row order in polars, so ``.sample(seed=42)`` after it draws a different 10k each run.

All data is read from the manuscript test_data at runtime (never copied into git: numbers only).
Heavy tools (tcrdist/ERGO/NetTCR) shell into the cmp-* conda envs and cache per-receptor scores under
bench/out/_noise_r50/; --assemble rebuilds the summary + figure data from those caches.

    .venv/bin/python bench/noise_recall50.py            # run all tools + assemble
    .venv/bin/python bench/noise_recall50.py --assemble # rebuild summary from cached scores
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench
import benchmark as B
from compare import TESTDATA
from tcrdist_samples import norm_gene
from vdjmatch import db
from vdjmatch.evalue import background

PANEL = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV",
         "YLQ": "YLQPRTFLL", "GLC": "GLCTLVAML"}
PEPTIDES = list(PANEL.values())
MHC = "HLA-A*02"
N_NOISE, SEED = 10000, 42
# threshold at the top-X percentile of a method's sample3 true-binder scores (X = recall level), then
# measure noise FP there. 50% = median (detect half the binders); 2.5% = only the top 2.5% (very strict).
RECALLS = (0.50, 0.25, 0.05, 0.025)
NOISE_FILE = {"TRB": "sample4", "TRA": "sample5"}                   # OLGA AIRR: 4=TRB, 5=TRA
ODS = TESTDATA / "immunedetect_results"
IMW_SAMPLE = {"TRB": "sample4", "TRA": "sample5"}                   # imw noise predictions

OUT = Path("bench/out/_noise_r50"); OUT.mkdir(parents=True, exist_ok=True)
WORK = OUT / "_work"; WORK.mkdir(parents=True, exist_ok=True)
ENV_TCRDIST = "cmp-tcrdist"
TCRDIST_BETA = Path(__file__).resolve().parent / "_tcrdist_compute_beta_nn.py"
TCRDIST_ALPHA = Path(__file__).resolve().parent / "_tcrdist_compute_alpha.py"
ERGO = Path(__file__).resolve().parent / "external" / "ERGO-II"
NET = Path(__file__).resolve().parent / "external" / "NetTCR-2.0"
RADIUS = 90                                                        # tcrdist NN search radius (as in detection)


# --------------------------------------------------------------------------- query sets (cdr3, v, j)
def true_binders(locus: str) -> pl.DataFrame:
    """sample3_vdjdb true binders for a locus: species=HomoSapiens, gene=locus, panel epitope, unique
    cdr3, with allele-stripped v/j. (sample3 = full VDJdb, not the vdjdb2026 shortlist -- a known
    limitation; future work recalibrates on the vdjdb2026 shortlist.)"""
    d = (pl.read_csv(TESTDATA / "sample3_vdjdb.txt", separator="\t")
         .filter((pl.col("species") == "HomoSapiens") & (pl.col("gene") == locus)
                 & pl.col("antigen.epitope").is_in(PEPTIDES))
         .select("cdr3", v="v.segm", j="j.segm").pipe(_bench.valid_cdr3).unique("cdr3"))
    return d


def noise(locus: str) -> pl.DataFrame:
    """10k OLGA noise receptors for a locus, sampled deterministically (seed=42). maintain_order=True
    is required: plain .unique does not fix row order in polars, so the seeded sample is irreproducible
    without it. v/j from the OLGA AIRR file (v_gene/j_gene)."""
    d = (pl.read_csv(TESTDATA / f"{NOISE_FILE[locus]}_olga_airr.txt", separator="\t")
         .select(cdr3="junction_aa", v="v_gene", j="j_gene")
         .drop_nulls("cdr3").unique("cdr3", maintain_order=True).sample(N_NOISE, seed=SEED)
         .pipe(_bench.valid_cdr3))
    return d


# --------------------------------------------------------------------------- vdjmatch (NED score)
def score_vdjmatch(locus: str, tb: pl.DataFrame, ns: pl.DataFrame) -> tuple[dict, dict]:
    """Per-receptor confidence = max NED score over the 5 epitopes' A*02 VDJdb2026 references,
    exact-CDR3 excluded (the committed detection score). Returns (tb_scores, noise_scores) cdr3->score.

    ONE combined reference over the 5 epitopes + ONE vdjmatch_classify call (it already loops epitopes
    internally and searches the control once): the per-epitope NED is identical to scoring each
    epitope separately (only rare cross-epitope shared CDR3s differ, and we take the max), but ~5x faster
    than re-searching the large OLGA control once per epitope."""
    allq = tb["cdr3"].to_list() + ns["cdr3"].to_list()
    qv = {c: B.vgene(v) for c, v in zip(tb["cdr3"], tb["v"])}
    qv.update({c: B.vgene(v) for c, v in zip(ns["cdr3"], ns["v"])})
    ref = pl.concat([B.epi_ref(estr, locus) for estr in PEPTIDES])
    tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v = B.ref_index(ref, locus)
    sc, _ = B.vdjmatch_classify(tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v, background(locus),
                                allq, [qv[q] for q in allq], PEPTIDES, 1e-3, True,
                                params=B.first_hit.scope(5, 2, 2))
    best = {q: max(sc[q][e][0] for e in PEPTIDES) for q in allq}
    return ({c: best[c] for c in tb["cdr3"].to_list()},
            {c: best[c] for c in ns["cdr3"].to_list()})


# --------------------------------------------------------------------------- tcrdist3 (-min NN dist)
def _combined_ref(locus: str) -> pl.DataFrame:
    """Union of all 5 panel epitopes' A*02 VDJdb2026 references for a locus (cdr3, v, j), unique cdr3."""
    v26 = db.load(_bench.source(), species="HomoSapiens")
    refs = []
    for estr in PEPTIDES:
        r = _bench.valid_cdr3(v26.filter((pl.col("epitope") == estr) & pl.col("mhc_a").str.contains(B.A02)
                                         & (pl.col("gene") == locus)))
        refs.append(r.select("cdr3", v="v", j="j"))
    return (pl.concat(refs).with_columns(v=norm_gene(pl.col("v")), j=norm_gene(pl.col("j")))
            .unique("cdr3"))


def score_tcrdist(locus: str, tb: pl.DataFrame, ns: pl.DataFrame) -> tuple[dict, dict]:
    """Per-receptor confidence = -min NN tcrdist over the combined 5-epitope reference (exact removed,
    radius=90). Queries with no neighbour within radius get NaN -> -inf. Beta uses _tcrdist_compute_beta_nn,
    alpha uses _tcrdist_compute_alpha; both emit (query_cdr3, nn_dist)."""
    ref = _combined_ref(locus)
    reff = WORK / f"tcrdist_{locus}_ref.tsv"; ref.write_csv(reff, separator="\t")
    q = (pl.concat([tb.select("cdr3", "v", "j"), ns.select("cdr3", "v", "j")])
         .with_columns(v=norm_gene(pl.col("v")), j=norm_gene(pl.col("j"))).unique("cdr3"))
    qf = WORK / f"tcrdist_{locus}_q.tsv"; q.write_csv(qf, separator="\t")
    outp = WORK / f"tcrdist_{locus}_nn.tsv"
    compute = TCRDIST_BETA if locus == "TRB" else TCRDIST_ALPHA
    subprocess.run(["conda", "run", "-n", ENV_TCRDIST, "python", str(compute),
                    "--ref", str(reff.resolve()), "--queries", str(qf.resolve()),
                    "--radius", str(RADIUS), "--out", str(outp.resolve())], check=True)
    nn = pl.read_csv(outp, separator="\t")
    nn_by_cdr3 = dict(zip(nn["query_cdr3"].to_list(), nn["nn_dist"].to_list()))

    def conf(c):
        d = nn_by_cdr3.get(c)
        return float("-inf") if d is None or d != d else -float(d)   # NaN/missing -> -inf
    return ({c: conf(c) for c in tb["cdr3"].to_list()},
            {c: conf(c) for c in ns["cdr3"].to_list()})


# --------------------------------------------------------------------------- TCRMatch (max kernel, TRB)
def _trim(cdr3: str) -> str:
    s = cdr3
    if s.startswith("C"):
        s = s[1:]
    if s.endswith(("F", "W")):
        s = s[:-1]
    return s


def score_tcrmatch(tb: pl.DataFrame, ns: pl.DataFrame) -> tuple[dict, dict]:
    """Per-receptor confidence = max TCRMatch BLOSUM62 kernel score over the combined 5-epitope TRB
    reference, exact self-hits dropped (mirrors exclude_exact). TRB only (beta kernel)."""
    import tempfile
    from collections import defaultdict
    from tcrmatch_samples import MAXLEN, run_tcrmatch
    ref = _combined_ref("TRB")
    with tempfile.TemporaryDirectory() as td:
        refp = Path(td) / "ref.tsv"
        with open(refp, "w") as f:
            f.write("trimmed_seq\toriginal_seq\treceptor_group\tepitopes\tsource_organisms\tsource_antigens\n")
            for cdr3 in ref["cdr3"].to_list():
                t = _trim(cdr3)
                if len(t) <= MAXLEN:
                    f.write(f"{t}\t{cdr3}\tpanel\tPANEL\t\t\n")
        allq = tb["cdr3"].to_list() + ns["cdr3"].to_list()
        t2o = defaultdict(list)
        for q in allq:
            if len(_trim(q)) <= MAXLEN:
                t2o[_trim(q)].append(q)
        qfile = Path(td) / "q.txt"; qfile.write_text("\n".join(t2o) + "\n")
        out = run_tcrmatch(qfile, refp, 16, 0.84)                  # low threshold -> full score gradient
        best: dict[str, float] = defaultdict(float)
        for line in out.splitlines():
            ff = line.split("\t")
            if len(ff) < 5:
                continue
            try:
                score = float(ff[2])
            except ValueError:
                continue
            ti, tm = ff[0], ff[1]
            if ti == tm:                                           # drop exact self
                continue
            for orig in t2o.get(ti, []):
                if score > best[orig]:
                    best[orig] = score
    return ({c: best.get(c, 0.0) for c in tb["cdr3"].to_list()},
            {c: best.get(c, 0.0) for c in ns["cdr3"].to_list()})


# --------------------------------------------------------------------------- imw-DETECT (.ods top-1 Score)
def _imw_top1(sample: str) -> dict:
    """cdr3 -> top-1 (max-over-epitopes) Score from a precomputed imw-DETECT .ods prediction file."""
    d = pd.read_excel(ODS / f"predictions_{sample}.tsv.ods", engine="odf").rename(
        columns={"junction_aa": "cdr3", "Score": "score"})
    return d.groupby("cdr3")["score"].max().astype(float).to_dict()


def score_imw(locus: str, tb: pl.DataFrame, ns: pl.DataFrame) -> tuple[dict, dict]:
    imw_tb = _imw_top1("sample3")                                   # sample3 holds both TRA+TRB true binders
    imw_ns = _imw_top1(IMW_SAMPLE[locus])
    return ({c: imw_tb[c] for c in tb["cdr3"].to_list() if c in imw_tb},
            {c: imw_ns[c] for c in ns["cdr3"].to_list() if c in imw_ns})


# --------------------------------------------------------------------------- ERGO-II (max P(bind))
def score_ergo(locus: str, tb: pl.DataFrame, ns: pl.DataFrame) -> tuple[dict, dict]:
    """Per-receptor confidence = max P(bind) over the 5 peptides (pretrained VDJdb model). TRA has no
    alpha-only path in ERGO, so the alpha CDR3 (+TRAV/TRAJ) goes through the TRB/TRBV/TRBJ encoder slots
    (per ergo_detect_alpha.py). cmp-ergo env."""
    q = pl.concat([tb.select("cdr3", "v", "j"), ns.select("cdr3", "v", "j")]).unique("cdr3")
    recs = []
    for cdr3, v, j in zip(q["cdr3"], q["v"], q["j"]):
        for pep in PEPTIDES:
            recs.append((cdr3, B.vgene(v), B.vgene(j), pep))
    ergo_in = pl.DataFrame(recs, schema=["TRB", "TRBV", "TRBJ", "Peptide"], orient="row").with_columns(
        TRA=pl.lit(None, dtype=pl.Utf8), TRAV=pl.lit(None, dtype=pl.Utf8),
        TRAJ=pl.lit(None, dtype=pl.Utf8)).with_columns(
        **{"T-Cell-Type": pl.lit("CD8")}, MHC=pl.lit(MHC)).select(
        "TRA", "TRB", "TRAV", "TRAJ", "TRBV", "TRBJ", "T-Cell-Type", "Peptide", "MHC")
    inp = WORK / f"ergo_{locus}_in.csv"; outp = WORK / f"ergo_{locus}_out.csv"
    ergo_in.write_csv(inp)
    subprocess.run(["conda", "run", "-n", "cmp-ergo", "python", "_run_predict.py",
                    "vdjdb", str(inp.resolve()), str(outp.resolve())],
                   check=True, cwd=ERGO, capture_output=True, text=True)
    out = pl.read_csv(outp)
    best = out.group_by("TRB").agg(pl.col("Score").cast(pl.Float64).max().alias("score"))
    by = dict(zip(best["TRB"].to_list(), best["score"].to_list()))
    return ({c: by[c] for c in tb["cdr3"].to_list() if c in by},
            {c: by[c] for c in ns["cdr3"].to_list() if c in by})


# --------------------------------------------------------------------------- NetTCR-2.0 (max score, TRB)
def score_nettcr(tb: pl.DataFrame, ns: pl.DataFrame) -> tuple[dict, dict]:
    """Per-receptor confidence = max NetTCR-2.0 beta score over the 5 peptides (model trained once on
    data/train_beta_90.csv). TRB only; the beta model takes no alpha chain. cmp-nettcr env."""
    q = pl.concat([tb.select("cdr3"), ns.select("cdr3")]).unique("cdr3")
    recs = []
    for cdr3 in q["cdr3"].to_list():
        if len(_trim(cdr3)) > 30:                 # NetTCR-2.0's encoder caps CDR3 at 30 aa; skip longer
            continue
        for pep in PEPTIDES:
            recs.append((_trim(cdr3), pep))
    inp = WORK / "nettcr_TRB_in.csv"; outp = WORK / "nettcr_TRB_out.csv"
    pl.DataFrame(recs, schema=["CDR3b", "peptide"], orient="row").write_csv(inp)
    subprocess.run(["conda", "run", "-n", "cmp-nettcr", "python", "_run_nettcr.py",
                    "--chain", "b", "--train", "data/train_beta_90.csv",
                    "--pred", f"TRB:{inp.resolve()}:{outp.resolve()}"],
                   check=True, cwd=NET, capture_output=True, text=True)
    out = pl.read_csv(outp)
    best = out.group_by("CDR3b").agg(pl.col("prediction").cast(pl.Float64).max().alias("score"))
    by_trim = dict(zip(best["CDR3b"].to_list(), best["score"].to_list()))
    return ({c: by_trim[_trim(c)] for c in tb["cdr3"].to_list() if _trim(c) in by_trim},
            {c: by_trim[_trim(c)] for c in ns["cdr3"].to_list() if _trim(c) in by_trim})


# --------------------------------------------------------------------------- cache + summary
def _cache(tool: str, locus: str, tb_s: dict, ns_s: dict):
    """Persist per-receptor scores so --assemble can rebuild without rerunning the heavy envs."""
    pl.DataFrame({"cdr3": list(tb_s), "score": [tb_s[c] for c in tb_s]}).write_csv(
        OUT / f"{tool}_{locus}_true.tsv", separator="\t")
    pl.DataFrame({"cdr3": list(ns_s), "score": [ns_s[c] for c in ns_s]}).write_csv(
        OUT / f"{tool}_{locus}_noise.tsv", separator="\t")


def _read_cache(tool: str, locus: str):
    ft = OUT / f"{tool}_{locus}_true.tsv"; fn = OUT / f"{tool}_{locus}_noise.tsv"
    if not (ft.exists() and fn.exists()):
        return None
    t = pl.read_csv(ft, separator="\t"); n = pl.read_csv(fn, separator="\t")
    return (np.asarray(t["score"].to_list(), float), np.asarray(n["score"].to_list(), float))


def recall_fp(tb_scores: np.ndarray, ns_scores: np.ndarray):
    """For each recall level r: s* = quantile(true, 1-r); noise FP = fraction noise >= s*."""
    out = []
    for r in RECALLS:
        sstar = float(np.quantile(tb_scores, 1.0 - r))
        fp = float((ns_scores >= sstar).mean() * 100.0)
        out.append((r, sstar, fp))
    return out


# tool -> (label, [loci]). Only the TRAINED classifiers are calibrated on sample3: sample3 is 100% inside
# vdjmatch's vdjdb2026 reference, so the matching tools (vdjmatch/tcrdist3/TCRMatch) are degenerate there
# (exact match -> trivial, exact-excluded -> NED 0). Matching-tool specificity is reported via their
# native thresholds (the E-value / radius) elsewhere; here we fix the trained classifiers' broken default.
TOOLS = [("imw-DETECT", ["TRB", "TRA"]), ("ERGO-II", ["TRB", "TRA"]), ("NetTCR-2.0", ["TRB"])]


def assemble():
    """Read cached scores, compute recall-50/25/10 thresholds + noise FP, write the TSV + noise.dat."""
    rows = []                                                       # tool,locus,recall_level,threshold,noise_fp_pct,n_true,n_noise
    fp50 = {}                                                       # (tool,locus) -> noise_fp_pct at recall 50
    for tool, loci in TOOLS:
        for locus in loci:
            c = _read_cache(tool, locus)
            if c is None:
                print(f"{tool}/{locus}: no cached scores; skipping")
                continue
            tb_s, ns_s = c
            for r, sstar, fp in recall_fp(tb_s, ns_s):
                rows.append((tool, locus, round(r * 100, 1), round(sstar, 6), round(fp, 3),
                             len(tb_s), len(ns_s)))
                if r == 0.50:
                    fp50[(tool, locus)] = round(fp, 3)
            print(f"{tool:11s} {locus}: n_true={len(tb_s)} n_noise={len(ns_s)} | "
                  f"FP @top-{{50,25,5,2.5}}% = " + ", ".join(f"{fp:.2f}" for _, _, fp in recall_fp(tb_s, ns_s)))

    res = pl.DataFrame(rows, schema=["tool", "locus", "top_pct", "threshold",
                                     "noise_fp_pct", "n_true", "n_noise"],
                       orient="row", schema_overrides={"top_pct": pl.Float64})
    res_path = Path("/Users/mikesh/vcs/manuscripts/2026-vdjmatch/benchmarks/results/noise_fp_recall50.tsv")
    res.write_csv(res_path, separator="\t")
    print(f"\nwrote {res_path} (trained classifiers; matching tools reported at native thresholds)")
    with pl.Config(tbl_rows=-1):
        print(res)


SCORERS = {
    "vdjmatch": lambda locus, tb, ns: score_vdjmatch(locus, tb, ns),
    "tcrdist3": lambda locus, tb, ns: score_tcrdist(locus, tb, ns),
    "TCRMatch": lambda locus, tb, ns: score_tcrmatch(tb, ns),
    "imw-DETECT": lambda locus, tb, ns: score_imw(locus, tb, ns),
    "ERGO-II": lambda locus, tb, ns: score_ergo(locus, tb, ns),
    "NetTCR-2.0": lambda locus, tb, ns: score_nettcr(tb, ns),
}


def run_tool(tool: str, loci: list[str]):
    for locus in loci:
        tb, ns = true_binders(locus), noise(locus)
        tb_s, ns_s = SCORERS[tool](locus, tb, ns)
        _cache(tool, locus, tb_s, ns_s)
        print(f"{tool}/{locus}: scored true={len(tb_s)}/{tb.height} noise={len(ns_s)}/{ns.height}")


def main():
    for tool, loci in TOOLS:
        run_tool(tool, loci)
    print("MixTCRpred: SKIPPED (paired-only model; no single-chain noise set exists)")
    assemble()


if __name__ == "__main__":
    if "--assemble" in sys.argv:
        assemble()
    else:
        main()
