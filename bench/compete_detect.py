"""Run the new competitors (pMTnet, TCRGP, TCR-BERT) on the held-out annotation set and score per-epitope
ROC + IMMREP AUC0.1, on the SAME pooled query set vdjmatch uses (pos = true==E, neg = all other queries).

All three are VDJdb-trained -> our held-out positives are likely in their training data: report with the
leakage caveat (per the user's decision: full naive re-run, leakage noted in text).

    .venv/bin/python bench/compete_detect.py pmtnet TRB
    .venv/bin/python bench/compete_detect.py tcrgp  TRB
    conda run -n cmp-tcrbert python bench/compete_detect.py tcrbert TRB
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
import holdout_eval as HE                                            # noqa: E402
from holdout_features import _jmap                                   # noqa: E402
from metrics import auc01, roc_auc                                   # noqa: E402

CODE = Path.home() / "vcs/code"
RES = Path.home() / "vcs/manuscripts/2026-vdjmatch/benchmarks/results"
PRED = Path(__file__).resolve().parent / "predictions"
A2 = "A*02:01"   # pMTnet matches HLA by key.startswith(input); IMGT keys are like A*02:01:01:01
EPI = HE.EPI                                                         # {sh: peptide}
TCRGP_SHIPPED = {"NLV", "LLW", "GLC", "YLQ"}                         # ELA/LLL have no shipped model


def load_pool(locus):
    """pooled held-out queries: (cdr3, v, j, true_epi_aa_or_None), valid CDR3 only."""
    jm = _jmap(locus)
    out = []
    for c, v, true, ds in HE.query_set(locus):
        if _bench.valid_cdr3(pl.DataFrame({"cdr3": [c]})).height:
            out.append((c, v, jm.get(c, ""), true))
    return out


def write_results(tool, locus, pool, scores):
    """scores: {epitope_aa: list aligned to pool}. Append per-epitope ROC/AUC0.1 to competitors_holdout.tsv."""
    rows = []
    for sh, e in EPI.items():
        sc = scores.get(e)
        if sc is None:
            continue
        y = [int(t == e) for _, _, _, t in pool]
        pairs = list(zip(y, sc))
        rows.append(dict(tool=tool, locus=locus, epitope=sh, n_pos=sum(y),
                         roc=round(roc_auc(pairs), 4), auc01=round(auc01(pairs), 4)))
    df = pl.DataFrame(rows)
    f = RES / "competitors_holdout.tsv"
    if f.exists():
        old = pl.read_csv(f, separator="\t").filter(~((pl.col("tool") == tool) & (pl.col("locus") == locus)))
        df = pl.concat([old, df], how="diagonal")
    df.write_csv(f, separator="\t")
    print(df.filter((pl.col("tool") == tool) & (pl.col("locus") == locus)))


# ---------------------------------------------------------------- pMTnet (TRB only) -------------------
def run_pmtnet(locus, pool):
    assert locus == "TRB", "pMTnet is beta-only"
    work = PRED / "pmtnet"; work.mkdir(parents=True, exist_ok=True)
    inp = work / "input.csv"
    rows = [(c, EPI[sh], A2) for (c, _, _, _) in pool for sh in EPI]
    pl.DataFrame(rows, schema=["CDR3", "Antigen", "HLA"], orient="row").write_csv(inp)
    repo = CODE / "pMTnet"
    subprocess.run(["conda", "run", "-n", "cmp-pmtnet", "python", "pMTnet.py", "-input", str(inp),
                    "-library", "library", "-output", str(work), "-output_log", str(work / "log.txt")],
                   cwd=repo, check=True)
    pred = pl.read_csv(work / "prediction.csv")
    rank = {(c, a): r for c, a, r in zip(pred["CDR3"], pred["Antigen"], pred["Rank"])}
    scores = {}
    for sh, e in EPI.items():
        scores[e] = [1.0 - rank.get((c, e), 1.0) for (c, _, _, _) in pool]   # lower Rank = binder
    return scores


# ---------------------------------------------------------------- TCRGP (beta, shipped epitopes) ------
def run_tcrgp(locus, pool):
    assert locus == "TRB", "TCRGP shipped models are beta-only"
    import numpy as np
    repo = CODE / "TCRGP"
    work = PRED / "tcrgp"; work.mkdir(parents=True, exist_ok=True)
    inp = work / "input.csv"
    pl.DataFrame({"cdr3b": [c for c, _, _, _ in pool]}).write_csv(inp)
    scores = {}
    for sh in TCRGP_SHIPPED:
        e = EPI[sh]
        out = work / f"{sh}.npy"
        subprocess.run(["conda", "run", "-n", "cmp-tcrgp", "python",
                        str(Path(__file__).resolve().parent / "_tcrgp_predict.py"), str(inp),
                        str(repo / "models" / "paper" / f"model_vdj_{e}_cdr3b"), str(out)], check=True)
        scores[e] = np.load(out).tolist()
    return scores


# ---------------------------------------------------------------- TCR-BERT (embeddings + kNN density) -
def _bert_embed(seqs, tag):
    """embed via the cmp-tcrbert env (torch); cache npy under predictions/tcrbert/."""
    import numpy as np
    work = PRED / "tcrbert"; work.mkdir(parents=True, exist_ok=True)
    lst, out = work / f"{tag}.txt", work / f"{tag}.npy"
    lst.write_text("\n".join(seqs))
    subprocess.run(["conda", "run", "-n", "cmp-tcrbert", "python", str(Path(__file__).resolve().parent /
                    "_tcrbert_embed.py"), str(lst), str(out)], check=True)
    X = np.load(out)
    return X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-9)


def run_tcrbert(locus, pool):
    import numpy as np
    from benchmark import A02, release                              # noqa: E402
    from holdout_controls import _airr
    qX = _bert_embed([c for c, _, _, _ in pool], f"query_{locus}")
    rdf = _bench.valid_cdr3(release("vdjdb2026").filter(pl.col("mhc_a").str.contains(A02) & (pl.col("gene") == locus)))
    bgX = _bert_embed(_airr("human", locus, 4000)["cdr3"].to_list(), f"bg_{locus}")
    scores = {}
    for sh, e in EPI.items():
        ref = rdf.filter(pl.col("epitope") == e).unique("cdr3")["cdr3"].to_list()
        if len(ref) < 20:
            continue
        rX = _bert_embed(ref[:4000], f"ref_{sh}_{locus}")
        # density: mean top-10 cosine sim to E binders minus to background (naive; reference leakage caveat)
        kpos = np.sort(qX @ rX.T, axis=1)[:, -10:].mean(1)
        kbg = np.sort(qX @ bgX.T, axis=1)[:, -10:].mean(1)
        scores[e] = (kpos - kbg).tolist()
    return scores


if __name__ == "__main__":
    tool, locus = sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "TRB"
    pool = load_pool(locus)
    fn = {"pmtnet": run_pmtnet, "tcrgp": run_tcrgp, "tcrbert": run_tcrbert}[tool]
    write_results(tool, locus, pool, fn(locus, pool))
