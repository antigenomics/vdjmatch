"""Do CDR3 length and J-gene carry leakage-safe signal beyond V?

J-gene usage (like V) is partly germline/protocol-driven; CDR3 length is a coarse but robust feature. Per
epitope we build length and J-gene log-odds vs airr_control and score the held-out queries, reporting ROC
full and cross-study (source study removed), next to V. Also per-(epitope,chain) dominant-J fraction and
length spread, to extend the chain-signal predictor (Step 3).

    .venv/bin/python bench/holdout_lenj.py TRB
"""
from __future__ import annotations

import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
import holdout_audit as HA                                          # noqa: E402
import holdout_eval as HE                                           # noqa: E402
from benchmark import A02, release, vgene                           # noqa: E402
from hardcase import vj_logodds                                     # noqa: E402
from holdout_controls import _airr                                  # noqa: E402
from metrics import roc_auc                                         # noqa: E402


def jmap(locus):
    """cdr3 -> J gene from the held-out TSVs (the query cache has no J)."""
    m = {}
    for name in HE.DATASETS[locus]:
        f = HE.HD / f"{name}.tsv"
        if f.exists():
            d = pl.read_csv(f, separator="\t", infer_schema_length=0)
            for c, j in zip(d["cdr3"], d["j"]):
                m.setdefault(c, vgene(j) if j else "")
    return m


def len_logodds(ref_len, bg_len):
    r, b = Counter(ref_len), Counter(bg_len)
    nr, nb = sum(r.values()) + len(r), sum(b.values()) + len(b)
    keys = set(r) | set(b)
    return {L: math.log(((r[L] + 1) / nr) / ((b[L] + 1) / nb)) for L in keys}


def main(locus):
    cache = HE.build(locus, "full")
    recs = cache["recs"]
    jm = jmap(locus)
    rdf = _bench.valid_cdr3(release("vdjdb2026").filter(pl.col("mhc_a").str.contains(A02)
                                                        & (pl.col("gene") == locus)))
    bgd = _airr("human", locus, 60000)
    bg_v, bg_j, bg_len = bgd["v"].to_list(), bgd["j"].to_list(), [len(c) for c in bgd["cdr3"]]
    print(f"\n=== {locus}: V / J / length channel ROC (full, and cross-study where in vdjdb) ===")
    print(f"{'epi':5}{'domJ':>6}{'Llen':>6}{'V':>7}{'J':>7}{'len':>7}{'J_xs':>7}{'len_xs':>8}")
    for sh in HE.EPI:
        e = HE.EPI[sh]
        sub = rdf.filter(pl.col("epitope") == e).group_by("cdr3").agg(pl.col("v").first(), pl.col("j").first())
        if sub.height < 20 or not any(r["true"] == e for r in recs):
            continue
        rj = [vgene(x) for x in sub["j"]]
        rl = [len(c) for c in sub["cdr3"]]
        Jlo = vj_logodds(rj, bg_j); Llo = len_logodds(rl, bg_len)
        y = [int(r["true"] == e) for r in recs]
        sJ = [Jlo.get(jm.get(r["cdr3"], ""), 0.0) for r in recs]
        sL = [Llo.get(len(r["cdr3"]), 0.0) for r in recs]
        sV = [vj_logodds([vgene(x) for x in sub["v"]], bg_v).get(r["v"], 0.0) for r in recs]
        domJ = max(Counter(rj).values()) / len(rj)
        # cross-study J / length
        jX = lX = float("nan")
        if sh in HA.SRC:
            sid = HA.SRC[sh].split("/")[-1]
            kept = rdf.filter((pl.col("epitope") == e)
                              & ~pl.col("reference_id").str.contains(sid, literal=True))
            kk = kept.group_by("cdr3").agg(pl.col("j").first())
            if kk.height >= 20:
                JloX = vj_logodds([vgene(x) for x in kk["j"]], bg_j)
                LloX = len_logodds([len(c) for c in kk["cdr3"]], bg_len)
                jX = roc_auc(list(zip(y, [JloX.get(jm.get(r["cdr3"], ""), 0.0) for r in recs])))
                lX = roc_auc(list(zip(y, [LloX.get(len(r["cdr3"]), 0.0) for r in recs])))
        print(f"{sh:5}{domJ:>6.2f}{np.std(rl):>6.1f}{roc_auc(list(zip(y, sV))):>7.3f}"
              f"{roc_auc(list(zip(y, sJ))):>7.3f}{roc_auc(list(zip(y, sL))):>7.3f}{jX:>7.3f}{lX:>8.3f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TRB")
