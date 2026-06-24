"""Feature-probe harness for the hard detection tasks (GLC, LLL, YLQ, NLV, LLW).

Loads enriched per-query features (cdr3, v, j, length, and for sample6 the paired alpha chain) from the
manuscript test_data AT RUNTIME (never copied into the repo) plus the epitope's A*02 VDJdb2026 reference
(cdr3, v, j), and the baseline vdjmatch caldens score. For GLC/LLL feature discovery only — does NOT
change any committed scoring. All TRB unless ``locus='TRA'`` (sample6 only).

    from _feat_probe import task_table, ref_table, baseline_scores
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench
import benchmark as B
from compare import TESTDATA
from vdjmatch import db
from vdjmatch.evalue import background

E = B.EPI                                                          # task -> epitope aa
A02 = B.A02
_REF = {}


def _strip(v):
    return B.vgene(v)


def task_table(task: str, locus: str = "TRB") -> pl.DataFrame:
    """Per-query feature frame for a detection task: columns query_id,label(1/0),cdr3,v,j,length and —
    for sample6 (YLQ/GLC) — paired a_cdr3,a_v,a_j. ``v``/``j`` are allele-stripped gene names."""
    if task == "NLV":
        s = (pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
             .filter(pl.col("gene") == "TRB").select("cdr3", lab="type", v="v.segm", j="j.segm"))
        s = _bench.valid_cdr3(s).unique("cdr3")
        d = s.with_columns(label=(pl.col("lab") == "cmv").cast(pl.Int8))
    elif task in ("LLW", "LLL"):
        s = (pl.read_csv(TESTDATA / "sample2_yf_bst2_5+reads.txt", separator="\t")
             .rename({"antigen.epitope": "lab"}).select("cdr3", "lab", v="v.segm", j="j.segm"))
        s = _bench.valid_cdr3(s).unique("cdr3")
        d = s.with_columns(label=(pl.col("lab") == E[task]).cast(pl.Int8))
    elif task in ("YLQ", "GLC"):
        t = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv").with_columns(pos=pl.col("padj") < 1e-5)
        t = t.filter(pl.col("epitope_aa") == E[task])
        cc, vc, jc = ("cdr3_beta_aa", "TRBV", "TRBJ") if locus == "TRB" else ("cdr3_alpha_aa", "TRAV", "TRAJ")
        ac, av, aj = ("cdr3_alpha_aa", "TRAV", "TRAJ") if locus == "TRB" else ("cdr3_beta_aa", "TRBV", "TRBJ")
        d = (t.select(cdr3=cc, v=vc, j=jc, a_cdr3=ac, a_v=av, a_j=aj, pos="pos")
             .pipe(_bench.valid_cdr3).unique("cdr3")
             .with_columns(label=pl.col("pos").cast(pl.Int8)))
    else:
        raise ValueError(task)
    d = d.with_columns(label=pl.col("label").fill_null(0),                # null padj -> negative
                       v=pl.col("v").map_elements(_strip, return_dtype=pl.Utf8),
                       j=pl.col("j").map_elements(_strip, return_dtype=pl.Utf8),
                       length=pl.col("cdr3").str.len_chars())
    if "a_v" in d.columns:
        d = d.with_columns(a_v=pl.col("a_v").map_elements(_strip, return_dtype=pl.Utf8),
                           a_j=pl.col("a_j").map_elements(_strip, return_dtype=pl.Utf8))
    return d.with_columns(query_id=pl.col("cdr3")).drop([c for c in ("lab", "pos") if c in d.columns])


def ref_table(task: str, locus: str = "TRB") -> pl.DataFrame:
    """Epitope's A*02 VDJdb2026 reference (cdr3, v, j allele-stripped), one row per unique cdr3."""
    key = (task, locus)
    if key not in _REF:
        v26 = db.load(_bench.source(), species="HomoSapiens")
        r = (_bench.valid_cdr3(v26.filter((pl.col("epitope") == E[task]) & pl.col("mhc_a").str.contains(A02)
                                          & (pl.col("gene") == locus)))
             .group_by("cdr3").agg(pl.col("v").first(), pl.col("j").first()))
        _REF[key] = r.with_columns(v=pl.col("v").map_elements(_strip, return_dtype=pl.Utf8),
                                   j=pl.col("j").map_elements(_strip, return_dtype=pl.Utf8))
    return _REF[key]


def baseline_scores(task: str, locus: str = "TRB") -> dict:
    """query cdr3 -> vdjmatch caldens hybrid score (the committed detection score), exact removed."""
    d = task_table(task, locus)
    ref = B.epi_ref(E[task], locus)
    tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v = B.ref_index(ref, locus)
    qv = {c: v for c, v in zip(d["cdr3"], d["v"])}
    allq = d["cdr3"].to_list()
    sc, _ = B.vdjmatch_classify(tgt, ref_epi, ref_v, n_epi, n_epi_v, n_v, background(locus),
                                allq, [qv[q] for q in allq], [E[task]], 1e-3, True,
                                params=B.first_hit.scope(5, 2, 2))
    return {q: sc[q][E[task]][0] for q in allq}


if __name__ == "__main__":
    for tk in ("GLC", "LLL", "YLQ"):
        d = task_table(tk)
        print(f"{tk}: n={d.height} pos={d['label'].sum()} neg={d.height - d['label'].sum()} "
              f"cols={d.columns} | ref={ref_table(tk).height}")
