"""Assemble the clean held-out benchmark datasets into ~/vcs/manuscripts/2026-vdjmatch/hold_out_data/.

Per (epitope, locus) a TSV with columns: cdr3, v, j, label (1=binder, 0=control), source.
Positives:
  NLV (TRA/TRB)         : evgeny-etal-2018-04 chunk (vdjdb-db #252), HLA-A*02:01 NLVPMVATV
  LLW/ELA (TRB), LLL (TRB+TRA): sewell-etal-2017-04-13 chunk (vdjdb-db #193), HLA-A*02
  GLC/YLQ (TRA/TRB)     : sample6 TCRvdb, padj < 1e-5
Negatives (here): NLV = sample1 tet-negative control; GLC/YLQ = sample6 padj >= 1e-5.
LLW/LLL/ELA negatives are Pgen-matched airr_control, added by mock_negatives_holdout.py.

    .venv/bin/python bench/build_holdout.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
from benchmark import vgene                                         # noqa: E402
from compare import TESTDATA                                        # noqa: E402

HD = Path.home() / "vcs/manuscripts/2026-vdjmatch/hold_out_data"
EPI = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV",
       "ELA": "ELAGIGILTV", "YLQ": "YLQPRTFLL", "GLC": "GLCTLVAML"}


def _clean(df, cd, v, j):
    """select cdr3/v/j (allele-stripped), drop invalid, unique by cdr3."""
    d = (df.select(cdr3=cd, v=v, j=j).filter(pl.col("cdr3").is_not_null() & (pl.col("cdr3") != ""))
         .pipe(_bench.valid_cdr3).unique("cdr3"))
    return d.with_columns(v=pl.col("v").map_elements(vgene, return_dtype=pl.Utf8),
                          j=pl.col("j").map_elements(vgene, return_dtype=pl.Utf8))


def chunk_pos(chunk, epi, locus):
    d = pl.read_csv(HD / chunk, separator="\t", infer_schema_length=0).filter(pl.col("antigen.epitope") == epi)
    cd, v, j = (("cdr3.beta", "v.beta", "j.beta") if locus == "TRB"
                else ("cdr3.alpha", "v.alpha", "j.alpha"))
    return _clean(d, cd, v, j).with_columns(label=pl.lit(1, dtype=pl.Int8), source=pl.lit("vdjdb_binder"))


def sample1_neg(locus):
    d = (pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
         .filter((pl.col("gene") == locus) & (pl.col("type") == "control")))
    return _clean(d, "cdr3", "v.segm", "j.segm").with_columns(label=pl.lit(0, dtype=pl.Int8),
                                                              source=pl.lit("sample1_tet_neg"))


def tcrvdb(epi, locus):
    t = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv").with_columns(pos=pl.col("padj") < 1e-5)
    t = t.filter(pl.col("epitope_aa") == epi)
    cd, v, j = (("cdr3_beta_aa", "TRBV", "TRBJ") if locus == "TRB" else ("cdr3_alpha_aa", "TRAV", "TRAJ"))
    pos = _clean(t.filter("pos"), cd, v, j).with_columns(label=pl.lit(1, dtype=pl.Int8), source=pl.lit("tcrvdb_pos"))
    neg = _clean(t.filter(~pl.col("pos")), cd, v, j).with_columns(label=pl.lit(0, dtype=pl.Int8), source=pl.lit("tcrvdb_neg"))
    return pl.concat([pos, neg]).unique("cdr3", keep="first")


def write(name, df):
    df = df.select("cdr3", "v", "j", "label", "source")
    df.write_csv(HD / f"{name}.tsv", separator="\t")
    p, n = int((df["label"] == 1).sum()), int((df["label"] == 0).sum())
    print(f"  {name:12} pos={p:5} neg={n:5} -> {name}.tsv")


def main():
    print("building hold_out_data (positives + NLV/GLC/YLQ negatives; LLW/LLL/ELA neg = Pgen, next step):")
    for loc in ("TRA", "TRB"):
        write(f"NLV_{loc}", pl.concat([chunk_pos("_chunk_evgeny.txt", EPI["NLV"], loc), sample1_neg(loc)]))
    write("LLW_TRB", chunk_pos("_chunk_sewell.txt", EPI["LLW"], "TRB"))      # neg added later
    write("ELA_TRB", chunk_pos("_chunk_sewell.txt", EPI["ELA"], "TRB"))
    write("LLL_TRB", chunk_pos("_chunk_sewell.txt", EPI["LLL"], "TRB"))
    write("LLL_TRA", chunk_pos("_chunk_sewell.txt", EPI["LLL"], "TRA"))
    for loc in ("TRA", "TRB"):
        write(f"GLC_{loc}", tcrvdb(EPI["GLC"], loc))
        write(f"YLQ_{loc}", tcrvdb(EPI["YLQ"], loc))


if __name__ == "__main__":
    main()
