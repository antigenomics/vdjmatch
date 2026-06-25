"""Collect the vdjmatch benchmark datasets into one unified schema under ~/hf/airr_benchmark/vdjmatch/.

Schema (tab-separated, gzipped):
  idx  species  cdr3_alpha  v_alpha  j_alpha  cdr3_beta  v_beta  j_beta  epitope  mhc_class  mhc_a  mhc_b  binder

dataset_1  NLV+/NLV- from Egorov et al. (evgeny chunk positives + sample1 tet-negative controls).
dataset_2  ELA/LLW/LLL from Sewell et al. (binder=1) + per-epitope, per-chain Pgen-matched airr_control (binder=0).
dataset_3  vdjdb2026 full, human+mouse, per (epitope,gene) >=30 clonotypes, spike study dropped, no cap (binder=1).
dataset_4  as dataset_3 but the >=2-reference shortlist (binder=1).
dataset_5  TCRvdb (sample6): binder=1 for padj<1e-5 else 0  [RESTRICTED -> gitignored].

    .venv/bin/python bench/build_airr_benchmark.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
import holdout_controls as HC                                       # noqa: E402
from benchmark import shortlist as _shortlist                       # noqa: E402
from compare import TESTDATA                                        # noqa: E402
from vdjmatch import db                                             # noqa: E402

HD = Path.home() / "vcs/manuscripts/2026-vdjmatch/hold_out_data"
OUT = Path.home() / "hf/airr_benchmark/vdjmatch"
OUT.mkdir(parents=True, exist_ok=True)
COLS = ["idx", "species", "cdr3_alpha", "v_alpha", "j_alpha", "cdr3_beta", "v_beta", "j_beta",
        "epitope", "mhc_class", "mhc_a", "mhc_b", "binder"]
EPI = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV", "ELA": "ELAGIGILTV"}
_VALID = r"^C[ACDEFGHIKLMNPQRSTVWY]{3,28}[FW]$"


def _sp(col):
    return (pl.col(col).str.to_lowercase().replace({"homosapiens": "human", "musmusculus": "mouse",
                                                    "macacamulatta": "macaque"}, default=pl.col(col)))


def _strip(col):  # allele/decoration -> gene-level (TRBV20-1*01 -> TRBV20-1)
    return pl.col(col).fill_null("").str.replace(r"[*/].*$", "")


def _chain(cdr3c, vc, jc, suffix):
    """three unified columns for one chain; cleared when the CDR3 is not a valid C...F/W."""
    m = pl.col(cdr3c).fill_null("").str.contains(_VALID)
    return [pl.when(m).then(pl.col(cdr3c)).otherwise(pl.lit("")).alias(f"cdr3_{suffix}"),
            pl.when(m).then(_strip(vc)).otherwise(pl.lit("")).alias(f"v_{suffix}"),
            pl.when(m).then(_strip(jc)).otherwise(pl.lit("")).alias(f"j_{suffix}")]


def finalize(df: pl.DataFrame) -> pl.DataFrame:
    """order columns, drop rows with neither chain, unique, add idx."""
    dc = [c for c in COLS if c != "idx"]
    df = (df.select(dc).filter((pl.col("cdr3_alpha") != "") | (pl.col("cdr3_beta") != ""))
          .unique(subset=dc, maintain_order=True))
    return df.with_row_index("idx").select(COLS)


def _per_gene(df, gene_col, cdr3_col, v_col, j_col):
    """single-chain rows: route (cdr3,v,j) into alpha/beta by gene; the other chain blank."""
    a = pl.col(gene_col) == "TRA"
    return df.with_columns(
        cdr3_alpha=pl.when(a & pl.col(cdr3_col).fill_null("").str.contains(_VALID)).then(pl.col(cdr3_col)).otherwise(pl.lit("")),
        v_alpha=pl.when(a & pl.col(cdr3_col).fill_null("").str.contains(_VALID)).then(_strip(v_col)).otherwise(pl.lit("")),
        j_alpha=pl.when(a & pl.col(cdr3_col).fill_null("").str.contains(_VALID)).then(_strip(j_col)).otherwise(pl.lit("")),
        cdr3_beta=pl.when(~a & pl.col(cdr3_col).fill_null("").str.contains(_VALID)).then(pl.col(cdr3_col)).otherwise(pl.lit("")),
        v_beta=pl.when(~a & pl.col(cdr3_col).fill_null("").str.contains(_VALID)).then(_strip(v_col)).otherwise(pl.lit("")),
        j_beta=pl.when(~a & pl.col(cdr3_col).fill_null("").str.contains(_VALID)).then(_strip(j_col)).otherwise(pl.lit("")),
    )


# ------------------------------------------------------------------ dataset 1: NLV Egorov ------------
def dataset_1():
    ch = pl.read_csv(HD / "_chunk_evgeny.txt", separator="\t", infer_schema_length=0) \
        .filter(pl.col("antigen.epitope") == EPI["NLV"])
    pos = ch.select(species=_sp("species"), *_chain("cdr3.alpha", "v.alpha", "j.alpha", "alpha"),
                    *_chain("cdr3.beta", "v.beta", "j.beta", "beta")) \
        .with_columns(epitope=pl.lit(EPI["NLV"]), mhc_class=pl.lit("MHCI"), mhc_a=pl.lit("HLA-A*02"),
                      mhc_b=pl.lit("B2M"), binder=pl.lit(1, dtype=pl.Int8))
    s1 = (pl.read_csv(TESTDATA / "sample1_cmv_5+reads.txt", separator="\t")
          .filter(pl.col("type") == "control"))
    neg = _per_gene(s1, "gene", "cdr3", "v.segm", "j.segm").select(
        species=pl.lit("human"), cdr3_alpha="cdr3_alpha", v_alpha="v_alpha", j_alpha="j_alpha",
        cdr3_beta="cdr3_beta", v_beta="v_beta", j_beta="j_beta", epitope=pl.lit(EPI["NLV"]),
        mhc_class=pl.lit("MHCI"), mhc_a=pl.lit("HLA-A*02"), mhc_b=pl.lit("B2M"),
        binder=pl.lit(0, dtype=pl.Int8))
    return finalize(pl.concat([pos, neg], how="diagonal_relaxed"))


# ------------------------------------------------------------------ dataset 2: Sewell + Pgen ---------
def dataset_2():
    sw = pl.read_csv(HD / "_chunk_sewell.txt", separator="\t", infer_schema_length=0)
    frames = []
    for sh in ("ELA", "LLW", "LLL"):
        e = EPI[sh]
        d = sw.filter(pl.col("antigen.epitope") == e)
        pos = d.select(species=pl.lit("human"), *_chain("cdr3.alpha", "v.alpha", "j.alpha", "alpha"),
                       *_chain("cdr3.beta", "v.beta", "j.beta", "beta")) \
            .with_columns(epitope=pl.lit(e), mhc_class=pl.lit("MHCI"), mhc_a=pl.lit("HLA-A*02"),
                          mhc_b=pl.lit("B2M"), binder=pl.lit(1, dtype=pl.Int8))
        frames.append(pos)
        for suffix, locus in (("alpha", "TRA"), ("beta", "TRB")):
            cd = pos.filter(pl.col(f"cdr3_{suffix}") != "").select(f"cdr3_{suffix}").unique() \
                [f"cdr3_{suffix}"].to_list()                        # unique positives -> matched neg count
            if not cd:
                continue
            neg = HC.matched_negatives(cd, "human", locus)              # [(cdr3, v, j)] Pgen-matched
            ndf = pl.DataFrame(neg, schema=["c", "v", "j"], orient="row").with_columns(
                species=pl.lit("human"),
                cdr3_alpha=pl.lit("") if suffix == "beta" else pl.col("c"),
                v_alpha=pl.lit("") if suffix == "beta" else _strip("v"),
                j_alpha=pl.lit("") if suffix == "beta" else _strip("j"),
                cdr3_beta=pl.col("c") if suffix == "beta" else pl.lit(""),
                v_beta=_strip("v") if suffix == "beta" else pl.lit(""),
                j_beta=_strip("j") if suffix == "beta" else pl.lit(""),
                epitope=pl.lit(e), mhc_class=pl.lit("MHCI"), mhc_a=pl.lit("HLA-A*02"),
                mhc_b=pl.lit("B2M"), binder=pl.lit(0, dtype=pl.Int8)).drop("c", "v", "j")
            frames.append(ndf)
    return finalize(pl.concat(frames, how="diagonal_relaxed"))


# ------------------------------------------------------------------ dataset 3/4: vdjdb full/shortlist
def _vdjdb(shortlist: bool):
    d = db.load(_bench.source()).filter(pl.col("species").is_in(["HomoSapiens", "MusMusculus"]))
    d = _bench.valid_cdr3(d)
    d = d.filter(~pl.col("reference_id").is_in(list(_bench.spike_studies(d))))   # drop len-14 spike study
    if shortlist:
        d = _shortlist(d, min_refs=2)                              # clonotype-pMHC in >=2 references
    rows = _per_gene(d, "gene", "cdr3", "v", "j").select(
        species=_sp("species"), cdr3_alpha="cdr3_alpha", v_alpha="v_alpha", j_alpha="j_alpha",
        cdr3_beta="cdr3_beta", v_beta="v_beta", j_beta="j_beta", epitope="epitope",
        mhc_class="mhc_class", mhc_a="mhc_a", mhc_b=pl.col("mhc_b").fill_null("B2M"),
        gene="gene", binder=pl.lit(1, dtype=pl.Int8)).unique(
        subset=[c for c in COLS if c != "idx"], maintain_order=True)
    keep = rows.group_by(["epitope", "gene"]).len().filter(pl.col("len") >= 30).select("epitope", "gene")
    return finalize(rows.join(keep, on=["epitope", "gene"], how="semi").drop("gene"))


# ------------------------------------------------------------------ dataset 5: TCRvdb ----------------
def dataset_5():
    t = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv", infer_schema_length=0)
    hla = pl.when(pl.col("hla_long").fill_null("") != "").then(pl.col("hla_long")).otherwise(pl.col("hla_short"))
    rows = t.with_columns(
        species=pl.lit("human"),
        *_chain("cdr3_alpha_aa", "TRAV", "TRAJ", "alpha"), *_chain("cdr3_beta_aa", "TRBV", "TRBJ", "beta"),
        epitope=pl.col("epitope_aa"), mhc_class=pl.lit("MHCI"), mhc_a=hla, mhc_b=pl.lit("B2M"),
        binder=(pl.col("padj").cast(pl.Float64, strict=False) < 1e-5).cast(pl.Int8).fill_null(0))
    return finalize(rows.select([c for c in COLS if c != "idx"]))


def write(name, df):
    import gzip
    import io
    f = OUT / f"{name}.tsv.gz"
    buf = io.BytesIO()
    df.write_csv(buf, separator="\t")
    with gzip.open(f, "wb") as gz:
        gz.write(buf.getvalue())
    print(f"  {name:10} n={df.height:>7}  pos={int((df['binder']==1).sum()):>7}  "
          f"neg={int((df['binder']==0).sum()):>7}  paired={(df.filter((pl.col('cdr3_alpha')!='') & (pl.col('cdr3_beta')!='')).height):>6}  -> {f.name}")


if __name__ == "__main__":
    print("building ~/hf/airr_benchmark/vdjmatch/ :")
    write("dataset_1", dataset_1())
    write("dataset_2", dataset_2())
    write("dataset_3", _vdjdb(shortlist=False))
    write("dataset_4", _vdjdb(shortlist=True))
    write("dataset_5", dataset_5())
