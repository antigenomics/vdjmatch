"""Collect the vdjmatch benchmark datasets into one unified schema under ~/hf/airr_benchmark/vdjmatch/.

Schema (tab-separated, gzipped):
  idx  species  cdr3_alpha  v_alpha  j_alpha  cdr3_beta  v_beta  j_beta  epitope  mhc_class  mhc_a  mhc_b  binder

dataset_1  NLV+/NLV- from Egorov et al. (evgeny chunk positives + sample1 tet-negative controls).
dataset_2  ELA/LLW/LLL from Sewell et al. (binder=1) + per-epitope, per-chain Pgen-matched negatives (binder=0).
dataset_3  vdjdb2026 full, human+mouse: single-chain (per epitope,gene >=30) and paired (complex_id, per
           epitope >=30) binders + Pgen-matched negatives. No cap; the length-14 spike study dropped.
dataset_4  as dataset_3 but the >=2-reference shortlist.

Negatives are REAL post-selection T-cell repertoire sequences (isalgo/airr_control, ~/hf/airr_control),
NOT generated -- OLGA is used only to compute the generation probability Pgen for matching. Per (epitope,
gene) the negatives match the positives' count and round(log2 Pgen) histogram; for paired records each chain
is matched on its own log2 Pgen (so log2 Pgen(alpha)+log2 Pgen(beta) is matched jointly). All CDR3 are the
20 standard amino acids, C...F/W, no * stops or _ frameshifts.
dataset_5  TCRvdb (sample6): binder=1 for padj<1e-5 else 0  [RESTRICTED -> gitignored].

    .venv/bin/python bench/build_airr_benchmark.py
"""
from __future__ import annotations

import multiprocessing as mp
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
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
    """clear any chain whose CDR3 is not 20-aa C...F/W (no */_), drop empty rows, unique, add idx."""
    dc = [c for c in COLS if c != "idx"]
    df = df.select(dc)
    for suf in ("alpha", "beta"):
        bad = (pl.col(f"cdr3_{suf}") != "") & ~pl.col(f"cdr3_{suf}").str.contains(_VALID)
        df = df.with_columns([pl.when(bad).then(pl.lit("")).otherwise(pl.col(f"{p}_{suf}")).alias(f"{p}_{suf}")
                              for p in ("cdr3", "v", "j")])
    df = (df.filter((pl.col("cdr3_alpha") != "") | (pl.col("cdr3_beta") != ""))
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


# ------------------------------------------------------------------ dataset 3/4: vdjdb + Pgen negatives
def _strip_s(x):
    return (x or "").split("*")[0].split("/")[0]


def _bins(cdr3s, species, locus):
    """parallel round(log2 Pgen) per CDR3 (OLGA model for species/locus) -> {cdr3: bin}. airr_control is
    REAL repertoire data; OLGA only computes Pgen here."""
    cdr3s = list(dict.fromkeys(cdr3s))
    if not cdr3s:
        return {}
    with mp.Pool(max(1, (os.cpu_count() or 2) - 1), initializer=HC._winit,
                 initargs=(HC._SUB[(species, locus)],)) as p:
        b = p.map(HC._wbin, cdr3s, chunksize=300)
    return {c: bb for c, bb in zip(cdr3s, b) if bb is not None}


def _bybin(species, locus):
    bb = defaultdict(list)
    for b, c, v, j in HC.pgen_pool(species, locus).iter_rows():
        bb[b].append((c, _strip_s(v), _strip_s(j)))
    return bb, sorted(bb)


def _draw(bb, avail, hist, rng):
    out = []
    for b, cnt in sorted(hist.items()):
        src = b if bb.get(b) else (min(avail, key=lambda x: abs(x - b)) if avail else None)
        if src is None:
            continue
        recs = bb[src]
        out += [recs[int(i)] for i in rng.choice(len(recs), cnt, replace=len(recs) < cnt)]
    return out


def _neg_df(records, species, epi, m, suffix):
    if not records:
        return None
    o = suffix == "alpha"
    return pl.DataFrame(records, schema=["c", "v", "j"], orient="row").with_columns(
        species=pl.lit(species),
        cdr3_alpha=pl.col("c") if o else pl.lit(""), v_alpha=pl.col("v") if o else pl.lit(""),
        j_alpha=pl.col("j") if o else pl.lit(""), cdr3_beta=pl.lit("") if o else pl.col("c"),
        v_beta=pl.lit("") if o else pl.col("v"), j_beta=pl.lit("") if o else pl.col("j"),
        epitope=pl.lit(epi), mhc_class=pl.lit(m["mhc_class"]), mhc_a=pl.lit(m["mhc_a"]),
        mhc_b=pl.lit(m["mhc_b"]), binder=pl.lit(0, dtype=pl.Int8)).drop("c", "v", "j")


def _negatives(sc, pair, seed=42):
    """single-chain neg per (species,gene,epitope) matching log2-Pgen; paired neg per (species,epitope)
    matching each chain's log2-Pgen (so log2Pgen_alpha+log2Pgen_beta is matched jointly)."""
    frames = []
    for species in sc["species"].unique().to_list():
        for suf, locus in (("alpha", "TRA"), ("beta", "TRB")):
            sub = sc.filter((pl.col("species") == species) & (pl.col(f"cdr3_{suf}") != ""))
            if not sub.height:
                continue
            bm = _bins(sub[f"cdr3_{suf}"].to_list(), species, locus)
            bb, avail = _bybin(species, locus)
            rng = np.random.default_rng(seed)
            for (epi,), g in sub.group_by("epitope"):
                cd = g.select(f"cdr3_{suf}").unique()[f"cdr3_{suf}"].to_list()
                frames.append(_neg_df(_draw(bb, avail, Counter(bm[c] for c in cd if c in bm), rng),
                                      species, epi, g.row(0, named=True), suf))
    for species in pair["species"].unique().to_list():
        sub = pair.filter(pl.col("species") == species)
        bmA = _bins(sub["cdr3_alpha"].to_list(), species, "TRA")
        bmB = _bins(sub["cdr3_beta"].to_list(), species, "TRB")
        bbA, avA = _bybin(species, "TRA"); bbB, avB = _bybin(species, "TRB")
        rng = np.random.default_rng(seed)
        for (epi,), g in sub.group_by("epitope"):
            grp = Counter((bmA.get(r["cdr3_alpha"]), bmB.get(r["cdr3_beta"])) for r in g.iter_rows(named=True))
            rows = []
            for (ba, bbn), cnt in grp.items():
                if ba is None or bbn is None:
                    continue
                for a, b in zip(_draw(bbA, avA, {ba: cnt}, rng), _draw(bbB, avB, {bbn: cnt}, rng)):
                    rows.append((a[0], a[1], a[2], b[0], b[1], b[2]))
            if rows:
                m = g.row(0, named=True)
                frames.append(pl.DataFrame(rows, orient="row", schema=["cdr3_alpha", "v_alpha", "j_alpha",
                              "cdr3_beta", "v_beta", "j_beta"]).with_columns(
                    species=pl.lit(species), epitope=pl.lit(epi), mhc_class=pl.lit(m["mhc_class"]),
                    mhc_a=pl.lit(m["mhc_a"]), mhc_b=pl.lit(m["mhc_b"]), binder=pl.lit(0, dtype=pl.Int8)))
    return pl.concat([f for f in frames if f is not None], how="diagonal_relaxed")


def _vdjdb(shortlist: bool):
    d = db.load(_bench.source()).filter(pl.col("species").is_in(["HomoSapiens", "MusMusculus"]))
    d = _bench.valid_cdr3(d).filter(~pl.col("reference_id").is_in(list(_bench.spike_studies(d))))
    if shortlist:
        d = _shortlist(d, min_refs=2)                              # clonotype-pMHC in >=2 references
    d = d.with_columns(species=_sp("species"), mhc_b=pl.col("mhc_b").fill_null("B2M"))
    sc = (_per_gene(d, "gene", "cdr3", "v", "j")
          .select("species", "cdr3_alpha", "v_alpha", "j_alpha", "cdr3_beta", "v_beta", "j_beta",
                  "epitope", "mhc_class", "mhc_a", "mhc_b", "gene").unique())
    sc = sc.join(sc.group_by(["epitope", "gene"]).len().filter(pl.col("len") >= 30).select("epitope", "gene"),
                 on=["epitope", "gene"], how="semi").drop("gene")
    pr = d.filter(pl.col("complex_id") != 0)
    a = pr.filter(pl.col("gene") == "TRA").select("complex_id", "species", "epitope", "mhc_class",
                                                  "mhc_a", "mhc_b", ca="cdr3", va="v", ja="j")
    b = pr.filter(pl.col("gene") == "TRB").select("complex_id", cb="cdr3", vb="v", jb="j")
    pair = a.join(b, on="complex_id").select(
        "species", "epitope", "mhc_class", "mhc_a", "mhc_b",
        cdr3_alpha="ca", v_alpha=_strip("va"), j_alpha=_strip("ja"),
        cdr3_beta="cb", v_beta=_strip("vb"), j_beta=_strip("jb")).unique()
    pair = pair.join(pair.group_by("epitope").len().filter(pl.col("len") >= 30).select("epitope"),
                     on="epitope", how="semi")
    pos = pl.concat([sc.with_columns(binder=pl.lit(1, pl.Int8)),
                     pair.with_columns(binder=pl.lit(1, pl.Int8))], how="diagonal_relaxed")
    return finalize(pl.concat([pos, _negatives(sc, pair)], how="diagonal_relaxed"))


# ------------------------------------------------------------------ dataset 5: TCRvdb ----------------
def dataset_5():
    t = (pl.read_csv(TESTDATA / "sample6_TCRvdb.csv", infer_schema_length=0)
         .filter(pl.col("epitope_aa").is_in(["YLQPRTFLL", "GLCTLVAML"])))   # YLQ + GLC only
    hla = pl.when(pl.col("hla_long").fill_null("") != "").then(pl.col("hla_long")).otherwise(pl.col("hla_short"))
    rows = t.with_columns(
        species=pl.lit("human"),
        *_chain("cdr3_alpha_aa", "TRAV", "TRAJ", "alpha"), *_chain("cdr3_beta_aa", "TRBV", "TRBJ", "beta"),
        epitope=pl.col("epitope_aa"), mhc_class=pl.lit("MHCI"), mhc_a=hla, mhc_b=pl.lit("B2M"),
        binder=(pl.col("padj").cast(pl.Float64, strict=False) < 1e-5).cast(pl.Int8).fill_null(0))
    return finalize(rows.select([c for c in COLS if c != "idx"]))


def _gz(df, path):
    """gzipped TSV; empty strings -> truly empty fields (no quotes)."""
    import gzip
    import io
    strc = [c for c in df.columns if df.schema[c] == pl.Utf8]
    dfw = df.with_columns([pl.when(pl.col(c) == "").then(None).otherwise(pl.col(c)).alias(c) for c in strc])
    buf = io.BytesIO()
    dfw.write_csv(buf, separator="\t")
    with gzip.open(path, "wb") as gz:
        gz.write(buf.getvalue())


def write(name, df):
    f = OUT / f"{name}.tsv.gz"
    _gz(df, f)
    print(f"  {name:10} n={df.height:>7}  pos={int((df['binder']==1).sum()):>7}  "
          f"neg={int((df['binder']==0).sum()):>7}  paired={(df.filter((pl.col('cdr3_alpha')!='') & (pl.col('cdr3_beta')!='')).height):>6}  -> {f.name}")


# AIRR-rearrangement reshape -> airr/. labels kept on every row so each file is benchmark-ready.
_LAB = ["species", "epitope", "mhc_class", "mhc_a", "mhc_b", "binder"]


def _airr_parts(df):
    """-> (tra, trb, al, be): single-chain TRA/TRB rows and the alpha/beta rows of paired clones (the
    latter carry clone_id = idx). All in AIRR junction_aa/v_call/j_call + labels schema."""
    a, b = pl.col("cdr3_alpha") != "", pl.col("cdr3_beta") != ""
    tra = (df.filter(a & ~b).select("cdr3_alpha", "v_alpha", "j_alpha", *_LAB)
           .rename({"cdr3_alpha": "junction_aa", "v_alpha": "v_call", "j_alpha": "j_call"}))
    trb = (df.filter(b & ~a).select("cdr3_beta", "v_beta", "j_beta", *_LAB)
           .rename({"cdr3_beta": "junction_aa", "v_beta": "v_call", "j_beta": "j_call"}))
    pr = df.filter(a & b)
    al = pr.select("idx", "cdr3_alpha", "v_alpha", "j_alpha", *_LAB).rename(
        {"idx": "clone_id", "cdr3_alpha": "junction_aa", "v_alpha": "v_call", "j_alpha": "j_call"})
    be = pr.select("idx", "cdr3_beta", "v_beta", "j_beta", *_LAB).rename(
        {"idx": "clone_id", "cdr3_beta": "junction_aa", "v_beta": "v_call", "j_beta": "j_call"})
    return tra, trb, al, be


def write_airr(name, df):
    """single-chain -> airr/{name}_TRA|TRB.tsv.gz; paired -> airr/{name}_paired.tsv.gz, one AIRR row per
    chain (alpha then beta) linked by clone_id (= the paired row's idx)."""
    out = OUT / "airr"
    out.mkdir(exist_ok=True)
    tra, trb, al, be = _airr_parts(df)
    if tra.height:
        _gz(tra, out / f"{name}_TRA.tsv.gz")
    if trb.height:
        _gz(trb, out / f"{name}_TRB.tsv.gz")
    if al.height:
        _gz(pl.concat([al, be]).sort("clone_id"), out / f"{name}_paired.tsv.gz")
    print(f"  airr/{name}: TRA={tra.height} TRB={trb.height} paired={al.height}(x2 rows)")


def write_airr2(name, df):
    """as airr/ but paired clones are split by chain -> airr2/{name}_paired_TRA.tsv.gz and
    {name}_paired_TRB.tsv.gz (rejoin on clone_id)."""
    out = OUT / "airr2"
    out.mkdir(exist_ok=True)
    tra, trb, al, be = _airr_parts(df)
    if tra.height:
        _gz(tra, out / f"{name}_TRA.tsv.gz")
    if trb.height:
        _gz(trb, out / f"{name}_TRB.tsv.gz")
    if al.height:
        _gz(al, out / f"{name}_paired_TRA.tsv.gz")
        _gz(be, out / f"{name}_paired_TRB.tsv.gz")
    print(f"  airr2/{name}: TRA={tra.height} TRB={trb.height} paired_TRA={al.height} paired_TRB={be.height}")


def breakdown(name, df):
    """flat rows: dataset, species, chain (TRA/TRB/paired), epitope, mhc_class, mhc_a, mhc_b, binder, n."""
    chain = (pl.when((pl.col("cdr3_alpha") != "") & (pl.col("cdr3_beta") != "")).then(pl.lit("paired"))
             .when(pl.col("cdr3_alpha") != "").then(pl.lit("TRA")).otherwise(pl.lit("TRB")))
    g = ["species", "chain", "epitope", "mhc_class", "mhc_a", "mhc_b", "binder"]
    return (df.with_columns(chain=chain).group_by(g).len().rename({"len": "n"})
            .with_columns(dataset=pl.lit(name)).select("dataset", *g, "n")
            .sort(["species", "chain", "epitope", "mhc_a", "mhc_b", "binder"],
                  descending=[False, False, False, False, False, True]))


def vdjdb_info(name, df):
    """flat rows: dataset, species, records, epitopes, mhc_a, MHCI, MHCII -- one row per (dataset, species)."""
    return (df.group_by("species").agg(records=pl.len(), epitopes=pl.col("epitope").n_unique(),
                                       mhc_a=pl.col("mhc_a").n_unique(),
                                       MHCI=(pl.col("mhc_class") == "MHCI").sum(),
                                       MHCII=(pl.col("mhc_class") == "MHCII").sum())
            .with_columns(dataset=pl.lit(name))
            .select("dataset", "species", "records", "epitopes", "mhc_a", "MHCI", "MHCII").sort("species"))


if __name__ == "__main__":
    print("building ~/hf/airr_benchmark/vdjmatch/ :")
    d1 = dataset_1(); write("dataset_1", d1)
    d2 = dataset_2(); write("dataset_2", d2)
    d3 = _vdjdb(shortlist=False); write("dataset_3", d3)
    d4 = _vdjdb(shortlist=True); write("dataset_4", d4)
    d5 = dataset_5(); write("dataset_5", d5)
    print("airr/ + airr2/ reshape:")
    for nm, d in (("dataset_1", d1), ("dataset_2", d2), ("dataset_3", d3), ("dataset_4", d4), ("dataset_5", d5)):
        write_airr(nm, d); write_airr2(nm, d)
    bd = pl.concat([breakdown(nm, d) for nm, d in (("dataset_1", d1), ("dataset_2", d2), ("dataset_5", d5))])
    bd.write_csv(OUT / "breakdown.tsv", separator="\t")
    bdv = pl.concat([breakdown(nm, d) for nm, d in (("dataset_3", d3), ("dataset_4", d4))])
    bdv.write_csv(OUT / "breakdown_vdjdb.tsv", separator="\t")
    rs = pl.concat([vdjdb_info(nm, d) for nm, d in (("dataset_3", d3), ("dataset_4", d4))])
    rs.write_csv(OUT / "reference_summary.tsv", separator="\t")
    print(f"\n== breakdown.tsv ({bd.height} rows, d1/d2/d5) =="); print(bd)
    print(f"\n== breakdown_vdjdb.tsv ({bdv.height} rows, d3/d4) ==")
    print("\n== reference_summary.tsv =="); print(rs)
