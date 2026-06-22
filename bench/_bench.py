"""Shared benchmark helpers: canonical VDJdb source + CDR3 sanitation.

The standing benchmark is the **VDJdb 2026-06-11-ZENODO** release, mirrored gzipped on the
``isalgo/airr_benchmark`` HuggingFace dataset and fetched via ``db.fetch_hf`` (reproducible; GitHub
release assets can change). Override with ``$VDJDB_SAMPLE`` (a local export path) or ``--pmhc`` for
ad-hoc runs. ``vdjdb.txt`` (asset ``default``) is the canonical long per-record table.
"""
from __future__ import annotations

import os

import polars as pl

from vdjmatch import db

TAG = "2026-06-11-ZENODO"
VALID = "^[ACDEFGHIKLMNPQRSTVWY]+$"  # standard 20 AA; drop entries with O/U/*/X/# etc.


def source(pmhc: str | None = None) -> str:
    """Resolve the benchmark table path: explicit ``pmhc`` > ``$VDJDB_SAMPLE`` > HF-pinned release."""
    return str(pmhc or os.environ.get("VDJDB_SAMPLE") or db.fetch_hf(tag=TAG, asset="default"))


def valid_cdr3(df: pl.DataFrame) -> pl.DataFrame:
    """Drop rows whose CDR3 has non-standard residues (else ``Index.build`` rejects them)."""
    return df.filter(pl.col("cdr3").str.contains(VALID))


def spectratype(df: pl.DataFrame, normal_min_n: int = 30) -> pl.DataFrame:
    """Per-epitope CDR3-length spectratype vs the normal repertoire. Columns: ``n`` clonotypes,
    ``modal_frac`` (largest single length), ``eff_len`` = exp(Shannon entropy) = effective #lengths
    (natural ~4-8, phage/clonal ~1), ``max_dev`` = largest per-length-bin deviation from the *normal*
    spectratype (macro-average of per-epitope length fractions over epitopes with n>=normal_min_n; so
    the 10x mega-epitopes don't define 'normal'), and ``dev_len`` = the length of that bin. A severe
    single-bin spike (high ``max_dev``) flags a study artefact even when overall entropy looks fine."""
    import math
    from collections import defaultdict
    u = valid_cdr3(df).select("cdr3", "epitope").unique().with_columns(
        pl.col("cdr3").str.len_chars().alias("L"))
    counts = u.group_by(["epitope", "L"]).len().rename({"len": "c"})
    by: dict[str, dict[int, int]] = defaultdict(dict)
    for epi, L, c in zip(counts["epitope"], counts["L"], counts["c"]):
        by[epi][L] = c
    lens = sorted(u["L"].unique().to_list())
    frac = {e: {L: d.get(L, 0) / sum(d.values()) for L in lens} for e, d in by.items()}
    big = [e for e, d in by.items() if sum(d.values()) >= normal_min_n]
    normal = {L: (sum(frac[e][L] for e in big) / len(big) if big else 0.0) for L in lens}
    out = []
    for epi, d in by.items():
        n = sum(d.values())
        ps = [c / n for c in d.values()]
        h = -sum(p * math.log(p) for p in ps)
        devs = {L: frac[epi][L] - normal[L] for L in lens}
        dl = max(devs, key=lambda L: abs(devs[L]))
        out.append((epi, n, max(d.values()) / n, math.exp(h), abs(devs[dl]), dl))
    return pl.DataFrame(out, schema=["epitope", "n", "modal_frac", "eff_len", "max_dev", "dev_len"],
                        orient="row")


def spectratype_anomalies(df: pl.DataFrame, min_n: int = 30, min_eff_len: float = 2.0,
                          max_modal_frac: float = 0.9, max_dev: float = 0.4) -> set[str]:
    """Epitopes with an anomalous CDR3 length spectratype at n >= ``min_n``: too few effective lengths
    (``eff_len < min_eff_len``), one length dominating (``modal_frac >= max_modal_frac``), or a single
    bin deviating severely from the normal spectratype (``max_dev >= max_dev``) -- phage-display /
    synthetic / clonal / single-study artefacts, quarantined from the natural-repertoire long list."""
    s = spectratype(df).filter(pl.col("n") >= min_n)
    return set(s.filter((pl.col("eff_len") < min_eff_len) | (pl.col("modal_frac") >= max_modal_frac)
                        | (pl.col("max_dev") >= max_dev))["epitope"].to_list())


def multiplex_studies(df: pl.DataFrame, max_epitopes: int = 100) -> pl.DataFrame:
    """References (studies) reporting > ``max_epitopes`` distinct epitopes -- multiplexed screens
    (e.g. 10x dextramer panels) whose bulk, low-confidence calls drive the score-0 mega-epitopes."""
    return (df.group_by("reference_id").agg(pl.col("epitope").n_unique().alias("n_epitopes"),
                                            pl.col("cdr3").n_unique().alias("n_cdr3"))
              .filter(pl.col("n_epitopes") > max_epitopes).sort("n_epitopes", descending=True))


def long_list(df: pl.DataFrame, cap: int = 3000, min_n: int = 30, seed: int = 0,
              drop_anomalous: bool = True) -> pl.DataFrame:
    """Composition-controlled long list of unique clonotypes ``(cdr3,v,j,epitope,mhc_class)``:
    keep epitopes with >= ``min_n`` clonotypes, cap each at a random ``cap`` (so a few 10x mega-studies,
    e.g. SLLMWITQV ~30k, don't dominate sampling or purity), and quarantine spectratype-anomalous
    (phage/clonal) epitopes. Deterministic via ``seed``. The full database stays available (don't pass
    through here) for sampling-depth / motif-coverage studies."""
    u = valid_cdr3(df).select("cdr3", "v", "j", "epitope", "mhc_class").unique()
    if drop_anomalous:
        u = u.filter(~pl.col("epitope").is_in(list(spectratype_anomalies(df, min_n=min_n))))
    keep = u.group_by("epitope").len().filter(pl.col("len") >= min_n)["epitope"]
    u = u.filter(pl.col("epitope").is_in(keep))
    rank = pl.int_range(pl.len()).shuffle(seed=seed).over("epitope")
    return u.with_columns(rank.alias("_r")).filter(pl.col("_r") < cap).drop("_r")
