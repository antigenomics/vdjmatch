"""Published naive precursor frequencies, and what this library predicts for the same epitopes.

`F(e)` estimates a **naive repertoire mass**: the fraction of a naive repertoire that can see an
epitope. That is a quantity other laboratories have measured directly, by tetramer or dextramer
enrichment of unexposed donors, umbilical cord blood and naive mice — so the estimator can be
checked against it rather than against a proxy. These two loaders fetch that comparison.

**The two tables are never merged.** :func:`load_compendium` is *experimental* — what a laboratory
counted. :func:`load_estimates` is *derived* — a mass under a generative recombination model.
Putting them in one column would conflate a measurement with a prediction of it, which is the whole
thing the comparison is for.

Both are mirrored at ``isalgo/airr_benchmark`` under ``vdjmatch/precursor_freq/`` and cached
alongside the VDJdb releases, so an installed user can reproduce the validation without checking out
the analysis repository.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl

from ..db.cache import cache_dir

#: HuggingFace mirror holding both tables.
HF_REPO = "isalgo/airr_benchmark"
HF_DIR = "vdjmatch/precursor_freq"

_FILES = {"compendium": "precursor_compendium.tsv", "estimates": "precursor_estimates.tsv"}

#: Columns of the estimates table that are numeric even where they are blank for a whole arm.
_NUMERIC = ("n_junctions", "observed_mass", "union", "retained", "n_ball", "occupancy",
            "event_ratio", "event_denominator", "measured_per_1e6", "measured_lo", "measured_hi",
            "n_measurements", "n_studies")


def _fetch(which: str, cache: str | os.PathLike | None = None, force: bool = False,
           repo: str = HF_REPO) -> Path:
    """Download (and cache) one of the two tables from the HuggingFace mirror."""
    name = _FILES[which]
    out = cache_dir(cache) / name
    if out.exists() and not force:
        return out
    from huggingface_hub import hf_hub_download  # optional dep (extras: bench/control)
    src = hf_hub_download(repo_id=repo, repo_type="dataset", filename=f"{HF_DIR}/{name}")
    out.write_bytes(Path(src).read_bytes())
    return out


def load_compendium(source: str | os.PathLike | None = None, *, measured_only: bool = True,
                    dedup_retabulations: bool = False, species: str | None = None,
                    cache: str | os.PathLike | None = None) -> pl.DataFrame:
    """Published measurements of naive precursor frequency — **experimental**, not computed.

    147 rows from 14 studies tracing to 31 originating laboratories (one is a review, whose rows
    carry the laboratory that made them in ``lab``). ``freq_per_1e6_naive`` is the harmonized
    column: cells per million of the *naive* compartment. The number as published is kept beside it
    in ``freq_per_1e6_published`` with the denominator it was published under, so every conversion
    is reversible.

    Args:
        source: path to an existing table; ``None`` fetches from the HuggingFace mirror.
        measured_only: keep only rows that *measured* an endogenous naive frequency — drops the
            adoptive-transfer titrations, where the frequency is an experimental input rather than
            an observation, and the memory arms. Almost always what you want.
        dedup_retabulations: drop the ten review rows that restate a primary already in the table
            (Obar 2008 on six mouse epitopes, Alanio 2010 on four human). The two copies differ only
            by the denominator the review renormalized to, so counting both turns one measurement
            into two. Off by default because a per-row median legitimately weights by how many
            donors stand behind a row; turn it on when you are counting *measurements*.
        species: ``"human"`` or ``"mouse"``.
        cache: cache directory (default ``~/.cache/vdjmatch`` or ``$VDJMATCH_CACHE``).

    Returns:
        One row per published record. See the mirror's ``README.md`` for the full schema.

    Note:
        Two denominator traps are already applied and must not be applied twice. Adult blood is
        ~50% naive-phenotype, so a per-total-CD8 figure is doubled to reach the naive scale; cord
        blood is essentially the naive compartment in full and is *not* doubled.
    """
    path = Path(source) if source is not None else _fetch("compendium", cache)
    df = pl.read_csv(path, separator="\t", infer_schema_length=0)
    df = df.with_columns(
        pl.col("freq_per_1e6_naive").cast(pl.Float64, strict=False),
        pl.col("freq_per_1e6_published").cast(pl.Float64, strict=False),
        pl.col("naive").cast(pl.Int64, strict=False),
    )
    if measured_only:
        df = df.filter((pl.col("design") == "measured") & (pl.col("naive") == 1)
                       & pl.col("freq_per_1e6_naive").is_not_null())
    if dedup_retabulations:
        held = set(zip(*[df.filter(pl.col("source_type") != "review")[c]
                         for c in ("species", "epitope_seq", "lab")]))
        df = df.filter(~pl.Series([r["source_type"] == "review" and r["epitope_seq"]
                                   and (r["species"], r["epitope_seq"], r["lab"]) in held
                                   for r in df.iter_rows(named=True)]))
    if species is not None:
        df = df.filter(pl.col("species") == species)
    return df


def load_estimates(source: str | os.PathLike | None = None, *, species: str | None = None,
                   chain: str | None = None, model: str | None = None,
                   cache: str | os.PathLike | None = None) -> pl.DataFrame:
    """What `F(e)` predicts for the epitopes the compendium measured — **derived**, not observed.

    One row per (species, model set, chain, epitope), at radius 1 with ``alpha = 0.1``. The
    estimator columns are repertoire *fractions*; multiply by ``1e6`` to compare against
    ``measured_per_1e6``, which is carried alongside with its across-study range.

    Args:
        source: path to an existing table; ``None`` fetches from the HuggingFace mirror.
        species: ``"human"`` or ``"mouse"``.
        chain: ``"TRA"`` or ``"TRB"``.
        model: recombination model set — ``"olga"`` or ``"arda"``.
        cache: cache directory.

    Note:
        Human is scored under **both** model sets on purpose, because mouse exists only under
        ``arda`` and a human-``olga`` against mouse-``arda`` comparison would confound species with
        model family. The sets are not interchangeable: ``olga`` scores 5.7% of human TRA junctions
        at exactly zero, where they contribute nothing to any mass and raise no error.
    """
    path = Path(source) if source is not None else _fetch("estimates", cache)
    # Pin the numeric columns rather than trusting inference: `n_ball` and `event_ratio` are empty
    # for whole arms (the ball census is not computed for every model set, and the event ratio only
    # exists where a pooled control repertoire does), so a leading run of blanks makes polars infer
    # them as strings and any arithmetic on them fails at the call site instead of here.
    df = pl.read_csv(path, separator="\t", schema_overrides={c: pl.Float64 for c in _NUMERIC})
    for col, val in (("species", species), ("chain", chain), ("source", model)):
        if val is not None:
            df = df.filter(pl.col(col) == val)
    return df
