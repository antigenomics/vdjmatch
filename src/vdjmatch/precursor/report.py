"""One row per TCR group: the `Pgen` spread, the union, the unseen set and the frequency.

This is what ``vdjmatch precursor`` writes. Every column is one of the estimators in this package;
nothing is computed here that the API does not also expose.
"""
from __future__ import annotations

import math

import polars as pl

from .frequency import N_T_CELLS, expected_cells, p_at_least
from .mass import ALPHA_PER_EDIT, shell_profile, union_mass
from .model import _JUNCTION_END, _JUNCTION_START, load_model
from .model import pgen as _pgen
from .unseen import unseen_junctions

#: Column order and dtype of the report, so the TSV is stable across releases. ``n_unseen`` is a
#: float on purpose: it is a Horvitz-Thompson estimate of a count, not a count.
SCHEMA: dict[str, object] = {
    "group": pl.String, "chain": pl.String,
    "n_records": pl.Int64, "n_junctions": pl.Int64, "n_dropped": pl.Int64,
    "pgen_median": pl.Float64, "pgen_q25": pl.Float64, "pgen_q75": pl.Float64,
    "pgen_log10_span": pl.Float64, "top1_share": pl.Float64, "top10_share": pl.Float64,
    "observed_mass": pl.Float64, "naive_sum": pl.Float64, "union": pl.Float64,
    "overlap": pl.Float64, "n_components": pl.Int64, "n_clustered": pl.Int64,
    "retained": pl.Float64, "F": pl.Float64, "q": pl.Float64, "alpha": pl.Float64,
    "radius": pl.Int64,
    "cells": pl.Float64, "lambda": pl.Float64, "p_ge_1": pl.Float64, "p_ge_10": pl.Float64,
    "n_ball": pl.Int64, "n_unseen_ball": pl.Int64, "unseen_ball_mass": pl.Float64,
    "n_unseen_ht": pl.Float64, "unseen_ht_mass": pl.Float64, "rarity_ratio": pl.Float64,
    "richness_reliable": pl.Boolean, "unseen_status": pl.String,
}


def _clean(junctions, multiplicity=None):
    """Deduplicate and junction-QC in one pass, keeping ``multiplicity`` aligned.

    Returns ``(seqs, mult, n_records, n_dropped)``. A repeated junction keeps its first
    multiplicity — the caller is expected to have aggregated already, and silently summing
    recaptures here would inflate exactly the quantity the coverage correction reads.
    """
    raw = list(junctions)
    mult = list(multiplicity) if multiplicity is not None else None
    if mult is not None and len(mult) != len(raw):
        raise ValueError(f"multiplicity has {len(mult)} entries for {len(raw)} junctions")
    seqs, keep, seen = [], [], set()
    dropped = 0
    for i, s in enumerate(raw):
        t = (s or "").strip().upper()
        if not (t.startswith(_JUNCTION_START) and t.endswith(_JUNCTION_END)):
            dropped += 1                    # an anchor-stripped IMGT CDR3 would score 0.0 silently
            continue
        if t in seen:
            continue
        seen.add(t)
        seqs.append(t)
        keep.append(i)
    return seqs, ([mult[i] for i in keep] if mult is not None else None), len(raw), dropped


def _spread(ps: list[float]) -> dict:
    """Order statistics of a `Pgen` vector, and how concentrated its sum is.

    ``top1_share``/``top10_share`` are the fraction of the total mass carried by the single largest
    and the ten largest junctions. They are the reason a mean is the wrong summary: a cognate set
    spanning many orders of magnitude has its sum set by a handful of public clonotypes.
    """
    nz = sorted((p for p in ps if p > 0), reverse=True)
    if not nz:
        return {}
    total, n, asc = sum(nz), len(nz), nz[::-1]
    return {"pgen_median": asc[n // 2], "pgen_q25": asc[n // 4], "pgen_q75": asc[(3 * n) // 4],
            "pgen_log10_span": math.log10(nz[0] / nz[-1]),
            "top1_share": nz[0] / total, "top10_share": sum(nz[:10]) / total}


def summarise_group(model, junctions, *, group: str = "", chain: str = "", r: int = 1,
                    alpha: float = ALPHA_PER_EDIT, q: float = 1.0, n_cells: float = N_T_CELLS,
                    compartment: float = 1.0, n_eff: float | None = None,
                    multiplicity=None, n_units: int | None = None, threads: int = 0) -> dict:
    """Every estimator for one set of junctions, as a flat row keyed by :data:`SCHEMA`.

    ``junctions`` are **junctions** (Cys104..Phe/Trp118 inclusive), not IMGT CDR3s; the ones that
    fail the anchor check are dropped and counted in ``n_dropped`` rather than scored as 0.
    ``multiplicity`` — one capture count per input junction, i.e. how many independent units
    re-reported it — switches on the unseen-species fields; without it they are null and
    ``unseen_status`` says why.
    """
    seqs, mult, n_records, dropped = _clean(junctions, multiplicity)
    row: dict = dict.fromkeys(SCHEMA)
    row.update(group=group, chain=chain, n_records=n_records, n_junctions=len(seqs),
               n_dropped=dropped, q=q, alpha=alpha, radius=r,
               unseen_status="no capture units given" if mult is None else "ok")
    if not seqs:
        return row

    ps = _pgen(model, seqs, threads=threads)
    row.update(_spread(ps))
    prof = shell_profile(model, seqs, r=r, alpha=alpha, threads=threads)
    top = union_mass(model, seqs, r=r, threads=threads)
    f = q * prof["retained"]
    observed = float(sum(ps))
    row.update(observed_mass=observed, naive_sum=top["naive_sum"], union=prof["union"],
               overlap=prof["overlap"], n_components=top["n_components"],
               n_clustered=top["n_clustered"], retained=prof["retained"], F=f,
               cells=expected_cells(f, n_cells=n_cells, compartment=compartment))

    # The finite census: every sequence in the union that no database has catalogued. Bounded and
    # exact given the ball, unlike the Horvitz-Thompson richness below.
    if top["n_union"] is not None:
        row.update(n_ball=top["n_union"], n_unseen_ball=top["n_union"] - len(seqs),
                   unseen_ball_mass=prof["union"] - observed)

    if n_eff is not None:
        row.update({"lambda": float(n_eff) * f, "p_ge_1": p_at_least(f, n_eff, 1),
                    "p_ge_10": p_at_least(f, n_eff, 10)})

    if mult is not None:
        u = unseen_junctions(model, seqs, mult, n_units=n_units, threads=threads)
        row["unseen_status"] = u["reason"] if u["degenerate"] else "ok"
        if not u["degenerate"]:
            row.update(n_unseen_ht=u["n_unseen"], unseen_ht_mass=u["unseen_mass"],
                       rarity_ratio=u["rarity_ratio"],
                       richness_reliable=u["richness_reliable"])
    return row


def summarise(frame: pl.DataFrame, *, junction_col: str = "cdr3", group_col: str | None = None,
              chain_col: str | None = None, capture_col: str | None = None,
              locus: str = "TRB", source: str = "olga", organism: str = "human",
              r: int = 1, alpha: float = ALPHA_PER_EDIT, q: float = 1.0,
              n_cells: float = N_T_CELLS, compartment: float = 1.0,
              n_eff: float | None = None, min_junctions: int = 1,
              threads: int = 0, progress: bool = False) -> pl.DataFrame:
    """Run :func:`summarise_group` over every group in ``frame``.

    ``group_col`` names the grouping column (``epitope`` on a normalised VDJdb table); without it
    the whole frame is one group. ``chain_col`` splits by locus and loads the matching recombination
    model per chain, so a mixed TRA/TRB table is scored correctly rather than with one model.
    ``capture_col`` (``reference_id`` on VDJdb) supplies the recapture multiplicities the
    unseen-species fields need — the number of distinct units that re-reported each junction.
    """
    from .._util import progress as _progress

    keys = [c for c in (chain_col, group_col) if c]
    groups = (list(frame.group_by(keys, maintain_order=True)) if keys
              else [((), frame)])
    models: dict[str, object] = {}
    rows = []
    for key, sub in _progress(groups, total=len(groups), desc="precursor", enable=progress):
        chain = str(key[0]) if chain_col else locus
        name = str(key[-1]) if group_col else ""
        if capture_col:
            agg = (sub.group_by(junction_col, maintain_order=True)
                      .agg(pl.col(capture_col).n_unique().alias("_m")))
            seqs, mult = agg[junction_col].to_list(), agg["_m"].to_list()
            n_units = int(sub[capture_col].n_unique())
        else:
            seqs, mult, n_units = sub[junction_col].unique(maintain_order=True).to_list(), None, None
        if len(seqs) < min_junctions:
            continue
        if chain not in models:
            models[chain] = load_model(chain, source=source, organism=organism)
        rows.append(summarise_group(
            models[chain], seqs, group=name, chain=chain, r=r, alpha=alpha, q=q,
            n_cells=n_cells, compartment=compartment, n_eff=n_eff,
            multiplicity=mult, n_units=n_units, threads=threads))
    return pl.DataFrame(rows, schema=SCHEMA)
