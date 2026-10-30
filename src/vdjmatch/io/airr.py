"""Read source-aware single-chain and linked-cell query repertoires."""

from __future__ import annotations

import os
import inspect
from pathlib import Path

import polars as pl

from . import columns

_SEP = {".csv": ",", ".tsv": "\t", ".txt": "\t"}


def _read_table(path: str | os.PathLike, *, columns=None) -> pl.DataFrame:
    p = Path(path)
    suffix = p.with_suffix("").suffix if p.suffix.lower() == ".gz" else p.suffix
    if columns is not None:
        scan = (pl.scan_parquet(p) if suffix.lower() in {".parquet", ".pq"} else
                pl.scan_csv(p, separator=_SEP.get(suffix.lower(), "\t"), quote_char='"', infer_schema_length=0))
        return scan.select(columns).collect(engine='streaming')
    if suffix.lower() in {".parquet", ".pq"}:
        return pl.read_parquet(p)
    return pl.read_csv(
        p,
        separator=_SEP.get(suffix.lower(), "\t"),
        quote_char='"',
        infer_schema_length=0,
    )


def _is_wide(raw):
    return bool({c.lower() for c in raw.columns} & {
        "junction_aa_alpha", "junction_aa_beta", "cdr3_aa_alpha", "cdr3_aa_beta"
    })


def _annotation_frame(raw, *, source="auto", sequence_convention=None):
    """Adapt supported raw rows without filtering or dropping source metadata."""
    names = {c.lower(): c for c in raw.columns}
    wide = _is_wide(raw)
    reserved = {"source_query_id", "source_sequence_id", "source_clone_id", "source_pair_id", "_population", "_service_query_id"}
    if reserved & names.keys():
        raise ValueError(f"input contains reserved identity fields: {sorted(reserved & names.keys())}")
    if wide:
        if "junction_aa" in names or "cdr3" in names or "cdr3_aa" in names:
            raise ValueError("ambiguous mixed long and wide sequence columns")
        if not {"junction_aa_alpha", "junction_aa_beta"} <= names.keys():
            raise ValueError("wide input requires junction_aa_alpha and junction_aa_beta; no anchors are fabricated")
        for identity in ("query_id", "sequence_id"):
            if identity in raw.columns:
                saved = "source_" + identity
                if saved in raw.columns:
                    raise ValueError(f"wide input reserves {saved} for the original {identity}")
                raw = raw.with_columns(pl.col(identity).alias(saved))
        if "query_id" in raw.columns:
            if raw["query_id"].null_count() or raw["query_id"].cast(pl.String).str.strip_chars().eq("").any():
                raise ValueError("wide query_id must be nonempty and non-null")
            if raw["query_id"].n_unique() != raw.height:
                raise ValueError("wide query_id must be unique")
    found = columns._resolve(raw)
    # A cell barcode and a clonotype ID can legitimately name different units.
    # Only alternate spellings of the requested clone ID must agree.
    links = [names[c] for c in ("clone_id", "cloneid") if c in names]
    if not links and "pair_id" in found:
        links = [found["pair_id"]]
    if links:
        chosen = pl.col(links[0]).cast(pl.String)
        for other in links[1:]:
            if raw.select((~chosen.eq_missing(pl.col(other).cast(pl.String))).any()).item():
                raise ValueError("conflicting pair/clone identity aliases")
        raw = raw.with_columns(chosen.alias("clone_id"))
    if "pair_id" in found:
        raw = raw.with_columns(pl.col(found["pair_id"]).cast(pl.String).alias("pair_id"))
    if wide:
        if "clone_id" not in raw.columns:
            raw = raw.with_columns(
                (pl.col("query_id").cast(pl.String) if "query_id" in raw.columns else
                 pl.int_range(pl.len()).cast(pl.String)).alias("clone_id")
            )
        if "pair_id" not in raw.columns:
            raw = raw.with_columns(pl.col("clone_id").alias("pair_id"))
        if raw["pair_id"].null_count() or raw["pair_id"].str.strip_chars().eq("").any():
            raise ValueError("wide pair identity must be nonempty and non-null")
        if raw["pair_id"].n_unique() != raw.height:
            raise ValueError("wide pair identity must be unique; multiple pairs cannot be selected arbitrarily")
        row_column = "_vdjmatch_input_row"
        if row_column in raw.columns:
            raise ValueError(f"input reserves {row_column}")
        raw = raw.with_row_index(row_column)
        lanes = []
        for offset, (lane, locus) in enumerate((("alpha", "TRA"), ("beta", "TRB"))):
            expressions = [pl.col(names["junction_aa_" + lane]).cast(pl.String).alias("junction_aa"),
                           pl.lit(locus).alias("locus"),
                           (pl.col(row_column) * 2 + offset).alias("query_id")]
            for gene in ("v", "j"):
                name = names.get(gene + "_" + lane)
                expressions.append((pl.col(name).cast(pl.String) if name else
                                    pl.lit(None, dtype=pl.String)).alias(gene + "_call"))
            if "sequence_id_" + lane in names:
                expressions.append(pl.col(names["sequence_id_" + lane]).cast(pl.String).alias("sequence_id"))
            lanes.append(raw.with_columns(expressions))
        raw = pl.concat(lanes, how="diagonal_relaxed").sort("query_id").drop(row_column)
    else:
        # Existing normalization owns source conventions and count validation.
        columns.normalize_query(raw, valid_aa=False, source=source,
                                sequence_convention=sequence_convention)
        expressions = [pl.col(found["cdr3"]).cast(pl.String).alias("junction_aa")]
        for gene in ("v", "j"):
            expressions.append((pl.col(found[gene]).cast(pl.String) if gene in found else
                                pl.lit(None, dtype=pl.String)).alias(gene + "_call"))
        locus = pl.col(found["locus"]).cast(pl.String) if "locus" in found else (
            pl.col(found["v"]).cast(pl.String).str.slice(0, 3) if "v" in found else
            pl.lit(None, dtype=pl.String))
        expressions.append(locus.str.to_uppercase().alias("locus"))
        raw = raw.with_columns(expressions)
        if "query_id" not in raw.columns:
            raw = raw.with_row_index("query_id")
        if raw["query_id"].null_count() or raw["query_id"].cast(pl.String).str.strip_chars().eq("").any():
            raise ValueError("query_id must be nonempty and non-null")
        if raw["query_id"].n_unique() != raw.height:
            raise ValueError("query_id must be unique")
    if wide:
        columns.normalize_query(raw, valid_aa=False, source=source,
                                sequence_convention=sequence_convention)
    return raw


def normalize_annotation_input(path: str | os.PathLike, *, source="auto", sequence_convention=None):
    """Return canonical raw long rows and ingestion diagnostics, without dropping rows.

    AIRR and simple junction tables retain source columns, identities and metadata.
    Wide ``junction_aa_alpha/beta`` rows expand to TRA/TRB rows in input order;
    original wide IDs remain in ``source_query_id``/``source_sequence_id`` and
    computational chain IDs are unique row offsets. Missing chains remain present.
    Explicit pair/cell links take priority over clone IDs and identify one wide pair
    each; omitted links use the
    supplied unique query ID or source row offset. No anchors are fabricated.
    ``source_rows`` counts physical input rows; ``input_rows`` and sequence
    diagnostics count expanded chain rows, and ``output_rows`` counts raw rows
    returned, including missing/invalid sequences.
    """
    original = _read_table(path)
    raw = _annotation_frame(original, source=source, sequence_convention=sequence_convention)
    _, report = columns.normalize_query(raw, source="airr", return_report=True)
    original_names = {c.lower() for c in original.columns}
    report.update(source_rows=original.height, output_rows=raw.height,
                  wide_input=_is_wide(original),
                  paired_input=_is_wide(original) or bool(original_names & {
                      "pair_id", "cell_id", "clone_id", "cloneid", "complex.id"}))
    return raw, report


def _read_query(path, *, valid_aa, source, sequence_convention, raw=None):
    raw = _read_table(path) if raw is None else raw
    if _is_wide(raw):
        raw = _annotation_frame(raw, source=source, sequence_convention=sequence_convention)
        return columns.normalize_query(raw, valid_aa=valid_aa, source="airr", return_report=True)
    if "query_id" not in raw.columns:
        raw = raw.with_row_index("query_id")
    # Validate original counts/convention before upstream readers can default counts.
    columns.validate_counts(raw)
    lower = {c.lower() for c in raw.columns}
    canonical = bool(lower & {"junction_aa", "cdr3_aa", "cdr3"})
    if canonical:
        columns.normalize_query(
            raw, valid_aa=False, source=source, sequence_convention=sequence_convention
        )
    # Custom CSV and short-form tables use our alias adapter. Standard repertoire
    # formats use the published readers, with AIRR collapsing disabled.
    from vdjtools import io as repertoire_io

    p = Path(path)
    suffix = p.with_suffix("").suffix if p.suffix.lower() == ".gz" else p.suffix
    try:
        fmt = repertoire_io.sniff_format(path) if suffix.lower() != ".csv" else "custom"
    except ValueError:
        fmt = "custom"
    if source in {"custom", "legacy"} or (
        canonical and not {"v_call", "j_call"} <= lower
    ):
        fmt = "custom"
    if fmt == "custom":
        frame = raw
    elif fmt == "airr_cell":
        from vdjtools.sc import read_airr_cell

        frame = read_airr_cell(path)
    elif fmt in {"airr", "arda"}:
        frame = repertoire_io.read_airr(path, collapse=False)
    elif fmt == "mitcr" and "collapse" in inspect.signature(repertoire_io.read_mitcr).parameters:
        frame = repertoire_io.read_mitcr(path, collapse=False)
    else:
        frame = repertoire_io.read(path, fmt=fmt)
    if frame.height != raw.height:
        raise ValueError(
            "upstream reader changed row count; original query identity cannot be preserved"
        )
    if fmt != "custom":
        # Some readers cast numeric-looking source IDs. Keep original identifiers
        # and linkage from the raw rows, whose order is guaranteed by these readers.
        extras = [c for c in ("query_id", "sequence_id", "clone_id", "cloneid", "clonotype_id") if c in raw.columns]
        targets = (
            ("pair_id", "count", "v", "j", "locus")
            if fmt in {"airr", "airr_cell", "arda", "vdjtools"}
            else ("pair_id", "count")
        )
        for target in targets:
            original = columns._resolve(raw).get(target)
            if target == "count" and original is None and fmt == "airr_cell":
                # AIRR Cell readers materialize absent optional counts as null.
                # An omitted source count means one observation, not a malformed count.
                frame = frame.drop(
                    [c for c in columns.ALIASES["count"] if c in frame.columns]
                )
                frame = frame.with_columns(pl.lit(1, dtype=pl.Int64).alias("count"))
            if original is not None:
                # Canonical aliases have priority; remove them before restoring raw.
                for alias in columns.ALIASES[target]:
                    if alias in frame.columns:
                        frame = frame.drop(alias)
                frame = frame.with_columns(raw[original].alias(target))
        # Link aliases were removed above to select the canonical pair ID; retain
        # their original spellings and values as source identities too.
        frame = frame.with_columns(raw[c] for c in extras)
    if 'locus' not in columns._resolve(raw) and 'locus' in frame.columns:
        frame = frame.drop('locus')  # Infer only recognised chains, not UNKNOWN -> UNK.
    return columns.normalize_query(
        frame,
        valid_aa=valid_aa,
        source=source if fmt == "custom" else "airr",
        sequence_convention=sequence_convention,
        return_report=True,
    )


def read_rearrangement(
    path: str | os.PathLike,
    dedup: bool = False,
    valid_aa: bool = True,
    *,
    source: str = "auto",
    sequence_convention: str | None = None,
    return_report: bool = False,
):
    """Read original query rows with stable IDs; explicitly request key deduplication.

    ``dedup=True`` sums counts per junction/V/J/locus/pair and retains original
    identities as ``query_ids`` and (when supplied) ``sequence_ids`` lists.
    ``return_report=True`` returns ``(frame, ingestion_counts)``.
    """
    df, report = _read_query(
        path, valid_aa=valid_aa, source=source, sequence_convention=sequence_convention
    )
    if dedup:
        ids = [c for c in ("query_id", "sequence_id") if c in df.columns]
        df = df.group_by("cdr3", "v", "j", "locus", "pair_id", maintain_order=True).agg(
            pl.col("count").sum(), *[pl.col(c).alias(c + "s") for c in ids]
        )
    report["output_rows"] = df.height
    report["deduplicated_rows"] = report["retained_rows"] - df.height
    return (df, report) if return_report else df


def read_cell(
    path: str | os.PathLike,
    link: str | None = None,
    valid_aa: bool = True,
    *,
    source: str = "auto",
    sequence_convention: str | None = None,
    ambiguity: str = "error",
    return_report: bool = False,
):
    """Pair TRA/TRB rows in first-cell order, preserving chain IDs and counts.

    ``link`` selects the original linkage column. Multiple rows for one cell/locus
    raise by default; ``ambiguity='first'`` explicitly selects the first. Missing
    chains have typed null fields and ``pair_status='missing_alpha'/'missing_beta'``.
    Wide ``junction_aa_alpha/beta`` input also retains wholly missing pairs and
    exposes original wide IDs as ``source_query_ida/b``/``source_sequence_ida/b``.
    """
    if ambiguity not in {"error", "first"}:
        raise ValueError("ambiguity must be 'error' or 'first'")
    raw = _read_table(path)
    wide_input = _is_wide(raw)
    if wide_input:
        raw = _annotation_frame(raw, source=source, sequence_convention=sequence_convention)
        df, report = columns.normalize_query(raw, valid_aa=valid_aa, source="airr", return_report=True)
    else:
        if link == "clone_id" and "clone_id" not in raw.columns:
            original = columns._resolve(raw).get("pair_id")
            if original is not None:
                raw = raw.with_columns(pl.col(original).cast(pl.String).alias("clone_id"))
        df, report = _read_query(path, valid_aa=valid_aa, source=source,
                                 sequence_convention=sequence_convention, raw=raw)
    if link is not None:
        if link not in raw.columns:
            raise ValueError(f"link column {link!r} is absent")
        # Reader-created IDs are row offsets, but user-supplied IDs need not be.
        if "query_id" not in raw.columns:
            raw = raw.with_row_index("query_id")
        if raw["query_id"].n_unique() != raw.height:
            raise ValueError(
                "query_id must be unique when selecting custom cell linkage"
            )
        df = df.drop("pair_id").join(
            raw.select("query_id", pl.col(link).cast(pl.String).alias("pair_id")),
            on="query_id",
            how="left",
            maintain_order="left",
        )
    if df["pair_id"].null_count() or df.filter(pl.col("pair_id") == "").height:
        raise ValueError("every retained cell row requires a nonempty cell/clone link")
    if df.filter(~pl.col("locus").is_in(["TRA", "TRB"]).fill_null(False)).height:
        raise ValueError(
            "paired TCR input requires a TRA or TRB locus for every retained row"
        )
    multiplicity = df.group_by("pair_id", "locus").len()
    n_ambiguous = multiplicity.filter(pl.col("len") > 1).height
    if n_ambiguous and ambiguity == "error":
        raise ValueError(
            f"ambiguous cell input: {n_ambiguous} cell/locus groups have multiple chains"
        )
    fields = ["cdr3", "v", "j", "count", "query_id"]
    fields.extend(c for c in ("sequence_id", "source_query_id", "source_sequence_id",
                             "clone_id", "cloneid", "clonotype_id") if c in df.columns)
    wide = (raw.select(pl.col(link or "pair_id").alias("pair_id")).unique(maintain_order=True)
            if wide_input else df.select("pair_id").unique(maintain_order=True))
    for locus, suffix in (("TRA", "a"), ("TRB", "b")):
        chain = df.filter(pl.col("locus") == locus).unique(
            "pair_id", maintain_order=True
        )
        chain = chain.select("pair_id", *[pl.col(c).alias(c + suffix) for c in fields])
        wide = wide.join(chain, on="pair_id", how="left", maintain_order="left")
    wide = wide.with_row_index("query_id").with_columns(
        pl.when(pl.col("cdr3a").is_null())
        .then(pl.lit("missing_alpha"))
        .when(pl.col("cdr3b").is_null())
        .then(pl.lit("missing_beta"))
        .otherwise(pl.lit("paired"))
        .alias("pair_status")
    )
    report.update(
        output_rows=wide.height,
        ambiguous_cell_locus_groups=n_ambiguous,
        selected_chain_rows=wide["cdr3a"].len()
        - wide["cdr3a"].null_count()
        + wide["cdr3b"].len()
        - wide["cdr3b"].null_count(),
    )
    return (wide, report) if return_report else wide


def read_tcrvdb(
    path: str | os.PathLike, valid_aa: bool = True, *, return_report: bool = False
):
    """TCRvdb wide-form CSV: α/β CDR3s in one row. → ``cdr3a, va, ja, cdr3b, vb, jb,
    epitope, mhc`` (epitope/mhc kept as ground-truth labels for benchmarking)."""
    df = _read_table(path)
    columns.validate_counts(df)
    if "query_id" not in df.columns:
        df = df.with_row_index("query_id")
    cols = {c.lower(): c for c in df.columns}
    g = lambda *names: next((cols[n] for n in names if n in cols), None)  # noqa: E731
    sel = {
        "cdr3a": g("cdr3_alpha_aa"),
        "va": g("trav"),
        "ja": g("traj"),
        "cdr3b": g("cdr3_beta_aa"),
        "vb": g("trbv"),
        "jb": g("trbj"),
        "epitope": g("epitope_aa", "antigen.epitope", "epitope"),
        "mhc": g("hla_short", "hla_long", "mhc.a"),
    }
    out = df.select(
        [pl.col(src).alias(dst) for dst, src in sel.items() if src is not None]
        + [pl.col(c) for c in ("query_id", "sequence_id", "pair_id") if c in df.columns]
    )
    if not {"cdr3a", "cdr3b"} <= set(out.columns):
        raise ValueError("TCRvdb input requires alpha and beta amino-acid columns")
    out = out.with_columns(
        pl.col("cdr3a").str.strip_chars().str.to_uppercase(),
        pl.col("cdr3b").str.strip_chars().str.to_uppercase(),
    )
    missing = (
        pl.col("cdr3a").is_null()
        | pl.col("cdr3b").is_null()
        | (pl.col("cdr3a") == "")
        | (pl.col("cdr3b") == "")
    )
    n_missing = out.select(missing.sum()).item()
    out = out.filter(~missing)
    n_nonmissing = out.height
    if valid_aa:
        rx = r"^[ACDEFGHIKLMNPQRSTVWY]+$"
        out = out.filter(
            pl.col("cdr3a").str.contains(rx) & pl.col("cdr3b").str.contains(rx)
        )
    report = {
        "input_rows": df.height,
        "missing_sequences": n_missing,
        "invalid_sequences": n_nonmissing - out.height,
        "retained_rows": out.height,
        "dropped_rows": df.height - out.height,
        "output_rows": out.height,
    }
    return (out, report) if return_report else out
