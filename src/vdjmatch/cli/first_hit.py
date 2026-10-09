"""Published global first-hit calibration CLI; searches fresh raw target/control sets.

This command reports whole-reference enrichment, not pMHC specificity confidence.
Parser registration stays stdlib-only so the parent can set native thread budgets.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.metadata
import importlib.resources
import json
from pathlib import Path
import resource
import sys
import time


def register(subparsers):
    parser = subparsers.add_parser(
        "first-hit", help="global first-hit control calibration"
    )
    parser.add_argument("sample")
    parser.add_argument("--vdjdb", required=True)
    parser.add_argument("--locus", choices=["TRA", "TRB"], default="TRB")
    parser.add_argument("--species", default="HomoSapiens")
    parser.add_argument("--scope", default="5,2,2,5")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-refs", type=int, default=2)
    parser.add_argument("--exclude-exact", action="store_true")
    parser.add_argument(
        "--control", help="raw junction repertoire; required for nonbundled controls"
    )
    parser.add_argument(
        "--input-format", choices=["auto", "airr", "legacy", "custom"], default="auto"
    )
    parser.add_argument("--sequence-convention", choices=["junction"])
    parser.add_argument("--output-prefix", required=True)
    parser.set_defaults(func=main)
    return parser


def _digest(path):
    path = Path(path)
    if path.is_dir():
        return {
            str(p.relative_to(path)): _digest(p)
            for p in sorted(path.rglob("*"))
            if p.is_file()
        }
    with path.open("rb") as handle:
        digest = hashlib.sha256()
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _control(args):
    import polars as pl
    from seqtree import Index
    from ..evalue.control import _organism

    if args.control:
        path = Path(args.control)
        suffix = path.with_suffix("").suffix if path.suffix == ".gz" else path.suffix
        separator = "," if suffix == ".csv" else "\t"
        headers = pl.read_csv(path, separator=separator, n_rows=0).columns
        column = next(
            (c for c in ("junction_aa", "cdr3_aa", "cdr3", "cdr3aa") if c in headers),
            None,
        )
        if column is None:
            raise ValueError("control requires a junction_aa or legacy junction column")
        if column != "junction_aa" and (
            "v_call" in headers or "sequence_id" in headers
        ):
            raise ValueError(
                "AIRR control requires junction_aa; bare CDR3 cannot be calibrated"
            )
        frame = pl.read_csv(
            path, separator=separator, columns=[column], infer_schema_length=0
        )
        frame = frame.with_columns(pl.col(column).str.strip_chars().str.to_uppercase())
        valid = frame.filter(pl.col(column).str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$"))
        sequences = sorted(set(valid[column]))
        provenance = {
            "kind": "explicit_raw_table",
            "path": str(path.resolve()),
            "sha256": _digest(path),
            "source_rows": frame.height,
            "productive_rows": valid.height,
            "excluded_rows": frame.height - valid.height,
            "declared_species": args.species,
            "declared_locus": args.locus,
        }
    else:
        if _organism(args.species) != "human" or args.locus != "TRB":
            raise ValueError(
                "this control is not bundled; supply --control with the matching raw species/locus repertoire"
            )
        asset = importlib.resources.files("seqtree").joinpath(
            "data/control_human_trb_aa.txt.gz"
        )
        payload = asset.read_bytes()
        raw_sequences = gzip.decompress(payload).decode().splitlines()
        sequences = sorted(set(raw_sequences))
        provenance = {
            "kind": "bundled_raw",
            "asset": "seqtree/data/control_human_trb_aa.txt.gz",
            "sha256": hashlib.sha256(payload).hexdigest(),
            "source_rows": len(raw_sequences),
        }
    if not sequences:
        raise ValueError("control has no valid junctions for the requested locus")
    return Index.build(sequences, "aa"), provenance


def main(args):
    started = time.perf_counter()
    if args.threads < 1:
        raise ValueError("threads must be a positive integer")
    if args.min_refs < 1:
        raise ValueError("min-refs must be a positive integer")
    import polars as pl
    from seqtree import Index
    from .. import db, io
    from ..evalue import first_hit
    from ..io.airr import _read_table
    from ..match.scope import search_params

    params = search_params(args.scope)
    raw = _read_table(args.sample).with_row_index("_input_row")
    query, ingestion = io.read_rearrangement(
        args.sample,
        valid_aa=False,
        source=args.input_format,
        sequence_convention=args.sequence_convention,
        return_report=True,
    )
    # The published reader preserves order and drops only empty sequences here.
    # Retain their positions separately so supplied duplicate IDs remain distinct.
    sequence_col = io.columns._resolve(raw).get("cdr3")
    if sequence_col is None:
        if raw.height != query.height:
            raise ValueError(
                "cannot preserve missing sequence positions for this input format"
            )
        positions = raw["_input_row"]
    else:
        positions = raw.filter(
            pl.col(sequence_col).is_not_null()
            & (
                pl.col(sequence_col).cast(pl.String).str.strip_chars().str.len_chars()
                > 0
            )
        )["_input_row"]
    if len(positions) != query.height:
        raise ValueError("query reader did not preserve original sequence positions")
    query = query.with_columns(positions)
    valid = pl.col("cdr3").str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
    same_locus = pl.col("locus").is_null() | (pl.col("locus") == args.locus)
    retained = query.filter(valid & same_locus)
    reference = db.load(args.vdjdb, species=args.species, gene=args.locus)
    source_records = reference.height
    if args.min_refs > 1:
        key = ["gene", "cdr3", "v", "j", "epitope"]
        study = pl.col("reference_id")
        keep = (
            reference.group_by(key)
            .agg(
                study.filter(study.is_not_null() & (study != ""))
                .n_unique()
                .alias("n_studies")
            )
            .filter(pl.col("n_studies") >= args.min_refs)
            .select(key)
        )
        reference = reference.join(keep, on=key, how="semi", nulls_equal=True)
    selected_records = reference.height
    reference = reference.filter(pl.col("reference_valid"))
    targets = sorted(set(reference["cdr3"]))
    if not targets:
        raise ValueError("reference selection has no valid junctions")
    target = Index.build(targets, "aa")
    control, control_source = _control(args)
    queries = retained["cdr3"].to_list()
    target_hits, control_costs = first_hit.scan(
        target,
        [""] * len(target),
        control,
        queries,
        params=params,
        threads=args.threads,
        exclude_exact=args.exclude_exact,
        chunk=max(1, len(queries)),
    )
    results = {
        row_id: first_hit.pvalue(hits, costs, len(target), len(control))
        for row_id, hits, costs in zip(
            retained["_input_row"], target_hits, control_costs
        )
    }
    normalized = {row["_input_row"]: row for row in query.iter_rows(named=True)}
    count_column = io.columns._resolve(raw).get("count")
    counts = raw.select(
        io.columns._integer_count(count_column)
        if count_column
        else pl.lit(1).repeat_by(pl.len()).explode()
    ).to_series()
    rows = []
    for source in raw.iter_rows(named=True):
        position = source["_input_row"]
        item = normalized.get(position)
        evidence = results.get(position)
        if item is None:
            status = "missing_sequence"
        elif not item["cdr3"] or any(aa not in io.columns.AA for aa in item["cdr3"]):
            status = "invalid_sequence"
        elif item["locus"] is not None and item["locus"] != args.locus:
            status = "locus_mismatch"
        else:
            status = "matched" if evidence["radius"] is not None else "no_hit"
        rows.append(
            {
                "query_id": source.get("query_id", position),
                "sequence_id": source.get("sequence_id"),
                "input_row": position,
                "query_cdr3": item["cdr3"] if item else source.get(sequence_col),
                "query_v": item["v"] if item else None,
                "locus": (item["locus"] or args.locus) if item else args.locus,
                "count": item["count"] if item else counts[position],
                "radius": evidence["radius"] if evidence else None,
                "n_target": evidence["n_target"] if evidence else 0,
                "n_control": evidence["n_control"] if evidence else 0,
                "n_reference": len(target),
                "control_size": len(control),
                "E": evidence["E"] if evidence else None,
                "p_enrichment": evidence["p_enrichment"] if evidence else None,
                "p_any": evidence.get("p_any") if evidence else None,
                "rule_of_three": evidence.get("rule_of_three", False)
                if evidence
                else False,
                "status": status,
            }
        )
    schema = {
        "query_id": raw.schema.get("query_id", pl.UInt32),
        "sequence_id": pl.String,
        "input_row": pl.UInt32,
        "query_cdr3": pl.String,
        "query_v": pl.String,
        "locus": pl.String,
        "count": pl.Int64,
        "radius": pl.Int64,
        "n_target": pl.Int64,
        "n_control": pl.Int64,
        "n_reference": pl.Int64,
        "control_size": pl.Int64,
        "E": pl.Float64,
        "p_enrichment": pl.Float64,
        "p_any": pl.Float64,
        "rule_of_three": pl.Boolean,
        "status": pl.String,
    }
    output = pl.DataFrame(rows, schema=schema)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    output.write_csv(str(prefix) + ".evidence.tsv", separator="\t")
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    manifest = {
        "calculation": "global-first-hit-v1",
        "interpretation": "whole-reference enrichment; not pMHC specificity confidence",
        "sample_sha256": _digest(args.sample),
        "reference_sha256": _digest(args.vdjdb),
        "parameters": {
            k: getattr(args, k)
            for k in (
                "locus",
                "species",
                "scope",
                "threads",
                "min_refs",
                "exclude_exact",
                "input_format",
                "sequence_convention",
            )
        },
        "software": {
            name: importlib.metadata.version(name)
            for name in ("vdjmatch", "seqtree", "vdjtools", "polars")
        },
        "queries": {
            "input_rows": raw.height,
            "searched_rows": retained.height,
            "excluded_rows": raw.height - retained.height,
            "status_counts": {
                s: n for s, n in output.group_by("status").len().iter_rows()
            },
            "ingestion": ingestion,
        },
        "reference": {
            "source_records": source_records,
            "selected_records": selected_records,
            "valid_records": reference.height,
            "unique_junctions": len(target),
        },
        "control": {**control_source, "unique_junctions": len(control)},
        "counting_unit": "unique target/control junctions; original query rows",
        "resources": {
            "wall_seconds": time.perf_counter() - started,
            "peak_rss_bytes": peak_bytes,
            "processes": 1,
            "rss_scope": "CLI process; no subprocesses created",
        },
    }
    Path(str(prefix) + ".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return 0
