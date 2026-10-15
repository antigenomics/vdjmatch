"""Cohort annotation through existing native population batches and shared ingestion."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import time


def _input_rows(inputs, metadata):
    from vdjtools.io import read_metadata
    import polars as pl

    if metadata is not None:
        if inputs is not None:
            raise ValueError("use inputs or metadata, not both")
        if isinstance(metadata, pl.DataFrame):
            frame, base = metadata, Path.cwd()
        else:
            metadata = Path(metadata).expanduser().resolve()
            frame, base = read_metadata(metadata), metadata.parent
        if "path" not in frame.columns:
            raise ValueError("metadata requires a path column")
        rows = []
        for row in frame.iter_rows(named=True):
            value = row["path"]
            if value is None or not str(value).strip():
                raise ValueError("metadata path must be nonempty")
            path = Path(value).expanduser()
            path = path if path.is_absolute() else base / path
            sid = row.get("sample_id", _sample_name(path))
            rows.append((sid, path, {c: v for c, v in row.items() if c not in {"path", "sample_id"}}))
    else:
        if inputs is None:
            raise ValueError("provide input files, a directory, or metadata")
        paths = [inputs] if isinstance(inputs, (str, Path)) else list(inputs)
        rows = []
        for value in paths:
            path = Path(value).expanduser()
            if path.is_dir():
                paths_in_dir = sorted(p for p in path.iterdir() if p.is_file() and
                    (p.suffix.lower() in {".csv", ".tsv", ".txt", ".parquet", ".pq"} or
                     (p.suffix.lower() == ".gz" and p.with_suffix("").suffix.lower() in {".csv", ".tsv", ".txt"})))
                rows.extend((_sample_name(p), p, {}) for p in paths_in_dir)
            else:
                rows.append((_sample_name(path), path, {}))
    if not rows:
        raise ValueError("no input repertoire files found")
    seen_ids, seen_paths, checked = set(), set(), []
    for sid, path, info in rows:
        if sid is None or not str(sid).strip():
            raise ValueError("sample_id must be nonempty")
        sid = str(sid)
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        if sid in seen_ids:
            raise ValueError(f"duplicate sample_id {sid!r}")
        if path in seen_paths:
            raise ValueError("duplicate input path")
        seen_ids.add(sid); seen_paths.add(path)
        checked.append((sid, path, info))
    return checked


def _validate_output_paths(prefix, protected):
    """Reject raw files in the output namespace, including existing symlinks."""
    prefix = Path(prefix).expanduser()
    if prefix.name in {"", ".."}:
        raise ValueError("output prefix must have a basename")
    # Writers append suffixes to the basename, not to a prefix symlink's target.
    prefix = prefix.parent.resolve() / prefix.name
    protected = {Path(path).resolve() for path in protected}
    for path in protected:
        if (path.is_dir() and prefix.parent.is_relative_to(path)) or (
                path.parent == prefix.parent and path.name.startswith(prefix.name + ".") and path.name.endswith((".tsv", ".json"))):
            raise ValueError("annotation output path conflicts with an input/reference/metadata file")
    if prefix.parent.exists():
        for output in prefix.parent.iterdir():
            if (output.name.startswith(prefix.name + ".") and output.name.endswith((".tsv", ".json"))
                    and any(output.resolve() == path or (path.is_dir() and output.resolve().is_relative_to(path)) or
                            (output.exists() and path.is_file() and output.samefile(path)) for path in protected)):
                raise ValueError("annotation output path conflicts with an input/reference/metadata file")
    return prefix


def _sample_name(path):
    path = Path(path)
    if path.suffix.lower() == ".gz":
        path = path.with_suffix("")
    return path.stem


def _core_defaults(reference, prefix):
    from ..cli import neighbours
    parser = argparse.ArgumentParser()
    neighbours.register(parser.add_subparsers(dest="command"), name="annotate")
    return parser.parse_args(["annotate", "unused", "--vdjdb", str(reference),
                              "--locus", "TRB", "--output-prefix", str(prefix)])


def _pair_mapping(raw, first_id):
    import polars as pl
    identity = [c for c in raw.columns if c not in {"query_id", "clone_id", "pair_id", "_population"}]
    fields = [pl.col("sample_id").first(), pl.col("source_clone_id").first(),
              pl.col("source_pair_id").first()]
    for col in identity:
        if col in {"sample_id", "source_pair_id"}:
            continue
        for locus, suffix in (("TRA", "a"), ("TRB", "b")):
            fields.append(pl.col(col).filter(pl.col("locus") == locus).first().alias(col + suffix))
    mapping = raw.group_by("clone_id", maintain_order=True).agg(fields).with_row_index("_service_query_id")
    return mapping.with_columns(pl.col("_service_query_id").cast(pl.String),
        (pl.int_range(pl.len()) + first_id).alias("query_id"),
        pl.col("source_pair_id").alias("pair_id"), pl.lit("paired").alias("locus"))


def _restore(frame, mapping):
    import polars as pl
    if "query_id" not in frame.columns:
        return frame
    frame = frame.with_columns(pl.col("query_id").cast(pl.String)).rename({"query_id": "_service_query_id"})
    if "clone_id" in frame.columns:
        frame = frame.drop("clone_id")
    if "pair_id" in frame.columns:
        frame = frame.drop("pair_id")
    # Generated score/status fields keep their public names; source metadata
    # with the same spelling is retained under an explicit input_ prefix.
    identities = {"_service_query_id", "query_id", "clone_id", "pair_id", "locus",
                  "sequence_id", "sequence_ida", "sequence_idb"}
    rename = {}
    occupied = set(mapping.columns) | set(frame.columns)
    for column in mapping.columns:
        if column in frame.columns and column not in identities:
            target = "input_" + column
            while target in occupied:
                target = "input_" + target
            rename[column] = target
            occupied.add(target)
    mapping = mapping.rename(rename)
    columns = ["_service_query_id"] + [c for c in mapping.columns if c not in frame.columns]
    result = frame.join(mapping.select(columns), on="_service_query_id", how="left", validate="m:1",
                        maintain_order="left").drop("_service_query_id")
    if "source_clone_id" in result.columns:
        result = result.with_columns(pl.col("source_clone_id").alias("clone_id"))
    return result


def annotate_repertoires(inputs=None, *, metadata=None, reference=None, output_prefix="annotation",
                         locus="auto", species="HomoSapiens", threads=1, selection_only=True, **options):
    """Annotate files, a directory, or path/sample_id metadata as population batches.

    Source rows are read once through shared ingestion. Query IDs are global numeric
    offsets; source IDs, metadata and sample-scoped pairing remain available in the
    returned tables. Defaults and advanced scorer options come from the existing
    targeted annotation parser. No per-file search pool or computed cache is used.
    Outputs are ``<prefix>.queries.tsv``, ``.selection.tsv``, optional detail tables,
    and ``.manifest.json``. The result contains frames, manifest and output paths.
    """
    if isinstance(threads, bool) or not isinstance(threads, int) or threads <= 0:
        raise ValueError("threads must be a positive integer")
    if locus not in {"auto", "TRA", "TRB", "paired"}:
        raise ValueError("locus must be auto, TRA, TRB or paired")
    table_budget = os.environ.get("POLARS_MAX_THREADS")
    if table_budget is not None and (not table_budget.isdecimal() or int(table_budget) < 1):
        raise ValueError("POLARS_MAX_THREADS must be a positive integer before import")
    os.environ.setdefault("POLARS_MAX_THREADS", str(threads))
    import polars as pl
    if pl.thread_pool_size() > threads:
        raise ValueError("initialized Polars pool exceeds threads; set POLARS_MAX_THREADS before import or start a new process")
    from vdjtools.cores import kernel_threads
    kernel_threads(threads)
    from .. import db, io
    from ..cli import neighbours
    from ..db.cache import sha256

    started = time.perf_counter()
    defaults = _core_defaults(reference or "unused", output_prefix)
    unknown = set(options) - set(vars(defaults))
    if unknown:
        raise ValueError(f"unknown scorer options: {sorted(unknown)}")
    if set(options) & {"sample", "vdjdb", "output_prefix", "func", "command"}:
        raise ValueError("input/reference/output options use the public runner arguments")
    if options.get('geometry_extra') is not None:
        vars(defaults).update(options)
        defaults.selection_only=selection_only
        defaults.locus=locus
        neighbours._apply_geometry_extra(defaults)
        options.update({key:getattr(defaults,key) for key in
            ('distance_matrix','junction_ends','stem_depth','stem_weight')})
    rows = _input_rows(inputs, metadata)
    prefix = Path(output_prefix).expanduser().resolve()
    protected = [path for _, path, _ in rows]
    protected.extend(Path(value).expanduser().resolve() for value in (metadata, reference, options.get("control"))
                     if isinstance(value, (str, Path)))
    _validate_output_paths(prefix, protected)
    frames, reports, offset = [], [], 0
    for sample_id, path, info in rows:
        raw, report = io.normalize_annotation_input(path)
        for col, value in info.items():
            if col in {"query_id", "sequence_id", "clone_id", "cloneid", "clonotype_id", "pair_id", "cell_id",
                       "source_query_id", "source_sequence_id", "source_clone_id", "source_pair_id", "_population", "_service_query_id"}:
                raise ValueError(f"sample metadata contains reserved identity field {col!r}")
            if col in raw.columns:
                raise ValueError(f"sample metadata column {col!r} conflicts with repertoire data")
            raw = raw.with_columns(pl.lit(value).alias(col))
        if "sample_id" in raw.columns:
            raise ValueError("sample_id is reserved for the input sample identity")
        source_query = pl.col("source_query_id") if "source_query_id" in raw.columns else pl.col("query_id")
        source_sequence = pl.col("source_sequence_id") if "source_sequence_id" in raw.columns else (
            pl.col("sequence_id") if "sequence_id" in raw.columns else pl.lit(None, dtype=pl.String))
        clone = pl.col("clone_id") if "clone_id" in raw.columns else pl.lit(None, dtype=pl.String)
        pair = pl.col("pair_id") if "pair_id" in raw.columns else clone
        linked = report["wide_input"] or (report["paired_input"] and "pair_id" in raw.columns and
            raw.group_by("pair_id").agg(pl.col("locus").n_unique()).filter(pl.col("locus") > 1).height > 0)
        population = ("paired" if linked else "single") if locus == "auto" else locus
        if population == "paired":
            if "pair_id" not in raw.columns or raw["pair_id"].null_count() or raw["pair_id"].str.strip_chars().eq("").any():
                raise ValueError("paired inputs require nonempty pair linkage")
            if raw.group_by("pair_id", "locus").len().filter(pl.col("len") > 1).height:
                raise ValueError("ambiguous paired input: multiple rows per pair/locus")
        raw = raw.with_columns(source_query.alias("source_query_id"), source_sequence.alias("source_sequence_id"),
            clone.alias("source_clone_id"), pair.alias("source_pair_id"), pl.lit(sample_id).alias("sample_id"),
            (pl.int_range(pl.len()) + offset).alias("query_id"), pl.lit(population).alias("_population"))
        raw = raw.with_columns(pl.concat_str([pl.col("sample_id").str.len_chars().cast(pl.String),
            pl.lit(":"), pl.col("sample_id"), pl.col("source_pair_id")]).alias("clone_id"))
        raw = raw.with_columns(pl.col("clone_id").alias("pair_id"))
        frames.append(raw); offset += raw.height
        reports.append({"sample_id": sample_id, "source_sha256": sha256(path), **report})
    cohort = pl.concat(frames, how="diagonal_relaxed")
    frames.clear()
    if cohort.filter(~pl.col("locus").is_in(["TRA", "TRB"]).fill_null(False)).height:
        raise ValueError("every query row requires TRA/TRB locus or a resolving V call")
    if reference is None:
        reference = db.fetch_latest()
        _validate_output_paths(prefix, [Path(reference).expanduser().resolve()])
    prefix.parent.mkdir(parents=True, exist_ok=True)
    tables, manifests, batches = {}, [], []
    paired_rows=cohort.filter(pl.col('_population')=='paired')
    mixed_union=(locus=='auto' and
        options.get('neighbour_weighting',defaults.neighbour_weighting)=='kernel' and
        paired_rows.height>0 and cohort.filter(pl.col('_population')=='single').height>0)
    marginal_tables={};service_batches=[];detail_columns={}
    from ..evalue.control import _organism
    from ..match.geometry import precompute_v_loop_distances
    v_loop_distances=precompute_v_loop_distances(species=_organism(species),
        alphabet=options.get('alphabet',defaults.alphabet),
        matrix=options.get('distance_matrix',defaults.distance_matrix),
        v_regions=options.get('v_regions',defaults.v_regions))
    with TemporaryDirectory(prefix="vdjmatch-annotation-") as temporary:
        for population in ("TRA", "TRB", "paired"):
            part = (cohort.filter(pl.col("_population") == "paired") if population == "paired" else
                cohort.filter((pl.col("_population") == "single") & (pl.col("locus") == population)) if locus == "auto" else
                cohort.filter(pl.col("_population") == population))
            service_part=part
            if mixed_union and population!='paired':
                # Both paired loci retain the existing nonselected-chain dispositions.
                service_part=pl.concat([part,paired_rows],how='vertical')
            if not service_part.height:
                continue
            sample = Path(temporary) / (population + ".tsv")
            local_prefix = Path(temporary) / population
            # Original wide sequence columns remain in identity/output metadata;
            # the scorer receives only the already-normalized long sequence axis.
            wide_sequences = [c for c in service_part.columns if c.lower() in {
                "junction_aa_alpha", "junction_aa_beta", "cdr3_aa_alpha", "cdr3_aa_beta"}]
            service_part.drop("_population", *wide_sequences).write_csv(sample, separator="\t")
            args = _core_defaults(reference, local_prefix)
            vars(args).update(options)
            args.sample, args.locus, args.species = str(sample), population, species
            args.threads, args.selection_only = threads, selection_only
            if mixed_union and population=='paired':
                neighbours.main(args,marginal_tables=marginal_tables,v_loop_distances=v_loop_distances)
            else:
                neighbours.main(args,v_loop_distances=v_loop_distances)
            service_batches.append({'locus':population,'service_chain_rows':service_part.height,
                'standalone_chain_rows':part.height if population!='paired' else 0,
                'paired_chain_rows':paired_rows.height if mixed_union else part.height if population=='paired' else 0})
            mapping = (_pair_mapping(part, offset) if population == "paired" else part.drop("_population").with_columns(
                pl.col("query_id").cast(pl.String).alias("_service_query_id"),
                pl.col("source_pair_id").alias("pair_id")))
            for output in sorted(Path(temporary).glob(population + ".*.tsv")):
                kind = output.name[len(population) + 1:-4]
                frame = pl.read_csv(output, separator="\t", schema_overrides={
                    "query_id": pl.String, "query_ida": pl.String, "query_idb": pl.String,
                    "sequence_id": pl.String, "sequence_ida": pl.String, "sequence_idb": pl.String,
                    "clone_id": pl.String, "pair_id": pl.String})
                if mixed_union and population!='paired':
                    marginal_tables.setdefault(population,{})[kind]=frame.filter(
                        pl.col('query_id').is_in(paired_rows['query_id'].cast(pl.String).implode()))
                    frame=frame.filter(pl.col('query_id').is_in(part['query_id'].cast(pl.String).implode()))
                lane_mapping = (part.drop("_population").with_columns(
                    pl.col("query_id").cast(pl.String).alias("_service_query_id"),
                    pl.col("source_pair_id").alias("pair_id"))
                    if population == "paired" and kind.startswith(("TRA.", "TRB.")) else mapping)
                restored=_restore(frame,lane_mapping)
                if mixed_union and not selection_only:
                    detail_columns.setdefault(population,{})[kind]=restored.columns
                tables.setdefault(kind, []).append(restored)
            manifest = json.loads(Path(str(local_prefix) + ".manifest.json").read_text())
            manifests.append(manifest)
            if part.height:
                batches.append({"locus": population, "query_chain_rows": part.height,
                                "sample_count": part["sample_id"].n_unique()})
    outputs = {}
    result = {}
    if not tables:
        empty = cohort.head(0).drop("_population").with_columns(pl.lit(None, dtype=pl.String).alias("status"))
        tables = {"queries": [empty], "selection": [empty.with_columns(
            pl.lit(None, dtype=pl.String).alias("selection_status"),
            pl.lit(None, dtype=pl.String).alias("selected_epitope"))]}
    for kind, pieces in tables.items():
        frame = pl.concat(pieces, how="diagonal_relaxed")
        if "query_id" in frame.columns:
            frame = frame.sort("query_id")
        path = Path(str(prefix) + "." + kind + ".tsv")
        frame.write_csv(path, separator="\t")
        result[kind], outputs[kind] = frame, str(path)
    manifest = {**(manifests[0] if len(manifests) == 1 else {}), "input_samples": len(rows), "input_chain_rows": cohort.height, "ingestion": reports,
        "populations": batches, "population_service_calls": len(manifests), "native_threads": threads,
        "table_threads": pl.thread_pool_size(), "species": species, "locus": locus,
        "selection_only": selection_only, "core_manifests": manifests, "outputs": outputs,
        **({'mixed_chain_union':True,'service_batches':service_batches} if mixed_union else {}),
        **({'population_table_columns':detail_columns} if detail_columns else {}),
        "wall_seconds": time.perf_counter() - started}
    manifest_path = Path(str(prefix) + ".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    result.update(manifest=manifest, manifest_path=str(manifest_path), output_paths=outputs)
    return result
