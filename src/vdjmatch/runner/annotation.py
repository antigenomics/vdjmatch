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


def _validate_output_paths(prefix, protected, *, suffixes=None):
    """Reject raw files in the output namespace, including existing symlinks."""
    prefix = Path(prefix).expanduser()
    if prefix.name in {"", ".."}:
        raise ValueError("output prefix must have a basename")
    # Writers append suffixes to the basename, not to a prefix symlink's target.
    prefix = prefix.parent.resolve() / prefix.name
    protected = {Path(path).resolve() for path in protected}
    def is_output(name):
        if suffixes is not None:
            return name.removeprefix(prefix.name) in suffixes and name.startswith(prefix.name)
        return name.startswith(prefix.name + ".") and name.endswith((".tsv", ".json", ".parquet", ".pq"))

    for path in protected:
        if (path.is_dir() and prefix.parent.is_relative_to(path)) or (
                path.parent == prefix.parent and is_output(path.name)):
            raise ValueError("annotation output path conflicts with an input/reference/metadata file")
    if prefix.parent.exists():
        for output in prefix.parent.iterdir():
            if (is_output(output.name)
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


def _vx_population(sample, part, reference, context, threads, reference_population, exclude_exact, locus, tie_policy):
    """Score one receptor population and keep its service row axis."""
    import polars as pl
    from ..cli.neighbours import _paired_input
    from ..cli.paired_vx_base import CHAINS, _all_valid, _reference_catalogue, _score_binding, _target_axis, _vx_maxima
    from ..db.vdjdb import _select_reference
    from ..evalue.equal_tail import _eligible, _prepare

    species = context["species"]
    organism = {"human": "HomoSapiens", "mouse": "MusMusculus"}[species]
    chains = CHAINS[locus]
    reference = _select_reference(reference, species=organism)
    if reference_population == "all-valid":
        reference = _all_valid(reference)
    fixed, targets, groups, backgrounds, sources, literal_keys = _reference_catalogue(
        reference, context["model"], species=species, locus=locus, allow_empty=True,
        literal_ties=tie_policy == "literal")
    if locus == "paired":
        pairs, _ = _paired_input(sample, part, context["model"], species, resolved_v_nearest=True)
        queries = _prepare(pairs.with_columns(pl.col("v_call_alpha").alias("va"),
                                             pl.col("v_call_beta").alias("vb")), chains, context["model"])
    else:
        aliases = {"human": ["human", "homosapiens", "homo sapiens"],
                   "mouse": ["mouse", "musmusculus", "mus musculus"]}
        species_ok = (pl.col("species").str.strip_chars().str.to_lowercase().is_in(aliases[species]).fill_null(False)
                      if "species" in part.columns else pl.lit(True))
        queries = _prepare(part.with_columns(
            pl.col("junction_aa").str.strip_chars().str.to_uppercase().alias("cdr3" + chains),
            pl.col("v_call").alias("v" + chains)), chains, context["model"])
        queries = queries.with_columns(
            pl.when(~species_ok).then(pl.lit("unselected_species"))
            .when((pl.col("locus") != locus) | ~_eligible(chains, "native", 5))
            .then(pl.lit("invalid_query_or_gene")).otherwise(pl.lit("searched")).alias("status"))
    first_chain = "native" + chains[0]
    queries = queries.with_columns(
        pl.when(pl.col("status") == "searched").then(pl.col(first_chain)).otherwise(None).alias(first_chain))
    axis = _target_axis(targets)
    labels = axis["target"].to_list()
    maxima = _vx_maxima(queries, fixed, groups, backgrounds, context, labels, None,
                        threads=threads, exclude_exact=exclude_exact, locus=locus, sources=sources,
                        literal_keys=literal_keys)
    selected = targets.with_columns(axis["target"])
    selected = selected.rename({name: "selected_" + name for name in selected.columns if name != "target"})
    disposition = maxima.with_columns(
        pl.when(pl.col("score").is_not_null()).then(pl.lit("ranked_evidence"))
        .when("tied").then(pl.lit("ambiguous")).otherwise(pl.lit("no_support")).alias("selection_status"),
        pl.lit("vx").alias("score_model"),
    )
    disposition = disposition.join(selected, on="target", how="left", validate="m:1", maintain_order="left")
    query = queries.select("query_id", "status").join(
        maxima.select("query_id", "maximum"), on="query_id", validate="1:1", maintain_order="left")
    query = query.with_columns(
        pl.when(pl.col("status") != "searched").then(pl.col("status"))
        .when(pl.lit(targets.is_empty())).then(pl.lit("no_reference"))
        .when(pl.col("maximum").is_not_null()).then(pl.lit("matched"))
        .otherwise(pl.lit("no_neighbours")).alias("status")
    ).drop("maximum")
    return query, disposition, _score_binding(
        fixed.height, exclude_exact, species=species, locus=locus,
        reference_population=reference_population, tie_policy=tie_policy), axis


def annotate_repertoires(inputs=None, *, metadata=None, reference=None, output_prefix="annotation",
                         locus="auto", species="HomoSapiens", threads=0, selection_only=True,
                         score_model="auto", vx_null_dir=None, vx_target_fpr=.01, vx_tie_policy=None, **options):
    """Annotate files, a directory, or path/sample_id metadata as population batches.

    Source rows are read once through shared ingestion. Query IDs are global numeric
    offsets; source IDs, metadata and sample-scoped pairing remain available in the
    returned tables. Defaults and advanced scorer options come from the existing
    targeted annotation parser. No per-file search pool or computed cache is used.
    ``threads=0`` uses all allocated physical cores; positive counts may use allocated logical CPUs.
    Outputs are ``<prefix>.queries.tsv``, ``.selection.tsv``, optional detail tables,
    and ``.manifest.json``. The result contains frames, manifest and output paths.
    For calibrated calls, supply one vX null directory per scored locus.
    A paired bank determines the tie policy unless explicitly selected.
    """
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 0:
        raise ValueError("threads must be a nonnegative integer (0 selects allocated cores)")
    from vdjtools.cores import kernel_threads
    threads = kernel_threads(threads)
    if locus not in {"auto", "TRA", "TRB", "paired"}:
        raise ValueError("locus must be auto, TRA, TRB or paired")
    if score_model not in {"auto", "vx", "kernel"}:
        raise ValueError("score_model must be auto, vx or kernel")
    if vx_tie_policy not in {None, "abstain", "literal"}:
        raise ValueError("vx_tie_policy must be abstain or literal")
    if score_model == "kernel" and vx_null_dir is not None:
        raise ValueError("vX null calibration requires the vX annotation score")
    table_budget = os.environ.get("POLARS_MAX_THREADS")
    if table_budget is not None and (not table_budget.isdecimal() or int(table_budget) < 1):
        raise ValueError("POLARS_MAX_THREADS must be a positive integer before import")
    os.environ.setdefault("POLARS_MAX_THREADS", str(threads))
    import polars as pl
    if pl.thread_pool_size() > threads:
        raise ValueError("initialized Polars pool exceeds threads; set POLARS_MAX_THREADS before import or start a new process")
    from .. import db, io
    from ..cli import neighbours
    from ..db.cache import sha256

    started = time.perf_counter()
    defaults = _core_defaults(reference or "unused", output_prefix)
    original_defaults = vars(defaults).copy()
    unknown = set(options) - set(vars(defaults))
    if unknown:
        raise ValueError(f"unknown scorer options: {sorted(unknown)}")
    if set(options) & {"sample", "vdjdb", "output_prefix", "func", "command"}:
        raise ValueError("input/reference/output options use the public runner arguments")
    reference_population = options.get("reference_population", defaults.reference_population)
    if reference_population not in {"curated-benchmark", "all-valid"}:
        raise ValueError("reference_population must be curated-benchmark or all-valid")
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
    frames, counts, reports, offset = [], [], [], 0
    for sample_id, path, info in rows:
        raw, report = io.normalize_annotation_input(path)
        # Resolve validated counts before the cohort introduces absent-alias nulls.
        count_column = io.columns._count_column(raw)
        count = (io.columns._integer_count(count_column) if count_column is not None else
                 pl.repeat(1, pl.len(), dtype=pl.Int64))
        counts.append(raw.select(count.alias("count"))["count"])
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
    canonical_counts = pl.concat(counts)
    counts.clear()
    frames.clear()
    if cohort.filter(~pl.col("locus").is_in(["TRA", "TRB"]).fill_null(False)).height:
        raise ValueError("every query row requires TRA/TRB locus or a resolving V call")
    if reference is None:
        reference = db.fetch_latest()
        _validate_output_paths(prefix, [Path(reference).expanduser().resolve()])
    prefix.parent.mkdir(parents=True, exist_ok=True)
    tables, manifests, batches = {}, [], []
    paired_rows=cohort.filter(pl.col('_population')=='paired')
    from ..evalue.control import _organism
    organism = _organism(species)
    vx_profile = selection_only and all(value == original_defaults[name] for name, value in options.items()
                                        if name not in {"reference_population", "exclude_exact", "verbose"})
    vx_enabled = score_model == "vx" or (score_model == "auto" and (
        vx_null_dir is not None or (organism == "mouse" and paired_rows.height == cohort.height and vx_profile)))
    if vx_enabled and not vx_profile:
        raise ValueError("vX scoring requires ordinary annotation options")
    if vx_tie_policy == "literal" and (not vx_enabled or paired_rows.is_empty()):
        raise ValueError("literal tie resolution requires paired vX annotation")
    banks = {}
    if vx_null_dir is not None:
        active_lanes = ([locus] if locus != "auto" else
                        cohort.filter(pl.col("_population") == "single")["locus"].unique().to_list()
                        + (["paired"] if paired_rows.height else []))
        archive = Path(reference)
        if not archive.is_file() or archive.suffix != ".zip":
            raise ValueError("vX null calibration requires the same VDJdb release ZIP used to score the null")
        reference_digest = sha256(archive)
        exclude_exact = options.get("exclude_exact", defaults.exclude_exact)
        from ..cli.paired_vx_evalue import _calibrate, _source
        directories = [vx_null_dir] if isinstance(vx_null_dir, (str, os.PathLike)) else vx_null_dir
        for directory in directories:
            null_directory = Path(directory).expanduser()
            bank = _source(null_directory)
            _, binding, _, _, input_paths = bank
            lane = binding["locus"]
            if (lane not in active_lanes or binding["species"] != organism
                    or binding["reference_population"] != reference_population):
                raise ValueError("vX null locus, species and reference population must match the annotation run")
            if lane == "paired" and vx_tie_policy is not None and binding["tie_policy"] != vx_tie_policy:
                raise ValueError("vX null tie policy must match the annotation run")
            if binding["reference_sha256"] != reference_digest or binding["exact_hits_excluded"] != exclude_exact:
                raise ValueError("vX null reference and exact-hit policy must match the annotation run")
            if lane in banks:
                raise ValueError(f"more than one vX null directory for {lane}")
            _validate_output_paths(prefix, [null_directory, *input_paths])
            banks[lane] = bank
        if banks.keys() != set(active_lanes):
            raise ValueError("vX annotation requires one matching null directory per scored locus")
    if vx_tie_policy is None:
        vx_tie_policy = banks["paired"][1]["tie_policy"] if "paired" in banks else "abstain"
    mixed_union=(locus=='auto' and
        not vx_enabled and
        options.get('neighbour_weighting',defaults.neighbour_weighting)=='kernel' and
        paired_rows.height>0 and cohort.filter(pl.col('_population')=='single').height>0)
    marginal_tables={};service_batches=[];detail_columns={}
    from ..match.geometry import precompute_v_loop_distances
    v_loop_distances = (
        precompute_v_loop_distances(species=organism) if vx_enabled else
        precompute_v_loop_distances(species=organism,
            alphabet=options.get("alphabet", defaults.alphabet),
            matrix=options.get("distance_matrix", defaults.distance_matrix),
            v_regions=options.get("v_regions", defaults.v_regions))
    )
    normalized_reference = None
    def reference_loader():
        nonlocal normalized_reference
        if normalized_reference is None:
            normalized_reference = db.load(reference)
        return normalized_reference
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
            count_aliases = [c for c in service_part.columns if c.lower().lstrip("#") in
                            (*io.columns.ALIASES["count"], *io.columns._COUNT_FALLBACKS)]
            # Original aliases remain in the restoration mapping as source metadata.
            scoring_part = service_part.drop("_population", *wide_sequences, *count_aliases).with_columns(
                canonical_counts.gather(service_part["query_id"]).alias("count")
            )
            if vx_enabled:
                if population == "paired":
                    scoring_part.write_csv(sample, separator="\t")
                query, selection, binding, axis = _vx_population(
                    sample, scoring_part, reference_loader(), v_loop_distances, threads,
                    reference_population, options.get("exclude_exact", defaults.exclude_exact), population,
                    vx_tie_policy if population == "paired" else "abstain")
                calibration = {}
                if banks:
                    from hashlib import sha256 as digest_bytes
                    nulls, null_binding, source, null_input, null_paths = banks[population]
                    binding["reference_sha256"] = reference_digest
                    binding["target_axis_sha256"] = digest_bytes(axis.write_csv(separator="\t").encode()).hexdigest()
                    selection, recipe = _calibrate(
                        selection, binding, nulls, null_binding, source, vx_target_fpr, .95)
                    calibration = {"paired_vx_binding": binding, "calibration": recipe,
                                   "calibration_bank_replayed": any(report["source_sha256"] == null_input
                                                                    for report in reports),
                                   "null_inputs_sha256": {p.name: digest for p, digest in null_paths.items()}}
                mapping = (_pair_mapping(part, offset) if population == "paired" else
                           part.drop("_population").with_columns(
                               pl.col("query_id").cast(pl.String).alias("_service_query_id"),
                               pl.col("source_pair_id").alias("pair_id")))
                for kind, table in (("queries", query), ("selection", selection)):
                    tables.setdefault(kind, []).append(_restore(table, mapping))
                manifests.append({"score_model": "vx", "locus": population,
                                  "tie_policy": binding["tie_policy"],
                                  "reference": db.provenance(reference),
                                  "reference_population": reference_population,
                                  "exact_hits_excluded": options.get("exclude_exact", defaults.exclude_exact),
                                  "queries": query.height, **calibration})
                batches.append({"locus": population, "query_chain_rows": part.height,
                                "sample_count": part["sample_id"].n_unique()})
                continue
            scoring_part.write_csv(sample, separator="\t")
            args = _core_defaults(reference, local_prefix)
            vars(args).update(options)
            args.sample, args.locus, args.species = str(sample), population, species
            args.threads, args.selection_only = threads, selection_only
            if mixed_union and population=='paired':
                neighbours.main(args,marginal_tables=marginal_tables,v_loop_distances=v_loop_distances,
                    reference_loader=reference_loader)
            else:
                neighbours.main(args,v_loop_distances=v_loop_distances,reference_loader=reference_loader)
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
    if banks and (sha256(archive) != reference_digest or any(
            sha256(path) != digest for _, _, _, _, paths in banks.values() for path, digest in paths.items())):
        raise ValueError("vX score inputs changed during calibration")
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
        "score_models": sorted({manifest.get("score_model", "kernel") for manifest in manifests}),
        **({'mixed_chain_union':True,'service_batches':service_batches} if mixed_union else {}),
        **({'population_table_columns':detail_columns} if detail_columns else {}),
        "wall_seconds": time.perf_counter() - started}
    if "selection_sha256" in manifest:
        manifest["selection_sha256"] = sha256(Path(outputs["selection"]))
    manifest_path = Path(str(prefix) + ".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    result.update(manifest=manifest, manifest_path=str(manifest_path), output_paths=outputs)
    return result
