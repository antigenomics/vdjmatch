"""Calibrate and annotate vX calls with a matched null source."""
from pathlib import Path


def register(subparsers):
    parser = subparsers.add_parser("vx-evalue", aliases=["paired-vx-evalue"], help=__doc__)
    parser.add_argument("queries", type=Path, help="Ordinary vx-base output directory")
    parser.add_argument("--null-dir", type=Path, required=True, help="Matching vx-base output with --null-source")
    parser.add_argument("--target-fpr", type=float, required=True)
    parser.add_argument("--confidence", type=float, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.set_defaults(func=main)


def _source(directory):
    from hashlib import sha256 as digest_bytes
    import json
    import polars as pl
    from .paired_vx_base import LITERAL_MIN_SCORE
    from ..db.cache import sha256

    manifest_path = directory / "manifest.json"
    maxima_path = directory / "maxima.tsv"
    targets_path = directory / "targets.tsv"
    source = manifest_path.read_bytes()
    manifest = json.loads(source)
    if not isinstance(manifest, dict) or not isinstance(manifest.get("inputs_sha256"), dict):
        raise ValueError("incomplete paired-vx-base manifest")
    if (manifest.get("species") not in ("mouse", "human") or "withholding_folds" in manifest
            or "withhold_peptides" in manifest.get("inputs_sha256", {})):
        raise ValueError("ordered vX calibration requires ordinary score populations")
    try:
        keys = ("species", "exact_hits_excluded", "reference_keys", "region_weights", "competitive_rule",
                "software", "v_loop_model_sha256", "source_sha256")
        binding = {key: manifest[key] for key in keys}
        binding["locus"] = manifest.get("locus", "paired")
        binding["reference_population"] = manifest.get("reference_population", "curated-benchmark")
        binding["tie_policy"] = manifest.get("tie_policy", "abstain")
        binding["literal_tie_min_score"] = manifest.get("literal_tie_min_score")
        if binding["tie_policy"] == "literal":
            if binding["locus"] != "paired" or binding["literal_tie_min_score"] != LITERAL_MIN_SCORE:
                raise ValueError("unknown literal tie gate for this vX locus")
        elif binding["tie_policy"] != "abstain" or binding["literal_tie_min_score"] is not None:
            raise ValueError("unknown vX tie policy")
        if binding["locus"] not in ("paired", "TRA", "TRB") or binding["reference_population"] not in (
                "curated-benchmark", "all-valid"):
            raise ValueError("unknown vX locus or reference population")
        binding["reference_sha256"] = manifest["inputs_sha256"]["reference"]
        query_input_sha256 = manifest["inputs_sha256"]["queries"]
        signed = {manifest_path: digest_bytes(source).hexdigest()}
        for path in (maxima_path, targets_path):
            signed[path] = sha256(path)
            if signed[path] != manifest["outputs_sha256"][path.name]:
                raise ValueError(f"{path.name} differs from its paired-vx-base manifest")
        binding["target_axis_sha256"] = manifest["outputs_sha256"][targets_path.name]
    except (KeyError, TypeError) as exc:
        raise ValueError("incomplete paired-vx-base manifest") from exc
    if not isinstance(binding["exact_hits_excluded"], bool):
        raise ValueError("paired-vx-base manifest must declare its exact-hit policy")
    columns = {"query_id": pl.String, "target": pl.String, "maximum": pl.Float64,
               "score": pl.Float64, "tied": pl.Boolean}
    if binding["tie_policy"] == "literal":
        columns.update(literal_target=pl.String, tie_resolved=pl.Boolean)
    frame = pl.read_csv(maxima_path, separator="\t", infer_schema_length=0, schema_overrides=columns)
    if set(columns) - set(frame.columns):
        raise ValueError("paired-vx-base maxima require the complete winner and tie columns")
    if frame.height != manifest.get("queries"):
        raise ValueError("maxima row count differs from its paired-vx-base manifest")
    score, maximum, tied = (frame[column] for column in ("score", "maximum", "tied"))
    supported = score.is_not_null()
    if binding["tie_policy"] == "literal":
        literal, resolved = frame["literal_target"], frame["tie_resolved"]
        eligible_tie = literal.is_not_null()
        expected = eligible_tie & maximum.ge(LITERAL_MIN_SCORE).fill_null(False)
        if ((eligible_tie & ~tied).any() or not resolved.equals(expected)
                or not frame["target"].filter(resolved).equals(literal.filter(resolved))):
            raise ValueError("maxima do not match the literal tie policy")
    else:
        resolved = pl.Series([False] * frame.height, dtype=pl.Boolean)
    if (tied.null_count() or not (tied & supported).equals(resolved)
            or not frame["target"].is_not_null().equals(supported)
            or not maximum.is_not_null().equals(tied | supported)
            or not maximum.drop_nulls().is_finite().all()
            or not score.filter(supported).equals(maximum.filter(supported))):
        raise ValueError("maxima must preserve unsupported queries and exact winner ties")
    return frame, binding, manifest.get("null_source"), query_input_sha256, signed


def _calibrate(queries, binding, nulls, null_binding, source, target_fpr, confidence):
    from .calibration import apply_calibration, null_tail_probabilities
    from .ordered_calls import fit_ordered_calls

    if binding != null_binding:
        raise ValueError("query and null vX inputs must match in locus, reference, model, scorer, tie and exact-hit policy")
    if not source or "null_source" not in nulls.columns or set(nulls["null_source"]) != {source}:
        raise ValueError("null population must match its paired-vx-base --null-source declaration")
    recipe = fit_ordered_calls(nulls, target_fpr, confidence)
    accepted = apply_calibration(queries, recipe)
    calls, tails = null_tail_probabilities(accepted, nulls)
    return calls.with_columns(tails["null_queries"], tails["null_tail_count"]), recipe


def main(args):
    import json
    from ..db.cache import sha256
    from ..runner.annotation import _validate_output_paths

    queries, binding, _, query_input, query_paths = _source(args.queries)
    nulls, null_binding, source, null_input, null_paths = _source(args.null_dir)
    calls, recipe = _calibrate(queries, binding, nulls, null_binding, source, args.target_fpr, args.confidence)
    prefix = _validate_output_paths(args.output_prefix, [*query_paths, *null_paths],
                                    suffixes=(".calls.tsv", ".manifest.json"))
    if any(sha256(path) != expected for signed in (query_paths, null_paths) for path, expected in signed.items()):
        raise ValueError("vX score inputs changed during calibration")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    calls_path = Path(str(prefix) + ".calls.tsv")
    calls.write_csv(calls_path, separator="\t")
    Path(str(prefix) + ".manifest.json").write_text(json.dumps({
        "paired_vx_binding": binding,
        "query_manifest_sha256": query_paths[args.queries / "manifest.json"],
        "null_manifest_sha256": null_paths[args.null_dir / "manifest.json"],
        "query_maxima_sha256": query_paths[args.queries / "maxima.tsv"],
        "null_maxima_sha256": null_paths[args.null_dir / "maxima.tsv"],
        "calibration_bank_replayed": query_input == null_input,
        "queries": calls.height,
        "accepted": int((calls["call_status"] == "accepted").sum()),
        "calibration": recipe, "calls_sha256": sha256(calls_path),
    }, indent=2, allow_nan=False) + "\n")
    return 0
