"""Replay the frozen mouse paired-chain base score from a VDJdb snapshot."""
from pathlib import Path

from . import cpu_budget


REGION_WEIGHTS = {"cdr1": 1 + 1 / 256, "cdr2": 3 / 4 + 1 / 256, "cdr25": 3 / 2 + 1 / 256}
PAIR_RULE = {"margin_per_chain": 6, "margin_penalty": 1 / 40, "second_hit_gap_weight": -1 / 10}
CHAINS = ("junction_aa_alpha", "junction_aa_beta")


def register(subparsers):
    parser = subparsers.add_parser("paired-vx-base", help=__doc__)
    parser.add_argument("queries", type=Path, help="Mouse paired AIRR query table")
    parser.add_argument("--vdjdb", type=Path, required=True, help="Local VDJdb release ZIP")
    parser.add_argument("--exclude-exact", action="store_true",
                        help="Benchmark only: exclude matching paired junctions from the hit lists")
    parser.add_argument("--null-source", help="Declare all input queries as one named calibration null")
    parser.add_argument("--withhold-peptides", type=Path, help="TSV of query_id, epitope to withhold across restrictions")
    parser.add_argument("--threads", type=cpu_budget, default="0")
    parser.add_argument("--output-dir", type=Path, required=True, help="Fresh directory for paired motif rank inputs")
    parser.set_defaults(func=main)


def _reference_catalogue(raw, model, *, background_memberships=False):
    import polars as pl
    from .historical import assayed_targets
    from ..aggregate.candidates import PMHC
    from ..db.schema import mhc_compatible
    from ..evalue.equal_tail import _eligible, _prepare
    from ..evalue.paired import build_paired_ref

    original = raw.columns[:raw.columns.index("complex_id")]
    if "complex.id" in original:
        raw = raw.unique(subset=original, maintain_order=True)
    pairs = build_paired_ref(raw).rename(
        {"alpha": "cdr3a", "beta": "cdr3b", "v_alpha": "va", "v_beta": "vb"})
    pairs = _prepare(pairs, "ab", model).filter(_eligible("ab", "native", 5))
    key = ["cdr3a", "cdr3b", "nativea", "nativeb"]
    fixed = pairs.select(key).unique(maintain_order=True).with_row_index("reference_key")
    pairs = pairs.join(fixed, on=key, validate="m:1", maintain_order="left")
    declared = (pairs.select(PMHC).unique().drop_nulls()
                .filter(pl.all_horizontal(pl.col(*PMHC).str.strip_chars() != ""))
                .sort(PMHC).with_row_index("query_id").with_columns(pl.lit("mouse").alias("species")))
    if declared.is_empty():
        raise ValueError("mouse reference has no eligible paired pMHC targets")
    groups, backgrounds = [], []
    for target, _, selected in assayed_targets(declared.select("query_id"), declared, pairs, "mouse"):
        members = selected["reference_key"].unique().sort().to_numpy()
        compatible = pairs.filter(
            (pl.col("epitope") != target["epitope"])
            & (pl.col("mhc_class") == target["mhc_class"])
            & mhc_compatible("mhc_a", target["mhc_a"])
            & mhc_compatible("mhc_b", target["mhc_b"])
        )
        groups.append(members)
        if background_memberships:
            backgrounds.append((compatible["reference_key"].to_numpy(), compatible["epitope"].to_numpy()))
        else:
            backgrounds.append(compatible["reference_key"].unique().to_numpy())
    return fixed, declared.select(PMHC), groups, backgrounds


def _withholding_folds(path, query_ids, reference_peptides):
    import polars as pl
    from ..io.airr import _read_table

    frame = _read_table(path)
    if {"query_id", "epitope"} - set(frame.columns):
        raise ValueError("withholding requires query_id and epitope columns")
    frame = frame.select(pl.col("query_id", "epitope").cast(pl.String))
    if set(frame["query_id"]) != set(query_ids) or not set(frame["epitope"]) <= reference_peptides:
        raise ValueError("withholding requires every query ID and only known reference peptides")
    by_query = dict(frame.group_by("query_id").agg(pl.col("epitope").unique().sort()).iter_rows())
    folds = {}
    for i, query_id in enumerate(query_ids):
        folds.setdefault(tuple(by_query[query_id]), []).append(i)
    return folds


def _withheld_scores(queries, fixed, targets, groups, backgrounds, context, *, folds, threads, exclude_exact):
    import numpy as np

    result = np.full((queries.height, targets.height), -np.inf)
    for excluded, rows in folds.items():
        columns = np.flatnonzero(~targets["epitope"].is_in(excluded).to_numpy())
        if not len(columns):
            continue
        remaining = []
        for column in columns:
            keys, peptide = backgrounds[column]
            remaining.append(np.unique(keys[~np.isin(peptide, excluded)]))
        result[np.ix_(rows, columns)] = _base_scores(
            queries[rows], fixed, [groups[column] for column in columns], remaining, context,
            threads=threads, exclude_exact=exclude_exact)
    return result


def _regional_corrections(context, locus):
    from ..match.geometry import _v_loop_groups

    alleles = context["tables"][locus][0]
    correction = 0
    for region, weight in REGION_WEIGHTS.items():
        partial = {gene: {name: entry[name] if name == region else "" for name in REGION_WEIGHTS}
                   for gene, entry in context["model"].items()}
        distance, indices = _v_loop_groups(alleles.tolist(), partial, context["costs"])
        correction = correction + (weight - 1) * distance[indices[:, None], indices]
    return alleles, correction


def _base_scores(queries, fixed, groups, backgrounds, context, *, threads, exclude_exact):
    import numpy as np
    from ..evalue.competitive import reduced_hits
    from ..evalue.equal_tail import _eligible
    from ..match.geometry import paired_distance_matrix

    result = np.full((queries.height, len(groups)), -np.inf)
    eligible = queries.with_row_index("_i").filter(_eligible("ab", "native", 5))
    tables = [_regional_corrections(context, locus) for locus in ("TRA", "TRB")]
    refs = [fixed["cdr3" + chain].to_list() for chain in "ab"]
    ref_alleles = [fixed["native" + chain].to_list() for chain in "ab"]
    ri = [np.searchsorted(axis, values) for (axis, _), values in zip(tables, ref_alleles, strict=True)]
    for part in eligible.iter_slices(max(1, 4_000_000 // fixed.height)):
        if not part.height:
            continue
        sequences = [part["cdr3" + chain].to_list() for chain in "ab"]
        alleles = [part["native" + chain].to_list() for chain in "ab"]
        native = paired_distance_matrix(*sequences, *refs, *alleles, *ref_alleles, species="mouse",
                                        threads=threads, v_loop_distances=context, allow_short=True)
        native = native.astype(float)
        if exclude_exact:
            exact = ((np.asarray(sequences[0])[:, None] == np.asarray(refs[0])[None, :])
                     & (np.asarray(sequences[1])[:, None] == np.asarray(refs[1])[None, :]))
            native[exact] = np.inf
        correction = np.zeros(native.shape)
        for (axis, table), values, j in zip(tables, alleles, ri, strict=True):
            correction += table[np.ix_(np.searchsorted(axis, values), j)]
        hits = native + correction
        rows = part["_i"].to_numpy()
        for column, (members, background) in enumerate(zip(groups, backgrounds, strict=True)):
            competitor = native[:, background].min(axis=1) if len(background) else np.full(len(part), np.inf)
            index = min(1, len(members) - 1)
            top = np.partition(hits[:, members], index, axis=1)[:, :index + 1]
            first = top[:, 0]
            supported = np.isfinite(first)
            second = top[:, 1] if index else first
            second = np.where(np.isfinite(second), second, first)
            result[rows[supported], column] = reduced_hits(
                first[supported], second[supported], competitor[supported], PAIR_RULE, 2)
    return result


def _maxima_table(query_ids, base, labels, null_source):
    import numpy as np
    import polars as pl

    maximum = base.max(axis=1)
    finite = np.isfinite(maximum)
    tied = finite & (np.count_nonzero(base == maximum[:, None], axis=1) > 1)
    unique = finite & ~tied
    target = np.asarray(labels, dtype=object)[base.argmax(axis=1)]
    table = pl.DataFrame({
        "query_id": query_ids,
        "target": pl.Series(np.where(unique, target, None).tolist(), dtype=pl.String),
        "maximum": pl.Series(np.where(finite, maximum, np.nan)).fill_nan(None),
        "score": pl.Series(np.where(unique, maximum, np.nan)).fill_nan(None),
        "tied": tied,
    })
    if null_source is not None:
        table = table.with_columns(pl.lit("null").alias("outcome"),
                                   pl.lit(null_source).alias("null_source"))
    return table


def main(args):
    import json
    from importlib.metadata import version
    from importlib.resources import files

    import numpy as np
    import polars as pl
    from .equal_tail import _queries
    from .. import db
    from ..db.cache import sha256
    from ..evalue.equal_tail import _prepare
    from ..match.geometry import precompute_v_loop_distances

    out = args.output_dir
    if out.exists() or out.is_symlink():
        raise ValueError("paired-vx-base requires a fresh output directory")
    if not args.vdjdb.is_file() or args.vdjdb.suffix != ".zip":
        raise ValueError("paired-vx-base requires a local VDJdb release ZIP")
    null_source = args.null_source.strip() if args.null_source is not None else None
    if null_source == "":
        raise ValueError("--null-source must name a nonempty null population")
    inputs = {"reference": args.vdjdb, "queries": args.queries}
    if args.withhold_peptides is not None:
        inputs["withhold_peptides"] = args.withhold_peptides
        if null_source is None:
            null_source = "withheld"
    digests = {name: sha256(path) for name, path in inputs.items()}
    context = precompute_v_loop_distances(species="mouse")
    raw = db.load(args.vdjdb, species="MusMusculus")
    fixed, targets, groups, backgrounds = _reference_catalogue(
        raw, context["model"], background_memberships=args.withhold_peptides is not None)
    queries, ingestion = _queries(args.queries, "paired", "mouse", declared=False)
    if queries.is_empty() or queries.height * 2 != ingestion["output_rows"] or not queries.select(pl.all_horizontal(
        pl.col("cdr3" + chain).str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$").fill_null(False) for chain in "ab"
    ).all()).item():
        raise ValueError("provide complete mouse pairs with both canonical junctions")
    queries = _prepare(queries, "ab", context["model"])
    folds = None
    if args.withhold_peptides is None:
        base = _base_scores(queries, fixed, groups, backgrounds, context,
                            threads=args.threads, exclude_exact=args.exclude_exact)
    else:
        folds = _withholding_folds(args.withhold_peptides, queries["query_id"],
                                   set(raw["epitope"].drop_nulls()) - {""})
        base = _withheld_scores(queries, fixed, targets, groups, backgrounds, context,
                                folds=folds, threads=args.threads, exclude_exact=args.exclude_exact)
    labels = [json.dumps(target, separators=(",", ":")) for target in targets.iter_rows()]
    maxima = _maxima_table(queries["query_id"], base, labels, null_source)
    if digests != {name: sha256(path) for name, path in inputs.items()}:
        raise ValueError("paired-vx-base inputs changed during replay")
    out.mkdir(parents=True, exist_ok=False)
    np.savez(out / "base.npz", base=base)
    queries.select("query_id", pl.col("cdr3a").alias(CHAINS[0]),
                   pl.col("cdr3b").alias(CHAINS[1])).write_csv(out / "queries.tsv", separator="\t")
    pl.DataFrame({"target": labels, "mhc_class": targets["mhc_class"]}).write_csv(
        out / "targets.tsv", separator="\t")
    if folds is None:
        member_rows = pl.DataFrame({
            "reference_key": np.concatenate(groups),
            "target": np.repeat(labels, [len(group) for group in groups]),
        }).join(fixed.select("reference_key", pl.col("cdr3a").alias(CHAINS[0]),
                            pl.col("cdr3b").alias(CHAINS[1])), on="reference_key",
                validate="m:1", maintain_order="left")
        member_rows.write_csv(out / "memberships.tsv", separator="\t")
        fold_status = {"motif_reference_keys": member_rows["reference_key"].n_unique()}
    else:
        fold_status = {"withholding_folds": len(folds)}
    maxima.write_csv(out / "maxima.tsv", separator="\t")
    (out / "manifest.json").write_text(json.dumps({
        "inputs_sha256": digests, "species": "mouse", "exact_hits_excluded": args.exclude_exact,
        "queries": queries.height, "targets": targets.height, "reference_keys": fixed.height,
        "supported_queries": maxima["maximum"].count(), "unique_winners": maxima["score"].count(),
        "top_ties": int(maxima["tied"].sum()), "null_source": null_source,
        **fold_status,
        "threads": args.threads,
        "software": {package: version(package) for package in ("vdjmatch", "seqtree", "vdjtools")},
        "v_loop_model_sha256": sha256(Path(str(files("vdjmatch.resources") / "vloops" / "mouse_v_loops.tsv"))),
        "source_sha256": sha256(Path(__file__)),
        "region_weights": REGION_WEIGHTS, "competitive_rule": PAIR_RULE,
        "outputs_sha256": {file.name: sha256(file) for file in sorted(out.iterdir())},
    }, indent=2) + "\n")
    return 0
