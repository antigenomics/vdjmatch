"""Fresh amino-acid batch search with original row identities and explicit counts.

seqtm applies per-edit caps and reports edit decomposition. seqtrie uses only the
total unit-cost radius and reports null decomposition; its exact hits have score0.
"""

from __future__ import annotations

import importlib.metadata
import json
from pathlib import Path
import resource
import sys
import time


def register(subparsers):
    """Register without initializing the numerical runtime or its thread pools."""
    parser = subparsers.add_parser("search", help="native amino-acid batch search")
    parser.add_argument("query")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--sequence-col", required=True)
    parser.add_argument("--reference-sequence-col")
    parser.add_argument("--same-key-col", help="keep hits with equal non-null keys before output/rescoring")
    parser.add_argument("--scope", default="1,0,0,1")
    parser.add_argument("--engine", choices=["seqtm", "seqtrie"], default="seqtm")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--matrix", choices=["none", "vdjam", "blosum62"], default="none"
    )
    parser.add_argument(
        "--position-significance",
        action="store_true",
        help="end-anchored positional matrix weighting (seqtm, no indels)",
    )
    parser.add_argument("--exclude-exact", action="store_true")
    parser.add_argument("--counts-only", action="store_true")
    parser.add_argument("--output-prefix", required=True)
    parser.set_defaults(func=main)
    parser.epilog = "seqtrie uses total unit-cost radius; per-edit caps are ignored and edit decomposition is unavailable."
    return parser


def _sequences(path, column, role):
    import polars as pl
    from ..io.airr import _read_table

    frame = _read_table(path)
    if column not in frame.columns:
        raise ValueError(f"{role} is missing sequence column {column!r}")
    sequences = frame[column].cast(pl.String)
    null = sequences.null_count()
    invalid = int((~sequences.str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$")).sum())
    if null or invalid:
        raise ValueError(
            f"{role} has {null + invalid} invalid amino-acid rows: "
            f"{null} null, {invalid} empty/nonstandard; {len(sequences)} total rows"
        )
    return sequences


def _positional_scores(queries, references, pairs, matrix):
    """Vectorized integer rescoring; loops build model LUTs, never score pairs."""
    import numpy as np
    from seqtree import PositionalMatrix, amino_acids
    from ..match.regions import load_significance, significance_weights

    qrows = pairs["query_row"].to_numpy()
    rrows = pairs["reference_row"].to_numpy()
    scores = np.zeros(len(pairs), dtype=np.int64)
    if not len(pairs):
        return scores
    width = max(queries.str.len_bytes().max(), references.str.len_bytes().max())
    alphabet = amino_acids()
    ascii_lut = np.zeros(128, dtype=np.uint8)
    ascii_lut[np.frombuffer(alphabet.encode("ascii"), dtype=np.uint8)] = np.arange(
        len(alphabet)
    )
    query_bytes = (
        np.asarray(queries.to_list(), dtype=f"S{width}")
        .view(np.uint8)
        .reshape(-1, width)
    )
    reference_bytes = (
        np.asarray(references.to_list(), dtype=f"S{width}")
        .view(np.uint8)
        .reshape(-1, width)
    )
    qcodes, rcodes = ascii_lut[query_bytes], ascii_lut[reference_bytes]
    lengths = queries.str.len_bytes().to_numpy()[qrows]
    sig = load_significance()
    for length in np.unique(lengths):
        length = int(length)
        selected = np.flatnonzero(lengths == length)
        weights = [max(1, round(100 * w)) for w in significance_weights(length, sig)]
        positional = PositionalMatrix.from_weights(matrix, weights)
        lut = np.asarray(
            [
                [
                    [positional.penalty(pos, a, b) for b in range(len(alphabet))]
                    for a in range(len(alphabet))
                ]
                for pos in range(length)
            ],
            dtype=np.int64,
        )
        # One vectorized operation per position avoids a pair-by-length temporary.
        for pos in range(length):
            scores[selected] += lut[
                pos, qcodes[qrows[selected], pos], rcodes[rrows[selected], pos]
            ]
    return scores


def main(args):
    started = time.perf_counter()
    if isinstance(args.threads, bool) or args.threads < 1:
        raise ValueError("threads must be a positive integer")
    if args.engine == "seqtrie" and args.matrix != "none":
        raise ValueError(
            "matrix scoring requires --engine seqtm; seqtrie uses unit edit costs"
        )
    import polars as pl
    from seqtree import Index, SubstitutionMatrix
    from ..db.cache import sha256
    from ..match.scope import parse_scope, search_params
    from ..match.scoring import DEFAULT_SCALE, load_vdjam

    reference_col = args.reference_sequence_col or args.sequence_col
    queries = _sequences(args.query, args.sequence_col, "query")
    references = _sequences(args.reference, reference_col, "reference")
    key_column = getattr(args, "same_key_col", None)
    query_keys = reference_keys = None
    if key_column:
        from ..io.airr import _read_table

        keys = []
        for path, role in ((args.query, "query"), (args.reference, "reference")):
            table = _read_table(path)
            if key_column not in table.columns or table[key_column].null_count():
                raise ValueError(f"{role} requires non-null {key_column!r} keys")
            keys.append(table[key_column].cast(pl.String).to_list())
        query_keys, reference_keys = keys
    positional = getattr(args, "position_significance", False)
    _, ins, dels, _ = parse_scope(args.scope)
    if positional and (args.engine != "seqtm" or ins or dels or args.matrix == "none"):
        raise ValueError(
            "position significance requires seqtm, no indels and an explicit matrix"
        )
    matrix = (
        load_vdjam()
        if args.matrix == "vdjam"
        else SubstitutionMatrix.blosum62()
        if args.matrix == "blosum62"
        else ""
    )
    params = search_params(
        args.scope,
        engine=args.engine,
        matrix="" if positional else matrix,
        gap_open=DEFAULT_SCALE if args.matrix == "vdjam" else 1,
        gap_extend=DEFAULT_SCALE if args.matrix == "vdjam" else 1,
    )
    prefix = Path(args.output_prefix)
    paths = {
        kind: Path(str(prefix) + suffix)
        for kind, suffix in (
            ("counts", ".counts.tsv"),
            ("pairs", ".pairs.tsv"),
            ("manifest", ".manifest.json"),
        )
    }
    if {Path(args.query).resolve(), Path(args.reference).resolve()} & {
        p.resolve() for p in paths.values()
    }:
        raise ValueError("output prefix would overwrite an input table")
    tick = time.perf_counter()
    index = Index.build(references.to_list(), "aa")
    build_seconds = time.perf_counter() - tick
    if len(index) != len(references):
        raise RuntimeError("native index changed reference row cardinality")
    tick = time.perf_counter()
    hits = index.search_batch(queries.to_list(), params, args.threads)
    search_seconds = time.perf_counter() - tick

    def accepted(row, h):
        return (query_keys is None or query_keys[row] == reference_keys[h.ref_id]) and (not args.exclude_exact or (
            h.score > 0
            if args.engine == "seqtrie"
            else h.n_subs + h.n_ins + h.n_dels > 0
        ))

    counts = pl.DataFrame(
        {
            "query_row": pl.Series(range(len(queries)), dtype=pl.UInt32),
            "n_hits": pl.Series(
                [
                    sum(accepted(row, h) for h in hs) if args.exclude_exact or key_column else len(hs)
                    for row, hs in enumerate(hits)
                ],
                dtype=pl.UInt64,
            ),
        }
    )
    prefix.parent.mkdir(parents=True, exist_ok=True)
    counts.write_csv(paths["counts"], separator="\t")
    rescore_seconds = 0.0
    if not args.counts_only:
        pairs = pl.DataFrame(
            [
                (
                    row,
                    h.ref_id,
                    h.score,
                    h.n_subs if args.engine == "seqtm" else None,
                    h.n_ins if args.engine == "seqtm" else None,
                    h.n_dels if args.engine == "seqtm" else None,
                )
                for row, hs in enumerate(hits)
                for h in hs
                if accepted(row, h)
            ],
            schema={
                "query_row": pl.UInt32,
                "reference_row": pl.UInt32,
                "score": pl.Int64,
                "n_subs": pl.UInt16,
                "n_ins": pl.UInt16,
                "n_dels": pl.UInt16,
            },
            orient="row",
        )
        if positional:
            tick = time.perf_counter()
            pairs = pairs.with_columns(
                pl.Series(
                    "score", _positional_scores(queries, references, pairs, matrix)
                )
            )
            rescore_seconds = time.perf_counter() - tick
        pairs.write_csv(paths["pairs"], separator="\t")
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    options = {
        name: getattr(args, name)
        for name in (
            "sequence_col",
            "scope",
            "engine",
            "threads",
            "matrix",
            "position_significance",
            "exclude_exact",
            "counts_only",
        )
    }
    options["reference_sequence_col"] = reference_col
    options["same_key_col"] = key_column
    from importlib import resources

    model_inputs = {}
    if positional:
        model_inputs["position_significance"] = sha256(
            Path(
                resources.files("vdjmatch.resources")
                / "trimming"
                / "position_significance.tsv"
            )
        )
    if args.matrix == "vdjam":
        model_inputs["vdjam"] = sha256(
            Path(resources.files("vdjmatch.resources") / "vdjam.txt")
        )
    manifest = {
        "operation": "native-aa-search",
        "model_input_sha256": model_inputs,
        "score_contract": "positional significance integer penalty"
        if positional
        else "native engine score",
        "scores_reported": not args.counts_only,
        "positional_weight_scale": 100 if positional else None,
        "native_scope_contract": "per-edit caps"
        if args.engine == "seqtm"
        else "total unit-cost radius; per-edit caps ignored",
        "edit_decomposition_available": args.engine == "seqtm",
        "parameters": options,
        "sequence_conventions": {
            role: "declared junction_aa"
            if column == "junction_aa"
            else "generic amino-acid sequence"
            for role, column in (
                ("query", args.sequence_col),
                ("reference", reference_col),
            )
        },
        "inputs": {
            role: {"sha256": sha256(Path(path))}
            for role, path in (("query", args.query), ("reference", args.reference))
        },
        "software": {
            name: importlib.metadata.version(name)
            for name in (
                ("vdjmatch", "seqtree", "polars", "numpy")
                if positional
                else ("vdjmatch", "seqtree", "polars")
            )
        },
        "query_rows": len(queries),
        "reference_rows": len(references),
        "reference_unique_sequences": references.n_unique(),
        "hit_pairs": int(counts["n_hits"].sum()),
        "native_candidate_pairs": sum(len(hs) for hs in hits),
        "key_filter_stage": "after native batch, before output/rescoring" if key_column else None,
        "resources": {
            "index_build_seconds": build_seconds,
            "search_seconds": search_seconds,
            "positional_rescore_seconds": rescore_seconds,
            "wall_seconds": time.perf_counter() - started,
            "process_peak_rss_bytes": peak if sys.platform == "darwin" else peak * 1024,
            "native_threads": args.threads,
            "polars_pool_threads": pl.thread_pool_size(),
            "rss_scope": "CLI process; no subprocesses created",
            "native_batches": 1,
        },
        "outputs": {
            kind: str(path)
            for kind, path in paths.items()
            if kind != "pairs" or not args.counts_only
        },
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return 0
