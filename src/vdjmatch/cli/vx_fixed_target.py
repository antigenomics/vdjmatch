"""Score declared target rectangles produced by a native comparator."""

import json
from pathlib import Path

CANONICAL = frozenset("ACDEFGHIKLMNPQRSTVWY")


def register(subparsers):
    parser = subparsers.add_parser("vx-fixed-target", help=__doc__)
    parser.add_argument("rectangle", type=Path, help="Aligned native distances, receptor axes and background in NPZ")
    parser.add_argument("--locus", required=True, choices=("TRA", "TRB", "paired"))
    parser.add_argument("--species", default="human", choices=("human", "mouse"))
    parser.add_argument("--exclude-exact", action="store_true",
                        help="Mask exact junction hits, requiring both junctions for paired cases")
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.set_defaults(func=main)


def main(args):
    import numpy as np
    import polars as pl

    from .paired_vx_base import CHAINS, RULES, _regional_corrections
    from ..db.cache import sha256
    from ..evalue.vx_fixed_target import fixed_target_scores
    from ..match.geometry import precompute_v_loop_distances, resolve_v_alleles
    from ..runner.annotation import _validate_output_paths

    inputs = args.rectangle.resolve()
    prefix = _validate_output_paths(args.output_prefix, [inputs], suffixes=(".scores.tsv", ".manifest.json"))
    scores_path = Path(str(prefix) + ".scores.tsv")
    manifest_path = Path(str(prefix) + ".manifest.json")
    if scores_path.exists() or manifest_path.exists():
        raise FileExistsError("vx-fixed-target requires a fresh output prefix")
    digest = sha256(inputs)
    with np.load(inputs, allow_pickle=False) as data:
        fields = ("query_id", "reference_key", "native", "query_junction", "reference_junction",
                  "query_v", "reference_v", "background", "source_log2")
        if set(data.files) != set(fields):
            raise ValueError("rectangle requires " + ", ".join(fields))
        ids, keys, native, q_seq, r_seq, q_v, r_v, background, source = (data[name] for name in fields)

    if native.ndim != 2:
        raise ValueError("native distances must be a matrix")
    queries, references = native.shape
    chains = CHAINS[args.locus]
    receptor_axes = ((q_seq, queries), (r_seq, references), (q_v, queries), (r_v, references))
    if (ids.shape != (queries,) or keys.shape != (references,) or
            any(axis.shape != (size, len(chains)) or axis.dtype.kind != "U" for axis, size in receptor_axes) or
            background.shape != (queries,) or source.shape != (references,)):
        raise ValueError("query, reference and distance rectangle axes must agree")
    if any(axis.dtype.kind != "U" or np.any(axis == "") or len(np.unique(axis)) != len(axis)
           for axis in (ids, keys)):
        raise ValueError("query and reference identifiers must be nonempty unique Unicode strings")
    if (native.dtype.kind not in "fiu" or background.dtype.kind not in "fiu" or source.dtype.kind not in "fiu" or
            np.isnan(native).any() or np.isnan(background).any() or not np.isfinite(source).all() or
            np.any(native < 0) or np.any(background < 0) or np.any(source < 0)):
        raise ValueError("invalid native distances, background or source support")

    context = precompute_v_loop_distances(species=args.species)
    valid = []
    indices = []
    for junctions, genes in ((q_seq, q_v), (r_seq, r_v)):
        keep = np.ones(len(junctions), dtype=bool)
        alleles = []
        for column, chain in enumerate(chains):
            locus = "TRA" if chain == "a" else "TRB"
            resolved = list(resolve_v_alleles(genes[:, column].tolist(), context["model"], locus=locus))
            keep &= [gene is not None and len(seq) >= 5 and set(seq) <= CANONICAL
                     for gene, seq in zip(resolved, junctions[:, column], strict=True)]
            alleles.append(resolved)
        valid.append(np.flatnonzero(keep))
        indices.append(alleles)
    correction = np.zeros(native.shape)
    for column, chain in enumerate(chains):
        axis, region_cost = _regional_corrections(context, "TRA" if chain == "a" else "TRB")
        qi, ri = (np.searchsorted(axis, [a[column][i] for i in rows])
                  for a, rows in zip(indices, valid, strict=True))
        correction[np.ix_(*valid)] += region_cost[np.ix_(qi, ri)]

    exact = np.full(native.shape, args.exclude_exact, dtype=bool)
    if args.exclude_exact:
        for column in range(len(chains)):
            exact &= np.equal.outer(q_seq[:, column], r_seq[:, column])
    native_scores, vx_scores = fixed_target_scores(
        native, correction, exact, background, source, RULES[args.locus], len(chains))
    if sha256(inputs) != digest:
        raise ValueError("fixed target rectangle changed during scoring")

    def finite(values):
        return pl.Series(values).replace(-np.inf, None)

    prefix.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame({"query_id": ids, "native_score": finite(native_scores), "score": finite(vx_scores)}).write_csv(
        scores_path, separator="\t")
    manifest_path.write_text(json.dumps({
        "score_model": "vx", "locus": args.locus, "species": args.species, "rectangle_sha256": digest,
        "queries": queries, "target_keys": references, "score_sha256": sha256(scores_path),
        "exact_hits_excluded": args.exclude_exact,
    }, indent=2) + "\n")
    return 0
