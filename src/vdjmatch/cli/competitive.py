"""Replay vX research scores and fixed-reference rankings."""

from pathlib import Path

from . import cpu_budget


def register(subparsers):
    motif = subparsers.add_parser("anchored-motif-score", help="Score TRB junctions against a fixed motif catalogue")
    motif.add_argument("--reference", type=Path, required=True,
                       help="TSV with target, junction_aa, v_call and j_call")
    motif.add_argument("--queries", type=Path, required=True, help="TSV with query_id and junction_aa")
    motif.add_argument("--motif-length", type=int, choices=(3, 4), default=3,
                       help="Anchored motif length (default: 3)")
    motif.add_argument("--output", type=Path, required=True)
    motif.set_defaults(func=motif_main)

    assay = subparsers.add_parser("paired-assay-score", help="Rank one annotated paired-chain assay group")
    assay.add_argument("--distances", type=Path, required=True, help="Square NumPy .npy distance matrix")
    assay.add_argument("--queries", type=Path, required=True, help="TSV with query_id in matrix order")
    assay.add_argument("--exact-hits", type=Path, help="Optional aligned boolean .npy mask")
    assay.add_argument("--output", type=Path, required=True)
    assay.set_defaults(func=assay_main)

    paired = subparsers.add_parser("paired-motif-rank", help="Refine paired ranks within a fixed target context")
    paired.add_argument("--scores", type=Path, required=True,
                        help="NumPy .npz with base and, without --reference, tra, trb and target_key_counts")
    paired.add_argument("--queries", type=Path, required=True, help="TSV with query_id in matrix row order")
    paired.add_argument("--targets", type=Path, required=True, help="TSV with full target ID in matrix column order")
    paired.add_argument("--reference", type=Path,
                        help="Fixed paired reference memberships for scoring motifs separately by MHC class")
    paired.add_argument("--support-weight", type=float, required=True)
    paired.add_argument("--conflict-weight", type=float, required=True)
    paired.add_argument("--output", type=Path, required=True)
    paired.set_defaults(func=paired_main)

    mouse = subparsers.add_parser("mouse-paired-rank",
                                  help="Rank paired mouse queries against a fixed VDJdb release")
    mouse.add_argument("queries", type=Path, help="Paired AIRR query table")
    mouse.add_argument("--vdjdb", type=Path, required=True, help="Local VDJdb release ZIP")
    mouse.add_argument("--exclude-exact", action="store_true",
                       help="Benchmark only: exclude matching paired junctions from score hit lists")
    mouse.add_argument("--threads", type=cpu_budget, default="0")
    mouse.add_argument("--output-dir", type=Path, required=True)
    mouse.set_defaults(func=mouse_paired_main)


def mouse_paired_main(args):
    import json
    from tempfile import TemporaryDirectory
    from types import SimpleNamespace

    from .paired_vx_base import main as score_base
    from ..db.cache import sha256
    from ..evalue import competitive as reducer

    out = args.output_dir
    if out.exists() or out.is_symlink():
        raise ValueError("mouse-paired-rank requires a fresh output directory")
    out.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=out.parent, prefix=".mouse-paired-rank-") as temporary:
        staged = Path(temporary) / "results"
        score_base(SimpleNamespace(queries=args.queries, vdjdb=args.vdjdb, output_dir=staged,
                                   exclude_exact=args.exclude_exact, null_source=None,
                                   withhold_peptides=None, threads=args.threads, species="mouse",
                                   locus="paired", reference_population="curated-benchmark",
                                   maxima_only=False, tie_policy="abstain"))
        rank_path = staged / "ranking.tsv"
        weights = {"support_weight": 0, "conflict_weight": 4}
        paired_main(SimpleNamespace(scores=staged / "base.npz", queries=staged / "queries.tsv",
                                    targets=staged / "targets.tsv", reference=staged / "memberships.tsv",
                                    output=rank_path, **weights))
        (staged / "ranking.manifest.json").write_text(json.dumps({
            "base_manifest_sha256": sha256(staged / "manifest.json"),
            "ranking_sha256": sha256(rank_path),
            "ranking_source_sha256": sha256(Path(__file__)),
            "motif_reducer_source_sha256": sha256(Path(reducer.__file__)),
            **weights,
        }, indent=2) + "\n")
        # Claim the destination exclusively and publish the signed manifests last.
        out.mkdir()
        for path in sorted(staged.iterdir(), key=lambda path: path.suffix == ".json"):
            path.rename(out / path.name)
    return 0


def _read(path, columns):
    import polars as pl

    frame = pl.read_csv(path, separator="\t", infer_schema_length=0).select(*columns)
    if not frame.height or frame.select(pl.any_horizontal(pl.all().is_null() | (pl.all() == "")).any()).item():
        raise ValueError(f"expected nonempty, nonmissing {', '.join(columns)}")
    return frame


def _output(a, inputs, result):
    from ..runner.annotation import _validate_output_paths

    path = _validate_output_paths(a.output, inputs, suffixes=("",))
    path.parent.mkdir(parents=True, exist_ok=True)
    result.write_csv(path, separator="\t")
    return 0


def motif_main(a):
    import numpy as np
    import polars as pl
    from ..evalue.competitive import anchored_motif_scores

    reference = _read(a.reference, ("target", "junction_aa", "v_call", "j_call")).unique()
    query = _read(a.queries, ("query_id", "junction_aa"))
    for frame in (reference, query):
        valid = frame["junction_aa"].str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
        if not valid.all():
            raise ValueError("junction_aa must contain canonical amino acids")
    if query["query_id"].n_unique() != query.height:
        raise ValueError("query identifiers must be unique")
    targets = reference["target"].unique().sort().to_list()
    values = anchored_motif_scores(reference.select("target", "junction_aa").iter_rows(),
                                   query["junction_aa"], targets, motif_length=a.motif_length)
    maxima = (values == values.max(axis=1, keepdims=True)) & np.isfinite(values)
    result = pl.DataFrame({
        "query_id": np.repeat(query["query_id"].to_numpy(), len(targets)),
        "target": np.tile(targets, query.height),
        "score": values.ravel(),
        "winner": (maxima & (maxima.sum(axis=1, keepdims=True) == 1)).ravel(),
        "tied_maximum": (maxima & (maxima.sum(axis=1, keepdims=True) > 1)).ravel(),
    })
    return _output(a, [a.reference, a.queries], result)


def assay_main(a):
    import numpy as np
    import polars as pl
    from ..evalue.competitive import paired_assay_local_score

    ids = _read(a.queries, ("query_id",))["query_id"]
    if ids.n_unique() != len(ids):
        raise ValueError("query identifiers must be unique")
    distances = np.load(a.distances, allow_pickle=False)
    exact = np.load(a.exact_hits, allow_pickle=False) if a.exact_hits else None
    if not isinstance(distances, np.ndarray) or distances.shape != (len(ids), len(ids)):
        raise ValueError("distance axes must match the query identifiers")
    scores = paired_assay_local_score(distances, exact_hits=exact)
    inputs = [a.queries, a.distances]
    if a.exact_hits:
        inputs.append(a.exact_hits)
    return _output(a, inputs, pl.DataFrame({"query_id": ids, "score": scores}))


def _motifs_from_paired_reference(reference_path, queries, targets):
    """Calculate motif planes and distinct key counts in the supplied axis order."""
    import numpy as np
    import polars as pl
    from ..evalue.competitive import anchored_motif_scores

    chains = ("junction_aa_alpha", "junction_aa_beta")
    memberships = _read(reference_path, ("reference_key", "target", *chains)).unique()
    for frame in (memberships, queries):
        if not frame.select(pl.all_horizontal(
            pl.col(chain).str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$") for chain in chains
        ).all()).item():
            raise ValueError("paired junctions must contain canonical amino acids")
    if memberships.select("reference_key", *chains).unique().height != memberships["reference_key"].n_unique():
        raise ValueError("each reference key must identify one paired junction")
    if set(memberships["target"]) != set(targets["target"]):
        raise ValueError("reference memberships must cover exactly the target axis")
    if not targets["mhc_class"].is_in(["MHCI", "MHCII"]).all():
        raise ValueError("mhc_class must be MHCI or MHCII")

    counts_by_target = dict(memberships.group_by("target").len().iter_rows())
    counts = np.array([counts_by_target[target] for target in targets["target"]])
    planes = [np.empty((queries.height, targets.height)) for _ in chains]
    for mhc_class in targets["mhc_class"].unique(maintain_order=True):
        columns = np.flatnonzero((targets["mhc_class"] == mhc_class).to_numpy())
        labels = targets["target"][columns].to_list()
        pool = memberships.filter(pl.col("target").is_in(labels))
        for chain, plane in zip(chains, planes, strict=True):
            junctions, inverse = np.unique(queries[chain].to_numpy(), return_inverse=True)
            motif = anchored_motif_scores(pool.select("target", chain).iter_rows(), junctions, labels)
            plane[:, columns] = motif[inverse]
    return planes[0], planes[1], counts


def paired_main(a):
    import numpy as np
    import polars as pl
    from ..evalue.competitive import paired_motif_ranking_scores

    reference = getattr(a, "reference", None)
    query_table = _read(a.queries, ("query_id", "junction_aa_alpha", "junction_aa_beta")
                        if reference else ("query_id",))
    target_table = _read(a.targets, ("target", "mhc_class") if reference else ("target",))
    queries, targets = query_table["query_id"], target_table["target"]
    if queries.n_unique() != len(queries) or targets.n_unique() != len(targets):
        raise ValueError("query and target identifiers must each be unique")
    with np.load(a.scores, allow_pickle=False) as arrays:
        base = arrays["base"]
        motif_keys = ("tra", "trb", "target_key_counts")
        if reference:
            if any(key in arrays for key in motif_keys):
                raise ValueError("--reference requires a score archive without precomputed motifs or key counts")
        else:
            tra, trb, counts = (arrays[key] for key in motif_keys)
        lower = arrays["lower_tier_base"] if "lower_tier_base" in arrays else None
    if base.shape != (len(queries), len(targets)):
        raise ValueError("score axes must match the query and target identifiers")
    if reference:
        tra, trb, counts = _motifs_from_paired_reference(reference, query_table, target_table)
    scores = paired_motif_ranking_scores(base, tra, trb, counts,
                                         support_weight=a.support_weight,
                                         conflict_weight=a.conflict_weight,
                                         lower_tier_base=lower)
    table = pl.DataFrame({
        "query_id": np.repeat(queries.to_numpy(), len(targets)),
        "target": np.tile(targets.to_numpy(), len(queries)),
        "score": scores.ravel(),
    })
    return _output(a, [a.scores, a.queries, a.targets] + ([reference] if reference else []), table)
