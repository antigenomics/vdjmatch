"""Replay opt-in vX research scores on explicitly prepared inputs."""

from pathlib import Path


def register(subparsers):
    motif = subparsers.add_parser("anchored-motif-score", help="Score TRB junctions against a fixed motif catalogue")
    motif.add_argument("--reference", type=Path, required=True,
                       help="TSV with target, junction_aa, v_call and j_call")
    motif.add_argument("--queries", type=Path, required=True, help="TSV with query_id and junction_aa")
    motif.add_argument("--output", type=Path, required=True)
    motif.set_defaults(func=motif_main)

    assay = subparsers.add_parser("paired-assay-score", help="Rank one annotated paired-chain assay group")
    assay.add_argument("--distances", type=Path, required=True, help="Square NumPy .npy distance matrix")
    assay.add_argument("--queries", type=Path, required=True, help="TSV with query_id in matrix order")
    assay.add_argument("--exact-hits", type=Path, help="Optional aligned boolean .npy mask")
    assay.add_argument("--output", type=Path, required=True)
    assay.set_defaults(func=assay_main)


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
                                   query["junction_aa"], targets)
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
