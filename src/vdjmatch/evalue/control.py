"""Background control repertoires for E-value calibration (thin wrapper over seqtree)."""

from __future__ import annotations

from seqtree import Index
from seqtree.control import load_control

# query locus -> seqtree control name (bundled human_trb_aa; others via HuggingFace)
_CONTROL = {
    "TRB": "human_trb_aa",
    "human:TRB": "human_trb_aa",
    "TRA": "human_tra_aa",
    "human:TRA": "human_tra_aa",
    "mouse:TRB": "mouse_trb_aa",
    "mouse:TRA": "mouse_tra_aa",
}


def _organism(species: str) -> str:
    aliases = {
        "human": "human",
        "homosapiens": "human",
        "homo sapiens": "human",
        "mouse": "mouse",
        "musmusculus": "mouse",
        "mus musculus": "mouse",
    }
    organism = aliases.get(str(species).strip().lower())
    if organism is None:
        raise ValueError(f"no control for species={species!r}")
    return organism


def background(
    locus: str = "TRB",
    species: str = "human",
    size: int | None = None,
    cache_dir: str | None = None,
) -> Index:
    """Load a deduplicated background repertoire ``Index`` for the given locus/species.

    Bundled: human TRB. Others (human TRA, mouse TRA/TRB) download via ``seqtree[control]``.
    """
    organism = _organism(species)
    name = _CONTROL.get(f"{organism}:{locus}")
    if name is None:
        raise ValueError(f"no control for locus={locus!r} species={species!r}")
    return load_control(name, size=size, cache_dir=cache_dir)


def raw_background(locus="TRB", species="human", path=None):
    """Build a fresh control index from a raw table or bundled human TRB junctions.

    Returns the index and raw-source provenance; never reads/writes a persisted index.
    Supplied tables explicitly declare junctions of the requested species and locus.
    """
    import gzip
    import hashlib
    import importlib.resources
    from pathlib import Path
    from ..db.cache import sha256

    import polars as pl
    from seqtree import Index

    species = _organism(species)
    if locus not in {"TRA", "TRB"}:
        raise ValueError("raw controls support TRA or TRB only")
    if path:
        path = Path(path)
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
            any(c in headers for c in ("v_call", "j_call", "sequence_id"))
        ):
            raise ValueError(
                "AIRR control requires junction_aa; bare CDR3 cannot be calibrated"
            )
        frame = pl.read_csv(
            path, separator=separator, columns=[column], infer_schema_length=0
        )
        frame = frame.with_columns(pl.col(column).str.strip_chars().str.to_uppercase())
        valid = frame.filter(pl.col(column).str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$"))
        sequences = valid[column].unique().sort().to_list()
        provenance = {
            "kind": "explicit_raw_table",
            "path": str(path.resolve()),
            "sha256": sha256(path),
            "source_rows": frame.height,
            "productive_rows": valid.height,
            "excluded_rows": frame.height - valid.height,
            "declared_species": species,
            "declared_locus": locus,
            "sequence_column": column,
        }
    else:
        if _organism(species) != "human" or locus != "TRB":
            raise ValueError(
                "this control is not bundled; supply --control with the matching raw species/locus repertoire"
            )
        asset = importlib.resources.files("seqtree").joinpath(
            "data/control_human_trb_aa.txt.gz"
        )
        payload = asset.read_bytes()
        raw_sequences = gzip.decompress(payload).decode().splitlines()
        sequences = pl.Series(raw_sequences).unique().sort().to_list()
        provenance = {
            "kind": "bundled_raw",
            "asset": "seqtree/data/control_human_trb_aa.txt.gz",
            "sha256": hashlib.sha256(payload).hexdigest(),
            "source_rows": len(raw_sequences),
        }
    if not sequences:
        raise ValueError("control has no valid junctions for the requested locus")
    provenance["unique_junctions"] = len(sequences)
    return Index.build(sequences, "aa"), provenance
