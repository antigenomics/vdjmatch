"""Background control repertoires for E-value calibration (thin wrapper over seqtree)."""
from __future__ import annotations

from seqtree import Index
from seqtree.control import load_control

# query locus -> seqtree control name (bundled human_trb_aa; others via HuggingFace)
_CONTROL = {
    "TRB": "human_trb_aa", "human:TRB": "human_trb_aa",
    "TRA": "human_tra_aa", "human:TRA": "human_tra_aa",
    "mouse:TRB": "mouse_trb_aa", "mouse:TRA": "mouse_tra_aa",
}


def background(locus: str = "TRB", species: str = "human", size: int | None = None,
               cache_dir: str | None = None) -> Index:
    """Load a deduplicated background repertoire ``Index`` for the given locus/species.

    Bundled: human TRB. Others (human TRA, mouse TRA/TRB) download via ``seqtree[control]``.
    """
    key = f"{species.lower()}:{locus}"
    name = _CONTROL.get(key) or _CONTROL.get(locus)
    if name is None:
        raise ValueError(f"no control for locus={locus!r} species={species!r}")
    return load_control(name, size=size, cache_dir=cache_dir)
