"""Reproducible ``match`` parameters as a JSON-serialisable dataclass.

``Params`` captures every knob the ``match`` subcommand exposes (search scope, scoring matrix,
gene-match constraints, E-value/alignment toggles, species, threads). It can be loaded from a JSON
config, merged with explicit CLI overrides, and written back next to a run so any annotation is
reproducible from its ``<prefix>.params.json``.
"""
from __future__ import annotations

import dataclasses
import json
from pathlib import Path


@dataclasses.dataclass
class Params:
    """Effective ``match`` parameters (one run's full configuration)."""

    scope: str = "1,0,0,1"          # subs,ins,dels,total-edits
    matrix: str = "vdjam"           # "vdjam" | "none"
    min_score: int = 0              # minimum VDJdb confidence score
    match_v: bool = False           # require V gene to match
    match_j: bool = False           # require J gene to match
    evalue: bool = True             # compute control-calibrated E-value
    align: bool = True              # emit per-hit CIGAR/alignment
    species: str = "HomoSapiens"    # species filter
    threads: int = 0                # worker threads (0 = all cores)

    @classmethod
    def defaults(cls) -> "Params":
        """Return a ``Params`` with all default values."""
        return cls()

    @classmethod
    def from_json(cls, path) -> "Params":
        """Load params from a JSON file, ignoring unknown keys."""
        data = json.loads(Path(path).read_text())
        fields = {f.name for f in dataclasses.fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in fields})

    def to_json(self, path) -> None:
        """Write params to ``path`` as pretty-printed JSON."""
        Path(path).write_text(json.dumps(dataclasses.asdict(self), indent=2) + "\n")
