"""TCR substitution scoring: load the bundled VDJAM matrix into a seqtree SubstitutionMatrix.

``vdjam.txt`` is a 20x20 amino-acid *similarity* table (log-odds; higher = more interchangeable
without changing specificity). seqtree's ``SubstitutionMatrix.from_similarity`` consumes an
integer similarity grid in ``amino_acids()`` order (24 symbols incl. B/Z/X/*) and derives the
matching penalty internally. We scale the float scores to integers and fill the non-standard
symbols neutrally (they are filtered out of real CDR3 input upstream).
"""
from __future__ import annotations

import os
from importlib import resources

from seqtree import SubstitutionMatrix, amino_acids

DEFAULT_SCALE = 100  # float similarity -> integer grid


def _vdjam_path() -> os.PathLike:
    return resources.files("vdjmatch.resources") / "vdjam.txt"


def load_vdjam(path: str | os.PathLike | None = None, scale: int = DEFAULT_SCALE) -> SubstitutionMatrix:
    """Build a seqtree ``SubstitutionMatrix`` from a VDJAM-format similarity table.

    Args:
        path: VDJAM ``aa.1 aa.2 score`` TSV (default: bundled ``vdjam.txt``).
        scale: integer scale factor applied to the float similarities.
    """
    src = path or _vdjam_path()
    sim: dict[tuple[str, str], float] = {}
    with open(src) as fh:
        next(fh)  # header: aa.1 aa.2 score
        for line in fh:
            a, b, s = line.split()
            sim[(a, b)] = float(s)
    diag = [v for (a, b), v in sim.items() if a == b]
    lo = round(scale * min(sim.values()))
    hi = round(scale * (sum(diag) / len(diag)))  # mean self-similarity for neutral fill

    order = amino_acids()
    n = len(order)
    # neutral default: dominant diagonal so derived substitution penalties stay non-negative
    grid = [[lo for _ in range(n)] for _ in range(n)]
    for i in range(n):
        grid[i][i] = hi
    for i, a in enumerate(order):
        for k, b in enumerate(order):
            if (a, b) in sim:
                grid[i][k] = round(scale * sim[(a, b)])
    return SubstitutionMatrix.from_similarity(grid)
