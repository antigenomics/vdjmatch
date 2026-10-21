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


def load_cost_table(path: str | os.PathLike) -> "list[list[int]]":
    """Read an explicit integer substitution COST table into an ``amino_acids()``-ordered grid.

    Format: ``aa.1 aa.2 cost`` with a header, one row per unordered pair; costs are non-negative
    integers, zero on the diagonal, and the table is symmetrised on read. Unlike the VDJAM
    similarity table this is consumed directly, without the Gram transform -- that transform is
    lossy (seqtree keeps ``similarity`` and ``penalty`` separately for exactly this reason), so a
    matrix derived AS costs would not survive a round trip through similarities.

    Non-standard symbols (``B Z X *``) are NOT inferred from the table. A 190-entry table carries
    no information about them, and they are not inert: the V-loop templates are padded with
    characters outside the standard twenty, so whatever those rows contain silently rewrites the
    V-loop distances that feed the kernel. They are therefore carried over from BLOSUM62's own Gram
    penalties, rescaled to this table's range, so that every candidate matrix handles padding
    identically and a comparison between two tables isolates the 20x20 that is actually being
    varied. A cost table written out from BLOSUM62 reads back as BLOSUM62 exactly.
    """
    cost: dict[tuple[str, str], int] = {}
    with open(path) as fh:
        header = next(fh).split()
        if len(header) < 3:
            raise ValueError("cost table needs an 'aa.1 aa.2 cost' header")
        for line in fh:
            if not line.strip():
                continue
            a, b, value = line.split()
            number = float(value)
            if number < 0 or number != int(number):
                raise ValueError(f"substitution costs must be non-negative integers, got {value}")
            cost[(a, b)] = cost[(b, a)] = int(number)
    if not cost:
        raise ValueError("cost table is empty")
    if any(cost.get((a, a), 0) for a, _ in cost):
        raise ValueError("a substitution cost table must be zero on the diagonal")

    standard = "ACDEFGHIKLMNPQRSTVWY"
    missing = [(a, b) for i, a in enumerate(standard) for b in standard[i + 1:]
               if (a, b) not in cost]
    if missing:
        raise ValueError(f"cost table is missing {len(missing)} standard pairs, e.g. {missing[:3]}")

    from seqtree import SubstitutionMatrix
    reference = SubstitutionMatrix.blosum62()
    reference_top = max(reference.penalty(a, b) for a in standard for b in standard)
    scale = (max(cost.values()) / reference_top) if reference_top else 1.0

    order = amino_acids()
    grid = [[0 for _ in order] for _ in order]
    for i, a in enumerate(order):
        for k, b in enumerate(order):
            if a == b:
                continue
            if a in standard and b in standard:
                grid[i][k] = cost[(a, b)]
            else:
                grid[i][k] = int(round(reference.penalty(a, b) * scale))
    return grid
