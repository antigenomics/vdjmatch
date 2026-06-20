"""Convert seqtree alignment ops into CIGAR + a human-readable match line.

``seqtree.Alignment.ops`` is one char per alignment column: ``M`` match, ``S`` substitution,
``I`` insertion (in query), ``D`` deletion (from query). We emit an extended CIGAR using ``=``
for matches and ``X`` for substitutions (``I``/``D`` unchanged), and a BLAST-style midline.
"""
from __future__ import annotations

_CIGAR = {"M": "=", "S": "X", "I": "I", "D": "D"}


def to_cigar(ops: str) -> str:
    """Run-length-encode ops into extended CIGAR, e.g. ``"5=1X3="``."""
    if not ops:
        return ""
    out, run, prev = [], 0, ops[0]
    for c in ops:
        if c == prev:
            run += 1
        else:
            out.append(f"{run}{_CIGAR[prev]}")
            run, prev = 1, c
    out.append(f"{run}{_CIGAR[prev]}")
    return "".join(out)


def match_line(ops: str) -> str:
    """BLAST-style midline: ``|`` under a match, space under a substitution/gap."""
    return "".join("|" if c == "M" else " " for c in ops)


def counts(ops: str) -> dict[str, int]:
    """Per-op counts: matches, subs, ins, dels."""
    return {"matches": ops.count("M"), "subs": ops.count("S"),
            "ins": ops.count("I"), "dels": ops.count("D")}
