"""Search scope/budget parsing → seqtree SearchParams.

Legacy vdjmatch scope syntax ``s,i,d,t`` = max substitutions, insertions, deletions, total
edits. ``s,id,t`` (3 fields) treats the middle as a symmetric indel budget. ``t`` defaults to
the sum when omitted.
"""
from __future__ import annotations

from seqtree import SearchParams


def parse_scope(spec: str) -> tuple[int, int, int, int]:
    """Parse ``"s,i,d,t"`` / ``"s,id,t"`` / ``"s"`` → (subs, ins, dels, total)."""
    parts = [int(x) for x in str(spec).split(",")]
    if len(parts) == 1:
        s = parts[0]
        return s, 0, 0, s
    if len(parts) == 2:
        s, t = parts
        return s, t, t, t
    if len(parts) == 3:
        s, idl, t = parts
        return s, idl, idl, t
    s, i, d, t = parts[:4]
    return s, i, d, t


def search_params(scope: str = "1,0,0,1", *, engine: str = "seqtm", matrix="",
                  pos_matrix=None, max_penalty: int = 0, mode: str = "all",
                  gap_open: int = 1, gap_extend: int = 1) -> SearchParams:
    """Build a seqtree ``SearchParams`` from a scope spec and scoring options.

    ``pos_matrix`` is a settable attribute (not a constructor arg) in seqtree; result-count
    limiting (top-k) is done in Python, as seqtree has no ``max_hits``. When a scaled matrix is
    used with indel scopes, set ``gap_open``/``gap_extend`` to the matrix scale (~a typical
    substitution penalty) so gaps aren't absurdly cheap relative to substitutions.
    """
    s, i, d, t = parse_scope(scope)
    p = SearchParams(max_subs=s, max_ins=i, max_dels=d, max_total_edits=t, engine=engine,
                     matrix=matrix, max_penalty=max_penalty, mode=mode,
                     gap_open=gap_open, gap_extend=gap_extend)
    if pos_matrix is not None:
        p.pos_matrix = pos_matrix
    return p
