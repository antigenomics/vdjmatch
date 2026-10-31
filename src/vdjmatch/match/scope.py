"""Search scope/budget parsing → seqtree SearchParams.

Legacy vdjmatch scope syntax ``s,i,d,t`` = max substitutions, insertions, deletions, total
edits. ``s,id,t`` (3 fields) treats the middle as a symmetric indel budget. ``t`` defaults to
the sum when omitted.
"""

from __future__ import annotations

from seqtree import SearchParams

from ..config import GRADED_FAMILIES


def parse_scope(spec: str) -> tuple[int, int, int, int]:
    """Parse ``"s,i,d,t"`` / ``"s,id,t"`` / ``"s"`` → (subs, ins, dels, total)."""
    parts = [int(x) for x in str(spec).split(",")]
    if len(parts) not in {1, 2, 3, 4} or any(p < 0 for p in parts):
        raise ValueError("scope needs 1 to 4 nonnegative integers")
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


def search_params(
    scope: str = "1,0,0,1",
    *,
    engine: str = "seqtm",
    matrix="",
    pos_matrix=None,
    max_penalty: int = 0,
    mode: str = "all",
    gap_open: int = 1,
    gap_extend: int = 1,
) -> SearchParams:
    """Build a seqtree ``SearchParams`` from a scope spec and scoring options.

    ``pos_matrix`` is a settable attribute (not a constructor arg) in seqtree; result-count
    limiting (top-k) is done in Python, as seqtree has no ``max_hits``. When a scaled matrix is
    used with indel scopes, set ``gap_open``/``gap_extend`` to the matrix scale (~a typical
    substitution penalty) so gaps aren't absurdly cheap relative to substitutions.
    """
    s, i, d, t = parse_scope(scope)
    p = SearchParams(
        max_subs=s,
        max_ins=i,
        max_dels=d,
        max_total_edits=t,
        engine=engine,
        matrix=matrix,
        max_penalty=max_penalty,
        mode=mode,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )
    if pos_matrix is not None:
        p.pos_matrix = pos_matrix
    return p


def graded_search_params(radius, matrix=None):
    """Existing edit-bounded family geometry shared by single and paired search.

    Weighted families change alignment costs, not the edit acceptance budget.
    Their NED score scale remains the candidate reducer's separate parameter.
    """
    from .scoring import DEFAULT_SCALE
    if isinstance(radius, bool) or not isinstance(radius, int) or not 1 <= radius <= 5:
        raise ValueError('radius must be an integer from 1 through 5')
    gap = DEFAULT_SCALE if matrix is not None else 1
    indels = min(radius, 2)
    return search_params(f'{radius},{indels},{indels},{radius}',
                         matrix=matrix if matrix is not None else '',
                         gap_open=gap, gap_extend=gap)


def graded_family_params(family, radius, radius_factor):
    """Shared M3 family validation, native geometry and dimensionless cost scale."""
    import math
    from .scoring import load_vdjam
    if family not in GRADED_FAMILIES:
        raise ValueError('unknown graded family')
    if radius_factor is None or isinstance(radius_factor, bool) or not math.isfinite(radius_factor) or radius_factor < 0:
        raise ValueError('graded family requires a finite nonnegative dimensionless radius factor')
    weighted = family.endswith('weighted')
    scale = 400. if weighted else 1.
    params = graded_search_params(radius, load_vdjam() if weighted else None)
    if not weighted and radius_factor < 1:
        # Every nonidentity unit edit costs at least one: use native exact search.
        params.max_subs = params.max_ins = params.max_dels = params.max_total_edits = 0
    budget = scale * radius_factor
    if 0 < budget < 2**31:
        # Round outward; the caller retains the exact normalized comparison.
        # Native zero means unbounded, so zero-cutoff searches stay post-filtered.
        params.max_penalty = math.ceil(budget)
    return params, scale


def blosum_radius_params(matcher, radius_factor, length):
    """Fixed BLOSUM recipes; radius is relative to the matrix substitution scale."""
    import math
    from seqtree import SubstitutionMatrix
    from .geometry import distance_units, _substitution_costs
    from .regions import significance_pssm
    if matcher not in {'blosum62','positional-blosum62'}:raise ValueError('unknown BLOSUM matcher')
    if radius_factor is None or isinstance(radius_factor,bool) or not math.isfinite(radius_factor) or radius_factor<=0:
        raise ValueError('explicit positive finite dimensionless radius factor required')
    scale,_=distance_units(_substitution_costs('full','blosum62'))
    budget=radius_factor*scale*(100 if matcher=='positional-blosum62' else 1)
    if not math.isfinite(budget) or budget>2**31-1:raise ValueError('radius factor exceeds native integer support')
    bound=math.floor(budget)
    # Native zero penalty disables the bound; zero edit scope is its exact equivalent.
    if not bound:return search_params('0,0,0,0')
    options=({'pos_matrix':significance_pssm(length)} if matcher=='positional-blosum62' else {'matrix':SubstitutionMatrix.blosum62()})
    return search_params('3,0,0,3',max_penalty=bound,**options)
