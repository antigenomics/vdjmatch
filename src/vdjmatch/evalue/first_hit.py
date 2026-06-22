"""First-hit (adaptive-scope) control-calibrated E-values.

Rather than fix a search ball, we widen to each query's **nearest** VDJdb hit (within ``max_edits``
total, ``max_ins``/``max_dels`` indels) and evaluate the E-value **at that first-hit radius**. Because
the background control's neighbour count grows with the radius, a distance-1 hit is significant while a
distance-5-only hit is not — the control calibrates the recall/specificity tradeoff per query, with no
fixed scope to pick. One wide search per side; the nested per-radius counts come from the per-hit edit
decomposition. ``E = (N/M)·n_control(r*)``; ``p_enrichment = Poisson P(X ≥ n_target(r*) | E)``.
"""
from __future__ import annotations

import polars as pl
from seqtree import Index, SearchParams
from seqtree.evalue import _poisson_sf


def scope(max_edits: int = 5, max_ins: int = 2, max_dels: int = 2, **kw) -> SearchParams:
    """First-hit search ball: up to ``max_edits`` total edits, at most ``max_ins`` ins and ``max_dels``
    del (default 5 / 2 / 2)."""
    return SearchParams(max_subs=max_edits, max_ins=max_ins, max_dels=max_dels,
                        max_total_edits=max_edits, engine="seqtm", **kw)


def _cost_lists(idx: Index, queries, params, threads, exclude_exact):
    out = []
    for hl in idx.search_batch(list(queries), params, threads):
        cs = [(h.n_subs + h.n_ins + h.n_dels, h.ref_id) for h in hl]
        if exclude_exact:
            cs = [(c, r) for c, r in cs if c > 0]
        out.append(sorted(cs, key=lambda x: x[0]))               # by total-edit cost, nearest first
    return out


def scan(target: Index, target_epi: list[str], control: Index, queries, *,
         params: SearchParams | None = None, threads: int = 0, exclude_exact: bool = False):
    """One wide search per side. Returns ``(target_hits, control_costs)`` where ``target_hits[q]`` is a
    cost-sorted list of ``(total_edits, epitope)`` and ``control_costs[q]`` a cost-sorted list of edits.
    ``target_epi`` maps each target ``ref_id`` to its epitope."""
    params = params or scope()
    th = _cost_lists(target, queries, params, threads, exclude_exact)
    ch = _cost_lists(control, queries, params, threads, exclude_exact)
    return ([[(c, target_epi[r]) for c, r in t] for t in th],
            [[c for c, _ in cc] for cc in ch])


def pvalue(target_hits, control_costs, N: int, M: int, epitope: str | None = None) -> dict:
    """First-hit E-value at the nearest (optionally ``epitope``-restricted) target hit:
    ``E = (N/M)·n_control(r*)``, ``p_enrichment = P(X ≥ n_target(r*) | E)``. ``N`` is the target size
    (use the epitope's size when ``epitope`` is set), ``M`` the control size."""
    th = target_hits if epitope is None else [(c, e) for c, e in target_hits if e == epitope]
    if not th:
        return {"radius": None, "n_target": 0, "n_control": 0, "E": 0.0, "p_enrichment": 1.0}
    r = th[0][0]
    n_t = sum(1 for c, _ in th if c <= r)
    n_c = sum(1 for c in control_costs if c <= r)
    E = (N / M) * n_c
    p = _poisson_sf(n_t, E) if E > 0 else (0.0 if n_t > 0 else 1.0)
    return {"radius": r, "n_target": n_t, "n_control": n_c, "E": E, "p_enrichment": p}


def query_first_hit(target: Index, target_epi: list[str], control: Index, queries, *,
                    N: int | None = None, M: int | None = None, params: SearchParams | None = None,
                    threads: int = 0, exclude_exact: bool = False) -> pl.DataFrame:
    """Per-query first-hit E-value + predicted epitope (the nearest target hit). Columns:
    ``query_cdr3, first_radius, n_target, n_control, E, p_enrichment, top_epitope``."""
    N = N if N is not None else len(target)
    M = M if M is not None else len(control)
    th, cc = scan(target, target_epi, control, queries, params=params, threads=threads,
                  exclude_exact=exclude_exact)
    rows = []
    for q, t, c in zip(queries, th, cc):
        r = pvalue(t, c, N, M)
        rows.append((q, r["radius"], r["n_target"], r["n_control"], r["E"], r["p_enrichment"],
                     t[0][1] if t else None))
    return pl.DataFrame(rows, orient="row", schema=["query_cdr3", "first_radius", "n_target",
                        "n_control", "E", "p_enrichment", "top_epitope"])
