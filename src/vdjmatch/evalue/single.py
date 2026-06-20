"""Single-chain control-calibrated E-values (thin polars wrapper over ``seqtree.evalues``).

For each query CDR3 we count neighbours in the target (VDJdb) and in a matched background control
within the search ball; ``E = (N/M)·n_control`` is the number of VDJdb neighbours expected by the
generative/background process, and ``p_enrichment`` is the Poisson-tail significance that the query
has *more* VDJdb neighbours than chance — the hallmark of antigen-driven convergence. Theory:
seqtree ``appendix/evalue.tex``.
"""
from __future__ import annotations

import polars as pl
from seqtree import Index, SearchParams
from seqtree.evalue import evalues


def query_evalues(target: Index, control: Index, queries: list[str], params: SearchParams,
                  threads: int = 0, exclude_exact: bool = False) -> pl.DataFrame:
    """Per-query E-values → polars frame with columns ``query_cdr3, n_target, n_control, E,
    p_any, p_enrichment, rule_of_three`` (one row per input query, input order preserved).

    ``exclude_exact=True`` punctures distance-0 (self/duplicate) hits on both sides — use when
    queries may be members of the target/control (e.g. VDJdb-vs-VDJdb verification)."""
    rows = evalues(target, control, queries, params, threads=threads, exclude_exact=exclude_exact)
    return pl.DataFrame({
        "query_cdr3": queries,
        "n_target": [r["n_target"] for r in rows],
        "n_control": [r["n_control"] for r in rows],
        "E": [r["E"] for r in rows],
        "p_any": [r["p_any"] for r in rows],
        "p_enrichment": [r["p_enrichment"] for r in rows],
        "rule_of_three": [r["rule_of_three"] for r in rows],
    })
