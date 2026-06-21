"""V-gene similarity from germline antigen-contacting loops (CDR1 + CDR2).

A V-gene match is itself a strong specificity prior (same-V CDR3 neighbours share an epitope far more
often than cross-V ones; see ``bench/vgene_strat.py``). To extend that prior to germline-similar genes
(e.g. ``TRBV5-1`` ≈ ``TRBV5-8``) we compare V genes by their germline CDR1 and CDR2 loops — the two
principal antigen/MHC-contacting regions outside CDR3 — bundled in ``resources/vgene/human_v_cdr12.tsv``
(first allele per family, from ``mirpy``'s IMGT region annotation). The DE-loop / "CDR2.5" sits inside
FWR3 and is not isolated here (it needs IMGT numbering of FWR3); it is a future refinement.

Similarity is the ``difflib`` ratio over the concatenated CDR1+CDR2 string (stdlib, length-robust);
single-linkage clustering at a similarity cut yields fuzzy V groups.
"""
from __future__ import annotations

from difflib import SequenceMatcher
from functools import lru_cache
from importlib import resources

from .regions import gene_family


def load_v_cdr12(chain: str | None = None, path=None) -> dict[str, str]:
    """Map V gene family -> germline ``CDR1+CDR2`` string (optionally filtered to one ``chain``)."""
    src = path or (resources.files("vdjmatch.resources") / "vgene" / "human_v_cdr12.tsv")
    out: dict[str, str] = {}
    with open(src) as fh:
        next(fh)
        for line in fh:
            ch, v, c1, c2 = line.rstrip("\n").split("\t")
            if chain is None or ch == chain:
                out[v] = c1 + c2
    return out


_LOOPS: dict[str, str] = {}


@lru_cache(maxsize=None)
def vsim(v1: str | None, v2: str | None) -> float:
    """Germline CDR1+CDR2 similarity in ``[0,1]`` between two V genes (1.0 = same family or identical
    loops; 0.0 if either loop is unknown). Allele/decoration is stripped first."""
    if not _LOOPS:
        _LOOPS.update(load_v_cdr12())
    f1, f2 = gene_family(v1), gene_family(v2)
    if f1 == f2:
        return 1.0
    l1, l2 = _LOOPS.get(f1), _LOOPS.get(f2)
    if not l1 or not l2:
        return 0.0
    return SequenceMatcher(None, l1, l2).ratio()


def v_clusters(chain: str, cut: float = 0.8) -> dict[str, int]:
    """Single-linkage cluster V families at a CDR1+CDR2 similarity ``cut``; returns ``v -> cluster_id``.
    Genes whose contacting loops are >= ``cut`` similar are merged (fuzzy V matching)."""
    loops = load_v_cdr12(chain)
    vs = sorted(loops)
    parent = {v: v for v in vs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, a in enumerate(vs):
        for b in vs[i + 1:]:
            if SequenceMatcher(None, loops[a], loops[b]).ratio() >= cut:
                parent[find(a)] = find(b)
    roots = {}
    return {v: roots.setdefault(find(v), len(roots)) for v in vs}


def _demo():
    assert vsim("TRBV5-1", "TRBV5-1*02") == 1.0          # same family / allele -> identical
    assert vsim("TRBV5-1", "NOPE") == 0.0                # unknown gene -> 0
    # contacting-loop similarity tracks germline relatedness but NOT family number: TRBV5-1 and
    # TRBV5-8 share CDR1 yet differ in CDR2, so they are only moderately similar.
    assert vsim("TRBV5-1", "TRBV5-8") > vsim("TRBV5-1", "TRBV19")
    cl = v_clusters("TRB", cut=0.8)
    assert len(set(cl.values())) < len(cl)               # some genes merge into fuzzy V groups
    merged = next((a, b) for a in cl for b in cl if a < b and cl[a] == cl[b])
    print("vsim(5-1,5-8)=%.3f  vsim(5-1,19)=%.3f  TRB groups=%d/%d (e.g. %s~%s)"
          % (vsim("TRBV5-1", "TRBV5-8"), vsim("TRBV5-1", "TRBV19"),
             len(set(cl.values())), len(cl), merged[0], merged[1]))


if __name__ == "__main__":
    _demo()
