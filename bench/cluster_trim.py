"""Clustering: full vs shortlist reference, untrimmed vs apex-trimmed CDR3 (germline ends removed).

Reuses the manuscript's exact clustering pipeline (single-linkage on the same-V neighbour graph: d1 edges
+ 2-sub PSSM bonus at the max-retention/baseline-purity operating point). The trim variant removes the
first 3 (germline V) and last 4 (germline J) residues from every CDR3 BEFORE the index/search, so two
same-V TCRs with an identical antigen-contacting apex but different germline ends become neighbours.
Reports purity / retention / #clusters for each (reference, chain, trim) cell.

    .venv/bin/python bench/cluster_trim.py [TRA|TRB|paired]
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path.home() / "vcs/manuscripts/2026-vdjmatch/benchmarks/scripts"))
import benchmark as B                                                # noqa: E402
import cluster_results as CR                                         # noqa: E402


def apex(c, ns=3, ne=4):
    """drop first ns germline-V and last ne germline-J residues; keep a central >=3-mer if too short."""
    k = c[ns:len(c) - ne]
    if len(k) >= 3:
        return k
    m = len(c) // 2
    return c[max(0, m - 2):m + 2]


def cluster_labels(cdr3, v, epi, trim):
    """vdjmatch cluster id per clonotype (-1 = singleton), apex-trimmed or not."""
    cd = [apex(c) for c in cdr3] if trim else cdr3
    n = len(cd)
    edges = CR._operating_edges(n, CR._d1_edges(cd, v), CR._pssm2_pairs(cd, v), epi)
    lab = [-1] * n
    for cid, m in enumerate(CR._components(n, edges).values()):
        if len(m) >= 2:
            for i in m:
                lab[i] = cid
    return lab


def cluster_pr(cdr3, v, epi, trim):
    cd = [apex(c) for c in cdr3] if trim else cdr3
    n = len(cd)
    edges = CR._operating_edges(n, CR._d1_edges(cd, v), CR._pssm2_pairs(cd, v), epi)
    pur, ret, nc = CR._purity_retention(CR._components(n, edges), epi, n)
    return pur, ret, nc, n


def cluster_pr_paired(ca, va, cb, vb, epi, trim):
    if trim:
        ca, cb = [apex(c) for c in ca], [apex(c) for c in cb]
    n = len(ca)
    costA, d1a = CR._chain_cost(ca, va, 5); costB, d1b = CR._chain_cost(cb, vb, 5)
    base = d1a & d1b
    bonus = [(i, j, costA[(i, j)] + costB[(i, j)])
             for (i, j) in (costA.keys() & costB.keys()) if (i, j) not in base]
    edges = CR._operating_edges(n, base, bonus, epi)
    pur, ret, nc = CR._purity_retention(CR._components(n, edges), epi, n)
    return pur, ret, nc, n


def main(which):
    d = B.release("vdjdb2026"); sl = B.shortlist(d)
    print(f"{'ref':10}{'chain':7}{'trim':7}{'n':>7}{'purity':>8}{'retention':>11}{'nclust':>8}")
    for ref_name, df in (("shortlist", sl), ("full", d)):
        if which in ("TRA", "TRB"):
            cdr3, epi, v = CR.single_clonotypes(df, which)
            for trim in (False, True):
                pur, ret, nc, n = cluster_pr(cdr3, v, epi, trim)
                print(f"{ref_name:10}{which:7}{'apex' if trim else 'none':7}{n:>7}{pur:>8.3f}{ret:>11.3f}{nc:>8}")
        else:
            ca, va, cb, vb, epi = CR.paired_clonotypes(d, df)   # sl=d -> all paired; sl=shortlist -> shortlist
            for trim in (False, True):
                pur, ret, nc, n = cluster_pr_paired(ca, va, cb, vb, epi, trim)
                print(f"{ref_name:10}{'paired':7}{'apex' if trim else 'none':7}{n:>7}{pur:>8.3f}{ret:>11.3f}{nc:>8}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TRB")
