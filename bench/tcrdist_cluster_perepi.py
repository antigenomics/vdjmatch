#!/usr/bin/env python3
"""tcrdist3 single-linkage CLUSTERING -> per-epitope purity/retention on the shared clonotype sets.

Runs under the vdjmatch venv. For each chain (TRA, TRB, paired) it:
  1. rebuilds the EXACT shared clonotype set (CR.single_clonotypes / CR.paired_clonotypes) and recovers
     the J gene per clonotype (cdr3->j map; paired (ca,cb)->(ja,jb) map) without re-ordering the set;
  2. ships a (cdr3,v,j) TSV per chain to bench/_tcrdist_cluster_compute.py in the `cmp-tcrdist` env,
     which returns the within-set sparse TCRdist edge list (i,j,dist) with i,j = shared-set row indices;
  3. forms single-linkage connected components: single-chain at R=12, paired requires BOTH chains with
     dist_a + dist_b <= 48 (additive paired scheme; each chain computed within radius 48 so any summing
     pair is captured);
  4. builds labels aligned 1:1 to the shared set (>=2-member component id, else -1 singleton), and
     reports CR.per_epitope_scores + the aggregate (sanity vs anchors TRB 0.989/0.452, TRA 0.975/0.447,
     paired 0.989/0.660).

Emits rows (method=tcrdist3) for the combined per-epitope TSV. Best-effort: a chain that errors or
overruns is skipped (reported), the rest still produced.

    ./.venv/bin/python bench/tcrdist_cluster_perepi.py
"""
from __future__ import annotations

import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import polars as pl

BENCH = Path(__file__).resolve().parent
sys.path.insert(0, "/Users/mikesh/vcs/manuscripts/2026-vdjmatch/benchmarks/scripts")
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(BENCH.parent / "src"))
import _bench  # noqa: E402
import benchmark as B  # noqa: E402
import cluster_results as CR  # noqa: E402

ENV = "cmp-tcrdist"
COMPUTE = BENCH / "_tcrdist_cluster_compute.py"
COMPUTE_PAIRED = BENCH / "_tcrdist_cluster_paired.py"
WORK = BENCH / "out" / "_tcrdist_cluster"
R_SINGLE = 12
R_PAIRED = 48


def _components(n, edges):
    uf = list(range(n))

    def find(x):
        while uf[x] != x:
            uf[x] = uf[uf[x]]
            x = uf[x]
        return x

    for i, j in edges:
        a, b = find(i), find(j)
        if a != b:
            uf[a] = b
    comp = defaultdict(list)
    for i in range(n):
        comp[find(i)].append(i)
    return comp


def _labels(n, edges):
    comp = _components(n, edges)
    labels = [-1] * n
    for cid, members in enumerate(comp.values()):
        if len(members) >= 2:
            for i in members:
                labels[i] = cid
    return labels


def _run_compute(cdr3, v, j, chain, radius, tag):
    """Write input TSV, run the conda worker, return edge list [(i, j, dist), ...]."""
    WORK.mkdir(parents=True, exist_ok=True)
    inp = WORK / f"{tag}_in.tsv"
    out = WORK / f"{tag}_edges.tsv"
    pl.DataFrame({"cdr3": cdr3, "v": [vv + "*01" for vv in v], "j": j}).write_csv(inp, separator="\t")
    subprocess.run(["conda", "run", "-n", ENV, "python", str(COMPUTE),
                    "--in", str(inp.resolve()), "--chain", chain, "--radius", str(radius),
                    "--out", str(out.resolve())], check=True)
    e = pl.read_csv(out, separator="\t")
    return list(zip(e["i"].to_list(), e["j"].to_list(), e["dist"].to_list()))


def _jmap_single(sl, locus):
    m = (_bench.valid_cdr3(sl.filter(pl.col("gene") == locus)).group_by("cdr3")
         .agg(pl.col("j").first()))
    return dict(zip(m["cdr3"], m["j"]))


def single(sl, locus):
    cdr3, epi, v = CR.single_clonotypes(sl, locus)
    chain = "alpha" if locus == "TRA" else "beta"
    jmap = _jmap_single(sl, locus)
    j = [jmap[c] for c in cdr3]
    edges3 = _run_compute(cdr3, v, j, chain, R_SINGLE, locus)
    edges = [(i, jj) for (i, jj, dd) in edges3 if dd <= R_SINGLE]
    labels = _labels(len(cdr3), edges)
    return labels, epi


def paired(d, sl):
    ca, va, cb, vb, epi = CR.paired_clonotypes(d, sl)
    # recover ja/jb aligned to the (ca,cb) order
    pr = d.filter(pl.col("complex_id") != 0)
    a = _bench.valid_cdr3(pr.filter(pl.col("gene") == "TRA")).select("complex_id", ca="cdr3", ja="j")
    b = _bench.valid_cdr3(pr.filter(pl.col("gene") == "TRB")).select("complex_id", cb="cdr3", jb="j")
    p = a.join(b, on="complex_id").unique(["ca", "cb"])
    jm = {(r[0], r[1]): (r[2], r[3]) for r in p.select("ca", "cb", "ja", "jb").iter_rows()}
    ja = [jm[(x, y)][0] for x, y in zip(ca, cb)]
    jb = [jm[(x, y)][1] for x, y in zip(ca, cb)]
    # NATIVE tcrdist3 paired distance on the UNIQUE (alpha,beta) clonotypes: the full paired key is unique
    # (9214 distinct pairs), so each clonotype is ONE node — we must NOT dedup per-chain and expand (a
    # single shared alpha is in up to ~3000 different paired clonotypes; expanding would weld them into a
    # cross-epitope blob). Edge = joint tcrdist (alpha+beta) <= R_PAIRED.
    WORK.mkdir(parents=True, exist_ok=True)
    inp = WORK / "paired_in.tsv"
    out = WORK / "paired_edges.tsv"
    pl.DataFrame({"cdr3_a": ca, "v_a": [x + "*01" for x in va], "j_a": ja,
                  "cdr3_b": cb, "v_b": [x + "*01" for x in vb], "j_b": jb}).write_csv(inp, separator="\t")
    subprocess.run(["conda", "run", "-n", ENV, "python", str(COMPUTE_PAIRED),
                    "--in", str(inp.resolve()), "--radius", str(R_PAIRED),
                    "--out", str(out.resolve())], check=True)
    e = pl.read_csv(out, separator="\t")
    edges = [(i, j) for (i, j, dd) in zip(e["i"].to_list(), e["j"].to_list(), e["dist"].to_list())
             if dd <= R_PAIRED]
    labels = _labels(len(ca), edges)
    return labels, epi


def produce():
    """Return {chain: (labels, epi)} for the chains that succeed (best-effort)."""
    d = B.release("vdjdb2026")
    sl = B.shortlist(d)
    out = {}
    for locus in ("TRA", "TRB"):
        try:
            out[locus] = single(sl, locus)
        except Exception as ex:  # noqa: BLE001
            print(f"[tcrdist3 {locus}] SKIPPED: {ex}", file=sys.stderr)
    try:
        la, ep = paired(d, sl)
        out["paired"] = (la, ep)
    except Exception as ex:  # noqa: BLE001
        print(f"[tcrdist3 paired] SKIPPED: {ex}", file=sys.stderr)
    return out


if __name__ == "__main__":
    res = produce()
    for chain, (labels, epi) in res.items():
        agg = CR.score_labels(labels, epi)
        print(f"tcrdist3 {chain:7}: aggregate purity={agg[0]} retention={agg[1]} clusters={agg[2]} "
              f"(n={len(labels)})")
