#!/usr/bin/env python3
"""Assemble the Fig 4A PER-EPITOPE clustering purity/retention TSV across five methods x three chains.

Runs under the vdjmatch venv. For each (method, chain) it obtains a label vector aligned 1:1 to the
shared clonotype set, then scores it with CR.per_epitope_scores (the canonical per-epitope metric), and
writes one row per (method, chain, epitope) to bench/predictions/_per_epitope/cluster_per_epitope.tsv.

Methods x chains (only where the method supports the chain; GLIPH2/GIANA/iSMART are beta-only):
  vdjmatch  : TRA, TRB, paired   (CR.vdjmatch_labels)
  tcrdist3  : TRA, TRB, paired   (single-linkage on TCRdist; bench/tcrdist_cluster_perepi.py)
  GLIPH2    : TRB                 (irtools/GLIPH2 CD48 background; bench/gliph2_cluster_perepi.py)
  GIANA     : TRB                 (cached GIANA output, parsed to labels)
  iSMART    : TRB                 (cached iSMART output, parsed to labels)

Best-effort: a (method, chain) that cannot produce labels is skipped (logged), the rest still written.

    ./.venv/bin/python bench/cluster_per_epitope_assemble.py
"""
from __future__ import annotations

import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parent
sys.path.insert(0, "/Users/mikesh/vcs/manuscripts/2026-vdjmatch/benchmarks/scripts")
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(BENCH.parent / "src"))
import _cluster_common as C  # noqa: E402
import benchmark as B  # noqa: E402
import cluster_results as CR  # noqa: E402

OUT = BENCH / "predictions" / "_per_epitope" / "cluster_per_epitope.tsv"
GIANA_OUT = BENCH / "out" / "_giana" / "trb_in--RotationEncodingBL62.txt"
ISMART_OUT = BENCH / "out" / "_ismart" / "trb_in_clustered_v3.txt"
PAIRED_PKL = BENCH / "out" / "_tcrdist_cluster" / "paired_labels.pkl"


def _parse_giana(path, cdr3):
    groups = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        if len(p) < 2:
            continue
        groups[p[0]] = int(p[1])
    return C.labels_from_groups(cdr3, groups)


def _parse_ismart(path, cdr3):
    groups = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        if not (p[0].startswith("C") and len(p) >= 2):
            continue
        try:
            cid = int(p[-1])
        except ValueError:
            continue
        groups[p[0]] = cid
    return C.labels_from_groups(cdr3, groups)


def collect():
    """Return {(method, chain): (labels, epi)} for every (method, chain) that succeeds."""
    d = B.release("vdjdb2026")
    sl = B.shortlist(d)
    s = C.sets()
    res = {}

    # vdjmatch — all three sets
    for chain in ("TRA", "TRB", "paired"):
        try:
            labels, epi = CR.vdjmatch_labels(chain, sl, d if chain == "paired" else None)
            res[("vdjmatch", chain)] = (labels, epi)
        except Exception as ex:  # noqa: BLE001
            print(f"[vdjmatch {chain}] SKIPPED: {ex}", file=sys.stderr)

    # tcrdist3 — TRA/TRB recomputed (fast), paired loaded from the background pickle if present
    try:
        import tcrdist_cluster_perepi as T  # noqa: PLC0415
        for chain in ("TRA", "TRB"):
            try:
                labels, epi = T.single(sl, chain)
                res[("tcrdist3", chain)] = (labels, epi)
            except Exception as ex:  # noqa: BLE001
                print(f"[tcrdist3 {chain}] SKIPPED: {ex}", file=sys.stderr)
    except Exception as ex:  # noqa: BLE001
        print(f"[tcrdist3 single] import/setup SKIPPED: {ex}", file=sys.stderr)
    if PAIRED_PKL.exists():
        import pickle  # noqa: PLC0415
        labels, epi = pickle.load(PAIRED_PKL.open("rb"))
        res[("tcrdist3", "paired")] = (labels, epi)
    else:
        print("[tcrdist3 paired] SKIPPED: paired_labels.pkl not present "
              "(background compute unfinished/over budget)", file=sys.stderr)

    # GLIPH2 — TRB only
    try:
        import gliph2_cluster_perepi as Gl  # noqa: PLC0415
        labels, epi = Gl.produce()["TRB"]
        res[("GLIPH2", "TRB")] = (labels, epi)
    except Exception as ex:  # noqa: BLE001
        print(f"[GLIPH2 TRB] SKIPPED: {ex}", file=sys.stderr)

    # GIANA / iSMART — TRB only, from cached tool outputs (keyed by CDR3 -> position-independent)
    trb = s["TRB"]
    if GIANA_OUT.exists():
        res[("GIANA", "TRB")] = (_parse_giana(GIANA_OUT, trb["cdr3"]), trb["epi"])
    else:
        print("[GIANA TRB] SKIPPED: cached output missing", file=sys.stderr)
    if ISMART_OUT.exists():
        res[("iSMART", "TRB")] = (_parse_ismart(ISMART_OUT, trb["cdr3"]), trb["epi"])
    else:
        print("[iSMART TRB] SKIPPED: cached output missing", file=sys.stderr)
    return res


def main():
    res = collect()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    rows = ["method\tchain\tepitope\tpurity\tretention\tn_clonotypes"]
    order = ["vdjmatch", "tcrdist3", "GLIPH2", "GIANA", "iSMART"]
    chain_order = {"TRA": 0, "TRB": 1, "paired": 2}
    print("\n=== aggregate sanity (per (method, chain)) ===")
    for (method, chain) in sorted(res, key=lambda k: (order.index(k[0]), chain_order[k[1]])):
        labels, epi = res[(method, chain)]
        agg = CR.score_labels(labels, epi)
        print(f"  {method:9} {chain:7}: purity={agg[0]:.3f} retention={agg[1]:.3f} clusters={agg[2]}")
        pe = CR.per_epitope_scores(labels, epi)
        for ep in sorted(pe):
            pur, ret, n = pe[ep]
            rows.append(f"{method}\t{chain}\t{ep}\t{pur}\t{ret}\t{n}")
    OUT.write_text("\n".join(rows) + "\n")
    print(f"\nwrote {OUT} ({len(rows) - 1} rows)")


if __name__ == "__main__":
    main()
