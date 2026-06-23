#!/usr/bin/env python3
"""Shared hooks for the external TCR-clustering comparators (clusTCR, GIANA, iSMART, DeepTCR).

Every producer imports this so the clonotype sets and scoring are IDENTICAL to the in-house
manuscript method. Runs under the vdjmatch venv. Tools run in their own conda env `cmp-<tool>`;
this module only prepares inputs, parses tool outputs back into `labels`, and writes the canonical
prediction TSV.

The contract (from cluster_results.score_labels): labels[i] = cluster id of clonotype i, aligned
1:1 to the clonotype list; None or a negative id means "unclustered / singleton".
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent  # repo root
# manuscript scripts expose the exact sets + scorer
sys.path.insert(0, "/Users/mikesh/vcs/manuscripts/2026-vdjmatch/benchmarks/scripts")
sys.path.insert(0, str(ROOT / "bench"))
sys.path.insert(0, str(ROOT / "src"))

import cluster_results as CR  # noqa: E402
import benchmark as B  # noqa: E402

PRED = ROOT / "bench" / "predictions"


def sets():
    """Return the three shared clonotype sets, exactly as the in-house method clusters them.

    Returns a dict keyed by set name. Each value is a dict with at least cdr3/v/epi aligned lists;
    paired additionally carries alpha (ca/va) and beta (cb/vb).
    """
    d = B.release("vdjdb2026")
    sl = B.shortlist(d)
    out = {}
    for locus in ("TRB", "TRA"):
        cdr3, epi, v = CR.single_clonotypes(sl, locus)
        out[locus] = {"cdr3": cdr3, "v": v, "epi": epi}
    ca, va, cb, vb, epi_p = CR.paired_clonotypes(d, sl)
    out["paired"] = {"ca": ca, "va": va, "cb": cb, "vb": vb, "epi": epi_p,
                     # for tools that take a single CDR3 column, beta is the primary chain
                     "cdr3": cb, "v": vb}
    return out


def score(labels, epi):
    """(macro_purity, retention, n_clusters>=2) via the canonical scorer."""
    return CR.score_labels(labels, epi)


def labels_from_groups(cdr3, groups):
    """Build labels aligned to `cdr3` from a {cdr3 -> cluster_id} mapping.

    CDR3s absent from `groups` (i.e. the tool left them unclustered) get a unique negative id so the
    scorer counts them as singletons. Identical-CDR3 collisions in the input share the tool's label,
    which is the intended behaviour (the in-house set is unique-CDR3, so collisions are rare/none).
    """
    labels = []
    neg = -1
    for c in cdr3:
        g = groups.get(c)
        if g is None:
            labels.append(neg)
            neg -= 1
        else:
            labels.append(g)
    return labels


def write_tsv(tool, rows):
    """Write bench/predictions/<tool>/clustering.tsv. rows: list of dicts with keys
    set, macro_purity, retention, n_clusters, n_clonotypes, note."""
    out = PRED / tool
    out.mkdir(parents=True, exist_ok=True)
    cols = ["set", "macro_purity", "retention", "n_clusters", "n_clonotypes", "note"]
    lines = ["\t".join(cols)]
    for r in rows:
        lines.append("\t".join(str(r[c]) for c in cols))
    (out / "clustering.tsv").write_text("\n".join(lines) + "\n")
    print(f"[{tool}] wrote {out / 'clustering.tsv'}")
    for r in rows:
        print(f"  {r['set']:7} purity={r['macro_purity']} retention={r['retention']} "
              f"clusters={r['n_clusters']} n={r['n_clonotypes']} ({r['note']})")
