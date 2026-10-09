"""Compare scoring formulations for the per-epitope ranking (ROC), from the cached search.

Per epitope E the reference size N_E is constant across queries, so a size-AWARE enrichment p-value is a
valid (and statistically grounded) ranking score. We compare:
  sim        : size-invariant control-calibrated similarity (baseline argmax score)
  fh_E       : per-epitope first-hit Poisson enrichment, -log10 p at the query's nearest E hit
  r1 / r2    : per-epitope Poisson enrichment at a FIXED radius (edit<=1 / <=2) -- the close-neighbour test
  r1V / r2V  : the same but restricted to references sharing the query's V gene (V+CDR3 enrichment)
For diverse epitopes the first-hit radius is uninformatively large; a fixed close radius should separate.

    .venv/bin/python bench/holdout_variants.py [TRB|TRA]
"""
from __future__ import annotations

import bisect
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import holdout_eval as HE                                            # noqa: E402
from metrics import pr_auc_balanced, roc_auc                         # noqa: E402
from seqtree.evalue import _poisson_sf                               # noqa: E402
from vdjmatch.match import vgene as _vg                              # noqa: E402

EPIS = ["ELA", "LLW", "LLL", "NLV", "GLC", "YLQ"]


def _enrich(nE, nc, NE, M):
    E = (NE / M) * nc
    if nE == 0:
        return 0.0
    p = _poisson_sf(nE, E) if E > 0 else 0.0
    return -math.log10(max(p, 1e-300))


def score_query(rec, e, M, NE, kind):
    cs = rec["ctrl"]
    if kind == "sim":
        s = 0.0
        for ps, ed, he, hv in rec["hits"]:
            if he != e:
                continue
            w = 1.0 if _vg.gene_family(rec["v"]) == _vg.gene_family(hv) else 0.25 * _vg.vsim(rec["v"], hv)
            if w > 0:
                s += w * math.exp(-ps / 400.0) / (bisect.bisect_right(cs, ed) + 1)
        return s
    if kind == "fh_E":
        eh = [c for c, ep in rec["first"] if ep == e]
        if not eh:
            return 0.0
        r = min(eh)
        nE = sum(1 for c, ep in rec["first"] if ep == e and c <= r)
        return _enrich(nE, bisect.bisect_right(cs, r), NE, M)
    r = int(kind[1])
    vonly = kind.endswith("V")
    if vonly:
        nE = sum(1 for ps, ed, he, hv in rec["hits"]
                 if he == e and ed <= r and _vg.gene_family(hv) == _vg.gene_family(rec["v"]))
    else:
        nE = sum(1 for ps, ed, he, hv in rec["hits"] if he == e and ed <= r)
    return _enrich(nE, bisect.bisect_right(cs, r), NE, M)


def main(locus):
    cache = HE.build(locus, "full")
    M, n_epi, recs = cache["M"], cache["n_epi"], cache["recs"]
    kinds = ["sim", "fh_E", "r1", "r2", "r1V", "r2V"]
    present = [sh for sh in EPIS if any(r["true"] == HE.EPI[sh] for r in recs)]
    print(f"\n=== {locus}: per-epitope ROC by scoring variant (n={len(recs)}) ===")
    print(f"{'epitope':8}" + "".join(f"{k:>8}" for k in kinds))
    means = {k: [] for k in kinds}
    for sh in present:
        e = HE.EPI[sh]
        NE = n_epi.get(e, 1)
        y = [int(r["true"] == e) for r in recs]
        line = f"{sh:8}"
        for k in kinds:
            sc = [score_query(r, e, M, NE, k) for r in recs]
            roc = roc_auc(list(zip(y, sc)))
            means[k].append(roc)
            line += f"{roc:>8.3f}"
        print(line)
    print(f"{'MEAN':8}" + "".join(f"{sum(means[k])/len(means[k]):>8.3f}" for k in kinds))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TRB")
