"""Diagnose why a held-out epitope's binders do/don't separate from controls, from the cached search.

For epitope E: split queries into E-binders (true==E), the E dataset's own controls, and other-epitope
binders; report for each group the median ELA-similarity, the number of in-E neighbours found, the
first-hit radius, and the fraction passing the significance gate. Tells us if failure is (a) private
binders (few in-E neighbours), (b) cross-reactivity (other binders also match E), or (c) the gate.

    .venv/bin/python bench/holdout_study.py ELA TRB
"""
from __future__ import annotations

import bisect
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import holdout_eval as HE                                            # noqa: E402
from benchmark import PSSM_SCALE, SOFTV_BETA                         # noqa: E402
from vdjmatch.evalue import first_hit                                # noqa: E402
from vdjmatch.match import vgene as _vg                              # noqa: E402


def per_query(rec, e, M, Ntot, alpha=1e-3):
    cs = rec["ctrl"]
    nE = sum(1 for ps, ed, he, hv in rec["hits"] if he == e)            # in-E neighbours found
    sim = 0.0
    for ps, ed, he, hv in rec["hits"]:
        if he != e:
            continue
        w = 1.0 if _vg.gene_family(rec["v"]) == _vg.gene_family(hv) else SOFTV_BETA * _vg.vsim(rec["v"], hv)
        if w > 0:
            sim += w * math.exp(-ps / PSSM_SCALE) / (bisect.bisect_right(cs, ed) + 1)
    t1 = [(c, ep) for c, ep in rec["first"] if c <= 1]
    sig = first_hit.pvalue(t1, cs, Ntot, M)["p_enrichment"] < alpha
    return nE, sim, sig


def study(sh, locus):
    e = HE.EPI[sh]
    cache = HE.build(locus, "full")
    M, n_epi, recs = cache["M"], cache["n_epi"], cache["recs"]
    Ntot = sum(n_epi.values())
    NE = n_epi.get(e, 0)
    groups = {"E-binders": [], f"{sh}-controls": [], "other-binders": []}
    for r in recs:
        nE, sim, sig = per_query(r, e, M, Ntot)
        # per-epitope significance vs E (size-aware Poisson)
        t1E = [(c, ep) for c, ep in r["first"] if ep == e and c <= 1]
        sigE = first_hit.pvalue(t1E, r["ctrl"], NE if NE else 1, M)["p_enrichment"] < 1e-3
        row = (nE, sim, sig, sigE)
        if r["true"] == e:
            groups["E-binders"].append(row)
        elif r["true"] is None and r["dataset"].startswith(sh):
            groups[f"{sh}-controls"].append(row)
        elif r["true"] is not None:
            groups["other-binders"].append(row)
    print(f"\n=== {sh}-{locus}  (E reference size N_E={NE}, candidates={len(n_epi)}) ===")
    print(f"{'group':16}{'n':>6}{'med n_E-nbr':>12}{'%≥1 nbr':>9}{'med sim':>10}{'%sig(any)':>10}{'%sig(E)':>9}")
    for g, rows in groups.items():
        if not rows:
            continue
        a = np.array([r[0] for r in rows]); s = np.array([r[1] for r in rows])
        sg = np.array([r[2] for r in rows]); sge = np.array([r[3] for r in rows])
        print(f"{g:16}{len(rows):>6}{np.median(a):>12.1f}{100*(a>=1).mean():>8.0f}%"
              f"{np.median(s):>10.4f}{100*sg.mean():>9.0f}%{100*sge.mean():>8.0f}%")


if __name__ == "__main__":
    if len(sys.argv) == 3:
        study(sys.argv[1], sys.argv[2])
    else:
        for sh in ("ELA", "LLW", "LLL", "NLV"):
            study(sh, "TRB")
