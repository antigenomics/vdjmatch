"""Annotation: does apex-trimming the CDR3 before the edit/PSSM search help per-epitope ROC?

Trim the first 3 (germline V) and last 4 (germline J) residues of every reference and query CDR3, build
the index over the trimmed strings, and score same-V neighbours by PSSM cost -> NED similarity. Compare
per-epitope ROC trimmed vs untrimmed on the held-out queries (full A*02 reference). Untrimmed uses the
central-weighting PSSM (significance_pssm); trimmed uses a flat matrix (the ends are already removed).

    .venv/bin/python bench/holdout_trim.py TRB
"""
from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

import polars as pl
from seqtree import Index, SearchParams

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
import holdout_eval as HE                                            # noqa: E402
from benchmark import A02, PSSM_SCALE, SOFTV_BETA, release, vgene    # noqa: E402
from cluster_trim import apex                                        # noqa: E402
from metrics import roc_auc                                          # noqa: E402
from vdjmatch.match import regions                                   # noqa: E402
from vdjmatch.match import vgene as _vg                              # noqa: E402


def search_ned(ref_cdr3, ref_epi, ref_v, q_cdr3, q_v, trim, max_edits=4):
    """index-keyed PSSM search; return per-query {epitope: NED-like sim} (same-V weighted)."""
    rc = [apex(c) for c in ref_cdr3] if trim else list(ref_cdr3)
    qc = [apex(c) for c in q_cdr3] if trim else list(q_cdr3)
    tgt = Index.build(rc, "aa")
    by_len = defaultdict(list)
    for i, q in enumerate(qc):
        by_len[len(q)].append(i)
    sims = [defaultdict(float) for _ in qc]
    for L, idxs in by_len.items():
        sp = SearchParams(max_subs=max_edits, max_ins=0, max_dels=0, max_total_edits=max_edits, engine="seqtm")
        if not trim:
            try:
                sp.pos_matrix = regions.significance_pssm(L)        # central-weighted (ends already low)
            except Exception:
                pass
        for i, hits in zip(idxs, tgt.search_batch([qc[k] for k in idxs], sp, 0)):
            for h in hits:
                ne = h.n_subs + h.n_ins + h.n_dels
                if ne == 0:
                    continue                                        # drop exact (leave-one-out)
                hv = ref_v[h.ref_id]
                w = 1.0 if _vg.gene_family(q_v[i]) == _vg.gene_family(hv) else SOFTV_BETA * _vg.vsim(q_v[i], hv)
                if w > 0:
                    sims[i][ref_epi[h.ref_id]] += w * math.exp(-h.score / PSSM_SCALE)
    return sims


def main(locus):
    rdf = release("vdjdb2026").filter(pl.col("mhc_a").str.contains(A02))
    r = (_bench.valid_cdr3(rdf.filter(pl.col("gene") == locus)).group_by("cdr3")
         .agg(pl.col("epitope").first(), pl.col("v").first()))
    ref_cdr3 = r["cdr3"].to_list(); ref_epi = r["epitope"].to_list()
    ref_v = [vgene(x) for x in r["v"]]
    pool = HE.query_set(locus)                                      # (cdr3, v, true, dataset)
    q_cdr3 = [p[0] for p in pool]; q_v = [p[1] for p in pool]; truth = [p[2] for p in pool]
    print(f"\n=== {locus}: annotation per-epitope ROC, untrimmed vs apex-trimmed search (n={len(pool)}) ===")
    res = {}
    for trim in (False, True):
        sims = search_ned(ref_cdr3, ref_epi, ref_v, q_cdr3, q_v, trim)
        res["apex" if trim else "none"] = sims
    print(f"{'epi':5}{'n+':>6}{'ROC_none':>10}{'ROC_apex':>10}{'delta':>8}")
    for sh in HE.EPI:
        e = HE.EPI[sh]
        if not any(t == e for t in truth):
            continue
        y = [int(t == e) for t in truth]
        rn = roc_auc(list(zip(y, [res["none"][i].get(e, 0.0) for i in range(len(pool))])))
        ra = roc_auc(list(zip(y, [res["apex"][i].get(e, 0.0) for i in range(len(pool))])))
        print(f"{sh:5}{sum(y):>6}{rn:>10.3f}{ra:>10.3f}{ra - rn:>+8.3f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TRB")
