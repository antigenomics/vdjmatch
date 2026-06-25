"""Paired alpha/beta scoring on the fully-paired held-out epitopes (GLC/YLQ from TCRvdb).

Only GLC/YLQ have cell-level paired alpha+beta with paired positives AND negatives (sample6 TCRvdb;
the sewell chunk has zero pairing -> ELA/LLW/LLL stay single-chain). Per paired clonotype we search each
chain against its locus reference (full vdjdb2026 A*02, exact removed), score per-chain NED + V-gene
log-odds, and fuse the four channels by calibrated min-p (non-diluting). We compare paired vs alpha-only
vs beta-only per-epitope ROC -- does pairing lift the already-strong convergent epitopes?

    .venv/bin/python bench/holdout_paired.py
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
import holdout_eval as HE                                            # noqa: E402
import holdout_features as HF                                        # noqa: E402
from benchmark import A02, _pssm_targets, ref_index, release, vgene  # noqa: E402
from compare import TESTDATA                                         # noqa: E402
from metrics import roc_auc                                          # noqa: E402
from vdjmatch.evalue import background, first_hit                    # noqa: E402

PAIR_EPI = {"GLC": "GLCTLVAML", "YLQ": "YLQPRTFLL"}
CACHE = HE.CACHE / "paired.pkl"


def paired_clonotypes():
    """paired GLC/YLQ clonotypes from sample6: (a_cdr3,a_v, b_cdr3,b_v, true_epi_or_None)."""
    t = pl.read_csv(TESTDATA / "sample6_TCRvdb.csv")
    out = []
    for sh, e in PAIR_EPI.items():
        d = t.filter(pl.col("epitope_aa") == e)
        for ac, av, bc, bv, padj in zip(d["cdr3_alpha_aa"], d["TRAV"], d["cdr3_beta_aa"], d["TRBV"], d["padj"]):
            if not ac or not bc or ac == "" or bc == "":
                continue
            if padj is None or \
               not _bench.valid_cdr3(pl.DataFrame({"cdr3": [ac]})).height or \
               not _bench.valid_cdr3(pl.DataFrame({"cdr3": [bc]})).height:
                continue
            out.append((ac, vgene(av), bc, vgene(bv), e if padj < 1e-5 else None))
    # dedup by (alpha,beta)
    seen, uniq = set(), []
    for r in out:
        if (r[0], r[2]) not in seen:
            seen.add((r[0], r[2])); uniq.append(r)
    return uniq


def _search(locus, queries):
    rdf = release("vdjdb2026").filter(pl.col("mhc_a").str.contains(A02))
    tgt, ref_epi, ref_v, n_epi, _, _ = ref_index(rdf, locus)
    ctrl = background(locus)
    params = first_hit.scope(5, 2, 2)
    th, cc = first_hit.scan(tgt, ref_epi, ctrl, queries, target_v=ref_v, params=params, exclude_exact=True)
    pt = _pssm_targets(tgt, ref_epi, ref_v, queries, 5)
    recs = {}
    for c, t, ctr in zip(queries, th, cc):
        recs[c] = {"hits": [(ps, ed, he, hv) for ps, ed, he, hv in pt[c]],
                   "first": [(cost, e) for cost, e, *_ in t], "ctrl": sorted(ctr)}
    return recs, len(ctrl), dict(n_epi)


def build():
    if CACHE.exists():
        return pickle.loads(CACHE.read_bytes())
    cl0 = paired_clonotypes()
    a_q = sorted({r[0] for r in cl0}); b_q = sorted({r[2] for r in cl0})
    a_recs, Ma, _ = _search("TRA", a_q)
    b_recs, Mb, _ = _search("TRB", b_q)
    cache = {"clones": cl0, "a_recs": a_recs, "b_recs": b_recs, "Ma": Ma, "Mb": Mb}
    CACHE.write_bytes(pickle.dumps(cache))
    print(f"[paired] {len(cl0)} clonotypes ({len(a_q)} alpha, {len(b_q)} beta)", file=sys.stderr)
    return cache


def main():
    cache = build()
    clones = cache["clones"]
    _, models_a = HF.build_models("TRA", "airr")
    _, models_b = HF.build_models("TRB", "airr")

    def chan_scores(rec, models, M, e):
        m = models.get(e)
        if m is None:
            return 0.0, 0.0
        ned = HF.ned_sim({"v": rec["v"], "hits": rec["hits"], "ctrl": rec["ctrl"]}).get(e, 0.0)
        v = m["V"].get(rec["v"], 0.0)
        return ned, v

    # assemble per-clone per-epitope 4 channel scores
    A = [dict(v=r[1], **cache["a_recs"][r[0]]) for r in clones]
    B = [dict(v=r[3], **cache["b_recs"][r[2]]) for r in clones]
    y_true = [r[4] for r in clones]
    print(f"\n=== paired GLC/YLQ: per-epitope ROC ===")
    print(f"{'epi':5}{'n+':>5}{'alpha':>8}{'beta':>8}{'pair_max':>9}{'pair_fish':>10}")
    for sh, e in PAIR_EPI.items():
        y = [int(t == e) for t in y_true]
        ned_a = np.array([chan_scores(A[i], models_a, cache["Ma"], e)[0] for i in range(len(clones))])
        v_a = np.array([chan_scores(A[i], models_a, cache["Ma"], e)[1] for i in range(len(clones))])
        ned_b = np.array([chan_scores(B[i], models_b, cache["Mb"], e)[0] for i in range(len(clones))])
        v_b = np.array([chan_scores(B[i], models_b, cache["Mb"], e)[1] for i in range(len(clones))])
        null = [i for i, t in enumerate(y_true) if t is None]

        def logps(chans):                                          # per-channel -log10 empirical p
            cols = []
            for arr in chans:
                n0 = sorted(arr[null])
                cols.append(np.array([-np.log10(max(HF._emp_p(arr[i], n0), 1e-12)) for i in range(len(clones))]))
            return np.vstack(cols)

        cmax = lambda chans: np.max(logps(chans), axis=0)           # min-p (non-diluting)
        cfish = lambda chans: np.sum(logps(chans), axis=0)          # Fisher (additive evidence)
        roc_a = roc_auc(list(zip(y, cmax([ned_a, v_a]).tolist())))
        roc_b = roc_auc(list(zip(y, cmax([ned_b, v_b]).tolist())))
        roc_pm = roc_auc(list(zip(y, cmax([ned_a, v_a, ned_b, v_b]).tolist())))
        roc_pf = roc_auc(list(zip(y, cfish([ned_a, v_a, ned_b, v_b]).tolist())))
        print(f"{sh:5}{sum(y):>5}{roc_a:>8.3f}{roc_b:>8.3f}{roc_pm:>9.3f}{roc_pf:>10.3f}")


if __name__ == "__main__":
    main()
