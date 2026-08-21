"""Per-method bootstrap 95% CI on clustering macro-purity + retention (TRB shortlist), so we can test
whether the high-purity clusterers are statistically TIED in purity and where vdjmatch sits on the
purity/retention frontier. Re-runs each tool (its labels are not persisted), builds labels aligned to the
shared clonotype set, and bootstraps (cluster_results.bootstrap_cluster_ci). Deterministic tools.

    .venv/bin/python bench/cluster_method_ci.py giana ismart clustcr
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path.home() / "vcs/manuscripts/2026-vdjmatch/benchmarks/scripts"))
import _cluster_common as C                                          # noqa: E402
import benchmark as B                                                # noqa: E402
import cluster_results as CR                                         # noqa: E402

def labels_for(method, sl, cdr3, v, epi):
    """Each tool's labels aligned to the shared TRB shortlist set (re-runs the tool in its conda env)."""
    import clustcr_cluster, giana_cluster, gliph2_cluster_perepi, ismart_cluster, tcrdist_cluster_perepi
    if method == "giana":
        return C.labels_from_groups(cdr3, giana_cluster.run_giana(cdr3, v))
    if method == "ismart":
        return C.labels_from_groups(cdr3, ismart_cluster.run_ismart(cdr3, v))
    if method == "clustcr":
        return C.labels_from_groups(cdr3, clustcr_cluster.run_clustcr(cdr3))
    if method == "gliph2":
        return C.labels_from_groups(cdr3, gliph2_cluster_perepi.run_gliph2(cdr3, v, gliph2_cluster_perepi._jvec(cdr3)))
    if method == "tcrdist":
        return tcrdist_cluster_perepi.single(sl, "TRB")[0]
    raise ValueError(method)


def main(methods):
    d = B.release("vdjdb2026"); sl = B.shortlist(d)
    cdr3, epi, v = CR.single_clonotypes(sl, "TRB")
    print(f"\n=== TRB shortlist: per-method clustering CI (n={len(cdr3)}) ===")
    print(f"{'method':9}{'purity (95% CI)':>26}{'retention (95% CI)':>26}", flush=True)
    for m in (["vdjmatch"] + list(methods)):
        labels = CR.vdjmatch_labels("TRB", sl)[0] if m == "vdjmatch" else labels_for(m, sl, cdr3, v, epi)
        (p, plo, phi), (r, rlo, rhi) = CR.bootstrap_cluster_ci(labels, epi, B=2000)
        print(f"{m:9}{f'{p:.3f} [{plo:.3f},{phi:.3f}]':>26}{f'{r:.3f} [{rlo:.3f},{rhi:.3f}]':>26}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or ["tcrdist", "gliph2", "giana", "ismart", "clustcr"])
