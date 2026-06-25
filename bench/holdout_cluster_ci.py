"""Bootstrap 95% CIs on clustering purity + retention (held-out clustering settings), and test whether
the differences (shortlist vs full reference; untrimmed vs apex-trim) are real or within noise.

Stratified bootstrap: resample clonotypes within epitope, cluster labels fixed (cluster_results.
bootstrap_cluster_ci). Two settings differ significantly if their 95% CIs do not overlap.

    .venv/bin/python bench/holdout_cluster_ci.py [TRA|TRB]
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path.home() / "vcs/manuscripts/2026-vdjmatch/benchmarks/scripts"))
import benchmark as B                                                # noqa: E402
import cluster_results as CR                                         # noqa: E402
import cluster_trim as CT                                            # noqa: E402


def setting_ci(cdr3, v, epi, trim, B_=2000):
    lab = CT.cluster_labels(cdr3, v, epi, trim)
    return CR.bootstrap_cluster_ci(lab, epi, B=B_)                   # ((p,plo,phi),(r,rlo,rhi))


def overlap(a, b):
    return not (a[2] < b[1] or b[2] < a[1])                          # CI (val,lo,hi) overlap test


def main(which):
    d = B.release("vdjdb2026"); sl = B.shortlist(d)
    print(f"\n=== {which}: clustering purity / retention with 95% CI ===")
    print(f"{'ref':10}{'trim':6}{'purity (95% CI)':>26}{'retention (95% CI)':>26}")
    res = {}
    for ref_name, df in (("shortlist", sl), ("full", d)):
        cdr3, epi, v = CR.single_clonotypes(df, which)
        for trim in (False, True):
            p, r = setting_ci(cdr3, v, epi, trim)
            res[(ref_name, "apex" if trim else "none")] = (p, r)
            print(f"{ref_name:10}{'apex' if trim else 'none':6}"
                  f"{f'{p[0]:.3f} [{p[1]:.3f},{p[2]:.3f}]':>26}{f'{r[0]:.3f} [{r[1]:.3f},{r[2]:.3f}]':>26}")
    print("\nSignificance (non-overlapping 95% CIs = real difference):")
    cmps = [(("shortlist", "none"), ("full", "none"), "shortlist vs full (purity)", 0),
            (("shortlist", "none"), ("full", "none"), "shortlist vs full (retention)", 1),
            (("shortlist", "none"), ("shortlist", "apex"), "none vs apex, shortlist (retention)", 1),
            (("shortlist", "none"), ("shortlist", "apex"), "none vs apex, shortlist (purity)", 0)]
    for a, b, lbl, m in cmps:
        ov = overlap(res[a][m], res[b][m])
        va, vb = res[a][m][0], res[b][m][0]
        print(f"  {lbl:42} {va:.3f} vs {vb:.3f}  -> {'overlap (NOT sig)' if ov else 'DISJOINT (significant)'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TRB")
