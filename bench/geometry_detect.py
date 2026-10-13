#!/usr/bin/env python3
"""vdjmatch junction-geometry detection producer: the hook for a geometry arm.

2026-10-13. The geometry sweep that selected the joint V-J stem ran in the analysis
workspace, so until now no script in ``bench/`` could exercise a geometry arm at all and
the selected geometry was not reproducible from this repository. This closes that gap.

It scores the Group-A detection tasks -- NLV+/tet-, LLW vs LLL, and YLQ/GLC on TCRvdb,
which is the leakage-free stratum where every external method is comparable -- with our
own :func:`vdjmatch.match.tcrdist.distance_matrix` under a *selectable* geometry, and
emits the same ``predictions/<method>/<task>_<locus>.tsv`` contract as every other
producer, so the result drops straight into the existing comparison.

    python bench/geometry_detect.py                      # shipped TCRdist3 geometry
    python bench/geometry_detect.py --stem-depth 4       # the selected arm
    python bench/geometry_detect.py --compare            # both, side by side

**What this is not.** The arms recorded in the analysis workspace score the kernel-ladder
statistic against a control population; this producer scores nearest-neighbour distance
and an in-radius kernel sum, matching ``tcrdist_detect.py`` so the two are like-for-like.
It therefore measures *the geometry*, with the scoring held fixed and simple. Its numbers
are not the ladder-statistic numbers and must not be quoted as them.

Exact junction identities are punctured from the reference, as the honest-benchmark
contract requires. Queries whose V call does not resolve are dropped and counted, never
silently.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench
import benchmark as B
from vdjmatch import db
from vdjmatch.match.tcrdist import (distance_matrix, distance_units, load_v_loops,
                                    resolve_v_alleles, _substitution_costs)

OUT = Path("bench/predictions")
CELL_BUDGET = 8 << 20  # stay well inside distance_matrix's 64 MiB output cap


def resolvable(cdr3s, vs, model, locus):
    """Keep entries whose V call resolves to a modelled allele of this locus."""
    if not cdr3s:
        return [], [], 0
    alleles = resolve_v_alleles(pl.Series(vs), model, locus=locus).to_list()
    keep = [i for i, (c, a) in enumerate(zip(cdr3s, alleles))
            if a in model and model[a]["locus"] == locus and len(c) >= 8]
    return ([cdr3s[i] for i in keep], [alleles[i] for i in keep],
            len(cdr3s) - len(keep))


def score_task(queries, qv, refs, rv, species, threads, radius, geom):
    """Nearest distance and in-radius kernel sum per query, exact identities punctured."""
    nearest = np.full(len(queries), np.inf)
    kernel = np.zeros(len(queries))
    ref_arr = np.asarray(refs, dtype=object)
    batch = max(1, CELL_BUDGET // max(1, 4 * len(refs)))
    for off in range(0, len(queries), batch):
        sl = slice(off, min(off + batch, len(queries)))
        d = distance_matrix(queries[sl], refs, qv[sl], rv, species=species,
                            threads=threads, **geom).astype(np.float64)
        # Honest puncture: remove the query's own junction, not its whole epitope.
        same = np.asarray(queries[sl], dtype=object)[:, None] == ref_arr[None, :]
        d[same] = np.inf
        nearest[sl] = d.min(axis=1)
        kernel[sl] = np.where(d <= radius, np.exp(-d / radius), 0.0).sum(axis=1)
    return nearest, kernel


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--stem-depth", type=int, default=0,
                   help="joint V-J hairpin stem levels; 0 is TCRdist3, 4 is the selected arm")
    p.add_argument("--stem-weight", type=float, default=1.0)
    p.add_argument("--junction-ends", choices=["tcrdist", "trim3"], default="tcrdist")
    p.add_argument("--v-loop-weight", type=int, choices=[0, 1], default=1)
    p.add_argument("--alphabet", choices=["full", "imgt"], default="full")
    p.add_argument("--loci", nargs="+", default=["TRA", "TRB"])
    p.add_argument("--radius-factor", type=float, default=22.5,
                   help="in-radius kernel scale, as a factor of the matrix unit")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--species", default="human")
    p.add_argument("--compare", action="store_true",
                   help="run stem depth 0 and 4 and print the per-task delta")
    p.add_argument("--tag", default=None, help="override the predictions/<method> directory")
    a = p.parse_args()

    # The supplied-sample tasks read $VDJMATCH_TESTDATA; fail with the fix rather than
    # with a raw FileNotFoundError from three frames down.
    from compare import TESTDATA
    if not (TESTDATA / "sample1_cmv_5+reads.txt").exists():
        raise SystemExit(
            f"supplied test data not found under {TESTDATA}\n"
            "Set VDJMATCH_TESTDATA to the directory holding sample1_cmv_5+reads.txt, "
            "sample2_yf_bst2_5+reads.txt and sample6_TCRvdb.csv "
            "(the manuscript repository's test_data/).")

    depths = [0, a.stem_depth or 4] if a.compare else [a.stem_depth]
    pmax = int(distance_units(_substitution_costs(a.alphabet))[0])
    radius = a.radius_factor * pmax
    model = load_v_loops(a.species)
    v26 = db.load(_bench.source(), species="HomoSapiens")
    print(f"matrix unit {pmax}, kernel radius {radius:g} distance units, "
          f"loci {' '.join(a.loci)}")

    results: dict[tuple, dict] = {}
    for task, estr, locus, pos_qv, neg_qv in B.detection_tasks(a.loci):
        ref = (_bench.valid_cdr3(
            v26.filter((pl.col("epitope") == estr)
                       & pl.col("mhc_a").str.contains(B.A02)
                       & (pl.col("gene") == locus)))
            .select("cdr3", "v").unique("cdr3"))
        refs, rv, ref_dropped = resolvable(ref["cdr3"].to_list(), ref["v"].to_list(),
                                           model, locus)
        allq = list(pos_qv) + list(neg_qv)
        qs, qvv, q_dropped = resolvable(allq, [pos_qv.get(q, neg_qv.get(q)) for q in allq],
                                        model, locus)
        if not refs or not qs:
            print(f"{task}/{locus}: no usable reference or query after V resolution")
            continue
        print(f"{task}/{locus}: {len(qs)} queries ({q_dropped} dropped), "
              f"{len(refs)} reference keys ({ref_dropped} dropped)")
        for depth in depths:
            geom = dict(ctrim=2 if a.junction_ends == "tcrdist" else 3,
                        v_loop_weight=a.v_loop_weight, alphabet=a.alphabet,
                        stem_depth=depth, stem_weight=a.stem_weight)
            nearest, kernel = score_task(qs, qvv, refs, rv, a.species, a.threads,
                                         radius, geom)
            method = a.tag or f"vdjmatch-geom-stem{depth}"
            (OUT / method).mkdir(parents=True, exist_ok=True)
            pl.DataFrame({"query_id": qs, "epitope": [estr] * len(qs),
                          # higher is better; nearest distance is inverted
                          "score": (-nearest).tolist(),
                          "kernel_sum": kernel.tolist(),
                          "significant": (nearest <= radius).tolist()}
                         ).write_csv(OUT / method / f"{task}_{locus}.tsv", separator="\t")
            sc = B.read_predictions(OUT / method / f"{task}_{locus}.tsv", qs, [estr])
            m = B.classify_metrics(sc, [q for q in pos_qv if q in qs],
                                   [q for q in neg_qv if q in qs], estr)
            results[(task, locus, depth)] = m
            print(f"  stem {depth}: ROC={m['roc_auc']:.4f} PR={m['pr_auc']:.4f} "
                  f"n_pos={m['n_pos']} n_neg={m['n_neg']} -> {method}/")

    if a.compare and results:
        base, other = depths
        print(f"\nper-task delta, stem {other} minus stem {base} (positive favours the stem)")
        print(f"{'task':>8} {'locus':>6} {'n_pos':>6} {'n_neg':>6} "
              f"{'dROC':>9} {'dPR':>9}")
        keys = sorted({(t, l) for t, l, _ in results})
        dr = [results[(t, l, other)]["roc_auc"] - results[(t, l, base)]["roc_auc"]
              for t, l in keys if (t, l, other) in results and (t, l, base) in results]
        dp = [results[(t, l, other)]["pr_auc"] - results[(t, l, base)]["pr_auc"]
              for t, l in keys if (t, l, other) in results and (t, l, base) in results]
        for (t, l), x, y in zip(keys, dr, dp):
            m = results[(t, l, base)]
            print(f"{t:>8} {l:>6} {m['n_pos']:>6} {m['n_neg']:>6} {x:>+9.4f} {y:>+9.4f}")
        print(f"{'mean':>8} {'':>6} {'':>6} {'':>6} "
              f"{float(np.mean(dr)):>+9.4f} {float(np.mean(dp)):>+9.4f}  "
              f"over {len(dr)} tasks")


if __name__ == "__main__":
    main()
