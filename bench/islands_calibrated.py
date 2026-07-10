#!/usr/bin/env python
"""Do epitope-specific TCRs form islands, once the neighbour test is calibrated?

The published island analysis joined two CDR3s whenever ``gapblock_score <= theta`` for a
fixed theta. That is not a significance test. The control repertoire is far denser near
germline than in the rare-junction region, so a fixed theta buys many more chance neighbours
for a common junction than for a rare one. Measured here: at theta = 60 a single ball holds a
median of 3,245 of the 250,000 control sequences (1.3%), and random control sequences form
islands *more* readily than real same-epitope ones.

This script replaces the fixed cutoff with the per-query, control-calibrated one from
``seqtree.threshold_for_evalue``: node ``a`` gets the largest cutoff ``theta(a)`` whose ball
holds at most ``e_target`` expected chance neighbours *within a group of this size*, and an
edge is drawn only when ``score(a,b) <= min(theta(a), theta(b))`` -- mutual significance.

By construction both arms then admit the same expected number of chance edges per node
(``e_target``), so the arms are directly comparable and the control arm doubles as a
calibration check: its realised edge rate must land near ``e_target``.

Verdict criterion: real epitope groups must form more, and larger, islands than
size-matched random control groups. If they do not, the island abstraction is dead.

Usage:
    python bench/islands_calibrated.py [--e-target 0.05] [--theta-max 50] [--cap 600]

2026-07-10
"""
from __future__ import annotations

import argparse
import collections
import gzip
import itertools
import random
import statistics as stt
import sys
import time
from importlib import resources

import polars as pl

import seqtree
from seqtree.evalue import thetas_from_scores
from seqtree.gapblock import GapBlockIndex, central_prior, gapblock_score

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from _bench import source, valid_cdr3  # noqa: E402

from vdjmatch import db  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
D_MAX = 3  # gapblock is undefined past this; pairs further apart in length are never edges


# ---------------------------------------------------------------- data

def real_groups(src, species, gene, min_refs, min_size, cap, seed):
    """Epitopes backed by >= ``min_refs`` publications -> their unique CDR3s."""
    vdj = valid_cdr3(db.load(src, species=species)).filter(pl.col("gene") == gene)
    keep = (
        vdj.group_by("epitope")
        .agg(pl.col("reference_id").n_unique().alias("n_refs"))
        .filter(pl.col("n_refs") >= min_refs)["epitope"]
    )
    # polars group_by does not guarantee iteration order, so a shared RNG would draw a
    # different subsample on every run. Sort the epitopes and seed each draw from its name.
    sub = vdj.filter(pl.col("epitope").is_in(keep))
    by_epitope = collections.defaultdict(set)
    for ep, cdr3 in zip(sub["epitope"].to_list(), sub["cdr3"].to_list()):
        by_epitope[ep].add(cdr3)
    out = {}
    for ep in sorted(by_epitope):
        seqs = sorted(by_epitope[ep])
        if len(seqs) < min_size:
            continue
        out[ep] = random.Random(f"{seed}:{ep}").sample(seqs, cap) if len(seqs) > cap else seqs
    return out


def control_groups(pool, sizes, seed):
    """Size-matched groups of random control junctions -- the null the real arm must beat."""
    rng = random.Random(seed)
    picks = rng.sample(pool, sum(sizes))
    out, at = {}, 0
    for k, n in enumerate(sizes):
        out[f"ctrl{k:03d}"] = picks[at:at + n]
        at += n
    return out


# ---------------------------------------------------------------- graph

def components(n, edges):
    """Union-find connected components; returns a component id per node."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        parent[find(a)] = find(b)
    roots = {}
    return [roots.setdefault(find(i), len(roots)) for i in range(n)]


def pair_scores(seqs, ceiling, pen, gap_open, prior):
    """Every within-group pair scoring at or below ``ceiling``, tagged with its block length.

    Computed once, thresholded many. Carrying ``d`` lets us report whether the surviving edges
    ever involve a gap -- if ``gap_open`` exceeds the calibrated cutoff, they structurally
    cannot, and any "islands are length-flat" conclusion would be an artifact of the gap cost.
    """
    out = []
    for i, j in itertools.combinations(range(len(seqs)), 2):
        a, b = seqs[i], seqs[j]
        d = abs(len(a) - len(b))
        if d > D_MAX:
            continue
        s, _ = gapblock_score(a, b, gap_open=gap_open, gap_extend=1, gap_prior=prior, _pen=pen)
        if s <= ceiling:
            out.append((i, j, s, d))
    return out


def calibrated_edges(scores, cutoffs):
    """Mutual significance: score <= min(theta(a), theta(b)). A -1 cutoff admits nothing."""
    return [(i, j, d) for i, j, s, d in scores
            if cutoffs[i] >= 0 and cutoffs[j] >= 0 and s <= min(cutoffs[i], cutoffs[j])]


def fixed_edges(scores, theta):
    """The published rule, kept for reference: one theta for everybody."""
    return [(i, j, d) for i, j, s, d in scores if s <= theta]


# ---------------------------------------------------------------- reporting

def summarise(sizes_by_group):
    """Island size histogram and the fraction of nodes in islands of >= k members."""
    sizes = [s for g in sizes_by_group for s in g]
    total = sum(sizes)
    bins = collections.Counter()
    for s in sizes:
        for lo, hi, name in ((1, 1, "1"), (2, 2, "2"), (3, 4, "3-4"), (5, 9, "5-9"),
                             (10, 19, "10-19"), (20, 49, "20-49"), (50, 10**9, "50+")):
            if lo <= s <= hi:
                bins[name] += 1
                break
    frac = {k: sum(s for s in sizes if s >= k) / total for k in (5, 10, 20)}
    return bins, frac, total


def islands_of(groups, scores_by_group, cutoffs_by_group, fixed=None):
    """Island size lists per group, mean edges per node, members of islands >= 5, and the
    block-length histogram of the surviving edges."""
    sizes, edge_rate, big, edge_d = [], [], [], collections.Counter()
    for name, seqs in groups.items():
        sc = scores_by_group[name]
        e = fixed_edges(sc, fixed) if fixed is not None else calibrated_edges(sc, cutoffs_by_group[name])
        edge_rate.append(2 * len(e) / len(seqs))
        for _, _, d in e:
            edge_d[d] += 1
        comp = components(len(seqs), [(i, j) for i, j, _ in e])
        members = collections.defaultdict(list)
        for node, cid in enumerate(comp):
            members[cid].append(seqs[node])
        sizes.append([len(v) for v in members.values()])
        big.extend(v for v in members.values() if len(v) >= 5)
    return sizes, edge_rate, big, edge_d


def length_structure(islands):
    """Does an island actually span several lengths? This is what a PWM frame would need."""
    if not islands:
        return None
    ranges = [max(map(len, m)) - min(map(len, m)) for m in islands]
    members = [len(s) for m in islands for s in m]
    widest = [max(map(len, m)) for m in islands for _ in m]
    gapped = sum(1 for w, ln in zip(widest, members) if w > ln)
    hist = collections.Counter(min(r, 6) for r in ranges)
    return {
        "n_islands": len(islands),
        "range_hist": hist,
        "flat": sum(1 for r in ranges if r == 0),
        "median_range": stt.median(ranges),
        "frac_members_needing_a_gap": gapped / len(members),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pmhc", default=None)
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--gene", default="TRB")
    ap.add_argument("--min-refs", type=int, default=2)
    ap.add_argument("--min-size", type=int, default=20)
    ap.add_argument("--cap", type=int, default=600)
    ap.add_argument("--e-target", type=float, default=0.05)
    ap.add_argument("--theta-max", type=int, default=50)
    ap.add_argument("--gap-open", type=int, default=None,
                    help="default 2*matrix.scale() = 28; a single indel then costs 28, "
                         "which exceeds a typical calibrated cutoff and makes gapped edges "
                         "unreachable. Pass 14 (= 1*scale) to let indels compete.")
    ap.add_argument("--fixed", type=int, nargs="*", default=[40, 60])
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    t_start = time.perf_counter()
    with resources.files("seqtree").joinpath("data/control_human_trb_aa.txt.gz").open("rb") as fh:
        pool = [ln.strip() for ln in gzip.open(fh, "rt") if ln.strip()]
    M = len(pool)

    real = real_groups(source(args.pmhc), args.species, args.gene,
                       args.min_refs, args.min_size, args.cap, args.seed)
    sizes = [len(v) for v in real.values()]
    ctrl = control_groups(pool, sizes, args.seed + 1)
    n_nodes = sum(sizes)
    print(f"control pool M = {M:,}")
    print(f"real: {len(real)} epitopes, {n_nodes:,} nodes "
          f"(sizes {min(sizes)}-{max(sizes)}, median {stt.median(sizes):.0f})")
    print(f"e_target = {args.e_target}, theta_max = {args.theta_max}, d_max = {D_MAX}")

    mat = seqtree.SubstitutionMatrix.blosum62()
    scale = mat.scale()
    gap_open = args.gap_open if args.gap_open is not None else 2 * scale
    prior = central_prior(int(1.5 * scale))
    pen = {(a, b): mat.penalty(a, b) for a in AA for b in AA}

    t0 = time.perf_counter()
    gbi = GapBlockIndex(pool, "aa", d_max=D_MAX)
    print(f"built GapBlockIndex(d_max={D_MAX}) over the control in {time.perf_counter()-t0:.1f}s")

    def cutoffs(groups, label):
        """One control search per node -> its calibrated cutoff."""
        t0, done, out, clipped = time.perf_counter(), 0, {}, 0
        for name, seqs in groups.items():
            scores = [
                [h[1] for h in gbi.search(s, args.theta_max, mat, gap_open=gap_open,
                                          gap_extend=1, gap_prior=prior)]
                for s in seqs
            ]
            th = thetas_from_scores(scores, len(seqs), M, args.e_target, args.theta_max,
                                    exclude_exact=True)
            clipped += sum(1 for t in th if t == args.theta_max)
            out[name] = th
            done += len(seqs)
            if done % 2000 < len(seqs):
                el = time.perf_counter() - t0
                print(f"  {label}: {done:,}/{sum(len(v) for v in groups.values()):,} nodes "
                      f"({el:.0f}s, {1e3*el/done:.1f} ms/node)", flush=True)
        flat = [t for v in out.values() for t in v]
        usable = [t for t in flat if t >= 0]
        print(f"  {label}: theta median {stt.median(usable) if usable else float('nan'):.0f}, "
              f"unreachable {len(flat)-len(usable):,}/{len(flat):,}, "
              f"clipped at theta_max {clipped:,} ({100*clipped/len(flat):.1f}%), "
              f"{time.perf_counter()-t0:.0f}s", flush=True)
        return out

    cut_real = cutoffs(real, "real")
    cut_ctrl = cutoffs(ctrl, "ctrl")

    ceiling = max([args.theta_max] + list(args.fixed))
    scores = {}
    for label, groups in (("real", real), ("ctrl", ctrl)):
        t0 = time.perf_counter()
        scores[label] = {n: pair_scores(s, ceiling, pen, gap_open, prior) for n, s in groups.items()}
        print(f"  {label}: pair scores <= {ceiling} in {time.perf_counter()-t0:.0f}s", flush=True)

    arms = (("real", real, cut_real), ("ctrl", ctrl, cut_ctrl))
    hdr = f"  {'arm':<6}{'nodes':>8}{'edges/node':>12}{'>=5':>9}{'>=10':>9}{'>=20':>9}"
    bin_names = ("1", "2", "3-4", "5-9", "10-19", "20-49", "50+")

    print(f"\n=== CALIBRATED EDGES (mutual significance, E* = {args.e_target}) ===")
    print(hdr)
    keep = {}
    for label, groups, cut in arms:
        isl, er, big, ed = islands_of(groups, scores[label], cut)
        bins, frac, total = summarise(isl)
        keep[label] = (bins, big, ed)
        print(f"  {label:<6}{total:>8,}{stt.mean(er):>12.3f}"
              f"{frac[5]:>9.3f}{frac[10]:>9.3f}{frac[20]:>9.3f}")
    print(f"\n  calibration check: the control arm's edges/node must sit near E* = {args.e_target}")
    ed = keep["real"][2]
    tot = sum(ed.values()) or 1
    print(f"\n  real calibrated edges by block length d (gap_open = {gap_open}, "
          f"gap_cost(1) = {gap_open}):")
    print("    " + "  ".join(f"d={k}: {ed[k]:,} ({100*ed[k]/tot:.1f}%)" for k in range(D_MAX + 1)))
    if ed[0] == tot:
        print("    !! every edge is ungapped. gap_cost(1) exceeds the calibrated cutoff, so a\n"
              "       length-different pair CANNOT be an edge. Any length-flatness below is an\n"
              "       artifact of gap_open, not biology. Re-run with a smaller --gap-open.")
    print("\n  island size histogram")
    print(f"  {'arm':<6}" + "".join(f"{k:>9}" for k in bin_names))
    for label, _, _ in arms:
        print(f"  {label:<6}" + "".join(f"{keep[label][0][k]:>9,}" for k in bin_names))

    print("\n=== LENGTH STRUCTURE of real calibrated islands >= 5 (does a PWM frame need gaps?) ===")
    ls = length_structure(keep["real"][1])
    if ls is None:
        print("  no islands >= 5")
    else:
        print(f"  islands >= 5: {ls['n_islands']}   length-flat (range 0): {ls['flat']} "
              f"({100*ls['flat']/ls['n_islands']:.0f}%)   median range {ls['median_range']:.0f}")
        print("  length range: " + "  ".join(f"{r}{'+' if r == 6 else ''}:{ls['range_hist'][r]}"
                                             for r in range(7)))
        print(f"  members needing >= 1 gap column: {100*ls['frac_members_needing_a_gap']:.1f}%")

    for theta in args.fixed:
        print(f"\n=== FIXED theta = {theta} (the published rule, for reference) ===")
        print(hdr)
        for label, groups, cut in arms:
            isl, er, _, _ = islands_of(groups, scores[label], cut, fixed=theta)
            bins, frac, total = summarise(isl)
            print(f"  {label:<6}{total:>8,}{stt.mean(er):>12.3f}"
                  f"{frac[5]:>9.3f}{frac[10]:>9.3f}{frac[20]:>9.3f}")

    print(f"\ntotal wall time {time.perf_counter()-t_start:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
