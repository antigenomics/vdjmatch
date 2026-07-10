#!/usr/bin/env python
"""A4: does a per-island PWM beat min-over-members? Held out, at a matched FPR.

The plan pre-registered this test and never ran it. It was killed on two arguments, one of
which turned out to be an artifact of the biased control ("only ~19 net real islands exist" --
under calibrated edges and a uniform control, real forms 251 islands of >= 5 and the control
forms none). This settles it on evidence instead.

Setup. Build calibrated islands exactly as ``islands_calibrated.py`` does: mutual-significance
E-value edges over VDJdb same-epitope junctions. Keep islands of at least ``--min-island``
members. Hold out a fraction of each island's members; both scorers see only the rest.

    min-over-members   s(q) = min_j gapblock_score(q, member_j)
    per-island PWM     embed the training members in a frame at column c, then
                       s(q) = sum_j pen(j, embed(q)[j]),
                       pen(j, a) = round(lam * log(p_max_j / p_j(a)))

The gap is a column symbol, so a column never gapped in training charges for a gap there. No
separate affine term is needed: the island's own members say where a gap is tolerated.

That penalty form is not incidental. A textbook PWM log-odds score is **signed**, which breaks
``s >= 0`` and ``s(q, q) = 0`` and so cannot flow through ``seqtree.evalues``. Scoring each column
against its own consensus instead keeps the score a non-negative penalty that is zero on the
consensus sequence -- a ball centred on the consensus rather than on any one member.

The frame column c is chosen per island as the argmin of summed column entropy over the training
members, so the PWM is given the best frame it can have. lam only sets resolution: both scorers are
compared at a matched 1% false-positive rate against the same control negatives, so any monotone
rescaling of either score cancels.

Positives are the held-out members. Negatives are random control junctions -- the same null the
cutoffs were calibrated against.

Run (from the repo root, venv active):
    python bench/island_pwm.py --min-island 10 --repeats 5
"""
from __future__ import annotations

import argparse
import bisect
import collections
import math
import random
import statistics as stt
import sys
import time
from pathlib import Path

import numpy as np

import seqtree
from seqtree.gapblock import central_prior, embed_in_frame, score_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
from islands_calibrated import (  # noqa: E402
    D_MAX, calibrated_edges, components, pair_scores, real_groups,
)
from _bench import source  # noqa: E402
from seqtree.evalue import thetas_from_scores  # noqa: E402
from seqtree.gapblock import GapBlockIndex  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
ALPHA = AA + "-"
LAM = 1000            # PWM resolution only; matched-FPR comparison is scale-free
PSEUDO = 0.5


# ---------------------------------------------------------------- the PWM

def column_entropy(members, width, c):
    """Summed Shannon entropy of the frame's columns. Gap columns count as one symbol."""
    cols = [collections.Counter() for _ in range(width)]
    for s in members:
        for j, ch in enumerate(embed_in_frame(s, width, c)):
            cols[j][ch] += 1
    total = 0.0
    for col in cols:
        n = sum(col.values())
        total += -sum((v / n) * math.log(v / n) for v in col.values() if v)
    return total


def best_frame_column(members, width):
    """Give the PWM the frame it would choose for itself."""
    c_max = min(members, key=len)
    return min(range(len(c_max) + 1), key=lambda c: column_entropy(members, width, c))


def build_pwm(members, width, c):
    """pen[j][a] = round(LAM * log(p_max_j / p_j(a))): non-negative, zero at the consensus.

    The gap symbol is a column symbol like any other, so a column that is never gapped in
    training charges heavily for a gap there. No separate affine gap term is needed -- or wanted:
    the island's own members say where a gap is tolerated.
    """
    pen = []
    for j in range(width):
        counts = collections.Counter()
        for s in members:
            counts[embed_in_frame(s, width, c)[j]] += 1
        n = sum(counts.values())
        denom = n + PSEUDO * len(ALPHA)
        probs = {a: (counts.get(a, 0) + PSEUDO) / denom for a in ALPHA}
        pmax = max(probs.values())
        pen.append({a: int(round(LAM * math.log(pmax / p))) for a, p in probs.items()})
    return pen


def pwm_score(pen, width, c, q):
    """Score q in the island's frame. Longer than the frame, or shorter than c => rejected."""
    if len(q) > width or c > len(q):
        return math.inf
    return sum(pen[j].get(ch, LAM * 6) for j, ch in enumerate(embed_in_frame(q, width, c)))


GAP = ord("-")


def encode_by_length(seqs):
    """{length: (row_indices, int16 char-code matrix)} -- built once, reused by every island."""
    out = collections.defaultdict(list)
    for i, s in enumerate(seqs):
        out[len(s)].append(i)
    return {L: (np.array(idx), np.array([[ord(ch) for ch in seqs[i]] for i in idx], dtype=np.int16))
            for L, idx in out.items()}


def pwm_scores_bulk(pen, width, c, n, groups):
    """Vectorised pwm_score over a pre-encoded corpus. The operational FPR needs 250k negatives,
    and a Python loop over them costs minutes per island-split."""
    tbl = np.full((width, 128), LAM * 6, dtype=np.int64)
    for j, col in enumerate(pen):
        for a, v in col.items():
            tbl[j, ord(a)] = v
    out = np.full(n, np.inf)
    for L, (idx, G) in groups.items():
        if L > width or c > L:
            continue                                   # unrepresentable: stays +inf, i.e. rejected
        d = width - L
        # frame column j reads source column j (before the block), the gap (inside), or j-d (after)
        sc = np.zeros(len(idx), dtype=np.int64)
        for j in range(width):
            if j < c:
                sc += tbl[j, G[:, j]]
            elif j < c + d:
                sc += tbl[j, GAP]
            else:
                sc += tbl[j, G[:, j - d]]
        out[idx] = sc
    return out


# ---------------------------------------------------------------- evaluation

def recall_at_fpr(pos, neg, fpr=0.01):
    """Recall at the largest cutoff that admits AT MOST ``fpr`` of the negatives.

    Not ``sorted(neg)[k]``. Gap-block scores are small integers and tie heavily, so the k-th
    smallest can admit hundreds of negatives -- handing min-over-members an effective FPR far
    looser than the one claimed, and an unearned recall advantage over a finely-graded PWM. The
    conservative rule below is also exactly what ``thetas_from_scores`` does: it picks the largest
    theta whose control-hit count stays within budget.

    A negative the scorer cannot represent scores +inf: correctly rejected, still in the
    denominator, so the PWM gets no free specificity either.

    Returns (recall, realised_fpr).
    """
    if not neg:
        return float("nan"), float("nan")
    s = sorted(neg)
    budget = int(fpr * len(s))
    best, admitted = None, 0
    for v in sorted(set(s)):
        n_le = bisect.bisect_right(s, v)
        if n_le > budget:
            break
        best, admitted = v, n_le
    if best is None:                       # even the smallest value blows the budget
        return sum(1 for x in pos if x < s[0]) / len(pos), 0.0
    return sum(1 for x in pos if x <= best) / len(pos), admitted / len(s)


# The middle two are the operating points, computed rather than guessed.
#
#   islands_calibrated passes N = the epitope group size (median 88 of 107 groups), so
#   k = floor(E* * M / N) has median 142 of M = 250k  ->  FPR 5.7e-4.
#
#   Annotating a whole repertoire instead puts N = 19,703. Then E* = 0.05 gives k = 0, which
#   thetas_from_scores reports as -1: the rule of three forbids certifying any E below 3N/M =
#   0.236. At that smallest certifiable E*, k = 3  ->  FPR 1.2e-5.
#
# 1% is a loose reference nobody uses; 4e-6 (one control neighbour) is the limit of this control.
FPRS = (0.01, 5.68e-4, 1.2e-5, 4e-6)
FPR_LABEL = {0.01: "loose ref", 5.68e-4: "per-epitope", 1.2e-5: "repertoire", 4e-6: "limit"}


def evaluate(island, outside, negatives, neg_groups, mat, gap_open, prior, rng, holdout=0.3):
    """Two positive sets.

    ``held-out``  members of this island the scorers never saw. Both scorers should ace this: the
                  island was built from calibrated edges, so every member is within theta of some
                  other member *by construction*. It is a floor, not a test.
    ``outside``   junctions specific to the same epitope that fell in a different island (or in
                  none). This is the question a profile is supposed to answer better than a
                  nearest-neighbour rule: does the island generalise to its epitope?
    """
    k = max(1, int(round(holdout * len(island))))
    shuffled = list(island)
    rng.shuffle(shuffled)
    test, train = shuffled[:k], shuffled[k:]
    if len(train) < 3:
        return None

    def mins(qs):
        if not qs:
            return []
        return list(np.asarray(score_matrix(qs, train, mat, gap_open=gap_open,
                                            gap_prior=prior)).min(1))

    mn_neg = mins(negatives)
    mn_held, mn_out = mins(test), mins(outside)

    width = max(len(s) for s in island)
    c = best_frame_column(train, width)
    pen = build_pwm(train, width, c)
    pw_neg = list(pwm_scores_bulk(pen, width, c, len(negatives), neg_groups))
    pw = lambda qs: [pwm_score(pen, width, c, q) for q in qs]  # noqa: E731
    pw_held, pw_out = pw(test), pw(outside)

    row = {
        "n": len(island),
        "train": len(train),
        "outside": len(outside),
        "lengths": len({len(s) for s in island}),
        "c": c,
        "width": width,
        "pwm_unscorable_held": sum(1 for x in pw_held if x == math.inf) / len(pw_held),
    }
    for f in FPRS:
        row[f"min_held@{f}"], row[f"min_fpr@{f}"] = recall_at_fpr(mn_held, mn_neg, f)
        row[f"pwm_held@{f}"], row[f"pwm_fpr@{f}"] = recall_at_fpr(pw_held, pw_neg, f)
        if outside:
            row[f"min_out@{f}"] = recall_at_fpr(mn_out, mn_neg, f)[0]
            row[f"pwm_out@{f}"] = recall_at_fpr(pw_out, pw_neg, f)[0]
    return row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pmhc", default=None)
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--gene", default="TRB")
    ap.add_argument("--control", default="human_trb_aa")
    ap.add_argument("--min-refs", type=int, default=2)
    ap.add_argument("--min-size", type=int, default=20)
    ap.add_argument("--cap", type=int, default=600)
    ap.add_argument("--e-target", type=float, default=0.05)
    ap.add_argument("--theta-max", type=int, default=50)
    ap.add_argument("--gap-open", type=int, default=None)
    ap.add_argument("--min-island", type=int, default=10)
    ap.add_argument("--negatives", type=int, default=10**9, help="default: the whole control")
    ap.add_argument("--repeats", type=int, default=3, help="held-out splits per island")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    t0 = time.perf_counter()
    ctrl = seqtree.load_control(args.control)
    pool = [ctrl.ref_seq(i) for i in range(len(ctrl))]
    M = len(pool)

    mat = seqtree.SubstitutionMatrix.blosum62()
    scale = mat.scale()
    gap_open = args.gap_open if args.gap_open is not None else 2 * scale
    prior = central_prior(int(1.5 * scale))
    syms = seqtree.alphabet_symbols("aa")
    pen_tbl = {(a, b): mat.penalty(a, b) for a in syms for b in syms}

    real = real_groups(source(args.pmhc), args.species, args.gene,
                       args.min_refs, args.min_size, args.cap, args.seed)
    print(f"control M = {M:,}; real: {len(real)} epitopes, "
          f"{sum(len(v) for v in real.values()):,} junctions")
    print(f"gap_open = {gap_open}, e_target = {args.e_target}, d_max = {D_MAX}")

    gbi = GapBlockIndex(pool, "aa", d_max=D_MAX)
    print(f"built GapBlockIndex over the control in {time.perf_counter()-t0:.1f}s")

    # ---- islands, exactly as islands_calibrated builds them
    islands = []                       # (members, same-epitope junctions outside this island)
    t1 = time.perf_counter()
    for name, seqs in real.items():
        scores = [[h[1] for h in gbi.search(s, args.theta_max, mat, gap_open=gap_open,
                                            gap_extend=1, gap_prior=prior)] for s in seqs]
        cut = thetas_from_scores(scores, len(seqs), M, args.e_target, args.theta_max,
                                 exclude_exact=True)
        ps = pair_scores(seqs, args.theta_max, pen_tbl, gap_open, prior)
        edges = calibrated_edges(ps, cut)
        comp = components(len(seqs), [(i, j) for i, j, _ in edges])
        members = collections.defaultdict(list)
        for node, cid in enumerate(comp):
            members[cid].append(seqs[node])
        for v in members.values():
            if len(v) >= args.min_island:
                inside = set(v)
                islands.append((v, [s for s in seqs if s not in inside]))
    print(f"calibrated islands >= {args.min_island}: {len(islands)} "
          f"({sum(len(m) for m, _ in islands):,} members) in {time.perf_counter()-t1:.1f}s")

    multi = sum(1 for m, _ in islands if len({len(s) for s in m}) > 1)
    print(f"  spanning more than one length: {multi} ({100*multi/len(islands):.0f}%)")
    print(f"  same-epitope junctions outside their island, median: "
          f"{stt.median([len(o) for _, o in islands]):.0f}\n")

    negatives = pool if args.negatives >= M else random.Random(args.seed + 7).sample(pool, args.negatives)
    neg_groups = encode_by_length(negatives)
    print(f"negatives: {len(negatives):,} control junctions "
          f"(the {int(FPRS[-1]*len(negatives))+1}th lowest sets the operational cutoff)")

    rows = []
    t2 = time.perf_counter()
    for r in range(args.repeats):
        rng = random.Random(args.seed + 1000 * r)
        for i, (isl, outside) in enumerate(islands):
            got = evaluate(isl, outside, negatives, neg_groups, mat, gap_open, prior, rng)
            if got:
                got["repeat"], got["island"] = r, i
                rows.append(got)
    print(f"evaluated {len(rows)} island-splits in {time.perf_counter()-t2:.1f}s\n")

    def report(kind, sel, label):
        sub = [r for r in rows if sel(r) and f"min_{kind}@{FPRS[0]}" in r]
        if not sub:
            print(f"  {label:<30}  (none)")
            return
        cells = ""
        for f in FPRS:
            mn = stt.mean(r[f"min_{kind}@{f}"] for r in sub)
            pw = stt.mean(r[f"pwm_{kind}@{f}"] for r in sub)
            cells += f"{100*mn:>9.1f}%{100*pw:>9.1f}%"
        wins = sum(1 for r in sub if r[f"pwm_{kind}@{FPRS[-1]}"] > r[f"min_{kind}@{FPRS[-1]}"])
        print(f"  {label:<30}{len(sub):>6}{cells}{100*wins/len(sub):>11.1f}%")

    hdr = ("".join(f"{'min':>9}{'PWM':>9}" for _ in FPRS))
    fprhdr = "".join(f"{f'-- {FPR_LABEL[f]} {100*f:.4g}% --':>18}" for f in FPRS)
    for kind, title in (("held", "held-out members of the island (a floor: both should ace it)"),
                        ("out", "same-epitope junctions OUTSIDE the island (generalisation)")):
        print(f"=== {title} ===")
        print(f"  {'stratum':<30}{'splits':>6}{fprhdr}{'PWM wins':>11}")
        print(f"  {'':<30}{'':>6}{hdr}{f'@{100*FPRS[-1]:g}%':>11}")
        report(kind, lambda r: True, "all islands")
        report(kind, lambda r: r["lengths"] == 1, "length-flat (frame is a no-op)")
        report(kind, lambda r: r["lengths"] > 1, "multi-length (frame matters)")
        for lo, hi in ((10, 19), (20, 49), (50, 10 ** 9)):
            report(kind, lambda r, lo=lo, hi=hi: lo <= r["n"] <= hi,
                   f"size {lo}-{hi if hi < 10**9 else '+'}")
        print()

    # ---- the decisive statistic: paired PWM - min, bootstrapped over ISLANDS ----------------
    # The 5 repeats of one island share its members, so they are not independent. Average within
    # an island first, resample islands.
    by_island = collections.defaultdict(list)
    for r in rows:
        by_island[r["island"]].append(r)
    keys = sorted(by_island)

    print("=== paired difference (PWM - min-over-members), bootstrapped over islands ===")
    print(f"  {len(keys)} islands, {args.repeats} splits each; positive means the PWM is better\n")
    print(f"  {'positives':<16}{'regime':>12}{'FPR':>9}{'mean diff':>12}{'95% CI':>20}"
          f"{'PWM better':>12}{'tie':>8}{'worse':>8}")
    boot = random.Random(args.seed + 99)
    for kind, label in (("held", "held-out"), ("out", "outside")):
        for f in FPRS:
            per = []
            for k in keys:
                sub = [r for r in by_island[k] if f"min_{kind}@{f}" in r]
                if sub:
                    per.append(stt.mean(r[f"pwm_{kind}@{f}"] - r[f"min_{kind}@{f}"] for r in sub))
            if not per:
                continue
            draws = sorted(stt.mean([per[boot.randrange(len(per))] for _ in per])
                           for _ in range(2000))
            lo, hi = draws[50], draws[1949]
            better = sum(1 for d in per if d > 1e-9)
            worse = sum(1 for d in per if d < -1e-9)
            tie = len(per) - better - worse
            print(f"  {label:<16}{FPR_LABEL[f]:>12} {100*f:>8.4g}%{100*stt.mean(per):>+11.2f}%"
                  f"{f'[{100*lo:+.2f}, {100*hi:+.2f}]':>20}"
                  f"{better:>12}{tie:>8}{worse:>8}")

    # ---- cost ------------------------------------------------------------------------------
    sizes = [r["train"] for r in rows]
    med_train, med_w = stt.median(sizes), stt.median([r["width"] for r in rows])
    pwm_bytes = med_w * len(ALPHA) * 4
    mem_bytes = med_train * med_w
    print(f"\n=== cost ===")
    print(f"  median island trains on {med_train:.0f} members of length {med_w:.0f}: "
          f"{mem_bytes:.0f} B of raw sequence.")
    print(f"  its PWM is {med_w:.0f} columns x {len(ALPHA)} symbols x 4 B = {pwm_bytes:.0f} B "
          f"({pwm_bytes/mem_bytes:.1f}x larger).")
    print(f"  break-even: an island needs > {len(ALPHA)*4} members before the PWM is the smaller")
    big = sum(1 for r in rows if r["train"] > len(ALPHA) * 4)
    print(f"  of the two. {big}/{len(rows)} splits ({100*big/len(rows):.1f}%) reach that.")

    frames = collections.Counter(r["c"] for r in rows if r["lengths"] > 1)
    if frames:
        print(f"\n  entropy-optimal frame column c on multi-length islands: "
              f"{dict(sorted(frames.items()))}")
        print(f"  modal c = {frames.most_common(1)[0][0]}; the structurally measured block sits at "
              f"Cys-offset 6.")
    unsc = stt.mean(r["pwm_unscorable_held"] for r in rows)
    print(f"  held-out members the PWM cannot represent (longer than the frame): {100*unsc:.1f}%")

    print("\n=== realised FPR (the cutoff admits at most the target; ties make min coarse) ===")
    print(f"  {'target FPR':>12}{'min realised':>16}{'PWM realised':>16}")
    for f in FPRS:
        mf = stt.mean(r[f"min_fpr@{f}"] for r in rows)
        pf = stt.mean(r[f"pwm_fpr@{f}"] for r in rows)
        print(f"  {100*f:>11.4f}%{100*mf:>15.4f}%{100*pf:>15.4f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
