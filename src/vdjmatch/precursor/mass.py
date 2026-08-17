"""`Pgen` masses of a cognate set: the observed bound, the closed ball, the union, motifs.

The one thing this module exists to get right is that **cognate junctions are near-duplicates**, so
their neighbourhoods overlap and their ball masses must not be added. :func:`union_mass` is the
exact union; :func:`ball_mass` is the enumeration route kept as its oracle.
"""
from __future__ import annotations

import csv
import gzip
import itertools
from collections import Counter
from dataclasses import dataclass
from math import comb
from pathlib import Path

from .model import _native, pgen

#: Cognacy retention per unit of edit distance -- the probability that a junction one substitution
#: away from a cognate TCR is itself cognate. Measured by Mayer & Callan, *PNAS* 2023;120:e2213264120
#: (PMID 36649423) from near-coincidence statistics within epitope-specific repertoires: binding
#: probability falls roughly **ten-fold per Levenshtein unit**, consistently across disparate
#: experiments. It is a parameter of :func:`shell_profile`, not a constant of the library -- pass
#: your own if you have measured one for your data.
ALPHA_PER_EDIT = 0.1

#: Default per-position residue-frequency cut for turning a VDJdb cluster PWM into an allowed-set
#: motif. Criterion: the largest round threshold at which a motif still matches **>= 95% of its own
#: cluster's member junctions**. Measured on the 2026-06 release (1,791 clusters, 53,757 members in
#: ``cluster_members.txt``) — member recall 0.9999 at 0.000, 0.996 at 0.002, 0.984 at 0.004, 0.960
#: at 0.008, 0.948 at 0.010, 0.267 at 0.10, 0.177 at 0.15. Recall falls off a cliff because a member
#: must clear the cut at *every* one of ~15 positions, so per-position losses compound; anything
#: above ~0.01 keeps only the consensus. Raise it for a tighter, lower-mass motif; set it to 0 for
#: exactly the residue alphabet the cluster was observed to use.
MOTIF_FREQ_THRESHOLD = 0.008

#: Ceiling on the number of enumerated neighbourhood members, for :func:`ball_mass`. The
#: enumeration is materialised as Python strings: 300 measured junctions at ``r=2`` produce ~9.9M
#: sequences and cost ~1.8 GB, i.e. ~190 bytes each. The default keeps a single call under ~0.4 GB.
#: :func:`union_mass` enumerates only *within a connected component*, and only the multiply-covered
#: members, so it does not come near this.
MAX_BALL_MEMBERS = 2_000_000


def observed_mass(model, junctions, v=None, j=None, threads: int = 0) -> float:
    """Sum of `Pgen` over the given junctions -- a strict lower bound on `F(e)`.

    ``v``/``j`` are per-junction allele-resolution call lists (``TRBV27*01``, not ``TRBV27``), or
    ``None`` to marginalise over V/J. The two are **different quantities** -- the marginal is larger
    -- so never mix them within one comparison.
    """
    return float(sum(pgen(model, junctions, v=v, j=j, threads=threads)))


# --------------------------------------------------------------------------- closed ball, per centre
def _wildcard_queries(seq: str, r: int):
    """The degenerate motifs of :func:`closed_ball_mass`: ``seq`` with ``<=r`` positions freed."""
    L = len(seq)
    base = list(seq)
    for k in range(r + 1):
        for positions in itertools.combinations(range(L), k):
            q = list(base)
            for p in positions:
                q[p] = ""                                # "" is vdjtools' wildcard
            yield k, q


def closed_ball_mass(model, junctions, r: int = 1, threads: int = 0) -> list[float]:
    """Mass of the closed Hamming-``r`` ball around each junction, in closed form.

    No enumeration at any radius. ``r=1`` is vdjtools' own ``pgen(..., mismatches=1)``. For
    ``r >= 2`` the same masked transfer-matrix DP gives the ball by an alternating sum over
    wildcarded motifs::

        m(B_r(a)) = sum_{k=0..r} (-1)^(r-k) * C(L-k-1, r-k) * sum_{|S|=k} m(W_S)

    where ``W_S`` is ``a`` with the positions in ``S`` freed. The coefficient counts how many
    ``W_S`` a sequence at distance ``d`` from ``a`` falls into, and is constructed so the alternating
    sum leaves 1 for every ``d <= r`` and 0 beyond. Cost is ``sum_k C(L,k)`` DP passes rather than
    ``|B_r|`` `Pgen` calls — 106 against 33,117 at ``L=14, r=2``.

    V/J are deliberately **not** accepted: a substituted neighbour need not keep the centre's V/J
    assignment, so conditioning the ball on the centre's call would be wrong. This marginalises.
    """
    seqs = [s for s in junctions if s]
    if not seqs:
        return []
    if r < 0:
        raise ValueError("r must be >= 0")
    if r == 0:
        return pgen(model, seqs, threads=threads)
    if r == 1:
        return pgen(model, seqs, mismatches=1, threads=threads)

    native = _native()
    flat, owner, weight = [], [], []
    for i, s in enumerate(seqs):
        L = len(s)
        if r > L:
            raise ValueError(f"r={r} exceeds the length of {s!r}; the ball is the whole space")
        for k, q in _wildcard_queries(s, r):
            flat.append(q)
            owner.append(i)
            weight.append((-1) ** (r - k) * comb(L - k - 1, r - k))

    out = [0.0] * len(seqs)
    step = 200_000                                       # bound peak memory on the query list
    for lo in range(0, len(flat), step):
        chunk = flat[lo:lo + step]
        vals = native.pgen_aa_degenerate_batch(model, chunk, threads=threads)
        for off, m in enumerate(vals):
            k = lo + off
            out[owner[k]] += weight[k] * float(m)
    return out


# --------------------------------------------------------------------------- the union
def _components(seqs: list[str], radius: int, threads: int = 0) -> list[list[int]]:
    """Connected components of the graph joining centres within ``radius`` substitutions."""
    from ..cluster import overlap

    parent = list(range(len(seqs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    pairs = overlap(seqs, scope=f"{radius},0,0,{radius}", threads=threads)
    for a, b in zip(pairs["a_idx"], pairs["b_idx"]):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    groups: dict[int, list[int]] = {}
    for i in range(len(seqs)):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


def union_mass(model, junctions, r: int = 1, threads: int = 0,
               max_members: int = MAX_BALL_MEMBERS, count_members: bool = True) -> dict:
    """Mass of the **union** of Hamming-``r`` balls — exact, and without enumerating the union.

    Returns ``{"union", "naive_sum", "overlap", "n_seqs", "n_union", "n_multiply_covered",
    "n_components", "n_clustered"}`` where ``overlap = 1 - union/naive_sum`` is the share of the
    naive per-sequence sum that double-counting would have invented, and ``n_union`` counts the
    distinct sequences the union holds — ``n_union - n_seqs`` of which are candidate cognate
    junctions no database has catalogued. ``n_union`` is ``None`` when counting them would exceed
    ``max_members``; the masses are unaffected, since they never enumerate.

    **How.** The naive sum counts every ``x`` in the union ``cov(x) = #{a : d(x,a) <= r}`` times, so

        m(union) = sum_a m(B_r(a))  -  sum_{x : cov(x) >= 2} (cov(x) - 1) * Pgen(x)

    exactly, with no inclusion–exclusion and hence no truncation error (truncating I–E at pairs and
    triples is *not* safe: a component of four mutually-close junctions has a non-empty four-way
    term). The first sum is closed-form (:func:`closed_ball_mass`). The multiply-covered set is
    ``union_{a != b} B_r(a) ∩ B_r(b)``, which is empty unless ``d(a,b) <= 2r``, so only centres
    inside one connected component of the ``2r`` graph can contribute — and within a component only
    the members covered twice need a `Pgen` call. Singleton components cost nothing at all, which
    matters: 41.8% of VDJdb human TRB junctions are singletons.

    V/J are not accepted, for the same reason as :func:`closed_ball_mass`: a substituted neighbour
    need not keep the centre's V/J.

    **Memory.** Finding the multiply-covered members counts neighbours **one connected component at
    a time**, so peak memory is the largest component's ball union rather than the whole set's — for
    a set whose junctions are all mutually distant it is zero. A single component above
    ``max_members`` raises rather than thrashing; split that epitope, or drop ``r``.
    """
    from seqtree.distance import neighbourhood, union_size

    seqs = list(dict.fromkeys(s for s in junctions if s))
    if not seqs:
        return {"union": 0.0, "naive_sum": 0.0, "overlap": 0.0, "n_seqs": 0, "n_union": 0,
                "n_multiply_covered": 0, "n_components": 0, "n_clustered": 0}

    def count_union(ss):
        """Distinct sequences in the union, or ``None`` when counting them is too expensive."""
        if not count_members:
            return None
        bound = sum(1 + 19 * len(s) if r == 1 else union_size([s], r=r) for s in ss)
        return union_size(ss, r=r) if bound <= max_members else None

    naive = float(sum(closed_ball_mass(model, seqs, r=r, threads=threads)))
    if r == 0:                                           # radius-0 balls are the points themselves
        return {"union": naive, "naive_sum": naive, "overlap": 0.0, "n_seqs": len(seqs),
                "n_union": len(seqs), "n_multiply_covered": 0, "n_components": len(seqs),
                "n_clustered": 0}
    if len(seqs) == 1:
        return {"union": naive, "naive_sum": naive, "overlap": 0.0, "n_seqs": 1,
                "n_union": count_union(seqs), "n_multiply_covered": 0, "n_components": 1,
                "n_clustered": 0}

    comps = _components(seqs, radius=2 * r, threads=threads)
    multi: list[str] = []
    excess: list[int] = []
    clustered = 0
    for comp in comps:
        if len(comp) < 2:
            continue
        clustered += len(comp)
        members = [seqs[i] for i in comp]
        size = union_size(members, r=r)
        if size > max_members:
            raise MemoryError(
                f"one connected component of {len(comp)} junctions has a radius-{r} union of "
                f"{size:,} sequences, above max_members={max_members:,}. Split this group or "
                f"lower r; components are independent, so their union masses add exactly.")
        counts: Counter[str] = Counter()
        for s in members:
            counts.update(neighbourhood(s, r=r))
        for x, c in counts.items():
            if c >= 2:
                multi.append(x)
                excess.append(c - 1)

    correction = 0.0
    if multi:
        ps = pgen(model, multi, threads=threads)
        correction = float(sum(e * p for e, p in zip(excess, ps)))

    union = naive - correction
    return {"union": union, "naive_sum": naive,
            "overlap": (1.0 - union / naive) if naive > 0 else 0.0,
            "n_seqs": len(seqs), "n_union": count_union(seqs), "n_multiply_covered": len(multi),
            "n_components": len(comps), "n_clustered": clustered}


def ball_mass(model, junctions, r: int = 1, threads: int = 0) -> dict:
    """Mass of the union of Hamming-``r`` balls, **by enumeration** — the oracle for
    :func:`union_mass`.

    Same return shape as :func:`union_mass` minus the component fields. Materialises the whole
    deduplicated union as Python strings and scores every member, so it is exact but scales as
    ``19L`` per centre at ``r=1`` and ``~180 L^2/2`` at ``r=2``. Prefer :func:`union_mass` for real
    work and keep this to regression-test it.

    V/J are deliberately **not** accepted: a substituted neighbour need not keep the centre's V/J
    assignment, so conditioning the ball on the centre's call would be wrong. This marginalises.
    """
    from seqtree.distance import neighbourhood_union
    native = _native()
    seqs = [s for s in junctions if s]
    if not seqs:
        return {"union": 0.0, "naive_sum": 0.0, "overlap": 0.0, "n_union": 0, "n_seqs": 0}

    members = neighbourhood_union(seqs, r=r)
    union = float(sum(native.pgen_aa_batch(model, list(members), threads=threads)))
    naive = float(sum(closed_ball_mass(model, seqs, r=r, threads=threads)))
    return {"union": union, "naive_sum": naive,
            "overlap": (1.0 - union / naive) if naive > 0 else 0.0,
            "n_union": len(members), "n_seqs": len(seqs)}


def shell_profile(model, junctions, r: int = 1, alpha: float = ALPHA_PER_EDIT,
                  threads: int = 0) -> dict:
    """:func:`union_mass` resolved by exact edit distance, with cognacy retention applied per shell.

    A ball at radius ``r`` treats a junction ``r`` substitutions from an observed cognate TCR as
    fully cognate, which it is not. Shell ``k`` is the set of sequences whose distance to the
    **nearest** observed junction is exactly ``k``, and the retained estimate is

    ``F ~= sum_k alpha**k * mass(shell k)``

    with ``alpha`` the per-edit cognacy retention, default :data:`ALPHA_PER_EDIT` (0.1, Mayer &
    Callan 2023). ``alpha=1`` reproduces the raw union; ``alpha=0`` collapses to
    :func:`observed_mass`.

    The shells are obtained by **differencing unions**, ``mass(shell k) = union(r=k) -
    union(r=k-1)``, which is exact because the min-distance shells partition the ball. So nothing is
    enumerated and there is no memory ceiling — the ``r=2`` profile that used to cost ~9.9M
    materialised strings for 300 junctions costs ``r+1`` calls to :func:`union_mass`.

    Returns ``{"shells": [{"r", "n", "mass", "alpha"}...], "retained", "union", "n_union",
    "n_seqs", "alpha", "overlap"}``. A shell's ``n`` is ``None`` when the union was too large to
    census; its ``mass`` never is, because the masses do not enumerate.
    """
    seqs = [s for s in junctions if s]
    if not seqs:
        return {"shells": [], "retained": 0.0, "union": 0.0, "n_seqs": 0, "alpha": alpha,
                "overlap": 0.0}
    if r < 0:
        raise ValueError("r must be >= 0")

    cumulative = [union_mass(model, seqs, r=k, threads=threads) for k in range(r + 1)]
    shells, retained, prev_mass, prev_n = [], 0.0, 0.0, 0
    for k, cum in enumerate(cumulative):
        mass = cum["union"] - prev_mass
        prev_mass = cum["union"]
        n = None if (cum["n_union"] is None or prev_n is None) else cum["n_union"] - prev_n
        prev_n = cum["n_union"]
        w = alpha ** k
        shells.append({"r": k, "n": n, "mass": mass, "alpha": w})
        retained += w * mass
    return {"shells": shells, "retained": retained, "union": cumulative[-1]["union"],
            "n_union": cumulative[-1]["n_union"], "n_seqs": len(seqs), "alpha": alpha,
            "overlap": cumulative[-1]["overlap"]}


# --------------------------------------------------------------------------- cluster PWM motifs
@dataclass(frozen=True)
class ClusterMotif:
    """One VDJdb cluster PWM as a degenerate motif ready for :func:`motif_mass`.

    ``allowed`` is one string of permitted residues per position (``""`` = wildcard). The cluster is
    V/J/length-pinned, so ``v``/``j`` are the conditioning to pass alongside it.
    """
    cid: str
    epitope: str
    gene: str
    species: str
    v: str
    j: str
    length: int
    size: int
    allowed: tuple[str, ...]


def load_cluster_motifs(path, threshold: float = MOTIF_FREQ_THRESHOLD, species: str | None = None,
                        gene: str | None = None, epitope: str | None = None,
                        min_size: int = 2) -> list[ClusterMotif]:
    """Read VDJdb's ``motif_pwms.txt`` into per-position allowed-residue sets.

    The file is one row per (cluster, position, residue) with ``freq`` the residue's frequency
    *within* the cluster. Because a cluster is pinned to one V, one J and one length, thresholding
    ``freq`` per position yields the ``allowed`` argument :func:`motif_mass` wants — no alignment,
    no register search, no enumeration.

    ``threshold`` is the per-position frequency cut, default :data:`MOTIF_FREQ_THRESHOLD`. A
    position never comes back empty: if no residue clears the cut the modal residue(s) are kept, so
    raising the threshold shrinks the motif monotonically towards the consensus sequence rather than
    zeroing its mass.

    **The file is not a complete count table**, and assuming it is silently drops real cluster
    members. Measured on the 2026-06 release: **2,326 of 24,036 listed positions (9.7%), spread over
    1,042 of 1,791 clusters (58%), have listed frequencies summing to less than 1**, and 172
    positions are missing outright — residues the release filtered out. In one worked case
    (``H.B.ATDALMTGY.5``) the position-12 row lists only ``F`` at ``freq=0.4`` while 6 of the 10
    member junctions carry ``Y`` there. So: whenever the listed frequencies fall short of 1 by more
    than ``threshold``, the unlisted residues could individually clear the cut and there is no way
    to know which they are, so **the position becomes a wildcard** rather than a set that excludes
    known members. Positions absent from the file are wildcards for the same reason.

    ``species``/``gene``/``epitope`` filter exactly; ``min_size`` drops clusters below that many
    members (``csz``). Accepts a plain or gzipped path.
    """
    p = Path(path)
    opener = (lambda: gzip.open(p, "rt")) if p.suffix == ".gz" else (lambda: open(p, "r"))
    rows: dict[str, dict] = {}
    with opener() as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            cid = rec["cid"]
            c = rows.get(cid)
            if c is None:
                c = rows[cid] = {"meta": rec, "pos": {}}
            c["pos"].setdefault(int(rec["pos"]), {})[rec["aa"]] = float(rec["freq"])

    out = []
    for cid, c in rows.items():
        meta = c["meta"]
        if species is not None and meta["species"] != species:
            continue
        if gene is not None and meta["gene"] != gene:
            continue
        if epitope is not None and meta["antigen.epitope"] != epitope:
            continue
        size = int(meta["csz"])
        if size < min_size:
            continue
        length = int(meta["len"])
        allowed = []
        for i in range(length):
            d = c["pos"].get(i)
            if not d or (1.0 - sum(d.values())) > threshold:
                allowed.append("")                          # unlisted residue mass -> wildcard
                continue
            keep = sorted(a for a, f in d.items() if f >= threshold)
            if not keep:                                    # never emit an empty set
                top = max(d.values())
                keep = sorted(a for a, f in d.items() if f >= top)
            allowed.append("".join(keep))
        out.append(ClusterMotif(
            cid=cid, epitope=meta["antigen.epitope"], gene=meta["gene"], species=meta["species"],
            v=meta["v.segm.repr"], j=meta["j.segm.repr"], length=length, size=size,
            allowed=tuple(allowed)))
    return out


def motif_mass(model, allowed, v=None, j=None) -> float:
    """`Pgen` of every junction matching a degenerate motif.

    ``allowed`` is one entry per position, each a string of permitted residues; ``""`` or ``"X"``
    means any residue. A VDJdb cluster PWM is V/J/length-pinned, so thresholding it per position
    gives exactly this — and one call returns the whole cluster's mass with no enumeration and no
    inclusion–exclusion.

    Unlike :func:`observed_mass` and :func:`union_mass` this scores a **set**, so it carries no
    observed-sample coverage bias. Pass the motif's own ``v``/``j`` — cluster motifs are V/J-pinned
    and the conditioned quantity is the right one for them.
    """
    return float(_native().pgen_aa_degenerate(model, list(allowed), v=v, j=j))


# --------------------------------------------------------------------------- the A-vs-B check
def cross_check(model, junctions, allowed, r: int = 1, alpha: float = ALPHA_PER_EDIT,
                threads: int = 0) -> dict:
    """Two independent estimates of the same `F(e)` — and their disagreement is the missing mass.

    Route **A** (:func:`motif_mass`) scores a *set*: it asks the recombination model for the total
    mass of every junction the motif admits, so it never touches the observed sample and carries
    none of its coverage bias. Route **B** (:func:`observed_mass` / :func:`shell_profile`) starts
    from the junctions VDJdb actually recorded, so it is bounded below by how deeply the epitope was
    sampled. Both are estimates of `F(e)` for the same epitope.

    Returns ``{"set_mass", "observed_mass", "ball_mass", "retained_mass", "ratio_observed",
    "ratio_retained", "missing_fraction", "n_seqs"}``.

    **Interpretation.** ``ratio_observed = set_mass / observed_mass`` is the factor by which the
    sample under-counts, and ``missing_fraction = 1 - observed_mass/set_mass`` the share of the
    cognate mass never observed. ``ratio_retained`` repeats it against the shell-weighted estimate:
    if the neighbourhood correction is doing its job, ``ratio_retained`` is closer to 1 than
    ``ratio_observed`` — that is the whole claim of the ``ball``/``shell`` route, tested rather than
    assumed. A ratio **below 1** is informative in the other direction: the motif is tighter than
    the sample it was built from, i.e. the threshold is too strict or the cluster is a proper subset
    of the epitope's cognate TCRs (which it usually is — one epitope has several clusters).

    Both routes are **marginalised over V/J** so the two numbers are the same quantity. Cluster
    motifs are V/J-pinned and :func:`motif_mass` will condition on request, but a conditioned A
    against a marginal B is a category error and is not offered here.
    """
    set_mass = motif_mass(model, allowed)
    obs = observed_mass(model, junctions, threads=threads)
    prof = shell_profile(model, junctions, r=r, alpha=alpha, threads=threads)
    return {"set_mass": set_mass, "observed_mass": obs, "ball_mass": prof["union"],
            "retained_mass": prof["retained"],
            "ratio_observed": (set_mass / obs) if obs > 0 else float("inf"),
            "ratio_retained": (set_mass / prof["retained"]) if prof["retained"] > 0 else float("inf"),
            "missing_fraction": (1.0 - obs / set_mass) if set_mass > 0 else 0.0,
            "n_seqs": prof["n_seqs"]}
