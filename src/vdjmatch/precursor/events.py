"""`F(e)` counted off repertoire data: independent recombination events over independent events.

No `Pgen`, no model, no coverage correction. This is the estimand itself, which makes it the check
everything in :mod:`vdjmatch.precursor.mass` is measured against.
"""
from __future__ import annotations

from dataclasses import dataclass

from .mass import MAX_BALL_MEMBERS

_NT = frozenset("ACGTNacgtn")


@dataclass(frozen=True)
class RecombinationEvent:
    """One **independent recombination event**: the key is ``(donor, v, j, junction_nt)``.

    Not a sequence, not a clonotype, not a person. The same ``junction_nt`` observed in two donors
    is **two** events — two rearrangements that happened to converge on the same nucleotide string —
    and must never be deduplicated across donors. The same nucleotide junction rearranged onto a
    different V is likewise a separate event.

    ``junction_aa`` is carried only for *matching* against a cognate set; it is deliberately not
    part of the key, because convergent recombination reaches one amino-acid junction by many
    nucleotide paths and collapsing them would undercount exactly the events this estimator exists
    to count.

    Construction validates the pair so a mis-typed call fails loudly: ``junction_nt`` must be
    ACGTN-only and exactly ``3x`` the length of ``junction_aa``, and ``donor`` must be non-empty.
    Passing an amino-acid junction as the nucleotide key, or pooling donors by leaving ``donor``
    blank, are the two mistakes that silently produce a wrong answer rather than an error.
    """
    donor: str
    v: str
    j: str
    junction_nt: str
    junction_aa: str

    def __post_init__(self):
        if not str(self.donor).strip():
            raise ValueError(
                "RecombinationEvent needs a non-empty donor: two donors carrying the same "
                "junction are two independent events, so pooling them undercounts the numerator "
                "and the denominator alike")
        if not self.junction_nt or not set(self.junction_nt) <= _NT:
            raise ValueError(
                f"junction_nt must be a nucleotide string (ACGTN), got {self.junction_nt!r} -- "
                "the event key is nucleotide-level; an amino-acid junction here collapses "
                "convergent rearrangements into one event")
        if len(self.junction_nt) != 3 * len(self.junction_aa):
            raise ValueError(
                f"junction_nt must be exactly 3x junction_aa in length, got "
                f"{len(self.junction_nt)} vs 3x{len(self.junction_aa)}")

    @property
    def key(self) -> tuple[str, str, str, str]:
        return (self.donor, self.v, self.j, self.junction_nt)


def event_ratio(cognate, events, r: int = 1, denominator: int | None = None,
                max_members: int = MAX_BALL_MEMBERS) -> dict:
    """`F(e)` counted directly off repertoire data — events over events, no `Pgen` at all.

    ::

        F_hat(e) = #{distinct (donor, V, J, junction_nt) matching C_e within r mismatches}
                   ----------------------------------------------------------------------
                   #{distinct (donor, V, J, junction_nt) in the whole dataset}

    This **is** the estimand — the probability that a random naive-repertoire rearrangement
    recognises `e` — not a proxy for it, so it is a direct empirical check on the `Pgen` route
    rather than a correlate of it.

    **Two samplings, and only one of them divides out.** *Repertoire* depth does: numerator and
    denominator are counted with the same key over the same donors, so sequencing more deeply
    changes both and nothing needs coverage-correcting. That is the advantage over
    :func:`coverage_corrected_mass`, which needs recaptures and is undefined when every junction is
    a singleton. *Database* depth does **not**: ``cognate`` is whatever VDJdb happens to hold for
    the epitope, and a bigger set mechanically matches more events. Measured on 138 epitopes, the
    set total ``f_hat`` is predicted at out-of-fold R^2 ~0.5 by log cognate-set size alone —
    :func:`observed_mass` and :func:`union_mass` carry exactly the same defect, for exactly the same
    reason, since all three are set totals.

    **The size-invariant form is a median over junctions, and it needs no new function**: pass a
    mapping of singletons, ``event_ratio({tau: [tau] for tau in cognate}, events)``, and take the
    median of the per-junction ``f_hat``. That estimates the per-junction event rate, which is a
    property of the epitope rather than of how many of its TCRs were catalogued, keeps the
    events-over-events property exactly, and throws away no data. It is also the empirical
    counterpart of the *median* rearrangement probability that Pogorelyy et al. correlated — the
    quantity the replication gate is built on — which is why the two agree. Prefer it to dividing
    the set total by ``len(cognate)``: cognate junctions are near-duplicates, their one-substitution
    balls overlap, and dividing a union by a count under-corrects tight groups and over-corrects
    diverse ones.

    ``cognate`` is either one epitope's junction list or a ``{epitope: junctions}`` mapping; the
    mapping form streams ``events`` once and shares the denominator across epitopes. ``events`` is
    an iterable of :class:`RecombinationEvent`. ``r`` is the match radius in **amino-acid**
    substitutions. ``denominator`` overrides the counted total, for callers that pre-filtered the
    event stream; it must have been counted with the same key.

    Returns ``{"denominator", "matched", "f_hat", "epitopes": {e: {"matched", "f_hat"}},
    "n_cognate", "r"}``.

    **`r = 1` is part of the estimator, not a tuning knob.** At `r = 0` the repertoire is far too
    sparse for this ratio to be estimable: a cognate set of a few hundred exact junctions matches
    almost nothing in a bulk repertoire, and the counts are then dominated by whether one public
    clonotype happened to be sequenced. Measured on the Pogorelyy replication the same correlation
    is 0.51–0.61 at `r = 0` and 0.76–0.86 at `r = 1`. Both this estimator and the `Pgen` route
    (``pgen(..., mismatches=1)``) therefore take the closed one-substitution ball as their unit.

    **Counting objects, not recaptures.** :func:`coverage_corrected_mass` takes a ``multiplicity``
    that means "how many independent units *re*-reported this same junction" — recaptures of one
    object, used to estimate what was missed. Here distinct donors are **distinct objects**, each
    contributing its own event. The two are not the same quantity and the same donor counts must
    not be fed to both: a donor count passed as ``multiplicity`` says "seen again", while the same
    number here says "happened again".
    """
    from seqtree.distance import neighbourhood_union, union_size

    single = not isinstance(cognate, dict)
    groups = {"": list(cognate)} if single else {k: list(v) for k, v in cognate.items()}
    groups = {k: [s for s in v if s] for k, v in groups.items()}

    total = sum(union_size(v, r=r) for v in groups.values() if v)
    if total > max_members:
        raise MemoryError(
            f"the radius-{r} neighbourhood of {sum(len(v) for v in groups.values()):,} cognate "
            f"junctions holds {total:,} sequences, above max_members={max_members:,}. Split the "
            "epitopes into batches over the same event stream and combine the counts -- the "
            "denominator is shared, so that is exact.")

    index: dict[str, set[str]] = {}
    for name, seqs in groups.items():
        if not seqs:
            continue
        for nb in neighbourhood_union(seqs, r=r):
            index.setdefault(nb, set()).add(name)

    seen: set[tuple[str, str, str, str]] = set()
    hits: dict[str, set[tuple[str, str, str, str]]] = {k: set() for k in groups}
    for e in events:
        k = e.key
        if k in seen:
            continue
        seen.add(k)
        for name in index.get(e.junction_aa, ()):
            hits[name].add(k)

    n = int(denominator) if denominator is not None else len(seen)
    per = {name: {"matched": len(v), "f_hat": (len(v) / n) if n else 0.0}
           for name, v in hits.items()}
    out = {"denominator": n, "epitopes": per, "r": r,
           "n_cognate": sum(len(v) for v in groups.values())}
    if single:
        out["matched"] = per[""]["matched"]
        out["f_hat"] = per[""]["f_hat"]
    else:
        out["matched"] = len(set().union(*hits.values())) if hits else 0
        out["f_hat"] = (out["matched"] / n) if n else 0.0
    return out
