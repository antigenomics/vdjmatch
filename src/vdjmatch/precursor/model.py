"""The recombination model, junction QC and per-junction `Pgen` — the vdjtools bridge.

Nothing here reimplements `Pgen`. vdjtools ships a C++ transfer-matrix DP with a codon-mask
parameterisation; this module is the thin, guarded surface vdjmatch calls it through.
"""
from __future__ import annotations

_JUNCTION_START = "C"
_JUNCTION_END = ("F", "W")


def _native():
    """The vdjtools Pgen bindings, or a clear error naming the extra."""
    try:
        from vdjtools.model import native
    except ImportError as e:                     # pragma: no cover - depends on the environment
        raise ImportError(
            "vdjmatch.precursor needs the optional 'vdjtools' dependency for the recombination "
            "model (Pgen). Install it with: pip install 'vdjmatch[precursor]'"
        ) from e
    return native


def load_model(locus: str = "TRB", source: str = "olga", organism: str = "human"):
    """Bundled recombination model.

    vdjtools ships three sets, and they differ in ways that matter for a mass over a *set* of
    junctions, where a junction scoring exactly zero is silently dropped from the total:

    ``"olga"``
        A bit-faithful import of OLGA's published models — native `Pgen` matches OLGA's own to
        machine precision. Human, seven loci. Faithful includes faithful to OLGA's defects: its
        per-locus deletion grid puts mass on trims an allele is too short to reach, and the DP
        never visits those, so `Pgen` through the affected alleles is an underestimate or exactly
        zero. On VDJdb this costs **5.7% of human TRA junctions** (3,000 of 52,363) and 41 of
        111,331 on TRB.
    ``"learned"``
        Refit from real 5'RACE reads on arda germline, so it does not inherit that grid. Human,
        seven loci. Loses 270 TRA junctions (0.5%) and 2 on TRB.
    ``"arda"``
        The same refit on the arda IMGT allele namespace, and the **only** set carrying a non-human
        organism: human for seven loci plus mouse TRA/TRB. Loses 19 TRA junctions and 2 on TRB;
        neither mouse locus loses any.

    ``"olga"`` remains the default because exact agreement with the reference implementation is the
    right property for a published null. **For a mass over a set — which is what everything in this
    package computes — prefer ``"learned"`` or ``"arda"``**, and for anything comparing human with
    mouse use ``"arda"`` on both so the model family is held fixed.

    Sparse gene usage is not a proxy for this: mouse ``arda`` TRA has `P(V) = 0` for 60% of its V
    genes and loses *no* junctions, because marginalising over V and J routes around genes the fit
    never saw, while ``olga`` TRA has one zero-usage V gene and loses 5.7%.
    """
    from vdjtools.model import load_bundled
    return load_bundled(locus, source=source, organism=organism)


def check_junctions(seqs):
    """Split ``seqs`` into (junctions, suspect) by the conserved anchors.

    **CDR3 is not junction.** VDJdb's column is *named* ``cdr3`` but holds junctions — Cys104 and
    Phe118/Trp included. An anchor-stripped IMGT CDR3 scores **exactly 0.0** with no error, so a
    silently mis-typed input reports a precursor frequency of zero rather than failing. Callers
    should check before scoring and report the dropped count rather than letting it vanish.
    """
    ok, bad = [], []
    for s in seqs:
        t = (s or "").strip().upper()
        (ok if t.startswith(_JUNCTION_START) and t.endswith(_JUNCTION_END) else bad).append(t)
    return ok, bad


def pgen(model, junctions, v=None, j=None, mismatches: int = 0, threads: int = 0) -> list[float]:
    """Per-junction `Pgen` — the vector behind every mass in this package.

    ``mismatches=1`` returns the **closed Hamming-1 ball** mass of each junction in closed form
    (inclusion–exclusion inside vdjtools, no enumeration), which is the frequency proxy Pogorelyy
    et al. used. ``v``/``j`` are per-junction allele-resolution call lists or ``None`` to
    marginalise; marginal and conditioned are different quantities, so never mix them in one
    comparison.
    """
    seqs = list(junctions)
    if not seqs:
        return []
    return [float(x) for x in _native().pgen_aa_batch(
        model, seqs, v=v, j=j, mismatches=mismatches, threads=threads)]
