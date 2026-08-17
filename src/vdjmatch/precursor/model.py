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

    ``source="olga"`` is the default for null work. Do **not** use ``"learned"``: those models were
    EM-fit on ~2k clonotypes without a gene-usage pseudocount, so 68 of 89 bundled TRB V alleles
    have `P(V) = 0` and any junction using them scores 0. Mouse is available only under
    ``source="arda"``, and only for TRA/TRB.
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
