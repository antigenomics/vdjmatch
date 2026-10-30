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

    This legacy interface retains its explicit OLGA default for reproduction.
    vdjtools owns the installed catalogue and parameters. Canonical ``isalgo``
    contains current gene-level fits; ``learned`` and ``arda`` are aliases of
    these current weights, not historical trained models. ``olga`` aliases
    frozen human ``olga-legacy`` models. Mouse isalgo supports TRA/TRB/IGH/IGK/IGL;
    missing contexts fail explicitly. Historical learned/ARDA results require
    their original pinned parameters and calibration identities.
    """
    from vdjtools.model import load_bundled
    return load_bundled(locus, source=source, organism=organism)


def bundled_model_path(locus: str, source: str = "olga", organism: str = "human"):
    """Locate installed model files for provenance using the public catalogue.

    Loading remains owned by ``vdjtools.model.load_bundled``. This resolver is
    only for hashing the same immutable files, never for selecting old weights.
    """
    from importlib.resources import files
    from vdjtools.model import list_bundled

    family = {"olga": "olga-legacy", "learned": "isalgo", "arda": "isalgo"}.get(source, source)
    catalogue = list_bundled()
    if family not in catalogue or organism.lower() not in ("human", "mouse"):
        raise ValueError(f"unsupported model source/organism: {source}/{organism}")
    key = locus.upper() if organism.lower() == "human" else "mouse_" + locus.upper()
    if key not in catalogue[family]:
        raise FileNotFoundError(f"no bundled {family} model for {organism}/{locus}")
    directory = "olga" if family == "olga-legacy" else family
    return files("vdjtools.model") / "_bundled" / directory / key


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
