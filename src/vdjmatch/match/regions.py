"""Region-aware CDR3 scoring: weight substitutions by where they fall in the rearrangement.

The germline-encoded V/J flanks of a CDR3 are near-invariant and their residues are determined by
gene choice, not antigen-driven selection; the non-template (NDN) core is where specificity-relevant
variation lives (see the scoring appendix, and Shcherbinin et al. 2023 on contacting regions). We
therefore weight each CDR3 position by ``1 - germline_retention`` — NDN positions get full weight,
germline flanks little — and score substitutions with the NDN-derived VDJAM matrix. The per-gene
germline-retention profiles are precomputed from the OLGA recombination model by ``mirpy``
(``mir.basic.trimming``) and shipped as ``resources/trimming/human_vj_retention.tsv``.
"""
from __future__ import annotations

import re
from importlib import resources

_ALLELE = re.compile(r"\*.*$")


def gene_family(g: str | None) -> str:
    """Strip the allele / IMGT decoration: ``TRBV19*01`` -> ``TRBV19``."""
    return _ALLELE.sub("", g.split("/")[0]) if g else ""


def load_retention(path=None) -> dict[tuple[str, str, str], list[float]]:
    """Load germline-retention profiles, keyed by ``(chain, segment, gene_family)`` ->
    ``[p_retain per offset]`` (V from the N-anchor, J from the C-anchor). First allele per family."""
    src = path or (resources.files("vdjmatch.resources") / "trimming" / "human_vj_retention.tsv")
    out: dict[tuple[str, str, str], list[float]] = {}
    with open(src) as fh:
        next(fh)
        for line in fh:
            chain, seg, gene, off, p, _aa = line.rstrip("\n").split("\t")
            key = (chain, seg, gene_family(gene))
            out.setdefault(key, []).append(float(p))
    return out


def position_weights(length: int, v: str | None, j: str | None, chain: str,
                     ret: dict) -> list[float]:
    """Per-CDR3-position substitution weight ``= 1 - max(V-side, J-side) germline retention``.
    NDN core -> ~1; germline V/J flanks -> ~0. Unknown genes contribute 0 retention (full weight)."""
    rv = ret.get((chain, "V", gene_family(v)), [])
    rj = ret.get((chain, "J", gene_family(j)), [])
    w = []
    for i in range(length):
        a = rv[i] if i < len(rv) else 0.0
        k = length - 1 - i
        b = rj[k] if k < len(rj) else 0.0
        w.append(1.0 - max(a, b))
    return w


def load_significance(path=None) -> tuple[list[float], list[float]]:
    """End-anchored positional informativeness profiles ``(p_same_V, p_same_J)`` indexed by offset
    from the V- and J-anchors, from ``bench/scoring_analysis.py``. ``p_same[d]`` = P(neighbour shares
    epitope | a single substitution at offset ``d`` from that anchor), Beta-Binomial-smoothed; offsets
    ``>=`` the last index are pooled (the NDN core). High near the conserved anchors, low in the core."""
    src = path or (resources.files("vdjmatch.resources") / "trimming" / "position_significance.tsv")
    v: dict[int, float] = {}
    j: dict[int, float] = {}
    with open(src) as fh:
        next(fh)
        for line in fh:
            side, off, _n, p = line.rstrip("\n").split("\t")
            (v if side == "V" else j)[int(off)] = float(p)
    pv = [v[d] for d in range(len(v))]
    pj = [j[d] for d in range(len(j))]
    return pv, pj


_SIG: tuple[list[float], list[float]] | None = None


def significance_weights(length: int, sig=None) -> list[float]:
    """Per-position informativeness weight for a length-``length`` CDR3, **anchored from both ends**:
    a position is uninformative if it sits near *either* the germline V or J anchor (fixed absolute
    offsets), informative in the NDN core. ``omega(p) = 1 - max(p_same_V[dV], p_same_J[dJ])`` with
    ``dV=p``, ``dJ=L-1-p`` (capped at the profile length), normalised to mean 1. Gene-agnostic (cf. the
    gene-specific germline-retention weights of :func:`position_weights`)."""
    global _SIG
    if length < 2:
        return [1.0] * length
    pv, pj = sig if sig is not None else _SIG or load_significance()
    if sig is None and _SIG is None:
        _SIG = (pv, pj)
    raw = []
    for p in range(length):
        dv = min(p, len(pv) - 1)
        dj = min(length - 1 - p, len(pj) - 1)
        raw.append(1.0 - max(pv[dv], pj[dj]))
    mean = sum(raw) / len(raw) or 1.0
    return [w / mean for w in raw]


def significance_pssm(length: int, base: str = "blosum62", scale: int = 100):
    """Experiment-(2) as a native seqtree PSSM: a fixed-width ``PositionalMatrix`` that scales the
    base substitution matrix by the per-position informativeness factor (centre > V/J borders), so
    the engine itself up-weights central mismatches. ``pen(pos,a,b) = weight[pos] · base(a,b)``."""
    from seqtree import PositionalMatrix, SubstitutionMatrix
    b = getattr(SubstitutionMatrix, base)()
    w = [max(1, round(scale * x)) for x in significance_weights(length)]
    return PositionalMatrix.from_weights(b, w)


def vdjam_penalties(path=None) -> dict[tuple[str, str], float]:
    """Squared-distance penalty ``pen(a,b)=s_aa+s_bb-2*s_ab`` from the VDJAM similarity table
    (seqtree's ``from_similarity`` convention), computed in Python for region-aware rescoring."""
    src = path or (resources.files("vdjmatch.resources") / "vdjam.txt")
    sim: dict[tuple[str, str], float] = {}
    with open(src) as fh:
        next(fh)
        for line in fh:
            a, b, s = line.split()
            sim[(a, b)] = float(s)
    aas = sorted({a for a, _ in sim})
    return {(a, b): sim[(a, a)] + sim[(b, b)] - 2 * sim[(a, b)] for a in aas for b in aas}


def substitution_score(qseq: str, rseq: str, v: str | None, j: str | None, chain: str,
                       ret: dict, pen: dict) -> float:
    """Region-weighted substitution score between two equal-length CDR3s (the no-indel case, e.g.
    scope ``s,0,0,s``). Sum over mismatched positions of ``weight(pos) * VDJAM_penalty``; lower =
    more similar. For gapped alignments use :func:`aligned_score`."""
    w = position_weights(len(qseq), v, j, chain, ret)
    total = 0.0
    for i, (a, b) in enumerate(zip(qseq, rseq)):
        if a != b:
            total += w[i] * pen.get((a, b), 0.0)
    return total


def aligned_score(aligned_query: str, aligned_ref: str, ops: str, v: str | None, j: str | None,
                  chain: str, ret: dict, pen: dict, gap: float = 1.0) -> float:
    """Region-weighted score from a seqtree alignment (gapped query/ref strings + per-column ops
    M/S/I/D). Substitutions are weighted by the query position's NDN weight and the VDJAM penalty;
    indels carry a flat ``gap`` penalty. Lower = more similar."""
    w = position_weights(len(aligned_query.replace("-", "")), v, j, chain, ret)
    qpos, total = 0, 0.0
    for q, r, op in zip(aligned_query, aligned_ref, ops):
        if op == "M":
            qpos += 1
        elif op == "S":
            total += (w[qpos] if qpos < len(w) else 1.0) * pen.get((q, r), 0.0)
            qpos += 1
        elif op == "I":          # insertion in query: consumes a query residue
            total += gap; qpos += 1
        elif op == "D":          # deletion from query: consumes a ref residue only
            total += gap
    return total
