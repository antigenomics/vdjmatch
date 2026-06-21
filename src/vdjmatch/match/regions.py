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


def load_significance(path=None) -> list[tuple[float, float]]:
    """Empirical positional informativeness factor (relpos, weight): how much a substitution at a
    relative CDR3 position changes specificity, from ``bench/scoring_analysis.py`` (weight = 1 −
    P(neighbour shares epitope | substitution here), normalised to mean 1; centre > V/J borders)."""
    src = path or (resources.files("vdjmatch.resources") / "trimming" / "position_significance.tsv")
    out = []
    with open(src) as fh:
        next(fh)
        for line in fh:
            relpos, _p, w = line.rstrip("\n").split("\t")
            out.append((float(relpos), float(w)))
    return out


def _interp(x: float, xs: list[float], ys: list[float]) -> float:
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for k in range(1, len(xs)):
        if x <= xs[k]:
            t = (x - xs[k - 1]) / (xs[k] - xs[k - 1])
            return ys[k - 1] + t * (ys[k] - ys[k - 1])
    return ys[-1]


_SIG: list[tuple[float, float]] | None = None


def significance_weights(length: int, sig=None) -> list[float]:
    """Per-position informativeness weight for a length-``length`` CDR3 (centre upweighted), by
    interpolating the empirical profile over relative position. Gene-agnostic (cf. the gene-specific
    germline-retention weights of :func:`position_weights`)."""
    global _SIG
    if length < 2:
        return [1.0] * length
    prof = sig if sig is not None else _SIG
    if prof is None:
        prof = _SIG = load_significance()
    xs = [r for r, _ in prof]
    ws = [w for _, w in prof]
    return [_interp(i / (length - 1), xs, ws) for i in range(length)]


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
