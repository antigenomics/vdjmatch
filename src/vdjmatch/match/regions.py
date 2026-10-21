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
    """Strip allele suffixes, retaining IMGT gene names such as TRAV23/DV6: ``TRBV19*01`` -> ``TRBV19``."""
    return _ALLELE.sub("", g) if g else ""


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


def annotate_ndn_boundaries(frame, *, species):
    """ARDA half-open boundaries on unchanged junction/V keys, with explicit fallback.

    Supplied J calls must agree for a key. Missing/ambiguous calls and repaired or
    unusable boundaries retain the row with uniform full substitution weight.
    One existing ARDA batch annotates unique supported keys; no sequence is repaired.
    """
    import polars as pl
    try:
        from arda.cdr3fix import markup_batch
    except ImportError as exc:
        raise RuntimeError('NDN annotation requires installed arda.cdr3fix.markup_batch') from exc
    keys=['cdr3','_allele']
    if 'j' not in frame.columns:
        frame=frame.with_columns(pl.lit(None,dtype=pl.String).alias('j'))
    calls=frame.with_columns(pl.col('j').str.strip_chars().str.to_uppercase().alias('_boundary_j'))
    calls=calls.with_columns(pl.when(pl.col('_boundary_j').str.contains(r'[,;+|*]') | (pl.col('_boundary_j')==''))
        .then(pl.col('_boundary_j')).otherwise(pl.col('_boundary_j')+'*01').alias('_boundary_j'))
    catalog=calls.group_by(keys,maintain_order=True).agg(pl.col('_boundary_j').unique().alias('_j_calls'))
    catalog=catalog.with_columns(
        pl.when(pl.col('_j_calls').list.len()==1).then(pl.col('_j_calls').list.first()).alias('_j'),
        pl.when(pl.col('_j_calls').list.len()!=1).then(pl.lit('conflicting_j_uniform'))
        .when(pl.col('_j_calls').list.first().is_null() | (pl.col('_j_calls').list.first()==''))
        .then(pl.lit('missing_call_uniform')).otherwise(pl.lit('pending')).alias('boundary_status'))
    catalog=catalog.with_columns(
        pl.when(pl.col('boundary_status')=='pending').then(
            pl.when(pl.col('_j').str.contains(r'[,;+|]') | pl.col('_allele').str.contains(r'[,;+|]'))
            .then(pl.lit('ambiguous_call_uniform')).otherwise(pl.lit('pending')))
        .otherwise(pl.col('boundary_status')).alias('boundary_status'))
    supported=catalog.filter(pl.col('boundary_status')=='pending')
    if supported.height:
        marks=markup_batch(supported,cdr3='cdr3',v='_allele',j='_j',organism=species,max_replace=0)
        if marks.height!=supported.height:
            raise RuntimeError('ARDA batch changed boundary row count')
        supported=pl.concat([supported,marks.select('cdr3_repaired','v_end','j_start',pl.col('v_call').alias('_markup_v'),pl.col('j_call').alias('_markup_j'))],how='horizontal')
        supported=supported.with_columns(
            pl.when(pl.col('cdr3')!=pl.col('cdr3_repaired')).then(pl.lit('annotation_changed_uniform'))
            .when((pl.col('_markup_v')!=pl.col('_allele')) |
                  (pl.col('_markup_j')!=pl.when(pl.col('_j').str.contains(r'\*')).then(pl.col('_j')).otherwise(pl.col('_j')+'*01')))
            .then(pl.lit('gene_changed_uniform'))
            .when(pl.col('v_end').is_null() | pl.col('j_start').is_null() |
                  (pl.col('v_end')<0) | (pl.col('j_start')<pl.col('v_end')) |
                  (pl.col('j_start')>pl.col('cdr3').str.len_chars()))
            .then(pl.lit('unavailable_uniform')).otherwise(pl.lit('annotated')).alias('boundary_status'))
        catalog=catalog.join(supported.select(*keys,'v_end','j_start',pl.col('boundary_status').alias('_status')),
            on=keys,how='left',validate='1:1',maintain_order='left').with_columns(
                pl.coalesce('_status','boundary_status').alias('boundary_status'))
    else:
        catalog=catalog.with_columns(pl.lit(None,dtype=pl.Int64).alias('v_end'),
                                     pl.lit(None,dtype=pl.Int64).alias('j_start'))
    result=catalog.select(*keys,'v_end','j_start','boundary_status').with_columns(
        *(pl.when(pl.col('boundary_status')=='annotated').then(pl.col(c)).otherwise(None).alias(c)
          for c in ('v_end','j_start')))
    return frame.join(result,on=keys,how='left',validate='m:1',maintain_order='left')


def load_position_profile(path) -> dict[tuple[str, int], list[float]]:
    """Explicit per-position substitution-cost multipliers, keyed by ``(locus, junction length)``.

    TSV with header ``locus length position multiplier``. ``length`` is the **full junction**
    length (conserved anchors included, the frame the search keys use) and ``position`` is 0-based
    in that frame. Every position of a declared ``(locus, length)`` must appear exactly once;
    multipliers are finite and nonnegative, and ``1`` leaves a position's cost unchanged.

    The multiplier scales substitution costs only. It does not scale gaps, and it does not change
    ``geometry.distance_units``, so the search radius, the cutoff ladder and the gap cost stay
    derived from the matrix maximum: a profile whose mean exceeds 1 retrieves fewer neighbours at
    the same radius. Normalise the file itself if an arm must stay radius-comparable.
    """
    rows: dict[tuple[str, int], dict[int, float]] = {}
    with open(path) as fh:
        header = next(fh, '').rstrip('\n').split('\t')
        if header != ['locus', 'length', 'position', 'multiplier']:
            raise ValueError('position profile needs a locus/length/position/multiplier header')
        for line in fh:
            if not line.strip():
                continue
            locus, length, position, multiplier = line.rstrip('\n').split('\t')
            length, position, multiplier = int(length), int(position), float(multiplier)
            if locus not in ('TRA', 'TRB'):
                raise ValueError('position profile locus must be TRA or TRB')
            if length < 1 or not 0 <= position < length:
                raise ValueError('position profile coordinates must lie inside their junction length')
            if not multiplier == multiplier or multiplier in (float('inf'), float('-inf')) or multiplier < 0:
                raise ValueError('position profile multipliers must be finite and nonnegative')
            frame = rows.setdefault((locus, length), {})
            if position in frame:
                raise ValueError(f'duplicate position profile row for {locus} length {length}')
            frame[position] = multiplier
    if not rows:
        raise ValueError('position profile is empty')
    profile = {}
    for (locus, length), frame in rows.items():
        if len(frame) != length:
            raise ValueError(f'position profile for {locus} length {length} covers {len(frame)}/{length} positions')
        profile[(locus, length)] = [frame[i] for i in range(length)]
    return profile


def profile_weights(profile, locus, lengths, ctrim, scale=100) -> dict[int, list[int]]:
    """Native ``position_weights_by_length`` table for one locus, from a loaded profile.

    The flanks the geometry trims (first 3, last ``ctrim``) carry weight 0 whatever the profile
    says there, exactly as the uniform table does, so an all-ones profile reproduces the
    unmodified weights bit-for-bit. Interior weights are ``round(scale*multiplier)``; unlike the
    significance route they are not floored to 1, so a 0 multiplier genuinely ignores a column.
    """
    missing = sorted(length for length in lengths if (locus, length) not in profile)
    if missing:
        raise ValueError(f'position profile does not cover {locus} junction lengths {missing}')
    weights = {}
    for length in lengths:
        row = profile[(locus, length)]
        weights[length] = [0]*3 + [round(scale*m) for m in row[3:length-ctrim]] + [0]*ctrim
    return weights


def position_profile_summary(profile) -> dict:
    """Per-locus multiplier range and covered lengths, for provenance in a run manifest."""
    summary = {}
    for locus in sorted({locus for locus, _ in profile}):
        values = [m for (l, _), row in profile.items() if l == locus for m in row]
        summary[locus] = dict(lengths=sorted(length for l, length in profile if l == locus),
                              positions=len(values), min=min(values), max=max(values),
                              mean=sum(values)/len(values))
    return summary
