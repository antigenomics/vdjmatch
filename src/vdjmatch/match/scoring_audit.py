"""Accepted legacy matrix/LOO calculations, migrated from bench; not annotation confidence.

Legacy contracts deliberately retain unique-junction eligibility, composition caps,
first-epitope assignment and balanced trapezoidal PR area rather than modern AP.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict

import numpy as np
import polars as pl
from seqtree import Index, SearchParams, SubstitutionMatrix

from . import regions

AA = "ACDEFGHIKLMNPQRSTVWY"
VALID = "^[ACDEFGHIKLMNPQRSTVWY]+$"



def valid_cdr3(df: pl.DataFrame) -> pl.DataFrame:
    """Drop rows whose CDR3 has non-standard residues (else ``Index.build`` rejects them)."""
    return df.filter(pl.col("cdr3").str.contains(VALID))

def spike_studies(df: pl.DataFrame, spike_min: int = 3000, spike_modal: float = 0.5) -> set[str]:
    """References that single-handedly create a spectratype spike: one study dumping >= ``spike_min``
    clonotypes into a single epitope with one CDR3 length dominating (modal fraction >= ``spike_modal``)
    -- e.g. PMID:40498839, which contributes ~30k SLLMWITQV clonotypes (68% at length 14) and reports
    no other epitope. We drop the offending *study's records* (not the epitope), so the epitope keeps
    its clonotypes from other studies and no other epitope is touched."""
    u = (valid_cdr3(df).select("reference_id", "epitope", "cdr3").unique()
         .with_columns(pl.col("cdr3").str.len_chars().alias("L")))
    per_len = u.group_by(["reference_id", "epitope", "L"]).len()
    tot = u.group_by(["reference_id", "epitope"]).len().rename({"len": "n"})
    top = per_len.group_by(["reference_id", "epitope"]).agg(pl.col("len").max().alias("top"))
    m = tot.join(top, on=["reference_id", "epitope"]).with_columns(
        (pl.col("top") / pl.col("n")).alias("modal"))
    return set(m.filter((pl.col("n") >= spike_min) & (pl.col("modal") >= spike_modal))
               ["reference_id"].to_list())

def long_list(df: pl.DataFrame, cap: int = 3000, min_n: int = 30, seed: int = 0,
              drop_spikes: bool = True) -> pl.DataFrame:
    """Composition-controlled long list of unique clonotypes ``(cdr3,v,j,epitope,mhc_class)``:
    drop the records of spike studies (one study dumping a length-skewed mega-set into one epitope;
    :func:`spike_studies`), keep epitopes with >= ``min_n`` clonotypes, and cap each at a random ``cap``
    (so legitimate multi-study mega-epitopes, e.g. NLV ~13k, don't dominate sampling or purity).
    Deterministic via ``seed``. The full database stays available (don't pass through here) for
    sampling-depth / motif-coverage studies."""
    d = valid_cdr3(df)
    if drop_spikes:
        d = d.filter(~pl.col("reference_id").is_in(list(spike_studies(df))))
    u = d.select("cdr3", "v", "j", "epitope", "mhc_class").unique()
    keep = u.group_by("epitope").len().filter(pl.col("len") >= min_n)["epitope"]
    # sort to a canonical order before the seeded shuffle: polars unique() does not preserve row order,
    # so without this the cap-`cap` draw varies run-to-run even at a fixed seed.
    u = u.filter(pl.col("epitope").is_in(keep.to_list())).sort(["epitope", "cdr3", "v", "j"])
    rank = pl.int_range(pl.len()).shuffle(seed=seed).over("epitope")
    return u.with_columns(rank.alias("_r")).filter(pl.col("_r") < cap).drop("_r")

def position_background(cdr3s: list[str]) -> dict[tuple[int, int], dict[str, float]]:
    """bg[(L, pos)][aa] = residue frequency at position pos among length-L CDR3s."""
    counts: dict[tuple[int, int], Counter] = defaultdict(Counter)
    for s in cdr3s:
        L = len(s)
        for p, a in enumerate(s):
            if a in AA:
                counts[(L, p)][a] += 1
    return {k: {a: c / sum(v.values()) for a, c in v.items()} for k, v in counts.items()}

def estimate_matrix(events, bg, alpha=1.0) -> dict[tuple[str, str], float]:
    """Position-debiased substitution log-odds from events (a,b,pos,L)."""
    f = Counter()
    slots = Counter()  # (L,pos) -> n events there
    for a, b, p, L in events:
        key = tuple(sorted((a, b)))
        f[key] += 1
        slots[(L, p)] += 1
    e = defaultdict(float)
    for (L, p), n in slots.items():
        freq = bg.get((L, p), {})
        items = list(freq.items())
        for ai in range(len(items)):
            a, fa = items[ai]
            for bj in range(ai, len(items)):
                b, fb = items[bj]
                key = tuple(sorted((a, b)))
                e[key] += n * (fa * fb if a == b else 2 * fa * fb)
    score = {}
    for a in AA:
        for b in AA:
            if a < b or a == b:
                key = (a, b)
                score[key] = math.log2((f.get(key, 0) + alpha) / (e.get(key, 0.0) + alpha))
    return score

def corr(m1: dict, m2: dict) -> float:
    """Pearson correlation of off-diagonal substitution scores shared by both matrices."""
    xs, ys = [], []
    for a in AA:
        for b in AA:
            if a < b:
                k = (a, b)
                if k in m1 and k in m2:
                    xs.append(m1[k]); ys.append(m2[k])
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return cov / (sx * sy) if sx and sy else float("nan")

def legacy_vdjam_offdiag(path) -> dict[tuple[str, str], float]:
    out = {}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            a, bb, s = ln.split()
            if a < bb and a in AA and bb in AA:
                out[(a, bb)] = float(s)
    return out

def region_of(pos: int, L: int, vlen: int, jlen: int) -> str:
    if pos < vlen:
        return "V"
    if pos >= L - jlen:
        return "J"
    return "NDN"

def dissim_from_scores(scores: dict[tuple[str, str], float]) -> dict[tuple[str, str], float]:
    """similarity log-odds -> non-negative dissimilarity (most-similar pair = 0)."""
    off = {k: v for k, v in scores.items() if k[0] != k[1]}
    hi = max(off.values())
    return {k: hi - v for k, v in off.items()}

_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

def codon_dissim() -> dict[tuple[str, str], float]:
    """VDJAMr: a substitution matrix from the genetic code alone -- how mutationally accessible one
    amino acid is from another. M(a,b) = P(a single random nucleotide substitution of a uniformly
    random codon of a yields a codon of b); symmetrised, then converted to off-diagonal dissimilarity
    (codon-adjacent pairs = cheap). This is the recombination/mutation null: it scores the substitution
    propensity built into the code, with no chemistry or selection."""
    from collections import defaultdict
    codons = defaultdict(list)
    for c, a in _CODE.items():
        if a != "*":
            codons[a].append(c)
    nts = "ACGT"
    m = {}
    for a in AA:
        for b in AA:
            tot = hit = 0
            for c in codons[a]:
                for i in range(3):
                    for nt in nts:
                        if nt == c[i]:
                            continue
                        tot += 1
                        hit += _CODE.get(c[:i] + nt + c[i + 1:]) == b
            m[(a, b)] = hit / tot if tot else 0.0
    sim = {(a, b): (m[(a, b)] + m[(b, a)]) / 2 for a in AA for b in AA}
    off = {k: v for k, v in sim.items() if k[0] != k[1]}
    hi = max(off.values())
    return {k: hi - v for k, v in off.items()}

def structural_dissim(inc=None) -> dict | None:
    """seqtree's TeXshade structural similarity matrix (sidechain volume + hydropathy) -> dissim.
    Parsed from the seqtree source .inc (24-symbol order ARNDCQEGHILKMFPSTWYVBZX*); None if absent.
    Path from arg or $STRUCTURAL_INC (e.g. <seqtree>/src/structural.inc)."""
    import re
    from pathlib import Path
    inc = inc or ""
    if not inc or not Path(inc).exists():
        return None
    order = "ARNDCQEGHILKMFPSTWYVBZX*"
    text = Path(inc).read_text()
    marker = "kStructural[24 * 24] = {"
    if marker not in text:
        raise ValueError("structural-inc lacks kStructural[24 * 24]")
    txt = text.split(marker, 1)[1].split("};", 1)[0]
    nums = [int(x) for x in re.findall(r"-?\d+", txt)]
    if len(nums) != 24 * 24:
        raise ValueError("structural-inc must contain exactly 576 matrix entries")
    sim = {(order[i], order[j]): nums[i * 24 + j] for i in range(24) for j in range(24)}
    vals = {(a, c): sim[(a, c)] for a in AA for c in AA if a != c}
    hi = max(vals.values())
    return {k: hi - v for k, v in vals.items()}

def _curve(pairs):
    """Cumulative (tp, fp) points after each tied-score GROUP (descending score), plus P, N. Grouping
    equal scores gives them 0.5 credit (Mann-Whitney-consistent) — without it, tied positives that
    happen to precede tied negatives in the input inflate every AUC."""
    s = sorted(pairs, key=lambda x: -x[1])
    P = sum(l for l, _ in s)
    N = len(s) - P
    pts = [(0, 0)]
    tp = fp = i = 0
    while i < len(s):
        sc = s[i][1]
        j = i
        while j < len(s) and s[j][1] == sc:
            tp += s[j][0]
            fp += 1 - s[j][0]
            j += 1
        pts.append((tp, fp))
        i = j
    return pts, P, N

def pr_auc_balanced(pairs, pi0: float = 0.5) -> float:
    """PR-AUC with precision re-normalized to reference prevalence ``pi0`` (balanced by default)."""
    pts, P, N = _curve(pairs)
    if P == 0 or N == 0:
        return float("nan")
    pr, pp, area = 0.0, 1.0, 0.0
    for tp, fp in pts[1:]:
        tpr, fpr = tp / P, fp / N
        denom = tpr * pi0 + fpr * (1 - pi0)
        p0 = (tpr * pi0) / denom if denom else 1.0
        area += (tpr - pr) * (p0 + pp) / 2
        pr, pp = tpr, p0
    return area


def named_scores(name):
    """Raw native similarities, preserving the legacy max-offdiagonal transform."""
    factory = getattr(SubstitutionMatrix, name.lower(), None)
    if factory is None:
        raise RuntimeError(f'seqtree lacks the native {name} substitution matrix')
    matrix = factory()
    if not callable(getattr(matrix, 'similarity', None)):
        raise RuntimeError('seqtree requires raw SubstitutionMatrix.similarity capability')
    return {(a, b): float(matrix.similarity(a, b)) for a in AA for b in AA}


def named_dissim(name):
    values = {k: v for k, v in named_scores(name).items() if k[0] != k[1]}
    hi = max(values.values())
    return {k: hi-v for k, v in values.items()}


def substitution_events(frame, *, threads=1):
    """One locus-wide Hamming-1 native batch, then same-epitope incidence joins.

    An unordered sequence pair contributes once per shared eligible epitope,
    exactly as the old separate epitope searches did.
    """
    membership = frame.select('cdr3', 'epitope').unique().sort('cdr3', 'epitope')
    seqs = membership['cdr3'].unique().sort().to_list()
    schema = {'epitope': pl.String, 'a': pl.String, 'b': pl.String, 'pos': pl.Int64, 'L': pl.Int64}
    if not seqs:
        return pl.DataFrame(schema=schema)
    result = Index.build(seqs, 'aa').search_batch(seqs,
        SearchParams(max_subs=1, max_total_edits=1, engine='seqtm'), threads)
    pairs = [(i, h.ref_id) for i, hits in enumerate(result) for h in hits
             if h.ref_id > i and h.n_subs == 1 and not h.n_ins and not h.n_dels]
    if not pairs:
        return pl.DataFrame(schema=schema)
    ids = np.asarray(pairs, dtype=np.int64)
    lengths = np.asarray([len(s) for s in seqs])
    parts = []
    for length in np.unique(lengths[ids[:, 0]]):
        pair = ids[lengths[ids[:, 0]] == length]
        local_ids = np.flatnonzero(lengths == length)
        lookup = np.full(len(seqs), -1, dtype=np.int64)
        lookup[local_ids] = np.arange(len(local_ids))
        encoded = np.frombuffer(''.join(seqs[i] for i in local_ids).encode('ascii'), dtype=np.uint8).reshape(-1, length)
        left, right = encoded[lookup[pair[:, 0]]], encoded[lookup[pair[:, 1]]]
        differences = left != right
        if not np.all(differences.sum(axis=1) == 1):
            raise RuntimeError('native substitution-only candidates violate Hamming-1 contract')
        positions = differences.argmax(axis=1)
        offsets = np.arange(len(pair))
        part = pl.DataFrame({'q': [seqs[i] for i in pair[:, 0]], 'r': [seqs[i] for i in pair[:, 1]],
            'a': left[offsets, positions].view('S1').astype(str),
            'b': right[offsets, positions].view('S1').astype(str),
            'pos': positions, 'L': np.full(len(pair), length)})
        parts.append(part)
    events = pl.concat(parts).join(membership.rename({'cdr3': 'q'}), on='q', how='inner', validate='m:m')
    return events.join(membership.rename({'cdr3': 'r'}), on=['r', 'epitope'], how='inner', validate='m:m').select(*schema)


def matrix_audit(frame, *, min_epi=30, vlen=4, jlen=6, legacy, threads=1):
    """Flat positional-matrix scores and offdiagonal Pearson correlations."""
    frame = valid_cdr3(frame)
    sizes = frame.group_by('epitope').agg(pl.col('cdr3').n_unique().alias('n'))
    selected = frame.filter(pl.col('epitope').is_in(sizes.filter(pl.col('n') >= min_epi)['epitope'].to_list()))
    background = position_background(frame['cdr3'].unique().to_list())
    events = substitution_events(selected, threads=threads)
    blosum = named_scores('BLOSUM62')
    matrix_rows, correlations = [], []
    for region in ('all', 'V', 'NDN', 'J'):
        selected_events = events if region == 'all' else events.filter(
            (pl.col('pos') < vlen) if region == 'V' else
            ((pl.col('pos') >= vlen) & (pl.col('pos') >= pl.col('L')-jlen)) if region == 'J' else
            ((pl.col('pos') >= vlen) & (pl.col('pos') < pl.col('L')-jlen)))
        context = {'region': region, 'events': selected_events.height,
                   'eligible_epitopes': sizes.filter(pl.col('n') >= min_epi).height,
                   'background_junctions': frame['cdr3'].n_unique()}
        if not selected_events.height:
            correlations.append(context | {'r_blosum62': None, 'r_legacy': None})
            continue
        scores = estimate_matrix(selected_events.select('a', 'b', 'pos', 'L').iter_rows(), background)
        correlations.append(context | {'r_blosum62': corr(scores, blosum), 'r_legacy': corr(scores, legacy)})
        matrix_rows.extend(context | {'a': a, 'b': b, 'log2_odds': value} for (a, b), value in scores.items())
    schema = {'region': pl.String, 'events': pl.Int64, 'eligible_epitopes': pl.Int64,
              'background_junctions': pl.Int64}
    return (pl.DataFrame(matrix_rows, schema=schema | {'a': pl.String, 'b': pl.String, 'log2_odds': pl.Float64}),
            pl.DataFrame(correlations, schema=schema | {'r_blosum62': pl.Float64, 'r_legacy': pl.Float64}))


def _penalties(dis):
    table = np.zeros((256, 256))
    for (a, b), value in dis.items():
        table[ord(a), ord(b)] = value
        table[ord(b), ord(a)] = value
    return table


def _score_candidates(queries, refs, candidates, *, locus, retention, significance, matrices):
    """Sparse candidate reductions in length batches, no per-hit Python scoring."""
    pairs = [(i, h.ref_id) for i, hits in enumerate(candidates) for h in hits
             if refs['cdr3'][h.ref_id] != queries['cdr3'][i]]
    if not pairs:
        return np.empty((0, 2), dtype=np.int64), {name: np.empty(0) for name in matrices}
    ids = np.asarray(pairs, dtype=np.int64)
    qs, rs = queries['cdr3'].to_list(), refs['cdr3'].to_list()
    qv, qj = queries['v'].to_list(), queries['j'].to_list()
    lengths = np.asarray([len(s) for s in qs])
    scores = {name: np.empty(len(ids)) for name in matrices}
    for length in np.unique(lengths[ids[:, 0]]):
        positions = np.flatnonzero(lengths[ids[:, 0]] == length)
        chunk = ids[positions]
        qids, rids = np.unique(chunk[:, 0]), np.unique(chunk[:, 1])
        if any(len(rs[i]) != length for i in rids):
            raise RuntimeError('substitution-only native candidates changed junction length')
        left = np.frombuffer(''.join(qs[i] for i in qids).encode('ascii'), dtype=np.uint8).reshape(-1, length)
        right = np.frombuffer(''.join(rs[i] for i in rids).encode('ascii'), dtype=np.uint8).reshape(-1, length)
        qlookup, rlookup = np.searchsorted(qids, chunk[:, 0]), np.searchsorted(rids, chunk[:, 1])
        left, right = left[qlookup], right[rlookup]
        weights = (np.asarray([regions.position_weights(int(length), qv[i], qj[i], locus, retention)
                              for i in qids])[qlookup]
                   if any(mode == 'region' for _, mode in matrices.values()) else None)
        positional = (np.asarray(regions.significance_weights(int(length), significance))
                      if any(mode == 'possig' for _, mode in matrices.values()) else None)
        for name, (dis, mode) in matrices.items():
            costs = _penalties(dis)[left, right] * (left != right)
            if mode == 'region':
                costs *= weights
            elif mode == 'possig':
                costs *= positional
            # Sequential accumulation preserves the old scorer's floating ties.
            np.add.accumulate(costs, axis=1, out=costs)
            scores[name][positions] = -costs[:, -1]
    return ids, scores


def loo_audit(frame, *, locus, retention, significance, structural, min_epi=30,
              top=8, max_queries=600, subs=2, seed=0, threads=1):
    """Legacy epitope-withheld matrix training and balanced PR candidate retrieval."""
    uc = long_list(frame, cap=3000, min_n=min_epi, seed=seed)
    if not uc.height:
        raise ValueError('no eligible legacy composition-controlled clonotypes')
    bg = position_background(uc['cdr3'].unique().to_list())
    sizes = (uc.group_by('epitope').agg(pl.col('cdr3').n_unique().alias('n'))
        .filter(pl.col('n') >= min_epi).sort(['n', 'epitope'], descending=[True, False]))
    held = sizes['epitope'].to_list()[:top]
    if not held:
        raise ValueError('no held-out epitopes meet the unique-junction minimum')
    # Preserve first-label semantics in the canonical long_list input order.
    refs = uc.group_by('cdr3', maintain_order=True).agg(pl.col('epitope').first())
    queries = pl.concat([uc.filter(pl.col('epitope') == epi).unique('cdr3', maintain_order=True)
        .sort('cdr3').head(max_queries) for epi in held]).with_row_index('query_id')
    candidates = Index.build(refs['cdr3'].to_list(), 'aa').search_batch(queries['cdr3'].to_list(),
        SearchParams(max_subs=subs, max_total_edits=subs, engine='seqtm'), threads)
    events = substitution_events(uc.filter(pl.col('epitope').is_in(sizes['epitope'].to_list())), threads=threads)
    events = events.filter((pl.col('pos') >= 4) & (pl.col('pos') < pl.col('L')-6))
    blo = named_dissim('BLOSUM62')
    matrices = {'unit': ({(a, b): 1. for a in AA for b in AA if a != b}, None),
        'blosum': (blo, None), 'pam': (named_dissim('PAM250'), None),
        'vdjamr': (codon_dissim(), None), 'bps': (blo, 'possig')}
    if structural is not None:
        matrices['struct'] = (structural, None)
    ids, scores = _score_candidates(queries, refs, candidates, locus=locus,
        retention=retention, significance=significance, matrices=matrices)
    rows, predictions = [], []
    qepi = np.asarray(queries['epitope'].to_list())
    refepi = np.asarray(refs['epitope'].to_list())
    for epi in held:
        query_ids = np.flatnonzero(qepi == epi)
        selected = np.isin(ids[:, 0], query_ids)
        evs = events.filter(pl.col('epitope') != epi)
        vdis = dissim_from_scores(estimate_matrix(evs.select('a', 'b', 'pos', 'L').iter_rows(), bg))
        # Candidates are already fixed by the one locus-wide retrieval batch.
        qpart = queries.filter(pl.col('epitope') == epi)
        candidate_part = [candidates[i] for i in query_ids]
        local_ids, local_scores = _score_candidates(qpart, refs, candidate_part, locus=locus,
            retention=retention, significance=significance,
            matrices={'vdjam': (vdis, None), 'region': (vdis, 'region')})
        values = {name: score[selected] for name, score in scores.items()} | local_scores
        pair_ids = ids[selected]
        if not np.array_equal(query_ids[local_ids[:, 0]], pair_ids[:, 0]) or not np.array_equal(local_ids[:, 1], pair_ids[:, 1]):
            raise RuntimeError('legacy held-out scoring changed candidate identities')
        labels = (refepi[pair_ids[:, 1]] == epi).astype(np.int64)
        context = {'epitope': epi, 'queries': len(query_ids), 'candidate_pairs': len(labels),
            'positive_pairs': int(labels.sum()), 'negative_pairs': int(len(labels)-labels.sum()),
            'training_events': evs.height, 'no_hit_queries': int(len(query_ids)-len(np.unique(pair_ids[:, 0])))}
        for name, score in values.items():
            rows.append(context | {'method': name, 'balanced_pr_auc': pr_auc_balanced(zip(labels, score))})
            predictions.append(pl.DataFrame({'query_id': pair_ids[:, 0], 'reference_key_id': pair_ids[:, 1],
                'label': labels, 'score': score}).with_columns(pl.lit(epi).alias('epitope'), pl.lit(name).alias('method')))
    results = pl.DataFrame(rows).sort('epitope', 'method')
    summary = results.group_by('method').agg(pl.len().alias('tasks'), pl.col('balanced_pr_auc').mean().alias('mean_balanced_pr_auc')).sort('method')
    return results, summary, pl.concat(predictions), queries, refs
