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


# Origin-aware block estimation; legacy matrix/LOO functions above remain unchanged.
def column_pair_counts(encoded, labels, mask, *, include_within=False):
    """BLOSUM between-cluster counts and the column-composition expectation.

    Boolean masks and fractional pure-N posterior weights use the same counter.
    Unknown origins must be accounted separately; NaNs are rejected.
    The existing cluster weights carry unit mass per cluster. Both returned
    matrices count ordered residue pairs, and therefore have equal total mass.
    """
    encoded = np.asarray(encoded)
    labels = np.asarray(labels)
    if (
        encoded.ndim != 2
        or labels.shape != (len(encoded),)
        or not np.issubdtype(labels.dtype, np.integer)
        or np.any(labels < 0)
    ):
        raise ValueError(
            "encoded sites need one nonnegative integer cluster label per row"
        )
    if (
        not np.issubdtype(encoded.dtype, np.integer)
        or np.any(encoded < 0)
        or np.any(encoded > 20)
    ):
        raise ValueError("encoded residues must be integers 0..20 (20 denotes a gap)")
    _, labels = np.unique(labels, return_inverse=True)
    mask = np.asarray(mask)
    if (
        mask.shape != encoded.shape
        or not np.isfinite(mask).all()
        or np.any(mask < 0)
        or np.any(mask > 1)
    ):
        raise ValueError(
            "origin weights must be finite probabilities matching encoded sites"
        )
    n_clusters = int(labels.max() + 1) if len(labels) else 0
    size = np.bincount(labels, minlength=n_clusters).astype(float)
    observed = np.zeros((20, 20))
    expected = np.zeros((20, 20))
    stats = dict(masked_columns=0, kept_columns=0, column_pairs=0.0)
    if n_clusters < 2 and not include_within:
        return observed, expected, stats
    for column in range(encoded.shape[1]):
        usable = (mask[:, column] > 0) & (encoded[:, column] < 20)
        if usable.sum() < 2:
            stats["masked_columns"] += 1
            continue
        c = np.zeros((n_clusters, 20))
        np.add.at(c, (labels[usable], encoded[usable, column]), mask[usable, column])
        c /= size[:, None]
        total = c.sum(axis=0)
        cross = np.outer(total, total)
        if not include_within:
            cross -= c.T @ c
        cross = np.maximum((cross + cross.T) / 2, 0.0)
        mass = cross.sum()
        if mass <= 0:
            continue
        observed += cross
        marginal = cross.sum(axis=1) / mass
        expected += mass * np.outer(marginal, marginal)
        stats["kept_columns"] += 1
        stats["column_pairs"] += float(mass / 2)
    return observed, expected, stats


def similarity_from_counts(
    observed, expected=None, lambda_bits=0.5, *, round_scores=True
):
    """Normalized ordered-pair log odds; shared by ungapped and gapped builders."""
    if not np.isfinite(observed).all() or observed.sum() <= 0:
        raise ValueError("no finite between-cluster column pairs")
    q = observed / observed.sum()
    p = q.sum(axis=1)
    if expected is None:
        expected = np.outer(p, p)
    else:
        if not np.isfinite(expected).all() or expected.sum() <= 0:
            raise ValueError("no finite expected column pairs")
        expected = expected / expected.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        score = np.where(
            (q > 0) & (expected > 0), np.log2(q / expected) / lambda_bits, np.nan
        )
    if not round_scores:
        return score, p
    floor = np.nanmin(score[np.isfinite(score)])
    score = np.where(np.isfinite(score), score, floor)
    return np.rint(score).astype(int), p


def costs_from_similarity(similarity):
    """The established nonnegative Gram transform used by VDJmatch matrix geometry."""
    cost = np.maximum(
        0, np.add.outer(np.diag(similarity), np.diag(similarity)) - 2 * similarity
    )
    np.fill_diagonal(cost, 0)
    return cost


def origin_block_matrices(clonotypes, origins, known_labels=None):
    """Same-epitope/V/J/length block counts with exact pure-N posterior weights.

    ``clonotypes`` has unique (epitope, query_id) rows and columns cdr3,v,j;
    ``origins`` has every zero-based codon of each query, including unsupported
    nulls. Unique clonotypes are unit clusters: this does not claim 62% clustering.
    Returns pooled/equal-epitope matrix rows and per-block support accounting.
    No pseudocount, reference matrix, unsupported-cell floor or model refit.
    """
    if clonotypes.select("epitope", "query_id").is_duplicated().any():
        raise ValueError("clonotype membership must be unique within each epitope")
    if clonotypes.select("epitope", "cdr3", "v", "j").is_duplicated().any():
        raise ValueError(
            "duplicate junction/V/J clonotypes must be pooled before counting"
        )
    if origins.select("query_id", "position").is_duplicated().any():
        raise ValueError("origin coordinates must be unique per query")
    queries = clonotypes.select("query_id", "cdr3").unique()
    if queries["query_id"].n_unique() != queries.height:
        raise ValueError("one query ID must identify exactly one junction")
    sizes = queries.select("query_id", pl.col("cdr3").str.len_chars().alias("L"))
    checked = origins.join(sizes, on="query_id", how="inner", validate="m:1")
    if (
        checked.height != origins.height
        or checked.filter(
            (pl.col("position") < 0) | (pl.col("position") >= pl.col("L"))
        ).height
    ):
        raise ValueError("origin coordinates do not match junctions")
    coverage = (
        checked.group_by("query_id").len().join(sizes, on="query_id", validate="1:1")
    )
    if (
        coverage.height != sizes.height
        or coverage.filter(pl.col("len") != pl.col("L")).height
    ):
        raise ValueError(
            "every query requires every codon, including unsupported sites"
        )
    values = origins["n_probability"].drop_nulls().to_numpy()
    if not np.isfinite(values).all() or np.any(values < 0) or np.any(values > 1):
        raise ValueError("pure-N posterior must be null or a finite probability")
    supported = origins.group_by("query_id").agg(
        pl.col("n_probability").is_not_null().all().alias("origin_supported")
    )
    selected = (
        clonotypes.join(supported, on="query_id", validate="m:1")
        .filter(pl.col("origin_supported"))
        .sort("epitope", "v", "j", "cdr3", "query_id")
    )
    rows = []
    counts = {}
    short_counts = {}
    contrast = None
    if selected.height:
        ids = queries["query_id"].to_list()
        lookup = {q: i for i, q in enumerate(ids)}
        width = int(sizes["L"].max())
        weights = np.zeros((len(ids), width))
        valid = origins.filter(pl.col("n_probability").is_not_null())
        weights[
            np.array([lookup[q] for q in valid["query_id"]]),
            valid["position"].to_numpy(),
        ] = valid["n_probability"].to_numpy()
        selected = selected.with_columns(pl.col("cdr3").str.len_chars().alias("L"))
        aa_lookup = np.full(128, 20, dtype=np.uint8)
        aa_lookup[np.frombuffer(AA.encode(), dtype=np.uint8)] = np.arange(20)
        for (epitope, v, j, length), block in selected.group_by(
            "epitope", "v", "j", "L", maintain_order=True
        ):
            if not block["cdr3"].str.contains(VALID).all():
                raise ValueError("block junctions require standard amino acids")
            encoded = aa_lookup[
                np.frombuffer("".join(block["cdr3"]).encode(), dtype=np.uint8)
            ].reshape(-1, length)
            mask = weights[[lookup[q] for q in block["query_id"]], :length]
            observed, _, stats = column_pair_counts(
                encoded, np.arange(block.height), mask
            )
            counts.setdefault(epitope, np.zeros((20, 20)))[:] += observed
            if length < 8:
                short_counts.setdefault(epitope, np.zeros((20, 20)))[:] += observed
            rows.append(
                dict(
                    epitope=epitope,
                    v=v,
                    j=j,
                    length=length,
                    clonotypes=block.height,
                    n_weighted_unordered_site_pairs=stats["column_pairs"],
                    columns_with_pairs=stats["kept_columns"],
                )
            )
        contrast = _matched_origin_background(
            selected, weights, lookup, aa_lookup, known_labels
        )
    matrix_rows = []
    contributing = [c for c in counts.values() if c.sum() > 0]
    pooled = sum(contributing, np.zeros((20, 20)))
    balanced = sum((c / c.sum() for c in contributing), np.zeros((20, 20)))
    for weighting, observed in [("pooled", pooled), ("equal_epitope", balanced)]:
        if observed.sum() > 0:
            score, p = similarity_from_counts(
                observed, lambda_bits=1.0, round_scores=False
            )
            q = observed / observed.sum()
            expected = np.outer(p, p)
            cost = costs_from_similarity(score)
            identified = (
                np.isfinite(score)
                & np.isfinite(np.diag(score))[:, None]
                & np.isfinite(np.diag(score))[None, :]
            )
            cost = np.where(identified, cost, np.nan)
        else:
            q = expected = np.zeros((20, 20))
            score = cost = np.full((20, 20), np.nan)
        for i, a in enumerate(AA):
            for k, b in enumerate(AA):
                matrix_rows.append(
                    dict(
                        weighting=weighting,
                        a=a,
                        b=b,
                        pair_weight=float(observed[i, k]),
                        q=float(q[i, k]),
                        background_q=float(expected[i, k]),
                        score_bits=float(score[i, k]),
                        distance_bits=float(cost[i, k]),
                        supported=bool(
                            np.isfinite(score[i, k])
                            and np.isfinite(score[i, i])
                            and np.isfinite(score[k, k])
                        ),
                        contributing_epitopes=len(contributing),
                    )
                )
    matrix = pl.DataFrame(matrix_rows).with_columns(
        pl.col("score_bits", "distance_bits").fill_nan(None)
    )
    blocks = pl.DataFrame(
        rows,
        schema={
            "epitope": pl.String,
            "v": pl.String,
            "j": pl.String,
            "length": pl.Int64,
            "clonotypes": pl.Int64,
            "n_weighted_unordered_site_pairs": pl.Float64,
            "columns_with_pairs": pl.Int64,
        },
    )
    epitope_pairs = pl.DataFrame(
        [
            dict(epitope=e, a=a, b=b, n_weighted_ordered_site_pairs=float(c[i, k]))
            for e, c in sorted(counts.items())
            for i, a in enumerate(AA)
            for k, b in enumerate(AA)
        ],
        schema={
            "epitope": pl.String,
            "a": pl.String,
            "b": pl.String,
            "n_weighted_ordered_site_pairs": pl.Float64,
        },
    )
    if contrast is None:
        contrast = _matched_origin_background(None, None, None, None, None)
    contrast["short-epitope-pairs"] = pl.DataFrame(
        [
            dict(
                epitope=e,
                a=AA[i],
                b=AA[j],
                n_weighted_ordered_site_pairs=float(c[i, j]),
            )
            for e, c in sorted(short_counts.items())
            for i, j in zip(*np.nonzero(c))
        ],
        schema={
            "epitope": pl.String,
            "a": pl.String,
            "b": pl.String,
            "n_weighted_ordered_site_pairs": pl.Float64,
        },
    )
    return (
        matrix,
        blocks.sort("epitope", "v", "j", "length"),
        supported.sort("query_id"),
        epitope_pairs,
        contrast,
    )


def _matched_origin_background(selected, weights, lookup, aa_lookup, known_labels):
    """Matched disjoint-label block background; no query-pair enumeration.

    Positive and background counts share epitope/V/J/length/position/N opportunity.
    Missing background excludes that site's positive mass from BOTH matrices.
    Complete known labels, including ineligible tails, define disjointness only.
    Different-labelled receptors are selected background, not assay negatives.
    """
    positive = {}
    background = {}
    directional_background = {}
    short_positive = {}
    short_background = {}
    short_directional = {}
    support_rows = []
    if selected is not None and selected.height:
        primary = selected.select("query_id", "epitope").unique()
        known_labels = (
            primary
            if known_labels is None
            else known_labels.select("query_id", "epitope").unique()
        )
        if primary.join(known_labels, on=["query_id", "epitope"], how="anti").height:
            raise ValueError("known label sets must retain every primary label")
        label_sets = {
            q: frozenset(part["epitope"].to_list())
            for (q,), part in known_labels.group_by("query_id")
        }
        for (v, j, length), membership in selected.group_by(
            "v", "j", "L", maintain_order=True
        ):
            block = membership.unique("query_id").sort("query_id")
            signatures = [label_sets[q] for q in block["query_id"]]
            groups = sorted(set(signatures), key=lambda x: tuple(sorted(x)))
            if len(groups) < 2:
                allowed = np.zeros((len(groups), len(groups)))
            else:
                allowed = np.array(
                    [[g.isdisjoint(h) for h in groups] for g in groups], dtype=float
                )
            group_id = {g: i for i, g in enumerate(groups)}
            ids = np.array([group_id[g] for g in signatures])
            encoded = aa_lookup[
                np.frombuffer("".join(block["cdr3"]).encode(), dtype=np.uint8)
            ].reshape(-1, length)
            mask = weights[[lookup[q] for q in block["query_id"]], :length]
            active = []
            for epitope in sorted(membership["epitope"].unique().to_list()):
                positive.setdefault(epitope, np.zeros((20, 20)))
                background.setdefault(epitope, np.zeros((20, 20)))
                directional_background.setdefault(epitope, np.zeros((20, 20)))
                short_positive.setdefault(epitope, np.zeros((20, 20)))
                short_background.setdefault(epitope, np.zeros((20, 20)))
                short_directional.setdefault(epitope, np.zeros((20, 20)))
                in_group = np.array([epitope in g for g in groups])
                in_query = in_group[ids]
                n = int(in_query.sum())
                if n >= 2:
                    active.append((epitope, in_group, in_query, n))
            for position in range(length):
                c = np.zeros((len(groups), 20))
                np.add.at(c, (ids, encoded[:, position]), mask[:, position])
                if c.sum() == 0:
                    continue
                others = allowed @ c
                for epitope, in_group, in_query, n in active:
                    pos, _, _ = column_pair_counts(
                        encoded[in_query, position : position + 1],
                        np.arange(n),
                        mask[in_query, position : position + 1],
                    )
                    pos_mass = float(pos.sum())
                    if pos_mass <= 0:
                        continue
                    directed = c[in_group].T @ others[in_group]
                    neg = directed + directed.T
                    neg_mass = float(neg.sum())
                    matched = neg_mass > 0
                    if matched:
                        positive[epitope] += pos
                        background[epitope] += neg * (pos_mass / neg_mass)
                        directional_background[epitope] += directed * (
                            2 * pos_mass / neg_mass
                        )
                        if length < 8:
                            short_positive[epitope] += pos
                            short_background[epitope] += neg * (pos_mass / neg_mass)
                            short_directional[epitope] += directed * (
                                2 * pos_mass / neg_mass
                            )
                    support_rows.append(
                        dict(
                            epitope=epitope,
                            v=v,
                            j=j,
                            length=length,
                            position=position,
                            positive_ordered_pair_weight=pos_mass,
                            background_ordered_pair_weight=neg_mass,
                            matched=matched,
                        )
                    )
    cells = []
    for epitope in sorted(positive):
        assert np.isclose(positive[epitope].sum(), background[epitope].sum())
        for i, a in enumerate(AA):
            for k, b in enumerate(AA):
                cells.append(
                    dict(
                        epitope=epitope,
                        a=a,
                        b=b,
                        positive_pair_weight=float(positive[epitope][i, k]),
                        background_pair_weight=float(background[epitope][i, k]),
                        directional_background_pair_weight=float(
                            directional_background[epitope][i, k]
                        ),
                        short_positive_pair_weight=float(short_positive[epitope][i, k]),
                        short_background_pair_weight=float(
                            short_background[epitope][i, k]
                        ),
                        short_directional_background_pair_weight=float(
                            short_directional[epitope][i, k]
                        ),
                    )
                )
    rows = []
    labels = [e for e in sorted(positive) if positive[e].sum() > 0]
    for weighting in ["pooled", "equal_epitope"]:
        scales = {
            e: 1.0 if weighting == "pooled" else 1.0 / positive[e].sum() for e in labels
        }
        pos = sum((scales[e] * positive[e] for e in labels), np.zeros((20, 20)))
        neg = sum((scales[e] * background[e] for e in labels), np.zeros((20, 20)))
        directed = sum(
            (scales[e] * directional_background[e] for e in labels), np.zeros((20, 20))
        )
        for estimator, bg in [
            ("unordered_density", neg),
            ("directional_log_mean", directed),
        ]:
            if pos.sum() > 0 and bg.sum() > 0:
                score, _ = similarity_from_counts(
                    pos, bg, lambda_bits=1.0, round_scores=False
                )
                if estimator == "directional_log_mean":
                    score = (score + score.T) / 2
                q = pos / pos.sum()
                q0 = bg / bg.sum()
                raw = np.add.outer(np.diag(score), np.diag(score)) - 2 * score
                identified = (
                    np.isfinite(score)
                    & np.isfinite(np.diag(score))[:, None]
                    & np.isfinite(np.diag(score))[None, :]
                )
                raw = np.where(identified, raw, np.nan)
                distance = np.where(identified, costs_from_similarity(score), np.nan)
            else:
                q = q0 = np.zeros((20, 20))
                score = raw = distance = np.full((20, 20), np.nan)
            for i, a in enumerate(AA):
                for k, b in enumerate(AA):
                    rows.append(
                        dict(
                            estimator=estimator,
                            weighting=weighting,
                            a=a,
                            b=b,
                            positive_pair_weight=float(pos[i, k]),
                            background_pair_weight=float(bg[i, k]),
                            q=float(q[i, k]),
                            background_q=float(q0[i, k]),
                            score_bits=float(score[i, k]),
                            raw_distance_bits=float(raw[i, k]),
                            distance_bits=float(distance[i, k]),
                            supported=bool(np.isfinite(distance[i, k])),
                            contributing_epitopes=len(labels),
                        )
                    )
    return {
        "recognition-matrix": pl.DataFrame(rows).with_columns(
            pl.col("score_bits", "raw_distance_bits", "distance_bits").fill_nan(None)
        ),
        "recognition-pairs": pl.DataFrame(
            cells,
            schema={
                "epitope": pl.String,
                "a": pl.String,
                "b": pl.String,
                "positive_pair_weight": pl.Float64,
                "background_pair_weight": pl.Float64,
                "directional_background_pair_weight": pl.Float64,
                "short_positive_pair_weight": pl.Float64,
                "short_background_pair_weight": pl.Float64,
                "short_directional_background_pair_weight": pl.Float64,
            },
        ),
        "recognition-support": pl.DataFrame(
            support_rows,
            schema={
                "epitope": pl.String,
                "v": pl.String,
                "j": pl.String,
                "length": pl.Int64,
                "position": pl.Int64,
                "positive_ordered_pair_weight": pl.Float64,
                "background_ordered_pair_weight": pl.Float64,
                "matched": pl.Boolean,
            },
        ).sort("epitope", "v", "j", "length", "position"),
    }


def n_region_residue_opportunities(model, junctions, positions, v, j, *, threads=1):
    """Exact residue probabilities conditional on unchanged flanks, V/J and pure N.

    Enumerate the 20 AA alternatives once per supplied context and reuse the two
    published native batches. Return context x AA probabilities in ``AA`` order;
    an all-NaN row means zero pure-N opportunity. Pgen values remain transient.
    Callers deduplicate flank/V/J/position contexts before this operation.
    """
    import re
    from vdjtools.model import pgen_aa_batch, n_region_posterior_batch

    n = len(junctions)
    if any(len(x) != n for x in (positions, v, j)) or threads < 1:
        raise ValueError('matching context/position/V/J lengths and positive threads required')
    if any(not isinstance(p, (int, np.integer)) or not 0 <= p < len(s)
           or re.fullmatch(VALID, s) is None for s, p in zip(junctions, positions)):
        raise ValueError('standard-AA junctions and in-range integer codon coordinates required')
    if not n:
        return np.empty((0, 20))
    variants = [s[:p] + a + s[p + 1:] for s, p in zip(junctions, positions) for a in AA]
    vs, js = np.repeat(v, 20).tolist(), np.repeat(j, 20).tolist()
    coordinates = [[int(p)] for p in positions for _ in AA]
    mass = np.asarray(pgen_aa_batch(model, variants, v=vs, j=js, threads=threads))
    posterior = n_region_posterior_batch(
        model, variants, positions=coordinates, v=vs, j=js, threads=threads)
    mass *= posterior['n_probability'].fill_null(0).to_numpy()
    mass = mass.reshape(n, 20)
    if not np.isfinite(mass).all() or (mass < 0).any():
        raise RuntimeError('native pure-N opportunity must be finite and nonnegative')
    total = mass.sum(axis=1, keepdims=True)
    return np.divide(mass, total, out=np.full_like(mass, np.nan), where=total > 0)
