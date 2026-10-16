"""Control-relative binary motif projections, conditional on V and length."""
import numpy as np
import polars as pl

from .structure import _gene_expression
from ..match.motifs import masked_motif_presence_counts


def motif_population(records, *, locus, mode='imgt_wild4'):
    """Transient fixed key/count populations for both declared length scopes.

    Input cdr3/v uses anchor-inclusive junctions. Alleles and J annotations do
    not multiply a full-junction/V-gene key; unknown reference genes stay null.
    Caller selects species/locus before this reduction. No event/read weighting.
    ``mode`` selects exact3/imgt3/wild4/imgt_wild4 in the existing tokenizer.
    """
    if locus not in {'TRA', 'TRB'}:
        raise ValueError('motif population requires TRA or TRB')
    keys = records.select('cdr3', _gene_expression(locus, 'V').alias('v_gene'),
        pl.col('cdr3').str.len_chars().alias('junction_length')).unique(['cdr3', 'v_gene'], maintain_order=True)
    strata = keys.group_by('v_gene', 'junction_length', maintain_order=True).len(name='keys').with_row_index('_stratum')
    frame = keys.join(strata.select('v_gene', 'junction_length', '_stratum'),
        on=['v_gene', 'junction_length'], nulls_equal=True, validate='m:1', maintain_order='left').select(
            pl.col('cdr3').alias('seq0'), pl.col('_stratum').alias('v0')).with_row_index('_id')
    counts, vocabulary = masked_motif_presence_counts(frame, mode=mode)
    pooled = strata.group_by('v_gene', maintain_order=True).agg(pl.col('keys').sum()).with_row_index('_pooled')
    indices = strata.select('v_gene').join(pooled, on='v_gene', nulls_equal=True, validate='m:1', maintain_order='left')['_pooled'].to_numpy()
    pooled_counts = np.zeros((pooled.height, vocabulary.height), dtype=np.int64)
    np.add.at(pooled_counts, indices, counts)
    pooled = pooled.rename({'_pooled': '_stratum'}).with_columns(pl.lit(None, dtype=pl.UInt32).alias('junction_length'))
    return {'exact_length': (strata, counts), 'all_lengths': (pooled, pooled_counts),
            'vocabulary': vocabulary, 'mode': mode}


def motif_selection_information(reference, control, *, locus, mode='imgt_wild4'):
    """Per-token presence likelihood ratios and Bernoulli KL within V/length.

    Each token is one binary projection. Overlapping tokens are dependent;
    their KL values are not an additive joint sequence information estimate.
    Emit reference-observed tokens, retaining missing control strata/support.
    ``control`` is this invocation's :func:`motif_population` from control keys.
    Whole-population counts also expose the matched stratum factor and joint
    stratum/token presence evidence. These are single-feature projections;
    overlapping tokens must not be added as independent evidence.
    Reference and control must use the same ``mode`` (default imgt_wild4).
    """
    if control['mode'] != mode:
        raise ValueError('reference and control motif mode must match')
    observed = motif_population(reference, locus=locus, mode=mode)
    parts = []
    for scope in ('all_lengths', 'exact_length'):
        strata, counts = observed[scope]
        background, background_counts = control[scope]
        matched = strata.join(background.select('v_gene', 'junction_length',
            pl.col('_stratum').alias('_control'), pl.col('keys').alias('control_keys')),
            on=['v_gene', 'junction_length'], how='left', nulls_equal=True, validate='1:1', maintain_order='left')
        row, token = np.nonzero(counts)
        n = counts[row, token]
        N = strata['keys'].to_numpy()[row]
        M = matched['control_keys'].fill_null(0).to_numpy()[row]
        ids = matched['_control'].fill_null(0).to_numpy()[row]
        m = np.zeros(len(row), dtype=np.int64)
        present = M > 0
        m[present] = background_counts[ids[present], token[present]]
        p = n / N
        q = np.divide(m, M, out=np.full(len(row), np.nan), where=present)
        bits = np.full(len(row), np.nan)
        valid_presence = present & (m > 0)
        bits[valid_presence] = np.log2(p[valid_presence] / q[valid_presence])
        valid = valid_presence & ((n == N) | (m < M))
        information = np.full(len(row), np.nan)
        information[valid] = p[valid] * bits[valid]
        absent = valid & (n < N)
        information[absent] += (1-p[absent]) * np.log2((1-p[absent]) / (1-q[absent]))
        status = np.full(len(row), 'available', dtype=object)
        status[~present] = 'no_control_stratum'
        status[present & (m == 0)] = 'motif_not_observed_in_control'
        status[valid_presence & (m == M) & (n < N)] = 'absence_not_observed_in_control'
        genes = strata['v_gene'].gather(row)
        status[genes.is_null().to_numpy()] = 'unresolved_reference_gene'
        parts.append(pl.DataFrame({'v_gene': genes,
            'junction_length': strata['junction_length'].gather(row), 'length_scope': pl.Series([scope] * len(row), dtype=pl.String),
            'motif': observed['vocabulary']['t0'].gather(token),
            'reference_keys': N, 'control_keys': M, 'reference_with_motif': n, 'control_with_motif': m,
            'reference_fraction': p, 'control_fraction': pl.Series(q, nan_to_null=True),
            'presence_bits': pl.Series(bits, nan_to_null=True),
            'motif_information_bits': pl.Series(information, nan_to_null=True), 'information_status': pl.Series(status, dtype=pl.String)}))
    reference_total = int(observed['exact_length'][0]['keys'].sum())
    control_total = int(control['exact_length'][0]['keys'].sum())
    return pl.concat(parts).with_columns(
        pl.lit(reference_total, dtype=pl.UInt64).alias('total_reference_keys'),
        pl.lit(control_total, dtype=pl.UInt64).alias('total_control_keys')).with_columns(
        pl.when((pl.col('control_keys') > 0) & (control_total > 0)).then(
            ((pl.col('reference_keys') / reference_total) /
             (pl.col('control_keys') / control_total)).log(2)).alias('stratum_bits'),
        pl.when((pl.col('control_with_motif') > 0) & (control_total > 0)).then(
            ((pl.col('reference_with_motif') / reference_total) /
             (pl.col('control_with_motif') / control_total)).log(2)).alias('joint_presence_bits')
    ).with_columns(pl.lit(mode).alias('motif_mode'),
        pl.lit('imgt11' if 'imgt' in mode else 'full20').alias('motif_alphabet')
    ).sort('length_scope', 'v_gene', 'junction_length', 'motif', nulls_last=True)
