"""Fixed label-shuffle evidence for V usage and V-conditional core motifs.

The default 999 shuffles and per-axis 0.005 cutoff are fixed before outcomes. Their
union bound is 0.01 per profile, not database-wide or final-score calibration.
"""
import numpy as np
import polars as pl
from scipy import sparse
from scipy.special import rel_entr

from .structure import _gene_expression
from ..evalue.epitope_rules import _tokens

NOISE_RECIPE = 'v3:label-independent-seed0-frequency-imputed-V;original-key-identities-retained;999-shuffles;seed0-independent-streams;global-V-KL;within-V-max-enriched-core-exact3-KL;axis-p<=.005'
_UNKNOWN = '__UNRESOLVED_V__'


def profile_noise(rows, *, locus, eligible_targets, motif_modes=('exact3',), shuffles=999):
    """Return selection-aware gene/motif statistics per eligible ``_target``.

    All supplied rows define the global distinct junction/V-key population.
    Shuffle complete binary target-membership vectors, retaining shared labels.
    G shuffles globally; M shuffles only within V, retaining target V counts.
    Each motif statistic is one best enriched presence projection, not a sum of
    dependent tokens. Missing V is imputed once per original junction/V key from
    the global resolved-key gene frequencies, independently of target labels.
    Original keys remain distinct and all denominators retain missing-V keys.
    An entirely unresolved population uses one neutral missing-gene category.
    ``motif_modes=('exact3','wild4')`` tests the union of full20 core 3-mers
    and all four single-X 4-mer masks. Maximum selection repeats over the
    complete union inside every shuffle, using the same V-conditional null.
    """
    if type(shuffles) is not int or shuffles < 1:
        raise ValueError('positive integer shuffles required')
    if locus not in {'TRA', 'TRB'} or not {'_target', 'cdr3', 'v'}.issubset(rows.columns):
        raise ValueError('profile noise requires TRA/TRB and target/junction/V rows')
    if tuple(motif_modes) not in {('exact3',), ('exact3', 'wild4')}:
        raise ValueError('profile motifs require exact3 or the exact3/wild4 union')
    targets = sorted(set(eligible_targets))
    if not targets:
        return {}
    source = rows.select('_target', 'cdr3', _gene_expression(locus, 'V').fill_null(_UNKNOWN).alias('_v'))
    if source['cdr3'].null_count() or not source['cdr3'].str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').all():
        raise ValueError('profile noise requires valid full amino-acid junctions')
    keys = source.select('cdr3', '_v').unique().sort('_v', 'cdr3').with_row_index('_key')
    if keys.is_empty():
        return {target: _empty() for target in targets}
    genes = keys['_v'].unique().sort().to_list()
    gene_ids = keys['_v'].replace_strict(dict(zip(genes, range(len(genes)))),return_dtype=pl.Int64).to_numpy().copy()
    if _UNKNOWN in genes and len(genes)>1:
        missing=gene_ids==genes.index(_UNKNOWN)
        counts=np.bincount(gene_ids[~missing],minlength=len(genes))
        imputation_rng=np.random.default_rng(np.random.SeedSequence(0).spawn(3)[2])
        gene_ids[missing]=imputation_rng.choice(len(genes),size=int(missing.sum()),p=counts/counts.sum())
    keys=keys.drop('_key').with_columns(pl.Series('_gene',gene_ids)).sort('_gene','cdr3','_v').with_row_index('_key')
    columns = pl.DataFrame({'_target': targets, '_column': range(len(targets))})
    membership = source.join(columns, on='_target', how='inner', validate='m:1').join(
        keys, on=['cdr3', '_v'], validate='m:1').select('_key', '_column').unique()
    B = sparse.coo_matrix((np.ones(membership.height, dtype=np.int64),
        (membership['_key'].to_numpy(), membership['_column'].to_numpy())),
        shape=(keys.height, len(targets))).tocsc()
    # Identical pooled/pMHC target memberships require only one null calculation.
    seen, representatives, inverse = {}, [], []
    for column in range(B.shape[1]):
        identity = B.indices[B.indptr[column]:B.indptr[column+1]].tobytes()
        if identity not in seen:
            seen[identity] = len(representatives)
            representatives.append(column)
        inverse.append(seen[identity])
    B = B[:, representatives].tocsr()
    n = np.asarray(B.sum(axis=0)).ravel()
    safe_n = np.maximum(n, 1)
    gene_ids = keys['_gene'].to_numpy()
    Z = sparse.coo_matrix((np.ones(keys.height, dtype=np.int64),
        (np.arange(keys.height), gene_ids)), shape=(keys.height, len(genes))).tocsr()
    global_v = np.bincount(gene_ids, minlength=len(genes))
    q_v = global_v/keys.height
    observed_v = (B.T@Z).toarray()
    lengths = keys['cdr3'].str.len_chars().to_numpy()
    _, length_ids = np.unique(lengths, return_inverse=True)
    L = sparse.coo_matrix((np.ones(keys.height, dtype=np.int64),
        (np.arange(keys.height), length_ids)), shape=(keys.height, int(length_ids.max())+1)).tocsr()
    q_length = np.asarray(L.sum(axis=0)).ravel()/keys.height
    p_length = (B.T@L).toarray()/safe_n[:, None]
    length_kl = rel_entr(p_length, q_length).sum(axis=1)/np.log(2.)
    mixture = (p_length+q_length)/2
    length_js = (rel_entr(p_length, mixture).sum(axis=1)+rel_entr(q_length, mixture).sum(axis=1))/(2*np.log(2.))

    def gene_stat(counts):
        return rel_entr(counts/safe_n[:, None], q_v[None, :]).sum(axis=1)/np.log(2.)

    token_keys = keys.select(pl.col('_key').alias('_id'), pl.col('cdr3').alias('seq0'), pl.col('_gene').alias('v0'))
    tokens = pl.concat([_tokens(token_keys, [0], mode, wildcard_positions=(0, 1, 2, 3)) for mode in motif_modes])
    vocabulary = tokens.select('v0', 't0').unique().sort('v0', 't0').with_row_index('_feature')
    present = tokens.join(vocabulary, on=['v0', 't0'], validate='m:1')
    X = sparse.coo_matrix((np.ones(present.height, dtype=np.int64),
        (present['_id'].to_numpy(), present['_feature'].to_numpy())),
        shape=(keys.height, vocabulary.height)).tocsr()
    feature_v = vocabulary['v0'].to_numpy().astype(np.int64)
    global_m = np.asarray(X.sum(axis=0)).ravel()
    q_m = global_m/global_v[feature_v]

    def motif_stat(counts, *, selected=False):
        counts = counts.tocoo()
        target, feature, k = counts.row, counts.col, counts.data
        nv = observed_v[target, feature_v[feature]]
        p, q = k/nv, q_m[feature]
        enriched = p > q
        target, feature, p, q, nv = [a[enriched] for a in (target, feature, p, q, nv)]
        values = (nv/safe_n[target])*(rel_entr(p, q)+rel_entr(1-p, 1-q))/np.log(2.)
        maximum = np.zeros(B.shape[1])
        np.maximum.at(maximum, target, values)
        if not selected:
            return maximum
        chosen = np.full(B.shape[1], vocabulary.height, dtype=np.int64)
        best = (values > 0) & (values == maximum[target])
        np.minimum.at(chosen, target[best], feature[best])
        return maximum, chosen

    gene_observed = gene_stat(observed_v)
    motif_counts = B.T@X
    motif_observed, chosen = motif_stat(motif_counts, selected=True)
    global_rng, motif_rng = [np.random.default_rng(seed) for seed in np.random.SeedSequence(0).spawn(2)]
    # Keys are sorted by V, so each block is contiguous and row shuffles stay sparse.
    boundaries = np.r_[0, np.cumsum(global_v)]
    groups = [np.arange(a, b) for a, b in zip(boundaries[:-1], boundaries[1:])]
    gene_sum, motif_sum = np.zeros(B.shape[1]), np.zeros(B.shape[1])
    gene_max, motif_max = np.zeros(B.shape[1]), np.zeros(B.shape[1])
    gene_ge, motif_ge = np.zeros(B.shape[1], dtype=int), np.zeros(B.shape[1], dtype=int)
    for _ in range(shuffles):
        g = gene_stat((B[global_rng.permutation(keys.height)].T@Z).toarray())
        permutation = np.concatenate([motif_rng.permutation(group) for group in groups])
        m = motif_stat(B[permutation].T@X)
        gene_sum += g; motif_sum += m
        gene_max = np.maximum(gene_max, g); motif_max = np.maximum(motif_max, m)
        gene_ge += g >= gene_observed; motif_ge += m >= motif_observed
    gene_mean, motif_mean = gene_sum/shuffles, motif_sum/shuffles
    gene_p, motif_p = (gene_ge+1)/(shuffles+1), (motif_ge+1)/(shuffles+1)
    unique_junctions = source.join(columns, on='_target', how='inner', validate='m:1').group_by('_column').agg(
        pl.col('cdr3').n_unique().alias('_junctions'))
    junction_counts = dict(unique_junctions.iter_rows())
    motif_lookup = motif_counts.tocsr()
    output = {}
    for column, target in enumerate(targets):
        i = inverse[column]
        status = ('empty_target' if not n[i] else 'insufficient_unique_junctions' if junction_counts.get(column, 0) < 2
            else 'whole_population_target' if n[i] == keys.height else 'observed')
        unresolved_axis = status != 'observed'
        selected = int(chosen[i])
        has_motif = selected < vocabulary.height
        v = int(feature_v[selected]) if has_motif else None
        output[target] = dict(gene_kl_bits=float(gene_observed[i]), gene_null_mean_bits=float(gene_mean[i]),
            length_database_kl_bits=float(length_kl[i]) if n[i] else None,
            length_database_js_bits=float(length_js[i]) if n[i] else None,
            gene_excess_bits=float(gene_observed[i]-gene_mean[i]), gene_null_max_bits=float(gene_max[i]),
            gene_p=float(gene_p[i]), gene_supported=None if unresolved_axis else bool(gene_p[i] <= .005 and gene_observed[i] > gene_mean[i]),
            motif_kl_bits=float(motif_observed[i]), motif_null_mean_bits=float(motif_mean[i]),
            motif_excess_bits=float(motif_observed[i]-motif_mean[i]), motif_null_max_bits=float(motif_max[i]),
            motif_p=float(motif_p[i]), motif_supported=None if unresolved_axis else bool(motif_p[i] <= .005 and motif_observed[i] > motif_mean[i]),
            motif_v_gene=genes[v] if has_motif else None, motif=vocabulary['t0'][selected] if has_motif else None,
            motif_reference_with_motif=int(motif_lookup[i, selected]) if has_motif else None,
            motif_reference_v_keys=int(observed_v[i, v]) if has_motif else None,
            motif_background_with_motif=int(global_m[selected]) if has_motif else None,
            motif_background_v_keys=int(global_v[v]) if has_motif else None,
            noise_shuffles=shuffles, noise_status=status)
    return output


def _empty():
    return dict(gene_kl_bits=0., gene_null_mean_bits=0., gene_excess_bits=0., gene_null_max_bits=0.,
        length_database_kl_bits=None, length_database_js_bits=None,
        gene_p=1., gene_supported=None, motif_kl_bits=0., motif_null_mean_bits=0., motif_excess_bits=0.,
        motif_null_max_bits=0., motif_p=1., motif_supported=None, motif_v_gene=None, motif=None,
        motif_reference_with_motif=None, motif_reference_v_keys=None, motif_background_with_motif=None,
        motif_background_v_keys=None, noise_shuffles=0, noise_status='empty_reference')
