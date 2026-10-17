"""All-label fixed rules using bounded native counts and sparse motif projections."""
import numpy as np
import polars as pl
from scipy import sparse
from scipy.special import digamma


def stream_reference_rules(queries, reference, *, paired=False, species='human', threads=1,
                           radius=90, consume):
    """Emit the scalar rule scores without retaining a detailed case grid.

    Inputs use resolved ``seq0/v0 [seq1/v1]`` and reference epitope. Query IDs
    retain their source identity. Full models, priors and motif counts are fixed;
    only exact hits are filtered. MHC pooling is the caller's explicit choice.
    Three query-by-label numeric arrays retain spatial masses/ranks in memory;
    one motif mode is prepared at a time. Nothing is persisted or resumed.
    """
    from .epitope_rules import _prepare_rule_frame, _tokens, MOTIFS
    from .local_rank import _harmonic_interval
    from ..match.geometry import (grouped_total_distance_count_batch, distance_units,
                                  distance_cutoff, precompute_v_loop_distances)

    axes = range(2 if paired else 1)
    seq, v = [f'seq{i}' for i in axes], [f'v{i}' for i in axes]
    keys = [*seq, *v]
    if (type(threads) is not int or threads < 1 or not callable(consume)
            or not 0 < radius < 2**28/200 or not distance_cutoff(radius)):
        raise ValueError('positive threads/radius and a score consumer required')
    if not {*keys, 'query_id'} <= set(queries.columns) or not {*keys, 'epitope'} <= set(reference.columns):
        raise ValueError('resolved query/reference keys and epitope required')
    if queries['query_id'].null_count() or queries['query_id'].n_unique() != queries.height:
        raise ValueError('unique non-null query identities required')
    if any(f[c].null_count() for f in (queries, reference) for c in keys):
        raise ValueError('complete resolved junction/V keys required')
    if reference.is_empty() or reference['epitope'].null_count():
        raise ValueError('nonempty labelled reference required')
    whole = reference.select(keys).unique(maintain_order=True).with_row_index('_ref')
    labels = reference.select('epitope').unique().sort('epitope').with_row_index('_target')
    membership = reference.select(*keys, 'epitope').unique().join(whole, on=keys, validate='m:1').join(
        labels, on='epitope', validate='m:1').select('_ref', '_target').sort('_ref', '_target')
    exclusive = membership.with_columns(pl.len().over('_ref').alias('_labels')).filter(pl.col('_labels') == 1).drop('_labels')
    Q, R, T = queries.height, whole.height, labels.height
    if not Q:
        return dict(queries=0, targets=T, cases=0, numeric_state_bytes=0,
                    native_batches=0, native_chunk_size=0)
    gene_keys = lambda f: f.select(*[pl.col(c).str.replace(r'\*.*$', '') for c in v])
    genes = gene_keys(whole).unique().sort(v).with_row_index('_gene')
    G = genes.height
    rg = gene_keys(whole).join(genes, on=v, validate='m:1', maintain_order='left')['_gene'].to_numpy()
    qg = gene_keys(queries).join(genes, on=v, how='left', validate='m:1', maintain_order='left')['_gene'].fill_null(G).to_numpy()
    required = 3*Q*T*8+2*(G+1)*T*8+256*1024**2
    if required > 12*1024**3:
        raise ValueError(f'rule numeric array budget exceeded: {required} bytes')
    def incidence(frame):
        return sparse.coo_matrix((np.ones(frame.height, dtype=np.int64),
            (frame['_ref'].to_numpy(), frame['_target'].to_numpy())), shape=(R, T)).tocsr()
    W, E = incidence(membership), incidence(exclusive)
    Z = sparse.coo_matrix((np.ones(R, dtype=np.int64), (rg, np.arange(R))), shape=(G+1, R)).tocsr()
    ng = (Z@W).toarray()
    global_g = np.asarray(Z.sum(axis=1)).ravel()
    mg = global_g[:, None]-(Z@E).toarray()
    N = np.asarray(W.sum(axis=0)).ravel()
    M = R-np.asarray(E.sum(axis=0)).ravel()
    q = queries.with_row_index('_query_row')
    exact = q.select('_query_row', *seq).join(whole.select('_ref', *seq), on=seq, validate='m:m')
    def exact_matrix(frame):
        return sparse.coo_matrix((np.ones(frame.height, dtype=np.int64),
            (frame['_query_row'].to_numpy(), frame['_ref'].to_numpy())), shape=(Q, R)).tocsr()
    H = exact_matrix(exact)
    same = exact.filter(pl.Series(qg[exact['_query_row'].to_numpy()] == rg[exact['_ref'].to_numpy()]))
    HV = exact_matrix(same)
    def removed(rows, conditional):
        h = (HV if conditional else H)[rows]
        return (h@W).toarray(), np.asarray(h.sum(axis=1)).ravel()[:, None]-(h@E).toarray()
    def populations(rows):
        n, m = ng[qg[rows]], mg[qg[rows]]
        pg = (n+m+.5)/(N+M+.5*(G+1))
        supported = (n > 0)&(m > 0)
        score = np.divide((n+pg)/(N+1), (m+pg)/(M+1), out=np.ones(n.shape), where=supported)
        return n, m, supported, score
    masses = np.zeros((Q, T, 2))
    rank = np.ones((Q, T))
    native = pl.concat([membership, whole.select('_ref').with_columns(pl.lit(T, dtype=pl.UInt32).alias('_target')),
        exclusive.with_columns((pl.col('_target')+T+1).alias('_target'))]).sort('_ref', '_target')
    models = 2*T+1
    offsets = np.r_[0, np.cumsum(np.bincount(native['_ref'].to_numpy(), minlength=R))]
    context = precompute_v_loop_distances(species=species)
    scale, quantum = distance_units()
    bound = 3*scale*sum(max((f[c].str.len_chars().max() or 0) for f in (queries, whole)) for c in seq)
    cuts = list(range(0, 200*bound+1, 600*quantum))
    weights = np.maximum(0, distance_cutoff(radius)-np.asarray(cuts, dtype=np.int64))
    # Keep native counts under128MiB; reduce labels in windows instead of
    # materializing every query/label/distance float scratch array together.
    per_query = 4*((1+2*models+1)*(len(cuts)+1)+models)
    per_query_total = per_query+min(T, 256)*(len(cuts)*80+256)
    chunk_size = min(65536, (128*1024**2)//per_query, (256*1024**2)//per_query_total)
    if chunk_size < 1:
        raise ValueError('one all-label rule query exceeds128MiB native count budget')
    fallback_gene = ~((ng > 0)&(mg > 0)).all(axis=1)
    fallback = fallback_gene[qg]
    batches = 0
    for conditional in (True, False):
        part = q if conditional else q.filter(pl.Series(fallback))
        if not part.height:
            continue
        batches += 1
        iterator = grouped_total_distance_count_batch(part.select(keys), whole.select(keys), cuts,
            paired=paired, species=species, threads=threads, ctrim=3, v_loop_weight=0,
            matrix='clipped', same_v_gene=conditional, exclude_exact=True, v_loop_distances=context,
            chunk_size=chunk_size, reference_target_offsets=offsets,
            reference_target_ids=native['_target'].to_numpy(), target_peptide_ids=[0]*models,
            reference_exclusive_peptide_ids=[-1]*R)
        for chunk in iterator:
            rows = part['_query_row'].slice(chunk.start, chunk.stop-chunk.start).to_numpy()
            n, m, supported, _ = populations(rows)
            if not conditional:
                n, m = np.broadcast_to(N, n.shape), np.broadcast_to(M, m.shape)
            er, eb = removed(rows, conditional)
            counts = np.asarray(chunk.counts)
            for start in range(0, T, 256):
                targets = np.arange(start, min(T, start+256))
                tc = counts[:, 1+targets, :].astype(np.int64)
                bc = counts[:, 1+T, None, :].astype(np.int64)-counts[:, 1+T+1+targets, :]
                if (not np.array_equal(tc[:, :, -1], (n-er)[:, targets]) or
                        not np.array_equal(bc[:, :, -1], (m-eb)[:, targets])):
                    raise ArithmeticError('all-label rank bound did not cover eligible hits')
                if conditional:
                    masses[np.ix_(rows, targets, [0])] = (np.diff(tc, axis=2, prepend=0)@weights/distance_cutoff(radius))[:, :, None]
                    masses[np.ix_(rows, targets, [1])] = (np.diff(bc, axis=2, prepend=0)@weights/distance_cutoff(radius))[:, :, None]
                shell = np.diff(tc, axis=2, prepend=0)
                a, b, c = np.nonzero(shell)
                before = np.where(c > 0, bc[a, b, c-1], 0)
                numerator = np.bincount(a*len(targets)+b, weights=shell[a, b, c]*
                    _harmonic_interval(before+1, bc[a, b, c]+1), minlength=len(rows)*len(targets)).reshape(len(rows), -1)
                nn, mm = n[:, targets], m[:, targets]
                result = np.ones(nn.shape)
                usable = (nn > 0)&(mm > 0)
                result[usable] = (mm[usable]+1)*numerator[usable]/(nn[usable]*(digamma(mm[usable]+2)+np.euler_gamma))
                old = rank[np.ix_(rows, targets)]
                rank[np.ix_(rows, targets)] = result if conditional else np.where(supported[:, targets], old, result)
    block_size = max(1, 65536//T)
    def emit(name, rows, values):
        valid = np.broadcast_to((N > 0)&(M > 0), values.shape)&np.isfinite(values)
        consume(name, pl.DataFrame({'query_id': queries['query_id'].gather(np.repeat(rows, T)),
            'epitope': labels['epitope'].gather(np.tile(np.arange(T), len(rows))),
            'score': pl.Series(np.where(valid, values, np.nan).ravel(), nan_to_null=True),
            'supported': valid.ravel()}))
    def kernel(vscore, n, m, rm, bm, exact_r, exact_b, available):
        c = (rm+bm+np.where(available[:, None], exact_r+exact_b, 0)+.5)/(n+m+1)
        return vscore*((rm+c)/(n+1))/((bm+c)/(m+1))
    for start in range(0, Q, block_size):
        rows = np.arange(start, min(Q, start+block_size))
        n, m, _, vscore = populations(rows)
        er, eb = removed(rows, True)
        emit('rule_v', rows, vscore)
        emit('rule_spatial', rows, kernel(vscore, n, m, masses[rows, :, 0], masses[rows, :, 1], er, eb, np.ones(len(rows), dtype=bool)))
        emit('rule_rank', rows, vscore*rank[rows])
    prepared_r = _prepare_rule_frame(whole.select(keys), axes)
    prepared_q = _prepare_rule_frame(queries, axes, True)
    for mode in MOTIFS:
        tokens = _tokens(prepared_r, axes, mode)
        token_keys = [*v, *[f't{i}' for i in axes]]
        vocabulary = tokens.select(token_keys).unique().sort(token_keys).with_row_index('_feature')
        if vocabulary.height == 0:
            X = sparse.csr_matrix((R, 0), dtype=np.int64)
        else:
            linked = tokens.join(vocabulary, on=token_keys, validate='m:1')
            X = sparse.coo_matrix((np.ones(linked.height, dtype=np.int64),
                (linked['_id'].to_numpy(), linked['_feature'].to_numpy())), shape=(R, vocabulary.height)).tocsr()
        A, B = (X.T@W).tocsr(), (X.T@E).tocsr()
        total = np.asarray(X.sum(axis=0)).ravel()
        qt = _tokens(prepared_q, axes, mode)
        nt = np.bincount(qt['_id'].to_numpy(), minlength=Q)
        linked = qt.join(vocabulary, on=token_keys, how='inner', validate='m:1')
        QT = sparse.coo_matrix((np.ones(linked.height, dtype=np.int64),
            (linked['_id'].to_numpy(), linked['_feature'].to_numpy())), shape=(Q, vocabulary.height)).tocsr()
        if int(ng.max(initial=0))*int(mg.max(initial=0))*int(nt.max(initial=0)) > 2**53-1:
            raise ValueError('conditional motif counts exceed exactly representable integer budget')
        fg = vocabulary.select(v).join(genes, on=v, validate='m:1', maintain_order='left')['_gene'].to_numpy()
        coo = A.tocoo()
        exclusive_values = np.asarray(B[coo.row, coo.col]).ravel() if coo.nnz else np.array([], dtype=np.int64)
        bg = total[coo.row]-exclusive_values
        delta = np.maximum(0, coo.data*mg[fg[coo.row], coo.col]-bg*ng[fg[coo.row], coo.col])
        D = sparse.coo_matrix((delta, (coo.row, coo.col)), shape=A.shape).tocsr()
        for start in range(0, Q, block_size):
            rows = np.arange(start, min(Q, start+block_size))
            n, m, vsupported, vscore = populations(rows)
            er, eb = removed(rows, True)
            present = nt[rows] > 0
            divisor = np.maximum(nt[rows], 1)[:, None]
            query = QT[rows]
            rm = (query@A).toarray()/divisor-np.where(present[:, None], er, 0)
            bm = (np.asarray(query@total).ravel()[:, None]-(query@B).toarray())/divisor-np.where(present[:, None], eb, 0)
            pure = kernel(vscore, n, m, rm, bm, er, eb, present)
            mixed_r = np.where(present[:, None], (rm+masses[rows, :, 0])/2, masses[rows, :, 0])
            mixed_b = np.where(present[:, None], (bm+masses[rows, :, 1])/2, masses[rows, :, 1])
            mixed = kernel(vscore, n, m, mixed_r, mixed_b, er, eb, np.ones(len(rows), dtype=bool))
            enrichment = (query@D).toarray()/(divisor*np.maximum(n, 1)*np.maximum(m, 1))
            enriched = vscore*rank[rows]*(1+np.where(vsupported, enrichment, 0))
            emit('rule_'+mode, rows, pure)
            emit('rule_mixed_'+mode, rows, mixed)
            emit('rule_enriched_'+mode, rows, enriched)
        del tokens, vocabulary, X, A, B, QT, D
    return dict(queries=Q, targets=T, cases=Q*T, numeric_state_bytes=required,
                native_batches=batches, native_chunk_size=chunk_size)
