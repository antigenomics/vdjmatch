"""All-label fixed rules using bounded native counts and sparse motif projections."""
import numpy as np
import polars as pl
from scipy import sparse
from scipy.special import digamma


def stream_reference_rules(queries, reference, *, paired=False, species='human', threads=1,
                           radius=90, consume, identity_keys=()):
    """Emit the scalar rule scores without retaining a detailed case grid.

    Inputs use resolved ``seq0/v0 [seq1/v1]`` and reference epitope. Query IDs
    retain their source identity. Full models, priors and motif counts are fixed;
    only exact hits are filtered. MHC pooling is the caller's explicit choice.
    Optional complete identity columns preserve original keys that share geometry.
    Three query-by-label numeric arrays retain spatial masses/ranks in memory;
    one motif mode is prepared at a time. Nothing is persisted or resumed.
    """
    from .epitope_rules import _prepare_rule_frame, _tokens, MOTIFS
    from .local_rank import _harmonic_interval
    from ..match.geometry import (grouped_total_distance_count_batch, distance_units,
                                  distance_cutoff, precompute_v_loop_distances)

    axes = range(2 if paired else 1)
    seq, v = [f'seq{i}' for i in axes], [f'v{i}' for i in axes]
    keys = [*seq, *v, *identity_keys]
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
        return dict(queries=0, reference_keys=R, targets=T, cases=0, numeric_state_bytes=0,
                    native_batches=0, native_chunk_size=0)
    gene_keys = lambda f: f.select(*[pl.col(c).str.replace(r'\*.*$', '') for c in v])
    genes = gene_keys(whole).unique().sort(v).with_row_index('_gene')
    G = genes.height
    rg = gene_keys(whole).join(genes, on=v, validate='m:1', maintain_order='left')['_gene'].to_numpy()
    qg = gene_keys(queries).join(genes, on=v, how='left', validate='m:1', maintain_order='left')['_gene'].fill_null(G).to_numpy()
    required = 3*Q*T*8+2*(G+1)*T*8+256*1024**2
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
    # One target per peptide: native whole/exclusive channels already supply
    # each fixed complement, without artificial duplicate target models.
    offsets = W.indptr
    exclusive_peptides = np.full(R, -1, dtype=np.int64)
    exclusive_peptides[exclusive['_ref'].to_numpy()] = exclusive['_target'].to_numpy()
    context = precompute_v_loop_distances(species=species)
    scale, quantum = distance_units()
    bound = 3*scale*sum(max((f[c].str.len_chars().max() or 0) for f in (queries, whole)) for c in seq)
    cuts = list(range(0, 200*bound+1, 600*quantum))
    native_radius = distance_cutoff(radius)
    weights = np.maximum(0, native_radius-np.asarray(cuts, dtype=np.int64))
    spatial_bins = int(np.count_nonzero(weights))
    # Sum CDF increments times weights by telescoping; zero-weight bins vanish.
    spatial_steps = weights[:spatial_bins]-np.r_[weights[1:spatial_bins], np.int64(0)]
    # Keep native counts under128MiB; reduce labels in windows instead of
    # materializing every query/label/distance float scratch array together.
    per_query = 4*((1+3*T)*(len(cuts)+1)+T)
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
            reference_target_ids=W.indices, target_peptide_ids=list(range(T)),
            reference_exclusive_peptide_ids=exclusive_peptides)
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
                bc = counts[:, 0, None, :].astype(np.int64)-counts[:, 1+2*T+targets, :]
                if (not np.array_equal(tc[:, :, -1], (n-er)[:, targets]) or
                        not np.array_equal(bc[:, :, -1], (m-eb)[:, targets])):
                    raise ArithmeticError('all-label rank bound did not cover eligible hits')
                if conditional:
                    masses[np.ix_(rows, targets, [0])] = (tc[:, :, :spatial_bins]@spatial_steps/native_radius)[:, :, None]
                    masses[np.ix_(rows, targets, [1])] = (bc[:, :, :spatial_bins]@spatial_steps/native_radius)[:, :, None]
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
    # Geometry may coincide after imputation while original source keys differ.
    # Preserve the already indexed reference axis through motif preparation too.
    prepared_r = (_prepare_rule_frame(whole.rename({'_ref': 'query_id'}), axes, True)
        if identity_keys else _prepare_rule_frame(whole.select(keys), axes))
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
    return dict(queries=Q, reference_keys=R, targets=T, cases=Q*T, numeric_state_bytes=required,
                native_batches=batches, native_chunk_size=chunk_size)


def stream_reference_kernels(queries,reference,*,species='human',threads=1,paired=False,consume,identity_keys=()):
    """Emit the six unchanged stem4 kernel alternatives from one native traversal.

    MHC is explicitly pooled by peptide by the owning reference classifier.
    Full reference populations remain fixed; only exact junction hits are masked.
    """
    from .kernel import grouped_kernel_selection,ladder_cutoffs,KERNEL_SCORERS
    from ..match.geometry import precompute_v_loop_distances
    query_names=({'seq0':'cdr3a','v0':'_allelea','seq1':'cdr3b','v1':'_alleleb'} if paired else {'seq0':'cdr3','v0':'_allele'})
    reference_names=({'seq0':'alpha','v0':'_allele_alpha','seq1':'beta','v1':'_allele_beta'} if paired else {'seq0':'cdr3','v0':'_allele'})
    q=queries.rename(query_names);r=reference.rename(reference_names)
    query_keys=[*query_names.values(),*identity_keys];search_keys=[*reference_names.values(),*identity_keys]
    context=precompute_v_loop_distances(species=species,matrix='blosum62')
    scale=context['distance_units'][0]
    groups=[(dict(epitope=epitope,mhc_a=None,mhc_b=None,mhc_class=None),q,r.filter(pl.col('epitope')==epitope)) for epitope in sorted(r['epitope'].unique())]
    rows=0
    def emit(frame):
        nonlocal rows
        rows+=frame.height
        for scorer in KERNEL_SCORERS:
            name='kernel_'+scorer.replace('-','_')
            consume(name,frame.select('query_id','epitope',pl.col(name).alias('score'),(pl.col(name).is_finite()&(pl.col(name)>0)).fill_null(False).alias('supported')))
    grouped_kernel_selection(q,r,groups,query_keys=query_keys,search_keys=search_keys,paired=paired,
        cutoffs=ladder_cutoffs(scale),radius=22.5*scale,kernel_scale=22.5*scale,scorer='v-bounded',
        _detailed=True,_consume=emit,progress=True,species=species,threads=threads,ctrim=3,matrix='blosum62',
        stem_depth=4,stem_weight=1.,v_loop_weight=1,exclude_exact=True,v_loop_distances=context)
    return dict(native_batches=1,query_keys=q.select(query_keys).unique().height,
        reference_keys=r.select(search_keys).unique().height,query_target_rows=rows,
        geometry='current stem4 BLOSUM62 trim3; fixed VDJdb-other background; six native kernel alternatives',MHC='pooled by peptide')


def stream_reference_policy(queries, reference, policy, *, species='human', locus='TRB',
                            threads=1, paired=False, consume, identity_keys=()):
    """Route fresh full-population rules/kernel scores through the frozen policy.

    Only routed V/exact3 columns are retained, indexed by original query identity.
    Every target, including missing-profile tails, reaches the policy and winner
    consumer. Default routes preserve the kernel's original support flag: finite
    zero kernel scores remain unsupported. No component probability is transferred.
    """
    from ..cli.score_policy import apply_reference_policy, FAST_ROUTES
    if not callable(consume) or not set(policy['reference_class']) <= set(FAST_ROUTES):
        raise ValueError('a consumer and frozen FAST_ROUTES policy classes required')
    organism = {'human':'HomoSapiens','mouse':'MusMusculus'}.get(species,species)
    lane = 'PAIRED' if paired else locus
    labels = reference.select('epitope').unique().sort('epitope')
    default = 'kernel_v_bounded'
    # The existing selector, rather than a second routing implementation, defines
    # which raw component columns are needed for this reference population.
    probe = labels.with_columns(pl.lit(organism).alias('species'),pl.lit(lane).alias('locus'),
        pl.lit(0.).alias(default),pl.lit(1.).alias('rule_v_score'),
        pl.lit(1.).alias('rule_exact3_score'),pl.lit(None,dtype=pl.Float64).alias('rule_rank_score'))
    plan = apply_reference_policy(probe,policy,original_column=default)
    query_rows = queries.select('query_id').with_row_index('_policy_row')
    target_index = dict(zip(labels['epitope'],range(labels.height)))
    arrays,seen,columns = {},{},{}
    for scorer,column in [('rule_v','rule_v_score'),('rule_exact3','rule_exact3_score')]:
        needed = plan.filter(pl.col('selected_score_column').is_in([column,'rule_exact3_score / rule_v_score']))['epitope']
        mapping = dict(zip(needed,range(len(needed))))
        columns[scorer] = np.array([mapping.get(label,-1) for label in labels['epitope']],dtype=np.int64)
        arrays[scorer] = np.full((queries.height,len(needed)),np.nan)
        seen[scorer] = np.zeros((queries.height,len(needed)),dtype=bool)
    def coordinates(frame):
        if frame.select('query_id','epitope').is_duplicated().any():
            raise ValueError('duplicate class-policy query/target component coordinates')
        frame = frame.join(query_rows,on='query_id',how='left',validate='m:1',maintain_order='left')
        if frame['_policy_row'].null_count():raise ValueError('unknown class-policy query identity')
        target = frame['epitope'].replace_strict(target_index,return_dtype=pl.UInt32).to_numpy()
        return frame['_policy_row'].to_numpy(),target
    def retain(scorer,frame):
        if scorer not in arrays or arrays[scorer].shape[1]==0:return
        row,target = coordinates(frame)
        col = columns[scorer][target];needed = col>=0
        row,col = row[needed],col[needed]
        if seen[scorer][row,col].any():raise ValueError('repeated class-policy component coordinates')
        values = frame['score'].to_numpy()
        arrays[scorer][row,col] = np.where(frame['supported'].to_numpy(),values,np.nan)[needed]
        seen[scorer][row,col] = True
    rule_info = (stream_reference_rules(queries,reference,species=species,paired=paired,
        threads=threads,consume=retain,identity_keys=identity_keys) if any(a.size for a in arrays.values()) else
        dict(native_batches=0,reason='all policy routes retain configured native default'))
    if any(not present.all() for present in seen.values()):
        raise ValueError('incomplete class-policy rule component coordinates')
    count = 0
    kernel_seen = np.zeros((queries.height,(labels.height+7)//8),dtype=np.uint8)
    def select(scorer,frame):
        nonlocal count
        if scorer!=default:return
        row,target = coordinates(frame)
        flags = np.left_shift(np.uint8(1),(target%8).astype(np.uint8))
        if np.any(kernel_seen[row,target//8]&flags):
            raise ValueError('repeated class-policy kernel coordinates')
        np.bitwise_or.at(kernel_seen,(row,target//8),flags)
        values = {}
        for name,column in [('rule_v','rule_v_score'),('rule_exact3','rule_exact3_score')]:
            col = columns[name][target];needed = col>=0
            value = np.full(frame.height,np.nan)
            value[needed] = arrays[name][row[needed],col[needed]]
            values[column] = pl.Series(column,value,nan_to_null=True)
        source = frame.rename({'score':default,'supported':'_native_supported'}).with_columns(
            pl.lit(organism).alias('species'),pl.lit(lane).alias('locus'),
            *values.values(),pl.lit(None,dtype=pl.Float64).alias('rule_rank_score'))
        selected = apply_reference_policy(source,policy,original_column=default)
        count += frame.height
        consume('reference_class_policy',selected.select('query_id','epitope',
            pl.col('policy_score').alias('score'),
            (pl.col('policy_score').is_finite() & pl.when(pl.col('selected_score_column')==default)
                .then(pl.col('_native_supported')).otherwise(pl.lit(True))).fill_null(False).alias('supported'),
            'reference_class','selected_score_column','policy_status'))
    kernel_info = stream_reference_kernels(queries,reference,species=species,paired=paired,
        threads=threads,consume=select,identity_keys=identity_keys)
    if count!=queries.height*labels.height:
        raise ValueError('incomplete class-policy kernel coordinates')
    return dict(rules=rule_info,kernel=kernel_info,query_target_rows=count,
        retained_component_bytes=sum(a.nbytes for a in arrays.values())+sum(a.nbytes for a in seen.values())+kernel_seen.nbytes,
        original_column=default,policy='unchanged FAST_ROUTES/apply_reference_policy',
        candidate_models='all valid fixed reference labels including tails',calibration='unavailable')
