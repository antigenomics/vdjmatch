"""Experimental continuous reference evidence; no assay labels or binder posterior."""
import numpy as np
import polars as pl
from scipy.special import ndtr


def reference_likelihood(queries, target, background, *, paired=False):
    """Common-covariance Gaussian likelihood of existing, correlated score features.

    Each query removes every full-junction match from both fitted populations.
    Joint moments capture correlations; missing features select a lower-dimensional
    joint model. Covariance is the pooled, unlabelled sample covariance, with its
    singular directions projected out. This uses reference annotations only.
    The normal tail is valid for the fitted Gaussian background model, not an
    empirical biological null. Full-pipeline null-call validation remains required.
    """
    fields = (['original_alpha_density_score', 'original_beta_density_score', 'original_germline_lr']
              if paired else ['original_score']) + ['conditional_score']
    keys = ['junction_alpha', 'junction_beta'] if paired else ['junction_aa']
    outputs = {'score', 'log_likelihood_ratio', 'p_gaussian_model', 'likelihood_status',
               'model_target_count', 'model_background_count', 'model_dimensions'}
    if outputs & set(queries.columns):
        raise ValueError('likelihood output columns already exist')
    tables = [('queries', queries), ('target', target), ('background', background)]
    for name, table in tables:
        if missing := {'query_id', *fields, *keys} - set(table.columns):
            raise ValueError(name+' missing columns: '+', '.join(sorted(missing)))
        if table['query_id'].null_count() or table['query_id'].n_unique() != table.height:
            raise ValueError(name+' query_id must be unique and non-null')
        if not table.select(pl.col(fields).cast(pl.Float64).is_finite() | pl.col(fields).is_null()).to_numpy().all():
            raise ValueError('scores must be finite or null')
    for field in ['species', 'locus', 'epitope', 'mhc_a', 'mhc_b', 'mhc_class', 'reference_sha256', 'scoring_sha256']:
        if any(field in t.columns for _, t in tables):
            if any(field not in t.columns for _, t in tables):
                raise ValueError('context column required in every table: '+field)
            values = [t[field].unique().to_list() for _, t in tables]
            if any(len(v) != 1 or v != values[0] for v in values):
                raise ValueError('one matching target/scorer context required: '+field)
    if not target.height or not background.height:
        raise ValueError('nonempty target and background reference populations required')

    def features(table):
        x = table.select(pl.col(fields).cast(pl.Float64)).to_numpy().copy()
        for i, field in enumerate(fields):
            if field in ('original_germline_lr',):
                continue
            if field == 'original_score':
                x[:, i] = np.arcsinh(x[:, i])
            else:
                if (x[:, i] < 0).any():
                    raise ValueError('density and conditional evidence must be nonnegative')
                x[:, i] = np.log1p(x[:, i])
        return x

    qx, tx, bx = [features(t) for _, t in tables]
    # Translation leaves the likelihood unchanged and avoids subtracting huge moments.
    for i in range(len(fields)):
        values = np.concatenate([tx[:, i], bx[:, i]])
        values = values[np.isfinite(values)]
        if len(values):
            center = values.mean()
            for x in [qx, tx, bx]:
                x[:, i] -= center
    # At most 16 availability patterns for four components; no per-query fitting loop.
    pattern = np.isfinite(qx) @ (1 << np.arange(len(fields)))
    valid_identity = queries.select(pl.all_horizontal(
        pl.col(keys).str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False))).to_series().to_numpy()
    score = np.full(queries.height, np.nan)
    p = np.ones(queries.height)
    counts = np.zeros((queries.height, 2), dtype=np.int64)
    dimensions = np.zeros(queries.height, dtype=np.int8)
    status = np.full(queries.height, 'no_evidence', dtype=object)
    for mask in np.unique(pattern[valid_identity]):
        if not mask:
            continue
        selected = np.flatnonzero((pattern == mask) & valid_identity)
        axes = [i for i in range(len(fields)) if mask & (1 << i)]
        d = len(axes)
        xy = qx[np.ix_(selected, axes)]
        fitted = []
        for table, matrix in [(target, tx), (background, bx)]:
            x = matrix[:, axes]
            valid = np.isfinite(x).all(axis=1)
            valid &= table.select(pl.all_horizontal(pl.col(keys).is_not_null())).to_series().to_numpy()
            x = x[valid]
            # Sufficient moments allow full-clonotype deletion without refitting each row.
            names = ['_x'+str(i) for i in range(d)]
            cross = ['_xx'+str(i)+'_'+str(j) for i in range(d) for j in range(d)]
            moments = table.filter(pl.Series(valid)).select(keys).with_columns(
                *[pl.Series(name, x[:, i]) for i, name in enumerate(names)],
                *[pl.Series(cross[i*d+j], x[:, i]*x[:, j]) for i in range(d) for j in range(d)])
            removed = moments.group_by(keys).agg(pl.len().alias('_n'), pl.col(names+cross).sum())
            removed = queries[selected.tolist()].select(keys).join(
                removed, on=keys, how='left', validate='m:1', maintain_order='left')
            n = len(x)-removed['_n'].fill_null(0).to_numpy()
            sums = x.sum(axis=0)-removed.select(pl.col(names).fill_null(0)).to_numpy()
            squares = (x.T @ x)[None, :, :]-removed.select(pl.col(cross).fill_null(0)).to_numpy().reshape(-1,d,d)
            fitted.append((n, sums, squares))
        (nt, st, sst), (nb, sb, ssb) = fitted
        counts[selected] = np.column_stack([nt, nb])
        usable = (nt > 0) & (nb > 0)
        status[selected] = np.where(nt == 0, 'no_target_model', np.where(nb == 0, 'no_background_model', 'available'))
        if not usable.any():
            continue
        ids = selected[usable]
        nt, nb, st, sb, sst, ssb, xy = [v[usable] for v in [nt, nb, st, sb, sst, ssb, xy]]
        mt, mb = st/nt[:, None], sb/nb[:, None]
        total, n = st+sb, nt+nb
        cov = (sst+ssb-total[:, :, None]*total[:, None, :]/n[:, None, None])/(n-1)[:, None, None]
        if not np.isfinite(cov).all():
            raise ValueError('reference moments exceed numerical range')
        precision = np.linalg.pinv(cov, hermitian=True)
        delta = mt-mb
        w = np.einsum('nij,nj->ni', precision, delta)
        lr = np.einsum('ni,ni->n', xy-(mt+mb)/2, w)
        separation = np.maximum(0., np.einsum('ni,ni->n', delta, w))
        if not np.isfinite(lr).all() or not np.isfinite(separation).all():
            raise ValueError('likelihood exceeds numerical range')
        score[ids] = lr
        informative = separation > 0
        p[ids[informative]] = ndtr(-(lr[informative]+separation[informative]/2)/np.sqrt(separation[informative]))
        dimensions[ids] = d
        status[ids[~informative]] = 'indistinguishable_models'
    return queries.with_columns(
        pl.Series('score', score, nan_to_null=True), pl.Series('log_likelihood_ratio', score, nan_to_null=True),
        pl.Series('p_gaussian_model', p), pl.Series('likelihood_status', status),
        pl.Series('model_target_count', counts[:, 0]), pl.Series('model_background_count', counts[:, 1]),
        pl.Series('model_dimensions', dimensions))
