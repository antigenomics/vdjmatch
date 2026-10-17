"""Fixed or fitted combinations of existing evidence, not binding confidence."""
import numpy as np

FEATURES=('likelihood_dense','likelihood_sparse','rank_dense','rank_sparse',
          'paired','likelihood_missing','rank_missing')


def reference_sparsity(frame):
    """Shared gate for fixed, fitted and null-normalized evidence."""
    import polars as pl
    structure=np.column_stack([frame[c].cast(pl.Float64).to_numpy() for c in
                               ('reference_coverage','knee_shape','knee_radius_ratio','original_support')])
    if not np.isfinite(structure).all() or np.any((structure<0)|(structure>1)):
        raise ValueError('reference structure and original support must lie in [0,1]')
    coverage,knee,radius,support=structure.T
    return 1-support*coverage*(1-(1-knee)*radius)


def evidence_features(frame):
    """One fixed formula; structure is reference-only and scale-independent.

    Coverage, normalized knee shape/radius and original support lie in [0,1]. Missing component scores
    have an explicit indicator; no row is discarded. An entirely unsupported
    receptor remains unavailable at prediction, rather than receiving an intercept.
    """
    import polars as pl
    values=[frame[c].cast(pl.Float64).to_numpy() for c in
        ('likelihood_score','conditional_rank_score')]
    missing=np.column_stack([~np.isfinite(x) for x in values])
    if np.any(values[1][np.isfinite(values[1])]<=0):
        raise ValueError('conditional rank evidence must be positive')
    likelihood=np.where(missing[:,0],0.,values[0])
    rank=np.log(np.where(missing[:,1],1.,values[1]))
    sparse=reference_sparsity(frame)
    paired=frame['paired'].cast(pl.Boolean).to_numpy()
    if frame['paired'].null_count():raise ValueError('paired mode must be explicit')
    return np.column_stack(((1-sparse)*likelihood,sparse*likelihood,
        (1-sparse)*rank,sparse*rank,paired,missing)),~missing.all(axis=1)


def fixed_evidence_combination(frame):
    """Full original log score refined by conditional rank using one unfitted gate.

    S=L+s log(B), s=1-u*c*[1-(1-k)*e], where u is original
    chain support, c reference radius coverage, k normalized knee shape, and
    e=epsilon/(epsilon+matrix-relative radius). Missing components use the
    available score; neither available returns NaN. This is an empirical ranking
    formula, not an independence claim or a transferred likelihood/P/E bound.
    """
    x,available=evidence_features(frame)
    result=x[:,0]+x[:,1]+x[:,3]
    missing_l,missing_r=x[:,-2:].T.astype(bool)
    result[missing_l]=x[missing_l,2]+x[missing_l,3]
    result[missing_r]=x[missing_r,0]+x[missing_r,1]
    result[~available]=np.nan
    return result


def fit_evidence_combination(frame, labels, *, weights=None):
    """Fit global coefficients with an analytic gradient; no parameter sweep.

    Caller owns disjoint training/evaluation identities and group weighting.
    Evidence coefficients are nonnegative, so stronger evidence cannot lower the
    score at fixed reference structure. A ridge penalty of one training-row unit
    stabilizes correlated columns; its strength is not tuned on benchmark labels.
    Returned linear margins require complete-pipeline null calibration for P/E.
    """
    from scipy.optimize import minimize
    from scipy.special import expit

    x,available=evidence_features(frame);y=np.asarray(labels,dtype=float)
    if y.shape!=(len(x),) or not np.isin(y,[0,1]).all() or len(np.unique(y))!=2:
        raise ValueError('aligned two-class binary training labels required')
    if not available.all():raise ValueError('unsupported training rows require an explicit disposition')
    w=np.ones(len(x)) if weights is None else np.asarray(weights,dtype=float)
    if w.shape!=y.shape or not np.isfinite(w).all() or np.any(w<=0):
        raise ValueError('aligned finite positive training weights required')
    w=w/w.sum();center=w@x;scale=np.sqrt(w@((x-center)**2));scale[scale<1e-12]=1
    design=np.column_stack((np.ones(len(x)),(x-center)/scale));ridge=1/len(x)
    def objective(beta):
        margin=design@beta
        loss=w@(np.logaddexp(0,margin)-y*margin)+ridge*(beta[1:]@beta[1:])/2
        gradient=design.T@(w*(expit(margin)-y))
        gradient[1:]+=ridge*beta[1:]
        return loss,gradient
    fitted=minimize(objective,np.zeros(design.shape[1]),jac=True,method='L-BFGS-B',
        bounds=[(None,None)]+[(0,None)]*4+[(None,None)]*3,
        options={'maxiter':500,'ftol':1e-12,'gtol':1e-7})
    if not fitted.success or not np.isfinite(fitted.x).all():
        raise RuntimeError('evidence combination did not converge: '+str(fitted.message))
    return {'schema':'evidence-combination-v2','features':list(FEATURES),
        'center':center.tolist(),'scale':scale.tolist(),'coefficients':fitted.x.tolist(),
        'training_rows':len(x),'positive_rows':int(y.sum()),'ridge':ridge,
        'objective':float(fitted.fun),'interpretation':'linear ranking margin; not calibrated binding probability'}


def combine_evidence(frame, model):
    """Apply the frozen formula; unsupported receptors return NaN."""
    if model.get('schema')!='evidence-combination-v2' or model.get('features')!=list(FEATURES):
        raise ValueError('incompatible evidence-combination model')
    x,available=evidence_features(frame)
    center=np.asarray(model['center']);scale=np.asarray(model['scale']);beta=np.asarray(model['coefficients'])
    if center.shape!=(len(FEATURES),) or scale.shape!=center.shape or beta.shape!=(len(FEATURES)+1,) or not np.isfinite(np.r_[center,scale,beta]).all() or np.any(scale<=0) or np.any(beta[1:5]<0):
        raise ValueError('invalid evidence-combination parameters')
    result=beta[0]+((x-center)/scale)@beta[1:];result[~available]=np.nan
    return result


THREE_COMPONENTS=('original_score','gapblock_score','stem4_score')
THREE_WEIGHTS=('original_weight','gapblock_weight','stem4_weight')


def weighted_component_support(frame):
    """Shared three-channel call support; missing flags preserve legacy behavior.

    A declared original_supported flag distinguishes actual original density
    support from its finite zero-density prior baseline. Positive weighted kernel
    evidence can independently support an axis. This never redistributes weights.
    """
    import polars as pl
    if 'original_supported' not in frame.columns:
        return pl.any_horizontal(*[pl.col(score).is_not_null() & (pl.col(weight) > 0)
            for score, weight in zip(THREE_COMPONENTS, THREE_WEIGHTS)])
    if frame.schema['original_supported'] not in (pl.Boolean, pl.Null):
        raise ValueError('original_supported must be Boolean or null')
    return ((pl.col('original_supported').fill_null(False) & pl.col('original_score').is_not_null() &
             (pl.col('original_weight') > 0)) |
            (pl.col('gapblock_score').cast(pl.Float64).gt(0).fill_null(False) & (pl.col('gapblock_weight') > 0)) |
            (pl.col('stem4_score').cast(pl.Float64).gt(0).fill_null(False) & (pl.col('stem4_weight') > 0)))


def null_normalized_evidence(frame, nulls, *, structure_conditioned=False, three_channel=False):
    """Put original/rank evidence on matching empirical null-tail scales.

    S=-log(min(p_L,p_R)) selects the stronger normalized component. It is a
    ranking statistic, not a final P-value: calibrate this complete selection.
    Ties, no-hit nulls and conservative source maxima use the existing helper.
    Each full-pMHC/species/locus group needs its own frozen null population.
    An optional non-null ``normalization_stratum`` in both tables additionally
    separates matching strata, without changing behavior when it is absent.
    With structure_conditioned, S=max(-log(p_L), -log(p_R)+log(s)), where
    s=1-c*[1-(1-k)*e] uses reference structure only. Query-local support
    does not switch routes; it already contributes to the component scores. The original tail
    retains unit weight; dense structure discounts the alternative route. A
    missing original component uses the available rank tail without discount.

    Three-channel mode instead returns sum(w_k * -log(p_k)), with explicit
    nonnegative weights summing to one. Missing components contribute zero,
    without redistributing their weights. Entirely unsupported weighted evidence
    remains unavailable. The caller supplies reference-derived weights; this
    function never fits them or reads assay labels. The final sum needs its own
    selected-label/paired calibration.
    """
    import polars as pl
    from ..aggregate.candidates import PMHC
    from ..evalue.directional import null_tail_probabilities, _null_group_axes
    keys=['species','locus',*PMHC]
    stratified='normalization_stratum' in frame.columns
    if stratified != ('normalization_stratum' in nulls.columns):
        raise ValueError('normalization_stratum must be present in both evidence and nulls')
    if stratified:
        if frame['normalization_stratum'].null_count() or nulls['normalization_stratum'].null_count():
            raise ValueError('normalization strata must be non-null')
        if frame.select('normalization_stratum').unique().join(
                nulls.select('normalization_stratum').unique(), on='normalization_stratum', how='anti').height:
            raise ValueError('missing matching normalization stratum')
        keys.append('normalization_stratum')
    if three_channel and structure_conditioned:
        raise ValueError('three-channel weights cannot use the two-channel structure gate')
    columns=list(THREE_COMPONENTS) if three_channel else ['likelihood_score','conditional_rank_score']
    tail_columns=['p_component_'+c.removesuffix('_score') for c in columns] if three_channel else ['p_component_likelihood','p_component_rank']
    required=set(keys+['query_id']+columns)
    if structure_conditioned:
        required.update(['reference_coverage','knee_shape','knee_radius_ratio'])
    if three_channel:
        required.update(THREE_WEIGHTS)
    for part in (frame,nulls):
        needed=required if part is frame else set(keys+['query_id']+columns)
        if missing:=needed-set(part.columns):
            raise ValueError('missing columns: '+', '.join(sorted(missing)))
        if part['query_id'].null_count():raise ValueError('query IDs must be non-null')
        for column in columns:
            values=part[column].cast(pl.Float64).drop_nulls()
            if not values.is_finite().all():raise ValueError('component scores must be finite or null')
        if not three_channel and (part['conditional_rank_score'].cast(pl.Float64).drop_nulls()<=0).any():
            raise ValueError('conditional rank evidence must be positive')
        if part.select(keys+['query_id']).n_unique()!=part.height:
            raise ValueError('query IDs must be unique within each calibration group')
    if not nulls.height:raise ValueError('nonempty null evidence required')
    if not {'outcome','null_source'}<=set(nulls.columns):raise ValueError('null evidence requires outcome and null_source')
    if nulls['outcome'].null_count() or not nulls['outcome'].eq('null').all():raise ValueError('only null outcomes define normalization')
    if nulls['null_source'].null_count() or nulls['null_source'].cast(pl.String).str.strip_chars().eq('').any():raise ValueError('null sources must be nonempty')
    channel_names=[c.removesuffix('_score') for c in columns]
    details=[prefix+name for prefix in ('log_evidence_','contribution_') for name in channel_names] if three_channel else []
    if {'_norm_row','normalized_score',*tail_columns,*details} & set(frame.columns):
        raise ValueError('normalization output columns already exist')
    if not frame.height:
        return frame.with_columns(*[pl.lit(None,dtype=pl.Float64).alias(c)
            for c in [*tail_columns,*details,'normalized_score']])
    values=frame.select(pl.col(columns).cast(pl.Float64)).to_numpy()
    available=np.isfinite(values).any(axis=1)
    if three_channel:
        weights=frame.select(pl.col(THREE_WEIGHTS).cast(pl.Float64)).to_numpy()
        if not np.isfinite(weights).all() or np.any(weights<0) or not np.allclose(weights.sum(axis=1),1.,rtol=0.,atol=1e-12):
            raise ValueError('three-channel weights must be finite, nonnegative and sum to one')
        available=(np.isfinite(values) & (weights>0)).any(axis=1)
    query_keys,null_keys=_null_group_axes(frame,nulls,keys)
    covered=query_keys.join(null_keys.unique().with_columns(pl.lit(True).alias('_has_null')),
        on=keys,how='left',nulls_equal=True,validate='m:1',maintain_order='left')['_has_null'].fill_null(False).to_numpy()
    if np.any(available & ~covered):
        raise ValueError('missing null calibration group')
    tails=[]
    for column in columns:
        observed=frame.select('query_id',*keys,
            pl.when(pl.Series(covered)).then(pl.col(column)).otherwise(None).alias(column))
        calibrated,_=null_tail_probabilities(observed,
            nulls.select('query_id',*keys,column,'outcome','null_source'),score_column=column,
            group_columns=keys,_details=False)
        tails.append(calibrated['p_null'].to_numpy())
    if three_channel:
        score=np.sum(weights * -np.log(np.column_stack(tails)),axis=1)
    elif structure_conditioned:
        sparse=reference_sparsity(frame.with_columns(pl.lit(1.).alias('original_support')))
        rank_tail=np.full(frame.height,-np.inf)
        np.log(sparse,out=rank_tail,where=sparse>0)
        rank_tail-=np.log(tails[1])
        missing_l=~np.isfinite(values[:,0])
        rank_tail[missing_l]=-np.log(tails[1][missing_l])
        score=np.maximum(-np.log(tails[0]),rank_tail)
    else:
        score=-np.log(np.minimum(*tails))
    score[~available]=np.nan
    extra=[]
    if three_channel:
        for index,name in enumerate(channel_names):
            evidence=-np.log(tails[index])
            extra.extend([pl.Series('log_evidence_'+name,evidence),
                pl.Series('contribution_'+name,weights[:,index]*evidence)])
    return frame.with_columns(*extra,*[pl.Series(name,tail) for name,tail in zip(tail_columns,tails)],
        pl.Series('normalized_score',score,nan_to_null=True))
