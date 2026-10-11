"""Small global logistic combination of existing evidence, not binding confidence."""
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

FEATURES=('likelihood_dense','likelihood_sparse','rank_dense','rank_sparse',
          'paired','likelihood_missing','rank_missing')


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
    structure=np.column_stack([frame[c].cast(pl.Float64).to_numpy() for c in
                               ('reference_coverage','knee_shape','knee_radius_ratio','original_support')])
    if not np.isfinite(structure).all() or np.any((structure<0)|(structure>1)):
        raise ValueError('reference structure and original support must lie in [0,1]')
    coverage,knee,radius,support=structure.T
    # Coverage handles flat curves; knee size and shape refine the covered part.
    sparse=1-support*coverage*(1-(1-knee)*radius)
    paired=frame['paired'].cast(pl.Boolean).to_numpy()
    if frame['paired'].null_count():raise ValueError('paired mode must be explicit')
    return np.column_stack(((1-sparse)*likelihood,sparse*likelihood,
        (1-sparse)*rank,sparse*rank,paired,missing)),~missing.all(axis=1)


def fixed_evidence_combination(frame):
    """Original log score completed by conditional rank using one unfitted gate.

    S=(1-s)L+s log(B), s=1-u*c*[1-(1-k)*e], where u is original
    chain support, c reference radius coverage, k normalized knee shape, and
    e=epsilon/(epsilon+matrix-relative radius). Missing components use the
    available score; neither available returns NaN. This is an empirical ranking
    formula, not an independence claim or a transferred likelihood/P/E bound.
    """
    x,available=evidence_features(frame)
    result=x[:,0]+x[:,3]
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
