"""Small global logistic combination of existing evidence, not binding confidence."""
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

FEATURES=('likelihood','rank','motif','sparse_rank','knee_rank','sparse_motif',
          'knee_motif','paired','likelihood_missing','rank_missing','motif_missing')


def evidence_features(frame):
    """One fixed formula; structure is reference-only and scale-independent.

    Coverage and normalized knee shape lie in [0,1]. Missing component scores
    have an explicit indicator; no row is discarded. An entirely unsupported
    receptor remains unavailable at prediction, rather than receiving an intercept.
    """
    import polars as pl
    values=[frame[c].cast(pl.Float64).to_numpy() for c in
        ('likelihood_score','conditional_rank_score','motif_score')]
    missing=np.column_stack([~np.isfinite(x) for x in values])
    if any(np.any(x[np.isfinite(x)]<0) for x in values[1:]):
        raise ValueError('rank and motif evidence must be nonnegative')
    likelihood,rank,motif=[np.where(np.isfinite(x),x,0.0) for x in values]
    rank=np.log1p(rank)
    structure=np.column_stack([frame[c].cast(pl.Float64).to_numpy() for c in
                               ('reference_coverage','knee_shape')])
    if not np.isfinite(structure).all() or np.any((structure<0)|(structure>1)):
        raise ValueError('reference coverage and knee shape must lie in [0,1]')
    sparse,knee=(1-structure).T
    paired=frame['paired'].cast(pl.Boolean).to_numpy()
    if frame['paired'].null_count():raise ValueError('paired mode must be explicit')
    return np.column_stack((likelihood,rank,motif,sparse*rank,knee*rank,
        sparse*motif,knee*motif,paired,missing)),~missing.all(axis=1)


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
        bounds=[(None,None)]+[(0,None)]*7+[(None,None)]*4,
        options={'maxiter':500,'ftol':1e-12,'gtol':1e-7})
    if not fitted.success or not np.isfinite(fitted.x).all():
        raise RuntimeError('evidence combination did not converge: '+str(fitted.message))
    return {'schema':'evidence-combination-v1','features':list(FEATURES),
        'center':center.tolist(),'scale':scale.tolist(),'coefficients':fitted.x.tolist(),
        'training_rows':len(x),'positive_rows':int(y.sum()),'ridge':ridge,
        'objective':float(fitted.fun),'interpretation':'linear ranking margin; not calibrated binding probability'}


def combine_evidence(frame, model):
    """Apply the frozen formula; unsupported receptors return NaN."""
    if model.get('schema')!='evidence-combination-v1' or model.get('features')!=list(FEATURES):
        raise ValueError('incompatible evidence-combination model')
    x,available=evidence_features(frame)
    center=np.asarray(model['center']);scale=np.asarray(model['scale']);beta=np.asarray(model['coefficients'])
    if center.shape!=(len(FEATURES),) or scale.shape!=center.shape or beta.shape!=(len(FEATURES)+1,) or not np.isfinite(np.r_[center,scale,beta]).all() or np.any(scale<=0) or np.any(beta[1:8]<0):
        raise ValueError('invalid evidence-combination parameters')
    result=beta[0]+((x-center)/scale)@beta[1:];result[~available]=np.nan
    return result
