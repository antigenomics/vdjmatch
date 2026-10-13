"""Default historical prior retains the accepted scalar operation order."""
import polars as pl
from vdjmatch.match.historical import germline_prior


def test_default_prior_preserves_accepted_scalar_denominator_bits(monkeypatch):
    # Synthetic count fixture: V reference10/background0, J reference3/background0;
    # populations345/300000. Column-valued denominators changed the final ULP.
    reference=pl.DataFrame({'cdr3':['junction'+str(i) for i in range(345)],
        'v':['TRBV19']*10+['TRBV7-9']*335,
        'j':['TRBJ1-1']*3+['TRBJ2-7']*342})
    background=pl.DataFrame({'cdr3':['background']*300000,
        'v':['TRBV7-9']*300000,'j':['TRBJ2-7']*300000})
    query=pl.DataFrame({'query_id':[0,1],'cdr3':['query']*2,'v':['TRBV19']*2,'j':['TRBJ1-1']*2})
    # Real multi-task joins materialize population columns. Reproduce that storage
    # condition; literal broadcast columns otherwise mask this arithmetic defect.
    join=pl.DataFrame.join
    def materialized_join(frame,*args,**kwargs):
        out=join(frame,*args,**kwargs)
        sizes=[pl.Series(c,out[c].to_list(),dtype=pl.UInt64)
               for c in ['n_reference','background_size'] if c in out.columns]
        return out.with_columns(sizes) if sizes else out
    monkeypatch.setattr(pl.DataFrame,'join',materialized_join)
    score=germline_prior(query,reference,background,include_length=False)
    assert score['germline_lr'].to_list()==[18.359856032115257]*2
    assert score['n_reference'][0]==345 and score['background_size'][0]==300000
