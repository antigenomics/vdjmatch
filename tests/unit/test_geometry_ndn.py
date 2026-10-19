"""Fixed boundary-weight hypothesis, independent of learned rival coefficients."""
from pathlib import Path
import polars as pl
import pytest
from vdjmatch.match.geometry import total_distance_count_batch


def test_fixed_ndn_weights_and_validation():
    from seqtree import gapblock
    if 'query_position_weights' not in gapblock._gapblock_count_batch.__doc__:
        with pytest.raises(RuntimeError,match='native per-input'):
            total_distance_count_batch(['CASSLGQAYEQYF'],['CASSAGQAYEQYF'],['TRBV19*01'],['TRBV19*01'],[[2400]],
                position_weighting='ndn',query_boundaries=[(4,9)],reference_boundaries=[(4,9)])
        return
    q=['CASSLGQAYEQYF']
    r=['AASSLGQAYEQYF','CASSAGQAYEQYF','CASSLGQAYEQAF']
    v=['TRBV19*01']
    kw=dict(position_weighting='ndn',query_boundaries=[(4,9)],reference_boundaries=[(4,9)]*3)
    # These equal-length substitutions each have clipped BLOSUM cost4:
    # flanks contribute240 native units, NDN2400, with unchanged V loops.
    assert total_distance_count_batch(q,r,v,v*3,[[239,240,2399,2400]],**kw)==[[0,2,2,3]]
    # Unavailable boundary uses100 throughout, never a fabricated flank/core split.
    assert total_distance_count_batch(q,r,v,v*3,[[240,2400]],position_weighting='ndn',
        query_boundaries=[None],reference_boundaries=[None]*3)==[[0,3]]
    for bounds in ([(9,4)],[(0,14)],[],[(True,9)]):
        with pytest.raises((ValueError,TypeError)):
            total_distance_count_batch(q,r,v,v*3,[[2400]],position_weighting='ndn',
                query_boundaries=bounds,reference_boundaries=[(4,9)]*3)


def test_boundary_annotation_preserves_identity_and_conflicts(monkeypatch):
    from vdjmatch.match.regions import annotate_ndn_boundaries
    arda=pytest.importorskip("arda.cdr3fix")
    def markup(frame,**kw):
        assert frame.height==1 and kw['max_replace']==0
        return pl.DataFrame({'cdr3_repaired':frame['cdr3'],'v_end':[4],'j_start':[9],
            'v_call':frame['_allele'],'j_call':frame['_j']})
    monkeypatch.setattr(arda,'markup_batch',markup)
    seq='CASSLGQAYEQYF'
    frame=pl.DataFrame({'cdr3':[seq]*4,'_allele':['TRBV19*01']*4,
        'j':['TRBJ2-7*01','TRBJ2-7*01','TRBJ1-1*01','TRBJ1-1*01']})
    # Contradictory J calls on one junction/V key never choose a first call.
    result=annotate_ndn_boundaries(frame,species='human')
    assert result.height==4 and result['boundary_status'].to_list()==['conflicting_j_uniform']*4
    clean=annotate_ndn_boundaries(frame.head(2),species='human')
    assert clean['v_end'].to_list()==[4,4] and clean['j_start'].to_list()==[9,9]
    assert clean['boundary_status'].to_list()==['annotated']*2
    equivalent=frame.head(2).with_columns(pl.Series('j',['TRBJ2-7','TRBJ2-7*01']))
    assert annotate_ndn_boundaries(equivalent,species='human')['boundary_status'].to_list()==['annotated']*2


@pytest.mark.parametrize('v_end,j_start', [(None, 9), (4, None), (None, None)])
def test_null_boundaries_use_uniform_fallback(monkeypatch, v_end, j_start):
    from vdjmatch.match.regions import annotate_ndn_boundaries
    arda = pytest.importorskip('arda.cdr3fix')
    frame = pl.DataFrame({'cdr3': ['CASSLGQAYEQYF'], '_allele': ['TRBV19*01'],
                          'j': ['TRBJ2-7*01']})
    def markup(rows, **kwargs):
        return rows.select(pl.col('cdr3').alias('cdr3_repaired'),
                           pl.col('_allele').alias('v_call'), pl.col('_j').alias('j_call')) \
            .with_columns(pl.lit(v_end, dtype=pl.Int64).alias('v_end'),
                          pl.lit(j_start, dtype=pl.Int64).alias('j_start'))
    monkeypatch.setattr(arda, 'markup_batch', markup)
    result = annotate_ndn_boundaries(frame, species='human')
    assert result['boundary_status'].to_list() == ['unavailable_uniform']
    assert result['v_end'].to_list() == result['j_start'].to_list() == [None]


def test_ndn_cli_exports_boundaries_and_retains_unavailable_queries(tmp_path):
    import json
    pytest.importorskip("arda.cdr3fix")
    from seqtree import gapblock
    if 'query_position_weights' not in gapblock._gapblock_count_batch.__doc__:pytest.skip('native per-input geometry unavailable')
    from vdjmatch.cli.__main__ import main
    q=tmp_path/'q.tsv';r=tmp_path/'r.tsv';prefix=tmp_path/'out'
    meta={'species':'HomoSapiens','epitope':'E','mhc_a':'HLA-A*02:01','mhc_b':'B2M','mhc_class':'MHCI'}
    pl.DataFrame([dict(junction_aa='CASSLGQAYEQYF',v_call='TRBV19*01',j_call='TRBJ2-7*01',locus='TRB',**meta),
                  dict(junction_aa='CASSLGQAYEQYF',v_call='unknown',j_call='TRBJ2-7*01',locus='TRB',**meta)]).write_csv(q,separator='\t')
    pl.DataFrame([dict(cdr3='CASSAGQAYEQYF',v='TRBV19*01',j='TRBJ2-7*01',gene='TRB',**meta),
                  dict(cdr3='CASSLGQAYEQAF',v='TRBV19*01',j='TRBJ2-7*01',gene='TRB',**(meta|{'epitope':'N'}))]).write_csv(r,separator='\t')
    assert main(['neighbours',str(q),'--vdjdb',str(r),'--locus','TRB',
        '--targets-from-sample','--background','vdjdb-other','--position-weighting','ndn',
        '--neighbour-weighting','local-rank','--exclude-exact','--output-prefix',str(prefix)])==0
    manifest=json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert manifest['background_enrichment']['boundaries']['weights']=={'V_J':0.1,'NDN':1.0}
    assert sum(manifest['background_enrichment']['boundaries']['status_counts'].values())==4
    queries=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t')
    assert queries.height==2 and queries['available'].all()
    assert queries['v_geometry_imputed'].to_list()==[False,True]
    assert pl.read_csv(str(prefix)+'.enrichment.tsv',separator='\t').height==2
