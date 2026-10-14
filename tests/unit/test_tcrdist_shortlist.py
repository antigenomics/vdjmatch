"""Actual sparse CDR3 retrieval must retain every full-distance radius neighbour."""
import numpy as np
from vdjmatch.match.tcrdist import distance_matrix,shortlist_audit_batch

def test_actual_shortlist_matches_dense_radius_and_puncture():
    q=['CASSLGQAYEQYF','CASSLGQAPAYEQYF'];r=['CASSLGQAYEQYF','CASSAGQAYEQYF','CASSLGQAPAYEQYF']
    qv=['TRBV19*01']*2;rv=['TRBV7-9*01','TRBV19*01','TRBV7-9*01']
    full,cdr3=distance_matrix(q,r,qv,rv,return_cdr3=True,ctrim=3)
    for excluded in [False,True]:
        for radius in [0,12,90]:
            a=shortlist_audit_batch(q,r,qv,rv,radius=radius,ctrim=3,exclude_exact=excluded)
            eligible=np.array(q)[:,None]!=np.array(r)[None,:] if excluded else np.ones(full.shape,dtype=bool)
            shortlisted=(cdr3<=radius)&eligible;retained=(full<=radius)&eligible
            assert a['cdr3_shortlist_keys']==shortlisted.sum(axis=1).tolist()
            assert a['post_v_radius_keys']==retained.sum(axis=1).tolist()
            for i,nearest in enumerate(a['nearest_within_shortlist_distance']):
                assert nearest==(float(full[i,shortlisted[i]].min()) if shortlisted[i].any() else None)

def test_cli_audit_preserves_exhaustive_output_and_unavailable_rows(tmp_path):
    import json
    import polars as pl
    from vdjmatch.cli.__main__ import main
    meta=dict(species='HomoSapiens',epitope='E',mhc_a='HLA-A*02:01',mhc_b='B2M',mhc_class='MHCI')
    q=tmp_path/'q.tsv';r=tmp_path/'r.tsv'
    pl.DataFrame([dict(junction_aa='CASSLGQAYEQYF',v_call=v,locus='TRB',**meta) for v in ['TRBV19*01','unknown','TRBV19*01']]).write_csv(q,separator='\t')
    pl.DataFrame([dict(cdr3='CASSAGQAYEQYF',v='TRBV7-9*01',gene='TRB',**meta)]).write_csv(r,separator='\t')
    for audit in [False,True]:
        prefix=tmp_path/str(audit)
        args=['tcrdist-neighbours',str(q),'--vdjdb',str(r),'--locus','TRB','--targets-from-sample','--unbounded-nearest','--output-prefix',str(prefix)]
        assert main(args+(['--shortlist-audit'] if audit else []))==0
    for suffix in ['queries','candidates']:
        assert pl.read_csv(tmp_path/f'False.{suffix}.tsv',separator='\t').equals(pl.read_csv(tmp_path/f'True.{suffix}.tsv',separator='\t'))
    diagnostic=pl.read_csv(tmp_path/'True.shortlist.tsv',separator='\t');assert diagnostic.height==3
    assert diagnostic['cdr3_shortlist_keys'].to_list()==[1,None,1]
    assert diagnostic['status'].to_list()==['matched','invalid_query_or_gene','matched']
    manifest=json.loads((tmp_path/'True.manifest.json').read_text())
    assert manifest['shortlist_audit']['stages'][0]['query_appearances']==2
    assert manifest['shortlist_audit']['stages'][0]['unique_computations']==1

def test_sparse_budget_fails_before_retrieval_and_parallel_results_equal(monkeypatch):
    import pytest
    from vdjmatch.match import tcrdist
    q=['CASSLGQAYEQYF']*2;r=['CASSAGQAYEQYF'];qv=['TRBV19*01']*2;rv=['TRBV7-9*01']
    serial=shortlist_audit_batch(q,r,qv,rv,threads=1)
    parallel=shortlist_audit_batch(q,r,qv,rv,threads=2)
    for key in ['cdr3_shortlist_keys','post_v_radius_keys','nearest_within_shortlist_distance']:
        assert serial[key]==parallel[key]
    monkeypatch.delattr(tcrdist.gapblock,'group_audit_batch',raising=False)
    from functools import wraps
    @wraps(tcrdist.gapblock.count_batch)
    def huge_counts(*a,**k):return [[100_000_000]]*len(q)
    monkeypatch.setattr(tcrdist.gapblock,'count_batch',huge_counts)
    @wraps(tcrdist.gapblock.ball_batch)
    def forbidden_retrieval(*a,**k):pytest.fail('retrieval must not allocate over budget')
    monkeypatch.setattr(tcrdist.gapblock,'ball_batch',forbidden_retrieval)
    with pytest.raises(MemoryError,match='no hits truncated'):
        shortlist_audit_batch(q,r,qv,rv)


def test_shortlist_requires_actual_native_geometry(monkeypatch):
    import pytest
    from vdjmatch.match import tcrdist
    name='_gapblock_group_audit_batch' if hasattr(tcrdist.gapblock,'group_audit_batch') else '_gapblock_ball_batch'
    monkeypatch.setattr(tcrdist.gapblock,name,object())
    with pytest.raises(RuntimeError,match='Python/native positional'):
        shortlist_audit_batch(['CASSLGQAYEQYF'],['CASSAGQAYEQYF'],['TRBV19*01'],['TRBV7-9*01'])


def test_streaming_and_bounded_legacy_audits_equal(monkeypatch):
    from vdjmatch.match import tcrdist
    if not hasattr(tcrdist.gapblock,'group_audit_batch'):
        import pytest
        pytest.skip('streaming native extension unavailable')
    args=(['CASSLGQAYEQYF','CASSLGQAPAYEQYF'],['CASSLGQAYEQYF','CASSAGQAYEQYF'],
          ['TRBV19*01']*2,['TRBV7-9*01','TRBV19*01'])
    streamed=shortlist_audit_batch(*args,exclude_exact=True)
    assert streamed['route']=='native_streaming' and streamed['projected_sparse_payload_bytes']==0
    monkeypatch.delattr(tcrdist.gapblock,'group_audit_batch')
    legacy=shortlist_audit_batch(*args,exclude_exact=True)
    assert legacy['route']=='native_ball_then_numpy'
    for key in ['cdr3_shortlist_keys','post_v_radius_keys','nearest_within_shortlist_distance']:
        assert streamed[key]==legacy[key]
