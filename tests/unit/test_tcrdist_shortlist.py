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
    pl.DataFrame([dict(junction_aa='CASSLGQAYEQYF',v_call=v,locus='TRB',**meta) for v in ['TRBV19*01','unknown']]).write_csv(q,separator='\t')
    pl.DataFrame([dict(cdr3='CASSAGQAYEQYF',v='TRBV7-9*01',gene='TRB',**meta)]).write_csv(r,separator='\t')
    for audit in [False,True]:
        prefix=tmp_path/str(audit)
        args=['tcrdist-neighbours',str(q),'--vdjdb',str(r),'--locus','TRB','--targets-from-sample','--unbounded-nearest','--output-prefix',str(prefix)]
        assert main(args+(['--shortlist-audit'] if audit else []))==0
    for suffix in ['queries','candidates']:
        assert pl.read_csv(tmp_path/f'False.{suffix}.tsv',separator='\t').equals(pl.read_csv(tmp_path/f'True.{suffix}.tsv',separator='\t'))
    diagnostic=pl.read_csv(tmp_path/'True.shortlist.tsv',separator='\t');assert diagnostic.height==2
    assert diagnostic['cdr3_shortlist_keys'].to_list()==[1,None]
    assert diagnostic['status'].to_list()==['matched','invalid_query_or_gene']
    assert 'shortlist_audit' in json.loads((tmp_path/'True.manifest.json').read_text())
