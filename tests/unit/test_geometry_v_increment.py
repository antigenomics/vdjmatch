"""Fixed CDR3 versus CDR3+V contrast, with identical eligibility and puncturing."""
import numpy as np
import pytest
from vdjmatch.match.geometry import distance_matrix,total_distance_count_batch,paired_distance_matrix,paired_total_distance_mass_batch


def test_v_increment_matches_dense_components_and_linked_counts():
    q=['CASSLGQAYEQYF'];r=['CASSAGQAYEQYF','CASSLGQAPAYEQYF'];qv=['TRBV19*01'];rv=['TRBV7-9*01']*2
    full,cdr3=distance_matrix(q,r,qv,rv,return_cdr3=True)
    assert np.all(full>cdr3)
    assert np.array_equal(distance_matrix(q,r,qv,rv,v_loop_weight=0),cdr3)
    limits=[0,*sorted(set((200*cdr3).ravel().tolist()+(200*full).ravel().tolist()))]
    for weight,expected in [(0,cdr3),(1,full)]:
        count,mass=total_distance_count_batch(q,r,qv,rv,[limits],v_loop_weight=weight,linear_mass=True)
        assert count==[[(200*expected[0]<=t).sum() for t in limits]]
        assert mass==[[np.maximum(0,t-200*expected[0]).sum() for t in limits]]
    alpha=['CAVRDSNYQLIW'];av=['TRAV12-2*01'];arv=['TRAV21*01']*2
    pair=paired_distance_matrix(alpha,q,alpha*2,r,av,qv,arv,rv,v_loop_weight=0)
    counts,mass=paired_total_distance_mass_batch(alpha,q,alpha*2,r,av,qv,arv,rv,[limits],v_loop_weight=0)
    assert counts==[[(200*pair[0]<=t).sum() for t in limits]]
    assert mass==[[np.maximum(0,t-200*pair[0]).sum() for t in limits]]
    # Native exact puncturing remains full junction based, even when V differs.
    assert total_distance_count_batch(q,q,qv,rv[:1],[[0]],v_loop_weight=0,exclude_exact=True)==[[0]]
    for invalid in [True,.5,2,-1]:
        with pytest.raises(ValueError,match='v_loop_weight'):
            total_distance_count_batch(q,r,qv,rv,[[0]],v_loop_weight=invalid)
    with pytest.raises(ValueError,match='unknown human V'):
        distance_matrix(q,r,['unknown'],rv,v_loop_weight=0)


def test_cli_v_increment_exports_same_population_and_explicit_geometry(tmp_path):
    import json
    import polars as pl
    from vdjmatch.cli.__main__ import main
    q=tmp_path/'query.tsv';r=tmp_path/'reference.tsv'
    meta=dict(species='HomoSapiens',epitope='E',mhc_a='HLA-A*02:01',mhc_b='B2M',mhc_class='MHCI')
    pl.DataFrame([dict(junction_aa='CASSLGQAYEQYF',v_call='TRBV19*01',locus='TRB',**meta),
                  dict(junction_aa='CASSLGQAYEQYF',v_call='unknown',locus='TRB',**meta)]).write_csv(q,separator='\t')
    pl.DataFrame([dict(cdr3='CASSAGQAYEQYF',v='TRBV7-9*01',gene='TRB',**meta),
                  dict(cdr3='CASSAGQAPAYEQYF',v='TRBV19*01',gene='TRB',**(meta|{'epitope':'N'}))]).write_csv(r,separator='\t')
    frames=[]
    for weight in (0,1):
        prefix=tmp_path/str(weight)
        assert main(['neighbours',str(q),'--vdjdb',str(r),'--locus','TRB',
            '--targets-from-sample','--background','vdjdb-other','--radius','12',
            '--v-loop-weight',str(weight),'--unbounded-nearest','--output-prefix',str(prefix)])==0
        manifest=json.loads((tmp_path/f'{weight}.manifest.json').read_text())
        assert manifest['parameters']['v_loop_weight']==weight
        assert manifest['background_enrichment']['v_loop_weight']==weight
        assert ('CDR3-only' in manifest['distance'])==(weight==0)
        queries=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t')
        assert queries.height==2 and queries['available'].all()
        assert queries['v_geometry_imputed'].to_list()==[False,True]
        frames.append(pl.read_csv(str(prefix)+'.enrichment.tsv',separator='\t'))
    assert frames[0]['n_reference'].to_list()==[1,1] and frames[1]['n_reference'].to_list()==[0,0]
    assert frames[0]['reference_population'].to_list()==frames[1]['reference_population'].to_list()
    assert frames[0]['control_population'].to_list()==frames[1]['control_population'].to_list()
