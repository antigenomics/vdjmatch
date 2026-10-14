"""Recorded targeted kernel arithmetic, independent of native distance implementation."""
import numpy as np
import pytest
from vdjmatch.evalue.kernel import kernel_ladder


@pytest.mark.parametrize('scorer', ['v-bounded','all-bounded','vonly-bounded',
                                   'v-rare','all-rare','vonly-rare'])
def test_advanced_scorer_selects_its_exact_retained_component(scorer):
    result=kernel_ladder([[1,3]],[[1,2]],[[0,2]],[3],[4],[0,26],
                         kernel_scale=585,scorer=scorer)
    np.testing.assert_array_equal(result['kernel_score'],
                                  result['kernel_'+scorer.replace('-','_')])
    compact=kernel_ladder([[1,3]],[[1,2]],[[0,2]],[3],[4],[0,26],
        kernel_scale=585,scorer=scorer,extras=False)
    np.testing.assert_array_equal(compact['kernel_score'],result['kernel_score'])
    assert len(compact)==2


def test_recorded_weighted_bounded_kernel():
    n = np.array([[1, 3, 4], [0, 0, 0]])
    same = np.array([[1, 2, 2], [0, 0, 0]])
    m = np.array([[0, 2, 5], [0, 1, 1]])
    t = np.array([0, 26, 52])
    result = kernel_ladder(n, same, m, [4, 4], [10, 10], t, kernel_scale=585)
    expected = (1 + (.25*2+.75)*np.exp(-26/585)/3 + .25*np.exp(-52/585)/6)/4
    assert result['kernel_score'][0] == pytest.approx(expected, abs=1e-15)
    assert result['kernel_score'][1] == 0
    assert np.allclose(result['kernel_all_bounded'], [
        (1 + 2*np.exp(-26/585)/3 + np.exp(-52/585)/6)/4, 0])
    assert result['kernel_v_rare'][0] > result['kernel_score'][0]


def test_kernel_empty_population_keeps_unavailable_rows():
    result = kernel_ladder([[0]], [[0]], [[0]], [0], [1], [0], kernel_scale=585)
    assert np.isnan(result['kernel_score'][0])


def test_unavailable_population_cannot_win_selection():
    import polars as pl
    from vdjmatch.evalue.kernel import select_kernel
    queries=pl.DataFrame({'query_id':[0,1]})
    scores=pl.DataFrame({'query_id':[0,0,1],'epitope':['EMPTY','SUPPORTED','EMPTY'],
        'kernel_score':[float('nan'),.5,float('nan')]})
    selected=select_kernel(queries,scores)
    assert selected['selected_epitope'].to_list()==['SUPPORTED',None]
    assert selected['selection_status'].to_list()==['ranked_evidence','no_support']


def test_default_annotation_route_selects_measured_geometry():
    import argparse
    from vdjmatch.cli.neighbours import register
    parser = argparse.ArgumentParser()
    register(parser.add_subparsers(), name='annotate')
    a = parser.parse_args(['annotate', 'sample.tsv', '--vdjdb', 'ref.tsv', '--locus', 'TRB',
                          '--epitope', 'E', '--output-prefix', 'result'])
    assert a.stem_depth == 4 and a.distance_matrix == 'blosum62'
    assert a.junction_ends == 'trim3' and a.neighbour_weighting == 'kernel'
    assert a.background == 'vdjdb-other'


def test_annotation_cli_reproduces_dense_kernel_and_retains_no_hits(tmp_path):
    import json
    import polars as pl
    from vdjmatch.cli.__main__ import main
    from vdjmatch.match.geometry import distance_matrix, distance_units, _substitution_costs
    from vdjmatch.evalue.kernel import ladder_cutoffs
    q = ['CASSLGQAYEQYF', 'WASSLGQAYEQYF']
    r = ['CASSAGQAYEQYF', 'CASSLGQAYEQWF', 'CASSLGQAPAYEQYF']
    rv = ['TRBV19', 'TRBV7-9', 'TRBV19']
    sample, reference = tmp_path/'q.tsv', tmp_path/'r.tsv'
    meta = dict(species='HomoSapiens', mhc_a='HLA-A*02:01', mhc_b='B2M', mhc_class='MHCI')
    pl.DataFrame([dict(sequence_id=str(i),junction_aa=s,v_call='TRBV19',locus='TRB')
                  for i,s in enumerate(q+['bad'])]).write_csv(sample,separator='\t')
    pl.DataFrame([dict(cdr3=s,v=v,gene='TRB',epitope='E' if i<2 else 'other',**meta)
                  for i,(s,v) in enumerate(zip(r,rv))]).write_csv(reference,separator='\t')
    d = distance_matrix(q,r,['TRBV19']*2,rv,ctrim=3,matrix='blosum62',stem_depth=4)
    scale = distance_units(_substitution_costs(matrix='blosum62'))[0]
    cuts=np.array(ladder_cutoffs(scale))
    rounded=cuts[np.searchsorted(cuts,d[:,:2])]
    m=np.sum(d[:,2,None,None]<=cuts[None,None,:],axis=1)
    control_at_shell = np.array([[int(d[i,2]<=t) for t in row] for i,row in enumerate(rounded)])
    weights=np.array([1,.25])
    expected=np.sum(weights[None,:]*np.exp(-rounded/(22.5*scale))/(1+control_at_shell),axis=1)/2
    for threads in (1,4):
        prefix=tmp_path/f'n{threads}'
        assert main(['annotate','--details',str(sample),'--vdjdb',str(reference),'--locus','TRB',
                     '--epitope','E','--threads',str(threads),'--output-prefix',str(prefix)])==0
        out=pl.read_csv(str(prefix)+'.enrichment.tsv',separator='\t')
        assert out['epitope'].to_list()==['E']*2
        np.testing.assert_allclose(out['kernel_score'].to_numpy(),expected,rtol=1e-14)
        status=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t')
        assert status.height==3 and not status['available'][2]
        manifest=json.loads((tmp_path/f'n{threads}.manifest.json').read_text())
        assert manifest['kernel']['matrix']=='blosum62'
        assert manifest['kernel']['cutoffs']==cuts.tolist()


def test_selection_pools_restrictions_and_preserves_ambiguity_and_no_hits():
    import polars as pl
    from vdjmatch.evalue.kernel import select_kernel
    scores=pl.DataFrame({'query_id':[0,0,0,1,1,2,4,4],
                         'epitope':['E','E','other','E','other','E','E','other'],
                         'kernel_score':[.8,.7,.2,.5,.5,0.,.9,.1],
                         'exact_reference_count':[0,0,0,0,0,0,1,1]})
    q=pl.DataFrame({'query_id':[0,1,2,3,4]})
    out=select_kernel(q,scores)
    assert out['kernel_support_share'].to_list()==[.8,.5,None,None,.9]
    assert out['selection_status'].to_list()==['ranked_evidence','ambiguous','no_support','no_support','ambiguous']
    assert out['kernel_score'].to_list()==[.8,.5,None,None,.9]


def test_paired_annotation_retains_unlinked_marginal_support(tmp_path):
    import polars as pl
    from vdjmatch.cli.__main__ import main
    sample,reference=tmp_path/'q.tsv',tmp_path/'r.tsv'
    pl.DataFrame({'clone_id':['pair','pair'],'sequence_id':['a','b'],'locus':['TRA','TRB'],
        'junction_aa':['CAVRDSNYQLIW','CASSLGQAYEQYF'],
        'v_call':['TRAV12-2','TRBV19']}).write_csv(sample,separator='\t')
    meta=dict(species='HomoSapiens',mhc_class='MHCI',mhc_a='HLA-A*02:01',mhc_b='B2M')
    pl.DataFrame([dict(cdr3=s,v=v,gene=g,epitope=e,**meta) for s,v,g,e in [
        ('CAVRDSNYQLIW','TRAV12-2','TRA','E'),
        ('CAVRDSNYQLVW','TRAV12-2','TRA','other'),
        ('CASSLGQAYEQYF','TRBV19','TRB','E'),
        ('CASSLGQAPYEQYF','TRBV19','TRB','other')]]).write_csv(reference,separator='\t')
    prefix=tmp_path/'out'
    assert main(['annotate','--details',str(sample),'--vdjdb',str(reference),'--locus','paired',
                 '--output-prefix',str(prefix)])==0
    evidence=pl.read_csv(str(prefix)+'.paired-evidence.tsv',separator='\t')
    assert evidence.height==2
    assert evidence['kernel_linked_score'].null_count()==2
    assert evidence['kernel_alpha_score'].min()>0 and evidence['kernel_beta_score'].min()>0
    selected=pl.read_csv(str(prefix)+'.selection.tsv',separator='\t')
    assert selected.height==1 and selected['has_support'][0]
    assert selected['linked_selection_status'][0]=='no_support'
    assert selected['selected_alpha_epitope'][0]=='E'


def test_kernel_rejects_missing_background_before_outputs(tmp_path):
    import argparse
    from vdjmatch.cli.neighbours import register,main
    parser=argparse.ArgumentParser()
    register(parser.add_subparsers())
    a=parser.parse_args(['neighbours','absent.tsv','--vdjdb','absent.zip',
        '--locus','TRB','--neighbour-weighting','kernel','--output-prefix',str(tmp_path/'out')])
    with pytest.raises(ValueError,match='requires.*background'):
        main(a)
    assert not list(tmp_path.iterdir())


def test_ordinary_annotation_validity_and_species_boundaries(tmp_path):
    import polars as pl
    from vdjmatch.cli.__main__ import main
    sample,reference=tmp_path/'q.tsv',tmp_path/'r.tsv'
    q='CASSLGQAYEQYF'
    pl.DataFrame({'query_id':[0,1],'junction_aa':[q,q],'v_call':['TRBV19']*2,
        'locus':['TRB']*2,'species':['HomoSapiens','MusMusculus']}).write_csv(sample,separator='\t')
    meta=dict(gene='TRB',v='TRBV19',species='HomoSapiens',mhc_class='MHCI',
              mhc_a='HLA-A*02:01',mhc_b='B2M')
    pl.DataFrame([dict(cdr3=s,epitope=e,**meta) for s,e in [
        (q,'E'),(q,'KLGGALQAK'),('CASSXXXAYEQYF','KLGGALQAK'),
        ('CASSLGQAPYEQYF','other')]]).write_csv(reference,separator='\t')
    for population,expected in [('all-valid',True),('curated-benchmark',False)]:
        prefix=tmp_path/population
        assert main(['annotate','--details',str(sample),'--vdjdb',str(reference),'--locus','TRB',
            '--epitope','E','--epitope','KLGGALQAK','--reference-population',population,
            '--output-prefix',str(prefix)])==0
        evidence=pl.read_csv(str(prefix)+'.enrichment.tsv',separator='\t')
        assert evidence['query_id'].unique().to_list()==[0]
        assert ('KLGGALQAK' in evidence['epitope'])==expected
        selected=pl.read_csv(str(prefix)+'.selection.tsv',separator='\t')
        assert selected['selection_status'][1]=='no_support'
        assert selected['selection_status'][0]==('ambiguous' if expected else 'ranked_evidence')
        single_prefix=tmp_path/(population+'-single')
        assert main(['annotate','--details',str(sample),'--vdjdb',str(reference),'--locus','TRB',
            '--epitope','E','--reference-population',population,'--output-prefix',str(single_prefix)])==0
        single=pl.read_csv(str(single_prefix)+'.enrichment.tsv',separator='\t')
        multi=evidence.filter(pl.col('epitope')=='E')
        assert multi['control_population'].to_list()==single['control_population'].to_list()
        np.testing.assert_allclose(multi['kernel_score'],single['kernel_score'],rtol=0,atol=0)
