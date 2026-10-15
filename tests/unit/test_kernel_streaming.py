"""Compact selection must equal detailed selection on synthetic evidence."""
import json
import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal
from vdjmatch.evalue.kernel import (select_kernel,kernel_peptides,
                                    update_kernel_selection,finish_kernel_selection)


@pytest.mark.parametrize('order',[['A','B','C'],['C','B','A']])
def test_streamed_selection_preserves_restrictions_ties_and_missing(order):
    queries=pl.DataFrame({'query_id':range(8),'available':[True]*7+[False]})
    scores=pl.DataFrame({'query_id':[0,0,0,1,1,2,2,3,3,4,4,5,6],
        'epitope':['A','A','B','A','B','A','B','A','B','A','B','C','B'],
        'kernel_score':[.8,.7,.2,.5,.5,0.,0.,float('nan'),.4,.9,.1,None,.3],
        'exact_reference_count':[0,0,0,0,0,0,0,0,0,1,1,0,0]})
    state=None
    for peptide in order:
        state=update_kernel_selection(state,kernel_peptides(scores.filter(pl.col('epitope')==peptide)))
        assert state.height<=queries.height
    assert_frame_equal(finish_kernel_selection(queries,state),select_kernel(queries,scores),check_dtypes=False)
    assert_frame_equal(finish_kernel_selection(queries,None),select_kernel(queries,scores.head(0)),check_dtypes=False)


def _annotation_main(args):
    import argparse
    from vdjmatch.cli.neighbours import register
    parser=argparse.ArgumentParser()
    register(parser.add_subparsers(),name='annotate')
    a=parser.parse_args(args)
    return a.func(a)


def _single_inputs(tmp_path):
    sample,reference=tmp_path/'q.tsv',tmp_path/'r.tsv'
    pl.DataFrame({'sequence_id':['x','dup','invalid','far'],'junction_aa':[
        'CASSLGQAYEQYF','CASSLGQAYEQYF','bad','CWWWWWWWWWWWF'],
        'v_call':['TRBV19']*4,'locus':['TRB']*4,
        'species':['HomoSapiens']*4,'epitope':['A']*4,'mhc_a':['HLA-A*02:01']*4,'mhc_b':['B2M']*4,'mhc_class':['MHCI']*4}).write_csv(sample,separator='\t')
    meta=dict(species='HomoSapiens',gene='TRB',v='TRBV19',mhc_class='MHCI',mhc_b='B2M')
    pl.DataFrame([dict(cdr3=s,epitope=e,mhc_a=a,**meta) for s,e,a in [
        ('CASSLGQAYEQYF','A','HLA-A*02:01'),('CASSLGQAYEQYF','B','HLA-A*02:01'),
        ('CASSLGQAPAYEQYF','A','HLA-A*03:01'),('CASSLGQAYEQWF','C','HLA-A*02:01')]]).write_csv(reference,separator='\t')
    return sample,reference


@pytest.mark.parametrize('selectors',[[],['--epitope','A','--epitope','B'],['--targets-from-sample'],['--withhold-epitope','A'],['--withhold-epitope','A','--withhold-epitope','B','--withhold-epitope','C']])
@pytest.mark.parametrize('radius',[None,'20.125'])
def test_cli_streamed_selection_matches_detailed(tmp_path,selectors,radius):
    main=_annotation_main
    sample,reference=_single_inputs(tmp_path)
    for compact in [False,True]:
        prefix=tmp_path/str(compact)
        args=['annotate',str(sample),'--vdjdb',str(reference),'--locus','TRB','--output-prefix',str(prefix),*selectors]
        if radius:args+=['--radius-factor',str(float(radius)/26)]
        if compact:args+=['--selection-only']
        assert main(args)==0
    for table in ['queries','selection']:
        full=pl.read_csv(str(tmp_path/'False')+'.'+table+'.tsv',separator='\t')
        stream=pl.read_csv(str(tmp_path/'True')+'.'+table+'.tsv',separator='\t')
        assert_frame_equal(full,stream,check_dtypes=False,check_exact=True)
    assert not (tmp_path/'True.enrichment.tsv').exists()
    assert not (tmp_path/'True.candidates.tsv').exists()
    assert json.loads((tmp_path/'True.manifest.json').read_text())['selection_only']


def test_native_stem_rounding_at_dense_membership_boundaries():
    from vdjmatch.match.geometry import distance_matrix,total_distance_count_batch
    q=['CASSLGQAYEQYF'];r=['WASSLGQAYEQYF','CASALGQAYEQWF','CASSLGQAPAYEQYF']
    for weight in [1,.25,1.5]:
        options=dict(ctrim=3,matrix='blosum62',stem_depth=4,stem_weight=weight)
        d=distance_matrix(q,r,['TRBV19'],['TRBV19']*3,**options)
        cuts=sorted(set([max(0,int(x*200)-1) for x in d[0]]+[int(x*200) for x in d[0]]))
        actual=total_distance_count_batch(q,r,['TRBV19'],['TRBV19']*3,cuts,shared_thresholds=True,**options)
        np.testing.assert_array_equal(actual,(d[:,:,None]*200<=np.array(cuts)).sum(axis=1))


def test_paired_streaming_preserves_unlinked_marginal_selections(tmp_path):
    main=_annotation_main
    sample,reference=tmp_path/'p.tsv',tmp_path/'r.tsv'
    pl.DataFrame({'clone_id':['pair','pair'],'sequence_id':['a','b'],'locus':['TRA','TRB'],
        'junction_aa':['CAVRDSNYQLIW','CASSLGQAYEQYF'],
        'v_call':['TRAV12-2','TRBV19']}).write_csv(sample,separator='\t')
    meta=dict(species='HomoSapiens',mhc_class='MHCI',mhc_a='HLA-A*02:01',mhc_b='B2M')
    pl.DataFrame([dict(cdr3=s,v=v,gene=g,epitope=e,**meta) for s,v,g,e in [
        ('CAVRDSNYQLIW','TRAV12-2','TRA','E'),('CAVRDSNYQLVW','TRAV12-2','TRA','other'),
        ('CASSLGQAYEQYF','TRBV19','TRB','E'),('CASSLGQAPYEQYF','TRBV19','TRB','other')]]).write_csv(reference,separator='\t')
    for compact in [False,True]:
        args=['annotate',str(sample),'--vdjdb',str(reference),'--locus','paired','--output-prefix',str(tmp_path/str(compact))]
        if compact:args+=['--selection-only']
        assert main(args)==0
    for suffix in ['','.TRA','.TRB']:
        assert_frame_equal(pl.read_csv(str(tmp_path/'False')+suffix+'.selection.tsv',separator='\t'),
                           pl.read_csv(str(tmp_path/'True')+suffix+'.selection.tsv',separator='\t'),check_dtypes=False,check_exact=True)
    assert not (tmp_path/'True.paired-evidence.tsv').exists()
