"""Gene-conditioned native counts retain the same fixed sequence-search geometry."""
import polars as pl
import pytest
from seqtree import Index
from vdjmatch.evalue.control import conditioned_counts,impute_gene_calls
from vdjmatch.match.scope import search_params

@pytest.mark.parametrize('exclude_exact',[False,True])
def test_conditioned_histogram_matches_explicit_fixed_populations(exclude_exact):
    q=pl.DataFrame({'query_id':[0,1,2,3],'cdr3':['CASSF','CASSF','CASRF','CAAAF'],'v':['TRBV19','TRBV5-1','TRBV19','TRBV9']})
    control=pl.DataFrame({'cdr3':['CASSF','CASSF','CASRF','CASSF'],'v':['TRBV19','TRBV19','TRBV19','TRBV5-1']})
    params=search_params('2,1,1,2',gap_open=1,gap_extend=1)
    cdf,totals,sizes=conditioned_counts(q,control,params,match_v=True,exclude_exact=exclude_exact)
    assert sizes.sort('query_id')['_control_size'].to_list()==[2,1,2,0]
    for row in q.iter_rows(named=True):
        sequences=control.filter(pl.col('v')==row['v'])['cdr3'].unique().to_list()
        hits=Index.build(sequences,'aa').search_batch([row['cdr3']],params,1)[0] if sequences else []
        expected=[sum(h.n_subs+h.n_ins+h.n_dels<=r and (not exclude_exact or sequences[h.ref_id]!=row['cdr3']) for h in hits) for r in range(3)]
        assert cdf.filter(pl.col('query_id')==row['query_id']).sort('edits')['_nc'].to_list()==expected
        assert totals.filter(pl.col('query_id')==row['query_id'])['n_control'][0]==expected[-1]

def test_missing_calls_are_retained_and_flagged_without_labels():
    q=pl.DataFrame({'query_id':[0,1,2],'v':['TRBV19*01',None,'unknown'],'count':[2,7,11]})
    result=impute_gene_calls(q,species='human',locus='TRB',kinds=('v',))
    assert result['query_id'].to_list()==q['query_id'].to_list()
    assert result['count'].to_list()==q['count'].to_list()
    assert result['original_v'].to_list()==q['v'].to_list()
    assert result['v'].to_list()==['TRBV19']*3
    assert result['v_imputed'].to_list()==[False,True,True]


def test_api_uses_same_gene_conditioned_populations_and_retains_missing_query():
    from vdjmatch.api import Annotator
    from seqtree.evalue import evalue_result
    ref=pl.DataFrame({'gene':['TRB']*3,'species':['HomoSapiens']*3,'cdr3':['CASSF','CASRF','CASTF'],'v':['TRBV19','TRBV19','TRBV5-1'],'epitope':['A']*3})
    ann=Annotator.from_frame(ref)
    control=pl.DataFrame({'cdr3':['CASSF','CASRF','CASTF'],'v':['TRBV19','TRBV19','TRBV5-1']})
    query=pl.DataFrame({'query_id':[0,1],'cdr3':['CASSF','CASSF'],'v':['TRBV19*01',None]})
    result=ann.graded_candidates(query,match_v=True,control=control,exclude_exact=True)
    assert result['query_id'].n_unique()==2
    one=result.filter(pl.col('radius')==1).sort('query_id')
    assert one['n_control'].to_list()==[1,1] and one['n_reference'].to_list()==[2,2]
    assert one['query_v_imputed'].to_list()==[False,True]
    assert one['query_original_v'].to_list()==['TRBV19*01',None]
    expected=evalue_result(1,1,2,2)
    assert one['E'].to_list()==[expected['E']]*2


def test_cli_gene_labelled_control_is_used_for_the_same_predicate(tmp_path):
    import json
    from vdjmatch.cli import __main__ as cli
    ref=tmp_path/'reference.tsv';query=tmp_path/'query.tsv';ctrl=tmp_path/'control.tsv'
    pl.DataFrame({'cdr3':['CASSF','CASRF','CASTF'],'gene':['TRB']*3,'species':['HomoSapiens']*3,'v.segm':['TRBV19','TRBV19','TRBV5-1'],'antigen.epitope':['E']*3}).write_csv(ref,separator='\t')
    pl.DataFrame({'sequence_id':['a','b'],'junction_aa':['CASSF','CASSF'],'v_call':['TRBV19*01',None],'locus':['TRB']*2}).write_csv(query,separator='\t')
    pl.DataFrame({'cdr3aa':['CASSF','CASRF','CASTF'],'v':['TRBV19','TRBV19','TRBV5-1']}).write_csv(ctrl,separator='\t')
    args=['match',str(query),'--vdjdb',str(ref),'--fresh-reference','--control','TRB='+str(ctrl),'--match-v','--search-mode','ball','--matrix','none','--exclude-exact','--no-align','--threads','1','-o',str(tmp_path/'out')]
    assert cli.main(args)==0
    got=pl.read_csv(tmp_path/'out.query.candidates.txt',separator='\t').filter(pl.col('radius')==1)
    assert got['query_id'].n_unique()==2 and got['n_control'].to_list()==[1,1]
    manifest=json.loads((tmp_path/'out.query.manifest.json').read_text())
    assert manifest['controls']['TRB']['gene_columns']['v']=='v'
    assert cli.main(args+['--no-hit-table','-o',str(tmp_path/'compact')])==0
    compact=pl.read_csv(tmp_path/'compact.query.candidates.txt',separator='\t').filter(pl.col('radius')==1)
    assert compact.sort('query_id').equals(got.sort('query_id'))
    assert not (tmp_path/'compact.query.hits.txt').exists()


def test_ambiguous_allele_calls_do_not_select_the_first_gene():
    q=pl.DataFrame({'v':['TRBV19*01','TRBV2*01,TRBV3*01','TRBV2,TRBV3','TRBV19*01,TRBV19*02']})
    got=impute_gene_calls(q,species='human',locus='TRB',kinds=('v',))
    assert got['v'].to_list()==['TRBV19']*4
    assert got['v_imputed'].to_list()==[False,True,True,False]
    again=impute_gene_calls(got,species='human',locus='TRB',kinds=('v',))
    assert again['original_v'].equals(q['v'].rename('original_v'))
    assert again['v_imputed'].equals(got['v_imputed'])


def test_mixed_control_populations_are_selected_before_imputation():
    from vdjmatch.evalue.control import control_population
    raw=pl.DataFrame({'cdr3':['CASSF','CAVRF','CASRF','CASTF'],'v':['TRBV19','TRAV1',None,'TRBV5-1'],'species':['human','human','mouse','human'],'locus':['TRB','TRA',None,None]})
    selected=control_population(raw,species='human',locus='TRB')
    assert selected['cdr3'].to_list()==['CASSF','CASTF']


def test_annotated_nohit_query_keeps_original_gene_and_flag():
    from vdjmatch.api import Annotator
    ref=pl.DataFrame({'gene':['TRB'],'species':['HomoSapiens'],'cdr3':['CASSF'],'v':['TRBV19'],'epitope':['A']})
    q=pl.DataFrame({'cdr3':['CYYYYYYYYF'],'v':[None]},schema_overrides={'v':pl.String})
    result=Annotator.from_frame(ref).annotate(q,match_v=True,scope='1')
    assert result['vdjmatch_status'].to_list()==['no_hit']
    assert result['vdjmatch_v_imputed'].to_list()==[True]
    assert result['vdjmatch_original_v'].to_list()==[None]


def test_graded_preimputed_metadata_survives_normalization():
    from vdjmatch.api import Annotator
    ref=pl.DataFrame({'gene':['TRB'],'species':['HomoSapiens'],'cdr3':['CASSF'],'v':['TRBV19'],'epitope':['A']})
    q=pl.DataFrame({'cdr3':['CASSF'],'v':['TRBV19'],'original_v':[None],'v_imputed':[True]},schema_overrides={'original_v':pl.String})
    got=Annotator.from_frame(ref).graded_candidates(q,match_v=True)
    assert got['query_original_v'].to_list()==[None]*5
    assert got['query_v_imputed'].to_list()==[True]*5


def test_mixed_invalid_locus_keeps_its_query_disposition():
    from vdjmatch.api import Annotator
    ref=pl.DataFrame({'gene':['TRB'],'species':['HomoSapiens'],'cdr3':['CASSF'],'v':['TRBV19'],'epitope':['A']})
    q=pl.DataFrame({'cdr3':['CASSF','CASSF'],'v':['TRBV19',None],'locus':['TRB',None]})
    got=Annotator.from_frame(ref).annotate(q,match_v=True)
    assert got.height==2
    assert got['vdjmatch_status'].to_list()==['matched','invalid_query']


def test_conditioned_counts_rejects_an_uncapped_nonexact_search():
    from seqtree import SearchParams
    q=pl.DataFrame({'query_id':[0],'cdr3':['CASSF'],'v':['TRBV19']})
    with pytest.raises(ValueError,match='explicit total-edit cap'):
        conditioned_counts(q,q.drop('query_id'),SearchParams(max_subs=2),match_v=True)


def test_missing_reference_genes_remain_in_the_fixed_conditioned_population():
    from vdjmatch.api import Annotator
    from vdjmatch.api import _prepare
    ref=pl.DataFrame({'gene':['TRB']*2,'species':['HomoSapiens']*2,'cdr3':['CASSF','CASRF'],'v':['TRBV19',None],'epitope':['A']*2})
    ann=Annotator.from_frame(ref)
    before=ann._index.reference.clone()
    _,q=_prepare(pl.DataFrame({'cdr3':['CASSF'],'v':['TRBV19']}))
    hits,c=ann._evidence(q,search_params('1',gap_open=1,gap_extend=1),match_v=True)
    assert hits.height==2 and int(hits['v_imputed'].sum())==1
    assert c['n_reference'].to_list()==[2]
    assert before.equals(ann._index.reference)
