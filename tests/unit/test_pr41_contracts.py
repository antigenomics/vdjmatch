"""Review regressions: fixed identities/populations and bounded CLI behavior."""
import numpy as np
import polars as pl
import pytest


def test_callback_measurements_follow_supplied_target_order():
    from vdjmatch.aggregate.repertoire_enrichment import matched_enrichment
    obs=pl.DataFrame(dict(sample_id=['s'],query_id=[0],count=[7],pgen_bin=[-8]))
    ctrl=pl.DataFrame(dict(query_id=[1],event_weight=[1],pgen_bin=[-8]))
    targets=pl.DataFrame(dict(epitope=['Y','X']))
    summary,mocks,_=matched_enrichment(obs,ctrl,None,targets,
        measure=lambda rows,counts:(np.array([1,0]),np.array([7,0])))
    assert summary.select('epitope','matched_reads').rows()==[('Y',7),('X',0)]
    assert set(mocks.filter(pl.col('epitope')=='Y')['mock_reads'])=={7}


def test_explicit_all_valid_index_keeps_selected_artifacts_but_not_invalid_aa():
    from vdjmatch.match.engine import VdjdbIndex
    ref=pl.DataFrame(dict(cdr3=['CASSF','CASRF','C?F'],gene=['TRB']*3,
        epitope=['KLGGALQAK','E','E'],reference_id=['other','PMID:40498839','other']))
    assert VdjdbIndex.build(ref).genes==[]
    got=VdjdbIndex.build(ref,reference_population='all-valid').records_for('TRB')
    assert got['cdr3'].to_list()==['CASSF','CASRF']


def test_imputed_geometry_never_merges_original_reference_keys():
    from vdjmatch.aggregate.structure import composition
    from vdjmatch.match.geometry import impute_reference_v_alleles
    ref=pl.DataFrame(dict(cdr3=['CASSF']*4,v=['TRBV1',None,'UNKNOWN',None]))
    allele,flag=impute_reference_v_alleles(ref['v'],{'TRBV1*01':dict(locus='TRB',cdr1='A',cdr2='A',cdr25='A')},locus='TRB')
    ref=ref.with_columns(allele.alias('_allele'),flag)
    full,_=composition(ref,ref.select('cdr3').unique(),include_punctures=False)
    assert full['clonotype_keys']==3 and full['reference_records']==4


def test_missing_gene_fallback_is_allele_pooled_and_idempotent():
    from vdjmatch.evalue.control import impute_gene_calls,conditioned_counts
    from vdjmatch.match.scope import search_params
    raw=pl.DataFrame(dict(cdr3=['CASSF'],v=[None]))
    got=impute_gene_calls(raw,species='human',locus='TRB',kinds=('v',))
    assert '*' not in got['v'].item()
    assert impute_gene_calls(got,species='human',locus='TRB',kinds=('v',)).equals(got)
    query=got.select('cdr3','v').with_row_index('query_id')
    _,counts,sizes=conditioned_counts(query,got,search_params('0'),match_v=True)
    assert counts['n_control'].item()==sizes['_control_size'].item()==1


def test_unsupplied_locus_with_unknown_v_remains_missing(tmp_path):
    from vdjmatch import io
    p=tmp_path/'rows.tsv'
    p.write_text('junction_aa\tv_call\tj_call\tduplicate_count\nCASSF\tUNKNOWN\t\t7\n')
    got=io.read_rearrangement(p,valid_aa=False)
    assert got['locus'].item() is None and got['count'].item()==7
    p.write_text(p.read_text().replace('duplicate_count','duplicate_count\tlocus').replace('\t7\n','\t7\tTRA\n'))
    assert io.read_rearrangement(p,valid_aa=False)['locus'].item()=='TRA'


def test_fixed_no_hit_retention_preserves_summary_and_native_streaming(tmp_path,monkeypatch):
    from vdjmatch import Annotator
    from vdjmatch.runner.multisample import annotate_sample
    from seqtree import Index
    ref=pl.DataFrame(dict(cdr3=['CASSF','CASRF','CASRF'],gene=['TRB']*3,
        epitope=['E']*3,antigen_species=['one','two','two']))
    ann=Annotator.from_frame(ref)
    p=tmp_path/'rows.tsv'
    pl.DataFrame(dict(cdr3=['CASSF']*40+['CWWWF'],v=['TRBV1']*41,count=[2]*41)).write_csv(p,separator='\t')
    expected=annotate_sample(ann._index,p,matrix=None,with_evalue=False,threads=1)
    native=Index.search_batch_consume;calls=[]
    def consume(self,*a,**kw):
        calls.append(1);return native(self,*a,**kw)
    monkeypatch.setattr(Index,'search_batch_consume',consume)
    actual=annotate_sample(ann._index,p,matrix=None,with_evalue=False,threads=1,retain_hits=False)
    assert calls==[1] and actual['hits'].is_empty()
    for key in ('summary','candidates','calls'):
        cols=actual[key].columns
        assert actual[key].sort(cols).equals(expected[key].sort(cols)),key


def test_species_any_paired_marginals_rejected_before_search(tmp_path,monkeypatch):
    from vdjmatch import Annotator
    from vdjmatch.runner.multisample import annotate_sample
    from vdjmatch.match.engine import VdjdbIndex
    monkeypatch.setattr(Annotator,'paired_scope_candidates',lambda *a,**kw:pytest.fail('search started'))
    with pytest.raises(ValueError,match='explicit human or mouse species'):
        annotate_sample(VdjdbIndex({}),tmp_path/'absent.tsv',paired=True,
            paired_marginals=True,with_evalue=False,species=None)
