"""Columnar zero rows preserve every scalar field and source appearance."""
import polars as pl
import pytest
from polars.testing import assert_frame_equal
from seqtree import Index

from vdjmatch.match import historical


def _reference(sequences, genes):
    return pl.DataFrame({'cdr3':sequences,'v':genes,'j':['TRBJ1-1']*len(sequences),
        'gene':['TRB']*len(sequences),'species':['HomoSapiens']*len(sequences),
        'epitope':['ARTIFICIAL']*len(sequences),'mhc_a':['A']*len(sequences),
        'mhc_b':['B2M']*len(sequences),'mhc_class':['MHCI']*len(sequences)})


def _compare(query,reference,*,similarities=None):
    control=Index.build(reference['cdr3'].to_list()+query['cdr3'].to_list()+['CASSLGKAYEQYF'])
    seqs=query['cdr3'].unique(maintain_order=True).to_list()
    rs=reference.unique('cdr3',keep='first',maintain_order=True)['cdr3'].to_list()
    hits,unit,radii=historical._density_neighbours(seqs,rs,1)
    prepared=historical._pack_prepared(tuple(seqs),rs,hits,unit,radii)
    counts=historical.control_histograms(control,seqs,radii=radii)
    options=dict(pool_reference=True,control_counts=counts,v_similarities=similarities)
    expected=historical.density(query,reference,control,**options)
    actual=historical.density(query,reference,control,prepared_neighbours=prepared,**options)
    assert_frame_equal(actual,expected,check_exact=True)
    assert actual['query_id'].equals(query['query_id'])
    return actual,prepared,control,options


def test_supported_mixed_v_shared_first_representative_and_duplicate_ids():
    query=pl.DataFrame({'query_id':[7,7,8,8,8,9,9],
        'cdr3':['CASSLGQAYEQYF','CASSLGQAYEQYF','CASSLGQAYEQFF','CASSLGQAYEQFF',
                'CASSLAGQAYEQYF','CWWWWWWWWWWWF','CASSLGQAYEQYF'],
        'v':['TRBV19','TRBV7-9',None,'','TRBV19',None,'TRBV19'],
        '_density_row':[100]*7})
    reference=_reference(['CASSLGQAYEQFF','CASSLGQAYEQYF','CASSLGQAYEQYF','CASSLAGQAYEQFF'],
                         ['TRBV7-9','TRBV19','TRBV7-9',None])
    result,_,_,_=_compare(query,reference)
    assert result['score'][0]>0 and result['score'][-1]==result['score'][0]
    assert result['score'][5]==0 and result['p_enrichment'][5]==1
    assert result['n_reference'].unique().to_list()==[3]
    # Reversing source-order representatives must still match the scalar owner.
    _compare(query,reference.reverse())


@pytest.mark.parametrize('genes',[[None,''],['TRBV19','TRBV7-9']])
def test_zero_edges_return_complete_typed_rows_without_regex_or_row_iteration(monkeypatch,genes):
    query=pl.DataFrame({'query_id':['repeat','repeat'], 'cdr3':['CWWWWWWWWWWWF']*2,'v':genes})
    reference=_reference(['CASSLGQAYEQYF'],['TRBV19'])
    result,prepared,control,options=_compare(query,reference)
    assert not prepared.query_indices
    assert result['score'].to_list()==[0.,0.] and result['p_enrichment'].to_list()==[1.,1.]
    # The zero-edge return must not visit source rows or perform per-row V regex.
    families=[]
    family=historical.gene_family
    def observed_family(value):
        families.append(value)
        return family(value)
    monkeypatch.setattr(historical,'gene_family',observed_family)
    monkeypatch.setattr(pl.DataFrame,'iter_rows',lambda *a,**k:pytest.fail('zero rows must stay columnar'))
    historical.density(query,reference,control,prepared_neighbours=prepared,
        **(options|{'v_similarities':{}}))
    assert families==['TRBV19']  # Reference normalized once; no query-row regex.


def test_unit_only_indel_and_exact_only_queries_keep_significance_and_zero_density():
    query=pl.DataFrame({'query_id':[1,1,2],
        'cdr3':['CASSLAGQAYEQYF','CASSGQAYEQYF','CASSLGQAYEQYF'],
        'v':[None,'','TRBV19']})
    reference=_reference(['CASSLGQAYEQYF'],['TRBV19'])
    result,prepared,_,_=_compare(query,reference)
    assert result['score'].to_list()==[0.,0.,0.]
    assert all(n>0 for n in list(prepared.unit_counts)[:2])
    assert result['p_enrichment'][0]<1 and result['p_enrichment'][1]<1
    assert result['p_enrichment'][2]==1.


def test_missing_similarity_error_and_prepared_axis_validation_remain_strict():
    query=pl.DataFrame({'query_id':[1],'cdr3':['CASSLGQAYEQFF'],'v':['TRBV7-9']})
    reference=_reference(['CASSLGQAYEQYF'],['TRBV19'])
    _,prepared,control,options=_compare(query,reference)
    for p in (None,prepared):
        with pytest.raises(ValueError,match='precomputed V similarities'):
            historical.density(query,reference,control,prepared_neighbours=p,
                **(options|{'v_similarities':{}}))
    with pytest.raises(ValueError,match='prepared.*axes'):
        historical.density(query.with_columns(pl.lit('CASSLGKAYEQYF').alias('cdr3')),reference,
            control,prepared_neighbours=prepared,**options)
