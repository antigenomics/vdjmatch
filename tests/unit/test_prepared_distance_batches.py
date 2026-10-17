"""Prepared model inputs preserve batch distances and strict gene validation."""
import numpy as np
import pytest
from vdjmatch.match import geometry as td


@pytest.mark.parametrize('species',['human','mouse'])
def test_prepared_alias_lookup_equals_original_model_resolution(species):
    prepared=td.precompute_v_loop_distances(species=species)
    model=prepared['model']
    aliases,_=td._v_allele_aliases(model)
    calls=list(aliases)+list(model)+[' '+v+' ' for v in model]+[None,'BAD','TRBV19*99','TRAV99,TRAV1']
    for locus in [None,'TRA','TRB']:
        assert td._resolved_v_calls(calls,model,prepared,locus=locus)==td.resolve_v_alleles(calls,model,locus=locus).to_list()


@pytest.mark.parametrize('species',['human','mouse'])
@pytest.mark.parametrize('matrix',['clipped','blosum62'])
def test_prepared_single_and_paired_batches_match_without_rebuilding(monkeypatch,species,matrix):
    prepared=td.precompute_v_loop_distances(species=species,matrix=matrix)
    model=prepared['model']
    alpha=next(v for v in sorted(model) if model[v]['locus']=='TRA')
    beta=next(v for v in sorted(model) if model[v]['locus']=='TRB')
    q=['CASSLGQAYEQYF','CASSLGQAYEQWF'];r=q[::-1]
    options=dict(species=species,matrix=matrix,ctrim=3,stem_depth=4)
    dense=td.distance_matrix(q,r,[beta]*2,[beta]*2,**options)
    paired=td.paired_distance_matrix(q,q,r,r,[alpha]*2,[beta]*2,[alpha]*2,[beta]*2,**options)
    counts=td.total_distance_count_batch(q,r,[beta]*2,[beta]*2,[0,200,50000],shared_thresholds=True,**options)
    paired_counts=td.paired_total_distance_mass_batch(q,q,r,r,[alpha]*2,[beta]*2,[alpha]*2,[beta]*2,
        [0,200,50000],shared_thresholds=True,**options)
    def fail(*args,**kwargs):
        raise AssertionError('batch rebuilt the model resolver or substitution matrix')
    monkeypatch.setattr(td,'resolve_v_alleles',fail)
    monkeypatch.setattr(td,'_substitution_costs',fail)
    np.testing.assert_array_equal(dense,td.distance_matrix(q,r,[beta]*2,[beta]*2,v_loop_distances=prepared,**options))
    np.testing.assert_array_equal(paired,td.paired_distance_matrix(q,q,r,r,[alpha]*2,[beta]*2,[alpha]*2,[beta]*2,v_loop_distances=prepared,**options))
    assert counts==td.total_distance_count_batch(q,r,[beta]*2,[beta]*2,[0,200,50000],shared_thresholds=True,v_loop_distances=prepared,**options)
    assert paired_counts==td.paired_total_distance_mass_batch(q,q,r,r,[alpha]*2,[beta]*2,[alpha]*2,[beta]*2,[0,200,50000],shared_thresholds=True,v_loop_distances=prepared,**options)


@pytest.mark.parametrize('bad',[None,'unknown','TRBV19*99','TRBV19,TRBV7-9'])
def test_prepared_unknown_alleles_still_fail(bad):
    prepared=td.precompute_v_loop_distances()
    with pytest.raises(ValueError,match='unknown human V allele'):
        td.total_distance_count_batch(['CASSLGQAYEQYF'],['CASSLGQAYEQYF'],[bad],['TRBV19'],[[0]],v_loop_distances=prepared)


def test_prepared_input_validation_and_matrix_mismatch():
    prepared=td.precompute_v_loop_distances()
    with pytest.raises(ValueError,match='one locus'):
        td.distance_matrix(['CASSLGQAYEQYF'],['CASSLGQAYEQYF'],['TRBV19'],['TRAV12-2'],v_loop_distances=prepared)
    with pytest.raises(ValueError,match='canonical amino acids'):
        td.distance_matrix(['bad'],['CASSLGQAYEQYF'],['TRBV19'],['TRBV19'],v_loop_distances=prepared)
    with pytest.raises(ValueError,match='matching substitution costs'):
        td.total_distance_count_batch(['CASSLGQAYEQYF'],['CASSLGQAYEQYF'],['TRBV19'],['TRBV19'],[[0]],matrix='blosum62',v_loop_distances=prepared)
