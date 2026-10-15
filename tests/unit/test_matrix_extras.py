"""Recorded extra matrix Gram penalties and shared distance units, synthetic only."""
import numpy as np
import pytest
from seqtree import SubstitutionMatrix, amino_acids
from vdjmatch.match.scoring import load_vdjam
from vdjmatch.match import geometry


def matrix_source(name):
    return load_vdjam() if name=='vdjam' else (
        SubstitutionMatrix.structural() if name=='mj' else SubstitutionMatrix.pam250())


@pytest.mark.parametrize('name,maximum,ac,wd',[('pam250',45,18,35),('mj',20,8,14),('vdjam',1885,655,624)])
def test_recorded_unclipped_gram_penalties(name,maximum,ac,wd):
    matrix=matrix_source(name);symbols=amino_acids()
    expected=np.asarray([[max(0,matrix.similarity(a,a)+matrix.similarity(b,b)-2*matrix.similarity(a,b))
        for b in symbols] for a in symbols],dtype=np.int32)
    got=geometry._substitution_costs('full',name)
    np.testing.assert_array_equal(got,expected)
    assert got.dtype==np.int32 and got.max()==maximum
    assert got[symbols.index('A'),symbols.index('C')]==ac
    assert got[symbols.index('W'),symbols.index('D')]==wd
    np.testing.assert_array_equal(got,got.T)
    assert not np.diag(got).any()
    assert geometry.distance_units(got)[0]==maximum
    with pytest.raises(ValueError,match='full alphabet'):
        geometry._substitution_costs('imgt',name)


@pytest.mark.parametrize('name',['pam250','mj','vdjam'])
def test_gap_and_stem_keep_the_recorded_shared_matrix_units(name):
    costs=geometry._substitution_costs('full',name)
    q='CAA'+'AAAA'+'AAF';r='CAA'+'AAAAA'+'AAF'
    gap=geometry.distance_matrix([q],[r],['TRBV19'],['TRBV19'],ctrim=3,
        matrix=name,v_loop_weight=0,stem_depth=0)
    assert gap[0,0]==3*int(costs.max())
    query='CASSLGQAYEQYF';reference='WASSLGQAYEQYW';symbols=amino_acids();index={a:i for i,a in enumerate(symbols)}
    expected=round(sum(max(int(costs[index[query[i]],index[reference[i]]]),
        int(costs[index[query[-1-i]],index[reference[-1-i]]])) for i in range(4))/4)
    assert geometry.stem_penalty([query],[reference],depth=4,costs=costs)[0,0]==expected


@pytest.mark.parametrize('name',['pam250','mj','vdjam'])
def test_v_loop_cost_uses_the_same_unclipped_matrix(name):
    costs=geometry._substitution_costs('full',name);symbols=amino_acids();index={a:i for i,a in enumerate(symbols)}
    model=geometry.load_v_loops();q='CASSLGQAYEQYF'
    loops=lambda allele: model[allele]['cdr1']+model[allele]['cdr2']+model[allele]['cdr25']
    expected=sum(int(costs[index[a],index[b]]) for a,b in zip(loops('TRBV19*01'),loops('TRBV7-9*01'))
        if a in 'ACDEFGHIKLMNPQRSTVWY' and b in 'ACDEFGHIKLMNPQRSTVWY')
    got=geometry.distance_matrix([q],[q],['TRBV19'],['TRBV7-9'],matrix=name,ctrim=3,stem_depth=0)
    assert got[0,0]==expected
    prepared=geometry.precompute_v_loop_distances(matrix=name)
    np.testing.assert_array_equal(got,geometry.distance_matrix([q],[q],['TRBV19'],['TRBV7-9'],
        matrix=name,ctrim=3,stem_depth=0,v_loop_distances=prepared))


@pytest.mark.parametrize('name',['pam250','mj','vdjam'])
def test_dense_and_native_counts_preserve_high_costs_and_exact_exclusion(name):
    q=['CASSLGQAYEQYF'];r=q+['CASSLWQAYEQYF','WASSLGQAYEQYW','CASSLGQAAYEQYF']
    qv=['TRBV19'];rv=['TRBV7-9','TRBV19','TRBV7-9','TRBV19']
    distance=geometry.distance_matrix(q,r,qv,rv,ctrim=3,matrix=name,stem_depth=4)
    bounds=sorted({0,*[200*int(x)+delta for x in distance[0] for delta in (-1,0,1)]})
    for exclude in (False,True):
        counts=geometry.total_distance_count_batch(q,r,qv,rv,[bounds],ctrim=3,matrix=name,
            stem_depth=4,exclude_exact=exclude)
        expected=[[sum(200*int(distance[0,j])<=bound and (not exclude or q[0]!=r[j])
            for j in range(len(r))) for bound in bounds]]
        assert counts==expected


def test_original_and_selected_matrices_remain_exact():
    symbols=amino_acids();blosum=SubstitutionMatrix.blosum62()
    original=np.asarray([[0 if a==b else max(0,min(4,4-blosum.similarity(a,b)))
        for b in symbols] for a in symbols],dtype=np.int32)
    selected=np.asarray([[blosum.penalty(a,b) for b in symbols] for a in symbols],dtype=np.int32)
    np.testing.assert_array_equal(geometry._substitution_costs(),original)
    np.testing.assert_array_equal(geometry._substitution_costs('full','blosum62'),selected)
    assert geometry._substitution_costs('imgt').max()==4
