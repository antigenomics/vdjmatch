"""Full-junction total-distance counts, independent dense and weighted geometry oracles."""
import numpy as np
import pytest
from vdjmatch.match import tcrdist


@pytest.mark.parametrize('species,genes',[
 ('human',['TRBV19*01','TRBV7-9*01','TRBV5-1*01','TRBV20-1*01']),
 ('mouse',['TRAV1*01','TRAV2*01','TRAV4D-2*01','TRAV5-2*01']),
])
@pytest.mark.parametrize('ctrim',[2,3])
def test_counts_equal_dense_distance(species,genes,ctrim):
    q=['CTSKGMKF','CPNKNINGEEPQLCF','CGNKELKNMHITGNARHAGFF','CSTNRTCCVAFAATHEDQYMWMNQTSYRWF']
    r=q+[q[0],q[0][:-3]+'A'+q[0][-2:]];rv=genes+[genes[0],genes[0]]
    scores=tcrdist.distance_matrix(q,r,genes,rv,species=species,ctrim=ctrim)
    cuts=[0,12,25,50,100,400,25,-1];thresholds=[[v*200 for v in cuts] for _ in q]
    for threads in (1,4):
        for excluded in (False,True):
            out=tcrdist.total_distance_count_batch(q,r,genes,rv,thresholds,species=species,ctrim=ctrim,threads=threads,exclude_exact=excluded)
            assert out==[[sum(scores[i,j]<=v and not(excluded and q[i]==r[j]) for j in range(len(r))) for v in cuts] for i in range(len(q))]


def test_full_identity_same_cores_v_offsets_and_aliases():
    q='CASSLGQAYEQYF';r='A'+q[1:]
    opts=dict(query_v=['TRBV19*01'],reference_v=['TRBV19*01','TRBV19*01','TRBV7-9*01'])
    assert tcrdist.total_distance_count_batch([q],[q,r,r],thresholds=[[0,100000]],exclude_exact=True,**opts)==[[1,2]]
    assert tcrdist.total_distance_count_batch([q],[q],['TRBV19*01'],['TRBV7-9*01'],[[0]],exclude_exact=False)==[[0]]
    # Resolve physical aliases on both axes through the existing model resolver.
    assert tcrdist.total_distance_count_batch([q],[r],['TRAV14'],['TRAV14/DV4*01'],[[0]])==[[1]]


def weighted_pair(q,r,ctrim):
    from vdjmatch.match.regions import significance_weights
    from seqtree import SubstitutionMatrix
    b=SubstitutionMatrix.blosum62()
    longer,shorter=(q,r) if len(q)>=len(r) else (r,q)
    n=len(shorter);delta=len(longer)-n
    weights=[max(1,round(100*w)) for w in significance_weights(len(longer))]
    def sub(a,c):return 0 if a==c else max(0,min(4,4-b.similarity(a,c)))
    if not delta:return sum(6*sub(a,c)*w for a,c,w in zip(q,r,weights))
    low,high=5,n-5
    while low>high:low-=1;high+=1
    return 2400*delta+min(sum(6*sub(shorter[j],longer[j if j<i else j+delta])*weights[j if j<i else j+delta] for j in range(n)) for i in range(low,high+1))


def test_significance_weighted_exhaustive_oracle():
    ctrim=2
    q=['CASSLGQAYEQYF','CASSLGQAPAYEQYF','CTSKGMKF']
    r=q+['AASSLGQAYEQYF','CASSLGRAYEQYF'];qv=['TRBV19*01']*3;rv=qv+['TRBV19*01','TRBV7-9*01']
    loops=[int(tcrdist.distance_matrix([q[0]],[q[0]],[qv[0]],[v])[0,0])*200 for v in rv]
    scores=np.array([[weighted_pair(a,b,ctrim)+loops[j] for j,b in enumerate(r)] for a in q])
    thresholds=[[-1,0,1000,3000,10000,*scores[i].tolist()] for i in range(len(q))]
    for threads in (1,4):
        got=tcrdist.total_distance_count_batch(q,r,qv,rv,thresholds,ctrim=ctrim,threads=threads,position_weighting='significance',exclude_exact=True)
        assert got==[[sum(a!=b and scores[i,j]<=v for j,b in enumerate(r)) for v in thresholds[i]] for i,a in enumerate(q)]


@pytest.mark.parametrize('kwargs',[{'ctrim':4},{'position_weighting':'fitted'},{'threads':0},{'ctrim':3,'position_weighting':'significance'}])
def test_count_invalid_geometry(kwargs):
    with pytest.raises(ValueError):tcrdist.total_distance_count_batch(['CASSLGQAYEQYF'],['CASSLGQAYEQYF'],['TRBV19'],['TRBV19'],[[0]],**kwargs)


def test_count_empty_axes_and_input_disposition():
    q='CASSLGQAYEQYF';v='TRBV19*01'
    assert tcrdist.total_distance_count_batch([],[],[],[],[])==[]
    assert tcrdist.total_distance_count_batch([q],[],[v],[],[[0,100]])==[[0,0]]
    assert tcrdist.total_distance_count_batch([],[q],[],[v],[])==[]
    with pytest.raises(ValueError,match='unknown human V'):
        tcrdist.total_distance_count_batch([q],[q],[None],[v],[[0]])
    with pytest.raises(ValueError,match='one locus'):
        tcrdist.total_distance_count_batch([q],[q],[v],['TRAV1-1'],[[0]])
    with pytest.raises(ValueError,match='requires a V call'):
        tcrdist.total_distance_count_batch([q],[q],[],[v],[[0]])
    with pytest.raises(ValueError,match='threshold'):
        tcrdist.total_distance_count_batch([q],[q],[v],[v],[[0.5]])


def test_count_capability_gate(monkeypatch):
    monkeypatch.setattr(tcrdist.gapblock,'count_batch',lambda queries,refs,thresholds:None)
    with pytest.raises(RuntimeError,match='group_distances'):
        tcrdist.total_distance_count_batch([],[],[],[],[])
