"""Exact compact arithmetic against a standalone frozen expression oracle."""
import numpy as np
import pytest
from vdjmatch.evalue.kernel import KERNEL_SCORERS, _kernel_scores


@pytest.mark.parametrize('scorer',KERNEL_SCORERS)
def test_owned_compact_buffers_preserve_old_arithmetic(scorer):
    rng=np.random.default_rng(104)
    shells=rng.integers(0,2**30,size=(73,153)).astype(float)
    same=np.floor(shells*rng.random(shells.shape))
    control=np.cumsum(rng.integers(0,2**20,size=shells.shape),axis=1).astype(float)
    shells[1]=shells[0];same[1]=same[0];control[1]=control[0]
    shells[2]=0;same[2]=0;control[2]=0
    N=shells.sum(axis=1);M=control[:,-1];N[3]=0;M[4]=0
    t=np.arange(153,dtype=float)*17
    name,mode=scorer.split('-')
    mass=shells if name=='all' else same if name=='vonly' else .25*shells+.75*same
    denominator=1+control if mode=='bounded' else np.maximum(control,.01)
    expected=np.einsum('ij,ij->i',mass,np.exp(-t/585)[None,:]/denominator)/np.where(N>0,N,1)
    expected=np.where((N>0)&(M>0),expected,np.nan)
    preserved=[x.copy() for x in (shells,same,control)]
    actual=_kernel_scores(shells.copy(),same.copy(),control.copy(),N,M,t,585,.25,.01,
        scorer,False,_reuse_buffers=True)['kernel_score']
    np.testing.assert_array_equal(actual.view(np.uint64),expected.view(np.uint64))
    assert actual[0]==actual[1]
    public=_kernel_scores(shells,same,control,N,M,t,585,.25,.01,scorer,True)['kernel_score']
    np.testing.assert_array_equal(public.view(np.uint64),expected.view(np.uint64))
    for original,copy in zip((shells,same,control),preserved):np.testing.assert_array_equal(original,copy)
