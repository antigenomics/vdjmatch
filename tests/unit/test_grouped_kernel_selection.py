"""Independent synthetic edge oracles for the grouped compact kernel route."""
import argparse
from itertools import product

import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal

from vdjmatch.evalue.kernel import KERNEL_SCORERS, kernel_ladder


def run(args):
    from vdjmatch.cli.tcrdist import register
    parser=argparse.ArgumentParser()
    register(parser.add_subparsers(),name='annotate')
    a=parser.parse_args(args)
    return a.func(a)


def sequences():
    return ['CASSL'+a+b+'AYEQYF' for a,b in list(product('ACDEFGHIKLMNPQRSTVWY',repeat=2))[:70]]


def fixture(tmp_path,paired=False):
    sample,reference=tmp_path/'q.tsv',tmp_path/'r.tsv'
    seqs=sequences()+[sequences()[0]]*5
    meta=dict(species='HomoSapiens',mhc_class='MHCI',mhc_a='HLA-A*02:01',mhc_b='B2M')
    queries=[]
    for i,s in enumerate(seqs):
        task=dict(meta,epitope='missing' if i==69 else 'B' if i%2 else 'A')
        queries.append(dict(sequence_id='b'+str(i),clone_id='p'+str(i),junction_aa=s,
            v_call='TRBV7-9' if i%3 else 'TRBV19',locus='TRB',**task))
        if paired:queries.append(dict(sequence_id='a'+str(i),clone_id='p'+str(i),
            junction_aa='CAVRDSNYQLIW' if i%2 else 'CAVRDSNYQLVW',v_call='TRAV12-2',locus='TRA',**task))
    pl.DataFrame(queries).write_csv(sample,separator='\t')
    refs=[]
    # Shared junction/V search keys retain multiple labels and restrictions.
    # Cross-V exact junctions must puncture the full junction, irrespective of V.
    labels=[(seqs[0],'TRBV19','A','HLA-A*02:01'),
            (seqs[0],'TRBV7-9','B','HLA-A*02:01'),
            (seqs[0],'TRBV19','B','HLA-A*02:01'),
            (seqs[2],'TRBV19','A','HLA-A*03:01'),
            (seqs[4],'TRBV7-9','C','HLA-A*02:01'),
            ('CWWWWWWWWWWWF','TRBV19','C','HLA-A*02:01')]
    for i,(s,v,e,mhc) in enumerate(labels):
        task=dict(meta,epitope=e,mhc_a=mhc)
        refs.append(dict(cdr3=s,v=v,gene='TRB',complex_id=str(i+1),**task))
        if paired:refs.append(dict(cdr3='CAVRDSNYQLIW' if i%2 else 'CAVRDSNYQLVW',
            v='TRAV12-2',gene='TRA',complex_id=str(i+1),**task))
    pl.DataFrame(refs).write_csv(reference,separator='\t')
    return sample,reference


@pytest.mark.parametrize('scorer',KERNEL_SCORERS)
@pytest.mark.parametrize('exclude',[False,True])
@pytest.mark.parametrize('paired',[False,True])
def test_grouped_matches_detailed_for_all_scorers_linkage_and_puncture(tmp_path,scorer,exclude,paired):
    sample,reference=fixture(tmp_path,paired)
    args=['annotate',str(sample),'--vdjdb',str(reference),'--locus','paired' if paired else 'TRB',
          '--kernel-scorer',scorer]
    if exclude:args+=['--exclude-exact']
    for compact in [False,True]:
        assert run(args+(['--selection-only'] if compact else [])+['--output-prefix',str(tmp_path/str(compact))])==0
    suffixes=['','.TRA','.TRB'] if paired else ['']
    for suffix in suffixes:
        full=pl.read_csv(str(tmp_path/'False')+suffix+'.selection.tsv',separator='\t')
        grouped=pl.read_csv(str(tmp_path/'True')+suffix+'.selection.tsv',separator='\t')
        assert_frame_equal(full,grouped,check_dtypes=False,check_exact=True)
    assert_frame_equal(pl.read_csv(str(tmp_path/'False')+'.queries.tsv',separator='\t'),
                       pl.read_csv(str(tmp_path/'True')+'.queries.tsv',separator='\t'),check_dtypes=False,check_exact=True)


@pytest.mark.parametrize('selectors',[['--targets-from-sample'],['--epitope','missing'],['--epitope','A','--epitope','B']])
@pytest.mark.parametrize('exclude',[False,True])
def test_assayed_masks_absent_and_explicit_target_scopes(tmp_path,selectors,exclude):
    sample,reference=fixture(tmp_path)
    args=['annotate',str(sample),'--vdjdb',str(reference),'--locus','TRB',*selectors]
    if exclude:args+=['--exclude-exact']
    for compact in [False,True]:
        assert run(args+(['--selection-only'] if compact else [])+['--output-prefix',str(tmp_path/str(compact))])==0
    for name in ['queries','selection']:
        assert_frame_equal(pl.read_csv(str(tmp_path/'False')+'.'+name+'.tsv',separator='\t'),
                           pl.read_csv(str(tmp_path/'True')+'.'+name+'.tsv',separator='\t'),check_dtypes=False,check_exact=True)


@pytest.mark.parametrize('scorer',KERNEL_SCORERS)
def test_extras_false_retains_identical_score_on_strided_curves(scorer):
    rng=np.random.default_rng(2)
    n=np.cumsum(rng.integers(0,3,size=(97,11)),axis=1)
    same=np.cumsum(rng.integers(0,2,size=(97,11)),axis=1)
    n+=same
    m=np.cumsum(rng.integers(0,4,size=(97,11)),axis=1)
    # Deliberately retain strided columns like the grouped chunk reducer.
    n,same,m=(x[:,::2] for x in (n,same,m))
    cuts=np.arange(6)*4
    args=(n,same,m,np.full(97,60),np.full(97,60),cuts)
    full=kernel_ladder(*args,kernel_scale=13.,scorer=scorer)
    selected=kernel_ladder(*args,kernel_scale=13.,scorer=scorer,extras=False)
    assert set(selected)=={'kernel_score','kernel_'+scorer.replace('-','_')}
    np.testing.assert_array_equal(selected['kernel_score'],full['kernel_score'])


def test_repeated_source_observations_do_not_expand_all_labels_before_reduction(tmp_path,monkeypatch):
    import vdjmatch.evalue.kernel as kernel
    sample,reference=fixture(tmp_path)
    raw=pl.read_csv(sample,separator='\t').head(1)
    raw=pl.concat([raw]*1000).with_row_index('_source').with_columns(
        pl.col('_source').cast(pl.String).alias('sequence_id')).drop('_source')
    raw.write_csv(sample,separator='\t')
    original=kernel._compact_kernel_selection
    seen=[]
    def bounded(scores,exact,*layout):
        seen.append(scores.size)
        assert scores.shape[0]<=64 and scores.shape[1]<=4
        return original(scores,exact,*layout)
    monkeypatch.setattr(kernel,'_compact_kernel_selection',bounded)
    assert run(['annotate',str(sample),'--vdjdb',str(reference),'--locus','TRB',
        '--selection-only','--output-prefix',str(tmp_path/'out')])==0
    assert seen and max(seen)<=4
    out=pl.read_csv(str(tmp_path/'out')+'.selection.tsv',separator='\t')
    assert out.height==1000


def test_assayed_duplicate_masks_bound_expansion_and_preserve_exact_ties(monkeypatch):
    import vdjmatch.evalue.kernel as kernel
    T=128;Q=700
    queries=pl.DataFrame({'query_id':range(Q),'cdr3':['CASSLGQAYEQYF']*Q,'_allele':['TRBV19*01']*Q})
    refs=pl.DataFrame({'cdr3':['CASSLGQAYEQYF']*T,'_allele':['TRBV19*01']*T,
                       'epitope':['E'+str(i) for i in range(T)]})
    groups=[({'epitope':e},queries,refs.filter(pl.col('epitope')==e)) for e in refs['epitope']]
    original=kernel.select_kernel;seen=[]
    def bounded(q,s):
        seen.append(s.height)
        assert s.height<=65536
        return original(q,s)
    monkeypatch.setattr(kernel,'select_kernel',bounded)
    state,hits=kernel.grouped_kernel_selection(queries,refs,groups,query_keys=['cdr3','_allele'],
        search_keys=['cdr3','_allele'],paired=False,cutoffs=[0,4],radius=4,kernel_scale=22.5,
        scorer='v-bounded',assayed=True,matrix='blosum62',ctrim=3,stem_depth=4)
    assert seen==[65536,(Q-512)*T]
    assert state.height==Q and hits==set(range(Q))
    assert state['exact_peptides'].to_list()==[T]*Q
    assert state['tied_peptides'].to_list()==[T]*Q


def test_assayed_missing_peptide_or_mhc_retains_detailed_dispositions(tmp_path):
    sample,reference=fixture(tmp_path)
    raw=pl.read_csv(sample,separator='\t').head(4)
    raw=raw.with_row_index('_i').with_columns(
        pl.when(pl.col('_i')==0).then(None).when(pl.col('_i')==3).then(pl.lit('missing')).otherwise(pl.col('epitope')).alias('epitope'),
        pl.when(pl.col('_i')==1).then(None).otherwise(pl.col('mhc_a')).alias('mhc_a'),
        pl.when(pl.col('_i')==2).then(None).otherwise(pl.col('mhc_b')).alias('mhc_b')).drop('_i')
    raw.write_csv(sample,separator='\t')
    args=['annotate',str(sample),'--vdjdb',str(reference),'--locus','TRB','--targets-from-sample']
    for compact in [False,True]:
        assert run(args+(['--selection-only'] if compact else [])+['--output-prefix',str(tmp_path/str(compact))])==0
    for name in ['queries','selection']:
        full=pl.read_csv(str(tmp_path/'False')+'.'+name+'.tsv',separator='\t')
        compact=pl.read_csv(str(tmp_path/'True')+'.'+name+'.tsv',separator='\t')
        assert_frame_equal(full,compact,check_dtypes=False,check_exact=True)
        if name=='selection':
            assert full['tied_peptides'].to_list()==[0]*4
            assert full['exact_peptides'].to_list()==[0]*4
    dispositions=pl.read_csv(str(tmp_path/'True')+'.queries.tsv',separator='\t')
    assert dispositions['status'].to_list()==['invalid_target']*3+['no_reference']
