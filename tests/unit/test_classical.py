"""Frozen manuscript coefficients, source-role CSR and native marginal reductions."""
import math

import numpy as np
import polars as pl
import pytest
from seqtree import Index

from vdjmatch.match.classical import alignment_probability, fitted_probability, frozen_model, gene_key
from vdjmatch.runner import classical

native = pytest.mark.skipif(not hasattr(Index, 'alignment_max_batch'),
    reason='requires isolated seqtree frozen-alignment capability; never use rounded fallback')


def frames():
    reference = pl.DataFrame({'cdr3': ['CASSF', 'CASSF', 'CASSF', 'CATSF', 'CASRF', 'WWWWW'],
        'v': ['TRBV19']*6, 'j': ['TRBJ2-1']*6, 'gene': ['TRB']*6,
        'epitope': ['PEP', 'PEP', 'OTHER', 'PEP', 'OTHER', 'OTHER'],
        'mhc_a': ['HLA-A*02:01']*6, 'mhc_b': ['B2M']*6, 'mhc_class': ['MHCI']*6})
    query = pl.DataFrame({'query_id': [7, 9], 'cdr3': ['CASSF']*2, 'v': ['TRBV19*01', 'TRBV19*02'],
        'j': ['TRBJ2-1']*2, 'locus': ['TRB']*2, 'epitope': ['PEP']*2,
        'mhc_a': ['HLA-A*02']*2, 'mhc_b': ['B2M']*2, 'mhc_class': ['MHCI']*2})
    return query, reference


@pytest.mark.parametrize('species,locus,b,c3,gamma', [
    ('human','TRA',.856831490352903,.0529834682577915,-.386684570923303),
    ('human','TRB',.247525036258933,.101047469405546,-.569365478686561),
    ('mouse','TRA',.632369694392929,.0550246931306362,-.475828115592542),
    ('mouse','TRB',1.05348709660423,.107419554461305,-.581148183244143)])
def test_frozen_full_precision_model_and_cloglog(species,locus,b,c3,gamma):
    model = frozen_model(species,locus)
    assert model['similarity'].dtype == np.float64
    assert model['similarity'][0,0] == 1.79412369344457
    assert model['coefficients'] == pytest.approx([b,c3,gamma],abs=1e-15)
    assert model['selection_gap'] == np.max(np.abs(model['similarity']))
    # Selection gap is not included in the final raw-similarity feature.
    expected = -math.expm1(-math.exp(b + .3 - .2 + c3*(8.-12.) + gamma*2))
    assert alignment_probability(model,8.,2,10.,12.,.3,-.2) == pytest.approx(expected,abs=1e-15)
    assert np.diag(model['group_offsets'][0]) == pytest.approx(0.,abs=1e-14)


def test_segment_offsets_are_already_normalized_and_alleles_not_fabricated():
    m=frozen_model('human','TRB');v,j=m['group_names']
    expected=.00478303816812313*(-334)+.0135579063429808*(-15.982522726059)+.000452292355176897*(-19.1639721393585)
    assert m['group_offsets'][0][v.index('TRBV19'),v.index('TRBV5-1')] == pytest.approx(expected,abs=1e-14)
    assert m['group_offsets'][1][j.index('TRBJ2-1'),j.index('TRBJ1-1')] == pytest.approx(.000876685891800022*(-68))
    f=pl.DataFrame({'v':[' TRBV19*02 ','TRAV23/DV6*01','TRBV19,TRBV5-1',None]})
    assert f.select(gene_key('v'))['v'].to_list()==['TRBV19','TRAV23/DV6',None,None]
    assert fitted_probability(np.asarray([-1000.,0.,1000.])) == pytest.approx([0.,-math.expm1(-1),1.])
    with pytest.raises(ValueError,match='finite'):fitted_probability(np.nan)
    with pytest.raises(ValueError,match='gap residues'):alignment_probability(m,1,-1,1,1)


def test_missing_native_api_fails_without_rounded_matrix_fallback(monkeypatch):
    query,reference=frames()
    monkeypatch.setattr(classical,'Index',object)
    with pytest.raises(RuntimeError,match='alignment_max_batch'):
        classical.classical_scores(query,reference)


@pytest.mark.parametrize('options,match', [({'max_indels':65536},'65535'),
    ({'scope':'1,1,0,2','exclude_exact':True},'symmetric'),({'scope':'65536,0,0,65536'},'65535')])
def test_invalid_native_budgets_fail_before_any_model_or_counting(monkeypatch,options,match):
    q,r=frames()
    monkeypatch.setattr(classical,'frozen_model',lambda *a:pytest.fail('model preparation preceded budget validation'))
    with pytest.raises(ValueError,match=match):classical.classical_scores(q,r,**options)


@native
def test_native_csr_multiplicity_crossreactive_roles_and_exact_puncture(monkeypatch):
    query,reference=frames();calls=[]
    class Spy:
        @staticmethod
        def build(sequences,alphabet):return Spy(Index.build(sequences,alphabet))
        def __init__(self,index):self.index=index
        def alignment_max_batch(self,*args,**kwargs):
            calls.append((args,kwargs));return self.index.alignment_max_batch(*args,**kwargs)
        def edit_histogram_batch(self,*args,**kwargs):return self.index.edit_histogram_batch(*args,**kwargs)
    monkeypatch.setattr(classical,'Index',Spy)
    result,report=classical.classical_scores(query,reference,scope='1,0,0,1',max_indels=0,exclude_exact=True)
    assert result['query_id'].to_list()==[7,9] and result.schema['query_id']==pl.Int64
    assert result['v'].to_list()==query['v'].to_list()
    assert len(calls)==1 and calls[0][0][0]==['CASSF']
    _,kw=calls[0]
    assert kw['reference_target_offsets'].tolist()==[0,0,1,2,2]
    assert kw['reference_multiplicities'].tolist()==[2,1]
    assert kw['reference_other_counts'].tolist()==[2,1]
    assert kw['reference_class_offsets'][-1]==5
    assert result['classical_reference_epitope_roles'].to_list()==[3,3]
    assert result['classical_excluded_exact_records'].to_list()==[2,2]
    assert result['classical_full_reference_records'].to_list()==[3,3]
    assert result['classical_reference_records'].to_list()==[1,1]
    assert result['classical_full_reference_junctions'].to_list()==[2,2]
    assert result['classical_reference_junctions'].to_list()==[1,1]
    assert result['classical_hit_records'].to_list()==[1,1]
    assert result['classical_best_reference_junction'].to_list()==['CATSF']*2
    m=frozen_model('human','TRB');lookup={s:i for i,s in enumerate(m['similarity_symbols'])};grid=m['similarity']
    q,r='CASSF','CATSF';raw=sum(grid[lookup[a],lookup[b]] for a,b in zip(q,r))
    selfq=sum(grid[lookup[a],lookup[a]] for a in q);selfr=sum(grid[lookup[a],lookup[a]] for a in r)
    p=-math.expm1(-math.exp(.247525036258933+.101047469405546*(raw-max(selfq,selfr))))
    assert result['classical_pfit'].to_numpy()==pytest.approx([p,p],abs=1e-13)
    assert result['classical_informativeness'].to_numpy()==pytest.approx([math.log10(3)]*2)
    assert result['classical_score'].to_numpy()==pytest.approx([p*math.log10(3)]*2,abs=1e-13)
    assert report['populations'][0]['alignment_batches']==1
    assert report['populations'][0]['reference_epitope_roles']==5
    assert all(not isinstance(x,(list,dict)) for x in report['populations'][0].values())


@native
def test_unknown_genes_nohits_invalid_and_absent_targets_retain_source_rows():
    q,r=frames();base=q.head(1)
    changes=[{}, {'v':'TRBV999'}, {'v':'TRBV19,TRBV5-1'}, {'cdr3':'WWWWW'},
        {'cdr3':'CAS?F'}, {'epitope':'MISSING'}, {'species':'mouse'}, {'locus':'TRD'}]
    rows=[]
    for i,change in enumerate(changes):
        row=base.row(0,named=True)|change|{'query_id':i,'species':change.get('species','human')};rows.append(row)
    result,_=classical.classical_scores(pl.DataFrame(rows),r,scope='1,0,0,1',max_indels=0)
    assert result.height==8
    assert result['classical_status'].to_list()==['scored','unknown_v_gene','unknown_v_gene','no_target_hit',
        'invalid_junction','no_reference','unselected_species','unknown_locus']
    assert result['classical_score'][1] is None and result['classical_score'][2] is None
    assert result['classical_hit_records'][1]==3
    assert result['classical_score'][3]==result['classical_score'][5]==0
    assert result['classical_nw_pairs'][3]==0  # Off-target exact candidate is filtered before NW.
    bad=r.with_columns(pl.lit('TRBV999').alias('v'))
    result,_=classical.classical_scores(base,bad,scope='1,0,0,1',max_indels=0)
    assert result['classical_status'][0]=='unknown_reference_genes'
    assert result['classical_score'][0] is None and result['classical_hit_records'][0]==3


@native
def test_exact_only_scope_and_empty_inputs_have_stable_columns():
    q,r=frames()
    ordinary,p=classical.classical_scores(q,r,scope='0,0,0,0',max_indels=0)
    excluded,_=classical.classical_scores(q,r,scope='0,0,0,0',max_indels=0,exclude_exact=True)
    assert ordinary['classical_hit_records'].to_list()==[2,2]
    assert ordinary['classical_excluded_exact_records'].to_list()==[0,0]
    assert ordinary['classical_reference_records'].to_list()==[3,3]
    assert excluded['classical_status'].to_list()==['no_target_hit']*2
    assert excluded['classical_score'].to_list()==[0.,0.]
    assert p['populations'][0]['histogram_batches']==0
    assert p['populations'][0]['histogram_seconds']==0.
    empty,_=classical.classical_scores(q.head(0),r)
    assert empty.height==0 and set(ordinary.columns)==set(empty.columns)


@native
def test_ambiguous_and_missing_pairs_keep_both_trained_marginals():
    q,r=frames();rows=[]
    for i,(chain,pair) in enumerate([('TRA','p'),('TRA','p'),('TRB','p'),('TRB',None)]):
        row=q.row(0,named=True)|dict(query_id=i,locus=chain,pair_id=pair)
        if chain=='TRA':row.update(v='TRAV1-2',j='TRAJ33')
        rows.append(row)
    alpha=r.with_columns(pl.lit('TRA').alias('gene'),pl.lit('TRAV1-2').alias('v'),pl.lit('TRAJ33').alias('j'))
    result,report=classical.classical_scores(pl.DataFrame(rows),pl.concat([r,alpha]),paired=True,scope='1,0,0,1',max_indels=0)
    assert result.height==4 and result['classical_score'].is_not_null().all()
    assert result['classical_joint_status'].to_list()==['unavailable_ambiguous_pair']*3+['unavailable_missing_pair_id']
    assert len(report['populations'])==2 and all(p['alignment_batches']==1 for p in report['populations'])


@native
def test_reference_annotation_csr_is_sorted_after_metadata_joins(monkeypatch):
    q,r=frames()
    # Two requested target groups, overlapping annotations and differing genes.
    q=pl.concat([q.head(1),q.tail(1).with_columns(pl.lit('OTHER').alias('epitope'))])
    r=r.with_columns(pl.when(pl.col('cdr3')=='CASRF').then(pl.lit('TRBV5-1')).otherwise(pl.col('v')).alias('v'))
    expected,_=classical.classical_scores(q,r,scope='1,0,0,1',max_indels=0)
    join=pl.DataFrame.join
    def unordered(self,other,*args,**kwargs):
        out=join(self,other,*args,**kwargs)
        if set(other.columns)=={'_classical_target_id','epitope'} and '_classical_r' in self.columns:
            return out.reverse()  # A valid unordered join order must not affect CSR.
        return out
    monkeypatch.setattr(pl.DataFrame,'join',unordered)
    actual,_=classical.classical_scores(q,r,scope='1,0,0,1',max_indels=0)
    assert actual.equals(expected)


@native
def test_all_missing_species_and_exhausted_reference_are_retained():
    q,r=frames()
    unknown=q.head(1).with_columns(pl.lit(None).alias('species'))
    result,_=classical.classical_scores(unknown,r)
    assert result.height==1 and result['species'][0] is None
    assert result['classical_status'][0]=='unselected_species' and result['classical_score'][0] is None
    sole=r.filter((pl.col('cdr3')=='CASSF')&(pl.col('epitope')=='PEP'))
    result,_=classical.classical_scores(q.head(1),sole,scope='1,0,0,1',max_indels=0,exclude_exact=True)
    assert result['classical_full_reference_records'][0]==2
    assert result['classical_reference_records'][0]==result['classical_reference_epitope_roles'][0]==0
    assert result['classical_status'][0]=='no_reference' and result['classical_score'][0]==0
