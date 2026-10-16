"""Independent small dense/serial oracle for joint V-conditioned kernel evidence."""
import polars as pl
import pytest
from vdjmatch.evalue.epitope_rules import reference_epitope_rules,MOTIFS,IMGT_CLASSES


def tokens(s,mode):
    s=s[3:-3]
    if 'imgt' in mode:s=''.join(IMGT_CLASSES[c] for c in s)
    k=4 if 'wild' in mode else 3
    return {s[i:i+2]+'X'+s[i+3] for i in range(len(s)-3)}|{s[i]+'X'+s[i+2:i+4] for i in range(len(s)-3)} if k==4 else {s[i:i+3] for i in range(len(s)-2)}


@pytest.mark.parametrize('paired',[False,True])
def test_rules_dense_oracle_with_fixed_models_and_hit_exclusion(paired):
    from vdjmatch.match.geometry import distance_matrix,paired_distance_matrix
    r=pl.DataFrame({'seq0':['CASSLGQAYEQYF','CASSLGRAYEQYF','CASSLGRAYEQYF','CASSWGRAYEQYF'],'v0':['TRBV19*01','TRBV19*01','TRBV20-1*01','TRBV19*01']})
    b=r.with_columns(pl.Series('v0',['TRBV20-1*01','TRBV19*01','TRBV20-1*01','TRBV20-1*01']))
    if paired:
        r=r.with_columns(pl.col('seq0').alias('seq1'),pl.col('v0').alias('v1'),pl.lit('TRAV1-1*01').alias('v0'))
        b=b.with_columns(pl.col('seq0').alias('seq1'),pl.col('v0').alias('v1'),pl.lit('TRAV1-1*01').alias('v0'))
    q=r.head(1).with_row_index('query_id');axes=range(2 if paired else 1)
    out=reference_epitope_rules(q,r,b,paired=paired,exclude_exact=True)
    oracle={}
    for side,pop in [('reference',r),('control',b)]:
        keep=[row for row in pop.to_dicts() if any(row[f'seq{i}']!=q[f'seq{i}'][0] for i in axes)]
        model_same=[row for row in pop.to_dicts() if all(row[f'v{i}'].split('*')[0]==q[f'v{i}'][0].split('*')[0] for i in axes)]
        same=[row for row in keep if all(row[f'v{i}'].split('*')[0]==q[f'v{i}'][0].split('*')[0] for i in axes)]
        assert out['rule_'+side+'_population'][0]==pop.height
        assert out['rule_'+side+'_v_count'][0]==len(model_same)
        assert out['rule_'+side+'_excluded_hits'][0]==pop.height-len(keep)
        if paired:
            d=paired_distance_matrix(q['seq0'],q['seq1'],[row['seq0'] for row in same],[row['seq1'] for row in same],q['v0'],q['v1'],[row['v0'] for row in same],[row['v1'] for row in same],ctrim=3,v_loop_weight=0)[0]
        else:d=distance_matrix(q['seq0'],[row['seq0'] for row in same],q['v0'],[row['v0'] for row in same],ctrim=3,v_loop_weight=0)[0]
        mass=sum(max(0,1-x/90) for x in d)
        assert out['rule_'+side+'_spatial_mass'][0]==pytest.approx(mass)
        oracle[side]={'spatial':mass,'distances':d,'same':same,'model_same':model_same,
                      'excluded_v':len(model_same)-len(same)}
        for mode in MOTIFS:
            mass=0
            for row in same:
                kernel=1
                for i in axes:
                    qt=tokens(q[f'seq{i}'][0],mode);rt=tokens(row[f'seq{i}'],mode)
                    kernel*=len(qt&rt)/len(qt)
                mass+=kernel
            assert out['rule_'+side+'_'+mode+'_mass'][0]==pytest.approx(mass)
            oracle[side][mode]=mass
    n=out['rule_reference_v_count'][0];m=out['rule_control_v_count'][0]
    # Hand factorization: V probability is multiplied once, CDR3 is same-V only.
    for kind in ('spatial',*MOTIFS):
        rm=oracle['reference'][kind];bm=oracle['control'][kind]
        prior=(rm+oracle['reference']['excluded_v']+bm+oracle['control']['excluded_v']+.5)/(n+m+1)
        expected=out['rule_v_score'][0]*(rm+prior)/(n+1)/((bm+prior)/(m+1))
        assert out['rule_'+kind+'_score'][0]==pytest.approx(expected)
    neutral=reference_epitope_rules(q,r,r,paired=paired)
    assert all(neutral[c][0]==pytest.approx(1) for c in neutral.columns if c.endswith('_score') and c!='rule_rank_score' and not c.startswith('rule_enriched_'))
    rd=oracle['reference']['distances'];bd=oracle['control']['distances']
    reciprocal=[sum(1/k for k in range(1+sum(bd<x),2+sum(bd<=x)))/(1+sum(bd<=x)-sum(bd<x)) for x in rd]
    expected=sum(reciprocal)/n*(m+1)/sum(1/k for k in range(1,m+2))
    assert out['rule_cdr3_rank'][0]==pytest.approx(expected)
    assert out['rule_rank_score'][0]==pytest.approx(expected*out['rule_v_score'][0])
    for mode in MOTIFS:
        import itertools
        qt=set(itertools.product(*[tokens(q[f'seq{i}'][0],mode) for i in axes]))
        presence={side:[set(itertools.product(*[tokens(row[f'seq{i}'],mode) for i in axes])) for row in oracle[side]['model_same']] for side in ('reference','control')}
        enrichment=sum(max(0,sum(t in s for s in presence['reference'])/n-sum(t in s for s in presence['control'])/m) for t in qt)/len(qt)
        assert out['rule_'+mode+'_enrichment'][0]==pytest.approx(enrichment)
        assert out['rule_enriched_'+mode+'_score'][0]==pytest.approx(out['rule_rank_score'][0]*(1+enrichment))
        assert neutral['rule_'+mode+'_enrichment'][0]==0
    ordinary=reference_epitope_rules(q,r,b,paired=paired)
    fixed=[c for c in ordinary.columns if c.endswith(('_population','_v_count','_probability','_enrichment')) or c=='rule_v_score']
    assert out.select(fixed).equals(ordinary.select(fixed))
    extra=q.with_columns(pl.lit(9).cast(pl.UInt32).alias('query_id'),pl.lit('TRBV20-1*01').alias('v1' if paired else 'v0'))
    assert out.equals(reference_epitope_rules(pl.concat([q,extra]),r,b,paired=paired,exclude_exact=True).head(1))
    longer=extra.with_columns(*[pl.lit('CAS'+'W'*35+'QYF').alias(f'seq{i}') for i in axes])
    assert out.equals(reference_epitope_rules(pl.concat([q,longer]),r,b,paired=paired,exclude_exact=True).head(1))


@pytest.mark.parametrize('paired',[False,True])
def test_excluded_only_hits_never_change_models_or_return_as_censored_rank(paired):
    ref=pl.DataFrame({'seq0':['CASSLGQAYEQYF'],'v0':['TRBV19*01']})
    if paired:
        ref=ref.with_columns(pl.col('seq0').alias('seq1'),pl.col('v0').alias('v1'),pl.lit('TRAV1-1*01').alias('v0'))
    query=ref.with_row_index('query_id')
    masked=reference_epitope_rules(query,ref,ref,paired=paired,exclude_exact=True)
    assert masked['rule_reference_population'].item()==masked['rule_reference_v_count'].item()==1
    assert masked['rule_v_supported'].item() and masked['rule_v_score'].item()==1
    assert masked['rule_cdr3_rank'].item()==masked['rule_rank_score'].item()==0
    for mode in MOTIFS:
        assert masked['rule_reference_'+mode+'_mass'].item()==0
        assert masked['rule_'+mode+'_enrichment'].item()==0


def test_paired_exact_filter_retains_references_sharing_only_one_chain():
    reference=pl.DataFrame({'seq0':['CAVRDSNYQLIW','CAVRDSNYQLIW','CAVRDTNYQLIW'],
        'seq1':['CASSLGQAYEQYF','CASSLGRAYEQYF','CASSLGQAYEQYF'],
        'v0':['TRAV1-1*01']*3,'v1':['TRBV19*01']*3})
    query=reference.head(1).with_row_index('query_id')
    out=reference_epitope_rules(query,reference,reference,paired=True,exclude_exact=True)
    assert out['rule_reference_population'].item()==3
    assert out['rule_reference_excluded_hits'].item()==out['rule_reference_removed_v'].item()==1
    assert out['rule_reference_spatial_mass'].item()>0 and out['rule_cdr3_rank'].item()>0


def test_filtering_only_same_v_hit_does_not_switch_to_global_v_fallback():
    reference=pl.DataFrame({'seq0':['CASSLGQAYEQYF','CASSLGRAYEQYF'],
        'v0':['TRBV19*01','TRBV20-1*01']})
    control=reference.with_columns(pl.Series('seq0',['CASSWGRAYEQYF','CASSLGRAYEQYF']))
    query=reference.head(1).with_row_index('query_id')
    out=reference_epitope_rules(query,reference,control,exclude_exact=True)
    assert out['rule_reference_v_count'].item()==out['rule_control_v_count'].item()==1
    assert out['rule_v_supported'].item() and out['rule_cdr3_rank'].item()==0


def test_unsupported_v_is_neutral_and_keeps_global_multireference_rank():
    r=pl.DataFrame({'seq0':['CASSLGRAYEQYF','CASSWGRAYEQYF'],'v0':['TRBV19*01']*2})
    b=pl.concat([r,r.with_columns(pl.lit('TRBV20-1*01').alias('v0'))])
    q=r.head(1).with_columns(pl.lit('CASSLGQAYEQYF').alias('seq0'),pl.lit('TRBV7-9*01').alias('v0')).with_row_index('query_id')
    out=reference_epitope_rules(q,r,b)
    assert not out['rule_v_supported'][0] and out['rule_v_score'][0]==1
    supported=reference_epitope_rules(q.with_columns(pl.lit('TRBV19*01').alias('v0')),r,b)
    assert out['rule_cdr3_rank'][0]>0
    assert all(out['rule_enriched_'+mode+'_score'][0]==out['rule_rank_score'][0] for mode in MOTIFS)
    # Changing the unseen V label changes neither fallback nor motif evidence.
    changed=reference_epitope_rules(q.with_columns(pl.lit('TRBV6-5*01').alias('v0')),r,b)
    assert out.equals(changed)
    assert supported['rule_v_supported'][0]


def test_rank_keeps_multiple_references_beyond_local_radius():
    q=pl.DataFrame({'query_id':[0],'seq0':['CAS'+'A'*20+'QYF'],'v0':['TRBV19*01']})
    r=pl.DataFrame({'seq0':['CAS'+'W'*20+'QYF','CAS'+'Y'*20+'QYF'],'v0':['TRBV19*01']*2})
    b=pl.DataFrame({'seq0':['CAS'+'A'*19+'WQYF','CAS'+'A'*18+'WWQYF'],'v0':['TRBV19*01']*2})
    out=reference_epitope_rules(q,r,b)
    assert out['rule_reference_spatial_mass'][0]==0
    # Both distant references follow both controls, tied at background rank3.
    assert out['rule_cdr3_rank'][0]==pytest.approx(1/(1+1/2+1/3))


@pytest.mark.parametrize('paired',[False,True])
def test_rules_reuse_grouped_curves_and_preserve_chunk_identity(monkeypatch,paired):
    from vdjmatch.match import geometry
    r=pl.DataFrame({'seq0':['CASSLGQAYEQYF','CASSLGRAYEQYF','CASSWGRAYEQYF'],
        'v0':['TRBV19*01','TRBV19*02','TRBV20-1*01']})
    b=r.head(2)
    if paired:
        r=r.rename({'seq0':'seq1','v0':'v1'}).with_columns(pl.lit('CAVRDSNYQLIW').alias('seq0'),pl.lit('TRAV12-2*01').alias('v0'))
        b=r.head(2)
    q=r.with_row_index('query_id')
    expected=reference_epitope_rules(q,r,b,paired=paired,exclude_exact=True,radius=89.7)
    calls=[];native=geometry.grouped_total_distance_count_batch
    def one_row(*args,**kwargs):
        calls.append((args,kwargs))
        return native(*args,**(kwargs|{'chunk_size':1}))
    monkeypatch.setattr(geometry,'grouped_total_distance_count_batch',one_row)
    actual=reference_epitope_rules(q,r,b,paired=paired,exclude_exact=True,radius=89.7)
    assert actual.equals(expected)
    # Conditional curves supply mass and rank; only the unsupported gene needs
    # an unconditional fallback. Shared reference/control keys remain shared.
    assert len(calls)==2
    assert calls[0][1]['same_v_gene'] and not calls[1][1]['same_v_gene']
    assert calls[0][0][1].height==r.height
    assert calls[1][0][0].height==1


def test_motif_classes_wildcards_and_short_sequences():
    from vdjmatch.evalue.epitope_rules import _tokens
    q=pl.DataFrame({'_id':[0,1],'seq0':['CASAILVSTDEQYF','CASSQYF'],'v0':['TRBV19','TRBV19']})
    for mode in MOTIFS:
        t=_tokens(q,range(1),mode)
        assert set(t.filter(pl.col('_id')==0)['t0'])==tokens(q['seq0'][0],mode)
        assert t.filter(pl.col('_id')==1).height==0
    assert IMGT_CLASSES['L']==IMGT_CLASSES['I'] and IMGT_CLASSES['F']!=IMGT_CLASSES['Y']


@pytest.mark.parametrize('diagnostic',['epitope-rules','chemistry'])
def test_cli_rules_preserve_existing_fields_and_empty_schema(tmp_path,diagnostic):
    from vdjmatch.cli.__main__ import main
    meta=dict(species='human',epitope='E',mhc_a='HLA-A*02',mhc_b='B2M',mhc_class='MHCI')
    sample=tmp_path/'q.tsv';ref=tmp_path/'r.tsv'
    pl.DataFrame([dict(junction_aa='CASSLGQAYEQYF',v_call='TRBV19*01',locus='TRB',**meta)]).write_csv(sample,separator='\t')
    pl.DataFrame([dict(cdr3=s,v='TRBV19*01',gene='TRB',**(meta|dict(species='HomoSapiens',epitope=e))) for s,e in [('CASSLGRAYEQYF','E'),('CASSWGRAYEQYF','OTHER')]]).write_csv(ref,separator='\t')
    args=['neighbours',str(sample),'--vdjdb',str(ref),'--locus','TRB','--targets-from-sample','--unbounded-nearest','--background','vdjdb-other','--neighbour-weighting','linear']
    for name,extra in [('old',[]),('new',['--'+diagnostic]),('parallel',['--'+diagnostic,'--threads','2']),('empty',['--'+diagnostic,'--species','mouse'])]:
        assert main(args+extra+['--output-prefix',str(tmp_path/name)])==0
    for suffix in ['queries','candidates','enrichment']:
        old=pl.read_csv(tmp_path/f'old.{suffix}.tsv',separator='\t');new=pl.read_csv(tmp_path/f'new.{suffix}.tsv',separator='\t')
        assert old.equals(new.select(old.columns))
        assert new.equals(pl.read_csv(tmp_path/f'parallel.{suffix}.tsv',separator='\t'))
    empty=pl.read_csv(tmp_path/'empty.enrichment.tsv',separator='\t')
    assert empty.height==0 and ('rule_mixed_imgt_wild4_score' if diagnostic=='epitope-rules' else 'chemistry_log_score') in empty.columns
