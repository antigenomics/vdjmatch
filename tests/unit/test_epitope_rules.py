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
def test_rules_dense_oracle_and_puncture(paired):
    from vdjmatch.match.tcrdist import distance_matrix,paired_distance_matrix
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
        same=[row for row in keep if all(row[f'v{i}'].split('*')[0]==q[f'v{i}'][0].split('*')[0] for i in axes)]
        assert out['rule_'+side+'_population'][0]==len(keep)
        assert out['rule_'+side+'_v_count'][0]==len(same)
        if paired:
            d=paired_distance_matrix(q['seq0'],q['seq1'],[row['seq0'] for row in same],[row['seq1'] for row in same],q['v0'],q['v1'],[row['v0'] for row in same],[row['v1'] for row in same],ctrim=3,v_loop_weight=0)[0]
        else:d=distance_matrix(q['seq0'],[row['seq0'] for row in same],q['v0'],[row['v0'] for row in same],ctrim=3,v_loop_weight=0)[0]
        mass=sum(max(0,1-x/90) for x in d)
        assert out['rule_'+side+'_spatial_mass'][0]==pytest.approx(mass)
        oracle[side]={'spatial':mass}
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
        rm=oracle['reference'][kind];bm=oracle['control'][kind];prior=(rm+bm+.5)/(n+m+1)
        expected=out['rule_v_score'][0]*(rm+prior)/(n+1)/((bm+prior)/(m+1))
        assert out['rule_'+kind+'_score'][0]==pytest.approx(expected)
    neutral=reference_epitope_rules(q,r,r,paired=paired)
    assert all(neutral[c][0]==pytest.approx(1) for c in neutral.columns if c.endswith('_score'))
    extra=q.with_columns(pl.lit(9).cast(pl.UInt32).alias('query_id'))
    assert out.equals(reference_epitope_rules(pl.concat([q,extra]),r,b,paired=paired,exclude_exact=True).head(1))


def test_motif_classes_wildcards_and_short_sequences():
    from vdjmatch.evalue.epitope_rules import _tokens
    q=pl.DataFrame({'_id':[0,1],'seq0':['CASAILVSTDEQYF','CASSQYF'],'v0':['TRBV19','TRBV19']})
    for mode in MOTIFS:
        t=_tokens(q,range(1),mode)
        assert set(t.filter(pl.col('_id')==0)['t0'])==tokens(q['seq0'][0],mode)
        assert t.filter(pl.col('_id')==1).height==0
    assert IMGT_CLASSES['L']==IMGT_CLASSES['I'] and IMGT_CLASSES['F']!=IMGT_CLASSES['Y']


def test_cli_rules_preserve_existing_fields_and_empty_schema(tmp_path):
    from vdjmatch.cli.__main__ import main
    meta=dict(species='human',epitope='E',mhc_a='HLA-A*02',mhc_b='B2M',mhc_class='MHCI')
    sample=tmp_path/'q.tsv';ref=tmp_path/'r.tsv'
    pl.DataFrame([dict(junction_aa='CASSLGQAYEQYF',v_call='TRBV19*01',locus='TRB',**meta)]).write_csv(sample,separator='\t')
    pl.DataFrame([dict(cdr3=s,v='TRBV19*01',gene='TRB',**(meta|dict(species='HomoSapiens',epitope=e))) for s,e in [('CASSLGRAYEQYF','E'),('CASSWGRAYEQYF','OTHER')]]).write_csv(ref,separator='\t')
    args=['tcrdist-neighbours',str(sample),'--vdjdb',str(ref),'--locus','TRB','--targets-from-sample','--unbounded-nearest','--background','vdjdb-other','--neighbour-weighting','linear']
    for name,extra in [('old',[]),('new',['--epitope-rules']),('parallel',['--epitope-rules','--threads','2']),('empty',['--epitope-rules','--species','mouse'])]:
        assert main(args+extra+['--output-prefix',str(tmp_path/name)])==0
    for suffix in ['queries','candidates','enrichment']:
        old=pl.read_csv(tmp_path/f'old.{suffix}.tsv',separator='\t');new=pl.read_csv(tmp_path/f'new.{suffix}.tsv',separator='\t')
        assert old.equals(new.select(old.columns))
        assert new.equals(pl.read_csv(tmp_path/f'parallel.{suffix}.tsv',separator='\t'))
    empty=pl.read_csv(tmp_path/'empty.enrichment.tsv',separator='\t')
    assert empty.height==0 and 'rule_mixed_imgt_wild4_score' in empty.columns
