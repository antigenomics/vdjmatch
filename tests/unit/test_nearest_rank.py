"""Uncensored nearest rank uses the same linked reference and V geometry."""
import numpy as np
import polars as pl
import pytest
from scipy.stats import hypergeom
from vdjmatch.cli.__main__ import main

@pytest.mark.parametrize('locus',['TRB','paired'])
def test_cli_nearest_rank_recovers_outside_ball_without_changing_existing_fields(tmp_path,locus):
    sample,reference=[tmp_path/n for n in ['q.tsv','r.tsv']]
    a='CASSLGQAYEQYF';b='CASSLGRAYEQYF';c='CASSWGRAYEQYF';d='CASAWGRAYEQYF'
    meta=dict(species='human',epitope='E',mhc_a='HLA-A*02',mhc_b='B2M',mhc_class='MHCI')
    loci=['TRA','TRB'] if locus=='paired' else ['TRB'];genes=['TRAV1-1*01','TRBV19*01'] if locus=='paired' else ['TRBV19*01']
    pl.DataFrame([dict(clone_id='q',junction_aa=a,v_call=v,locus=l,**meta) for l,v in zip(loci,genes)]).write_csv(sample,separator='\t')
    rows=[dict(complex_id=str(i),cdr3=s,v=v,gene=l,**(meta|{'species':'HomoSapiens','epitope':e})) for i,(s,e) in enumerate([(a,'E'),(b,'E'),(b,'OTHER'),(c,'OTHER'),(d,'OTHER')]) for l,v in zip(loci,genes)]
    pl.DataFrame(rows).write_csv(reference,separator='\t')
    args=['tcrdist-neighbours',str(sample),'--vdjdb',str(reference),'--locus',locus,'--targets-from-sample','--unbounded-nearest','--background','vdjdb-other','--neighbour-weighting','local-rank','--radius','0','--exclude-exact']
    for extra in [[],['--nearest-rank']]:
        prefix=tmp_path/('new' if extra else 'old');assert main(args+extra+['--output-prefix',str(prefix)])==0
    for suffix in ['queries','candidates']:
        assert pl.read_csv(tmp_path/f'old.{suffix}.tsv',separator='\t').equals(pl.read_csv(tmp_path/f'new.{suffix}.tsv',separator='\t'))
    old=pl.read_csv(tmp_path/'old.enrichment.tsv',separator='\t');new=pl.read_csv(tmp_path/'new.enrichment.tsv',separator='\t')
    assert old.equals(new.select(old.columns));assert old['n_reference'][0]==0 and old['p_nearest_rank'][0]==1
    distance=pl.read_csv(tmp_path/'new.candidates.tsv',separator='\t')['distance'][0]
    from vdjmatch.match.tcrdist import distance_matrix,paired_distance_matrix
    controls=[b,c,d]
    if locus=='paired':d=paired_distance_matrix([a],[a],controls,controls,[genes[0]],[genes[1]],[genes[0]]*len(controls),[genes[1]]*len(controls))[0]
    else:d=distance_matrix([a],controls,[genes[0]],[genes[0]]*len(controls))[0]
    m=int(np.count_nonzero(d<=distance));p=float(hypergeom.sf(0,4,1,m+1))
    assert new['nearest_rank_control_count'][0]==m
    assert new['p_nearest_unbounded'][0]==pytest.approx(p)
    assert new['nearest_rank_score'][0]==pytest.approx(-p)
    assert new['conditional_rank_bayes_factor'][0]==pytest.approx(4/(sum(1/i for i in range(1,5))*(m+1)))
    assert new['nearest_rank_status'][0]=='conditional_rank_model_only' and p<1
    assert main(args+['--nearest-rank','--threads','2','--output-prefix',str(tmp_path/'parallel')])==0
    assert new.equals(pl.read_csv(tmp_path/'parallel.enrichment.tsv',separator='\t'))
    # Opt-in residual is neutral on missing J and leaves every previous field exact.
    assert main(args+['--nearest-rank','--conditional-j','--output-prefix',str(tmp_path/'j')])==0
    residual=pl.read_csv(tmp_path/'j.enrichment.tsv',separator='\t')
    assert new.equals(residual.select(new.columns))
    assert residual['conditional_j_usage_ratio'][0]==1
    assert residual['j_conditioned_rank_score'][0]==new['conditional_rank_bayes_factor'][0]
    assert residual['conditional_j_status'][0]=='missing_query_j_neutral'
    # Puncturing the sole target key leaves explicit empty evidence, not a call.
    pl.DataFrame([r for r in rows if r['cdr3']!=b or r['epitope']!='E']).write_csv(reference,separator='\t')
    assert main(args+['--nearest-rank','--output-prefix',str(tmp_path/'empty')])==0
    empty=pl.read_csv(tmp_path/'empty.enrichment.tsv',separator='\t')
    assert empty['nearest_rank_score'][0] is None and empty['nearest_rank_status'][0]=='empty_population'
    assert empty['conditional_rank_bayes_factor'][0] is None


def test_conditional_j_empty_output_keeps_count_schema(tmp_path):
    sample=tmp_path/'sample.tsv';reference=tmp_path/'reference.tsv'
    pl.DataFrame({'junction_aa':['CASSLGQAYEQYF'],'v_call':['TRBV19*01'],'locus':['TRB'],'species':['human'],'epitope':['E'],'mhc_a':['HLA-A*02'],'mhc_b':['B2M'],'mhc_class':['MHCI']}).write_csv(sample,separator='\t')
    pl.DataFrame({'cdr3':['CASSLGRAYEQYF'],'v':['TRBV19*01'],'gene':['TRB'],'species':['HomoSapiens'],'epitope':['E'],'mhc_a':['HLA-A*02'],'mhc_b':['B2M'],'mhc_class':['MHCI']}).write_csv(reference,separator='\t')
    prefix=tmp_path/'empty'
    assert main(['tcrdist-neighbours',str(sample),'--vdjdb',str(reference),'--species','mouse','--locus','TRB','--targets-from-sample','--unbounded-nearest','--nearest-rank','--background','vdjdb-other','--neighbour-weighting','local-rank','--conditional-j','--output-prefix',str(prefix)])==0
    out=pl.read_csv(str(prefix)+'.enrichment.tsv',separator='\t')
    assert out.height==0
    assert {'j_reference_population','j_control_stratum','j_reference_probability','conditional_j_status','j_conditioned_rank_score'}<=set(out.columns)
