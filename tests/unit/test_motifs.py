"""Independent scalar oracle for the author-supplied masked-motif motif recipe."""
from collections import Counter
import math

import polars as pl
import pytest

from vdjmatch.match.motifs import masked_motif_scores


@pytest.mark.parametrize('source', ['sample', 'reference', 'directory', 'symlink', 'hardlink', 'prefix_symlink'])
def test_motif_cli_rejects_output_input_overlap_before_reading(tmp_path, source):
    from argparse import Namespace
    from vdjmatch.cli.motifs import main
    sample=tmp_path/'sample.tsv'; reference=tmp_path/'reference.tsv'; prefix=tmp_path/'result'
    sample.write_text('raw sample'); reference.write_text('raw reference')
    if source=='sample':
        sample=tmp_path/'result.scores.tsv'; sample.write_text('raw sample')
    elif source=='reference':
        reference=tmp_path/'result.manifest.json'; reference.write_text('raw reference')
    elif source=='directory':
        reference=tmp_path/'reference'; reference.mkdir(); prefix=reference/'result'
    elif source=='symlink':
        (tmp_path/'result.scores.tsv').symlink_to(sample)
    elif source=='hardlink':
        import os
        os.link(sample, tmp_path/'result.scores.tsv')
    else:
        prefix.symlink_to(tmp_path/'elsewhere')
        sample=tmp_path/'result.scores.tsv'; sample.write_text('raw sample')
    before=sample.read_bytes()
    with pytest.raises(ValueError, match='output path conflicts'):
        main(Namespace(sample=str(sample),vdjdb=str(reference),locus='TRB',species='human',
                       exclude_exact=False,output_prefix=str(prefix)))
    assert sample.read_bytes()==before


def test_masked_motif_counts_masks_and_full_junction_puncture():
    groups=('AILV','RHK','CM','ST','DE','NQ','G','F','P','W','Y')
    classes={a:g[0] for g in groups for a in g}
    def tokens(s):
        s=''.join(classes[a] for a in s)
        return {s[i:i+p]+'X'+s[i+p+1:i+4] for i in range(len(s)-3) for p in range(4)}
    rows=[{'cdr3':'CASSAAAA'+a+'F','v':'TRBV1*01','epitope':'NLV'} for a in 'AILVDEN']
    rows += [{'cdr3':'CASSGGGG'+a+'F','v':'TRBV2*01','epitope':'TTD'} for a in 'AILVDEN']
    rows += [dict(rows[0],v='TRBV2*02',epitope='TTD'),dict(rows[0])]
    r=pl.DataFrame(rows).with_columns(pl.lit('A').alias('mhc_a'),pl.lit('B').alias('mhc_b'),pl.lit('I').alias('mhc_class'))
    q=pl.DataFrame({'query_id':['q'],'cdr3':[rows[0]['cdr3']],'v':['TRBV1*02']})
    for exclude in (False,True):
        pop={(x['cdr3'],x['v'].split('*')[0],x['epitope']) for x in rows
             if not exclude or x['cdr3']!=rows[0]['cdr3']}
        counts=Counter((len(s),v,e,t) for s,v,e in pop for t in tokens(s))
        va=Counter();antigen=Counter();alltokens=Counter()
        for (length,v,e,t),n in counts.items():
            va[v,e,t]+=n;antigen[e,t]+=n;alltokens[t]+=n
        expected={}
        for (length,v,e,t),n in counts.items():
            if v!='TRBV1' or t not in tokens(rows[0]['cdr3']):continue
            if n>=2 and va[v,e,t]>=3 and antigen[e,t]>=5:
                expected[e,t]=max(expected.get((e,t),0),va[v,e,t]/alltokens[t]*math.log(len(pop)/alltokens[t]))
        out=masked_motif_scores(q,r,exclude_exact=exclude)
        for row in out.iter_rows(named=True):
            values=[x for (e,t),x in expected.items() if e==row['epitope']]
            assert row['motif_score']==pytest.approx(sum(values)/len(rows[0]['cdr3']))
            assert row['supported_motifs']==len(values)
        assert out['epitope'].to_list()==sorted({e for e,t in expected})
        targeted=masked_motif_scores(q.with_columns(pl.lit('NLV').alias('epitope'),
            pl.lit('A').alias('mhc_a'),pl.lit('B').alias('mhc_b'),pl.lit('I').alias('mhc_class')),
            r,exclude_exact=exclude,target_reference=r.filter(pl.col('epitope')=='NLV'))
        assert targeted.sort('epitope').equals(out.filter(pl.col('epitope')=='NLV').sort('epitope'))
    # A target with five instances loses its support gate after exact puncture.
    assert masked_motif_scores(q,r.head(5),exclude_exact=True).height==0
    with pytest.raises(ValueError,match='unique'):
        masked_motif_scores(pl.concat([q,q]),r)


@pytest.mark.parametrize('directory',[False,True])
def test_motif_cli_compatible_union_and_repeated_queries(tmp_path,monkeypatch,directory):
    from argparse import Namespace
    from vdjmatch import db
    from vdjmatch.cli.motifs import main
    sequences=['CASSAAAA'+a+'F' for a in 'AILVDEN']
    r=pl.DataFrame({'cdr3':sequences,'v':['TRBV1']*7,'gene':['TRB']*7,
        'epitope':['NLV']*7,'mhc_a':['HLA-A*02:01']*3+['HLA-A*02:02']*4,
        'mhc_b':['B2M']*7,'mhc_class':['MHCI']*7,'reference_valid':[True]*7})
    monkeypatch.setattr(db,'load',lambda *a,**kw:r)
    sample=tmp_path/'queries.tsv';reference=tmp_path/'reference.tsv'
    if directory:
        reference.mkdir();(reference/'vdjdb.txt').write_text('fixture')
    else:reference.write_text('fixture')
    pl.DataFrame({'junction_aa':[sequences[0],sequences[0],None,'BADX'],'v_call':['TRBV1']*4,'j_call':['TRBJ1-1']*4,
        'locus':['TRB']*4,'species':['human']*4,'epitope':['NLV']*4,'mhc_a':['HLA-A*02']*4,
        'mhc_b':['B2M']*4,'mhc_class':['MHCI']*4,'query_id':['a','b','missing','invalid']}).write_csv(sample,separator='\t')
    prefix=tmp_path/'result'
    main(Namespace(sample=str(sample),vdjdb=str(reference),locus='TRB',species='human',exclude_exact=True,output_prefix=str(prefix)))
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t')
    assert set(out['query_id'])=={'a','b','missing','invalid'}
    valid=out.filter(pl.col('query_id').is_in(['a','b']))
    assert valid['status'].to_list()==['scored','scored']
    assert valid['supported_motifs'].min()>0
    assert valid['motif_score'][0]==valid['motif_score'][1]
    invalid=out.filter(pl.col('query_id').is_in(['missing','invalid']))
    assert invalid['status'].to_list()==['invalid_query','invalid_query']
    assert invalid['motif_score'].null_count()==2
