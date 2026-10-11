"""Independent scalar oracle for the author-supplied VDJbet motif recipe."""
from collections import Counter
import math

import polars as pl
import pytest

from vdjmatch.match.motifs import vdjbet_motif_scores


def test_vdjbet_counts_masks_and_full_junction_puncture():
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
        out=vdjbet_motif_scores(q,r,exclude_exact=exclude)
        for row in out.iter_rows(named=True):
            values=[x for (e,t),x in expected.items() if e==row['epitope']]
            assert row['motif_score']==pytest.approx(sum(values)/len(rows[0]['cdr3']))
            assert row['supported_motifs']==len(values)
        assert out['epitope'].to_list()==sorted({e for e,t in expected})
    # A target with five instances loses its support gate after exact puncture.
    assert vdjbet_motif_scores(q,r.head(5),exclude_exact=True).height==0
    with pytest.raises(ValueError,match='unique'):
        vdjbet_motif_scores(pl.concat([q,q]),r)
