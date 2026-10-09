"""Joint puncture and finite-control arithmetic, without cohort calculations."""
from itertools import product
import math

import polars as pl
import pytest
from scipy.stats import beta, binom, multinomial, poisson
from seqtree import Index

from vdjmatch import Annotator
from vdjmatch.evalue.paired import _joint_result, PAIRED_BOUND_DELTA
from vdjmatch.match import search_params


def result(nt, na, nb, n, ma, mb, za, zb):
    return _joint_result(nt, na, nb, n, ma, mb, exclude_exact=True,
                         n_exact_alpha=za, n_exact_beta=zb)


def reference():
    pairs = [("self", "CAVVF", "CASSF"), ("aexact", "CAVVF", "CASRF"),
             ("bexact", "CAVAF", "CASSF"), ("bothnear", "CAVAF", "CASRF"),
             ("replicate", "CAVAF", "CASRF")]
    return pl.DataFrame([{"complex_id": cid, "gene": g, "cdr3": s, "epitope": "PEP",
                          "record_id": cid, "reference_id": cid}
                         for cid,a,b in pairs for g,s in [("TRA",a),("TRB",b)]])


def test_cartesian_cardinality_swap_and_tiny_controls():
    for ma,mb in product(range(1,5), repeat=2):
        for na,nb in product(range(ma+1), range(mb+1)):
            for za,zb in product((0,1), repeat=2):
                if za>na or zb>nb: continue
                ca=['exact']*za+['near']*(na-za)+['out']*(ma-na)
                cb=['exact']*zb+['near']*(nb-zb)+['out']*(mb-nb)
                k=sum(a!='out' and b!='out' and (a,b)!=('exact','exact') for a,b in product(ca,cb))
                r=result(1,na,nb,3,ma,mb,za,zb)
                s=result(1,nb,na,3,mb,ma,zb,za)
                assert r['n_control_joint']==k==na*nb-za*zb
                assert r['E_raw']==pytest.approx(3*k/(ma*mb))
                for name in ['E','E_raw','p_upper','p_enrichment','p_poisson_raw']:
                    assert r[name]==pytest.approx(s[name])
                    assert math.isfinite(r[name])
                assert 0<=r['p_upper']<=1 and 0<=r['p_enrichment']<=1


def test_cp_bound_all_counts_and_raw_tail_oracle():
    for na,nb,za,zb in [(0,0,0,0),(1,1,1,1),(3,4,1,0),(3,4,0,1),(20,30,1,1)]:
        ma,mb,n,nt=20,30,11,3
        eps=PAIRED_BOUND_DELTA/4
        def u(k,m): return 1. if k==m else beta.ppf(1-eps,k+1,m-k)
        p=min(1,u(za,ma)*u(nb-zb,mb)+u(na-za,ma)*u(zb,mb)+u(na-za,ma)*u(nb-zb,mb))
        r=result(nt,na,nb,n,ma,mb,za,zb)
        assert r['p_upper']==pytest.approx(p)
        assert r['E']==pytest.approx(n*p)
        assert r['p_enrichment']==pytest.approx(min(1,PAIRED_BOUND_DELTA+binom.sf(nt-1,n,p)))
        assert r['p_poisson_raw']==pytest.approx(poisson.sf(nt-1,r['E_raw']))
        assert r['E']>=r['E_raw']
    r=result(0,1,1,3,1,1,1,1)
    assert r['n_control_joint']==0 and r['E_raw']==0 and r['E']>0
    assert r['p_enrichment']==1 and r['p_poisson_raw']==1


def test_model_conditional_tail_rejection_oracle():
    # Exhaustively integrate IID multinomial controls and independent binomial targets.
    # A small software oracle, not a biological calibration claim or scientific sweep.
    grid=[(0.,.01,.99),(.1,.02,.88),(.1,.2,.7),(.5,.5,0.)]
    ma,mb,n,alpha=3,4,20,.001
    for pa,pb in product(grid,repeat=2):
        true=pa[0]*pb[1]+pa[1]*pb[0]+pa[1]*pb[1]
        rejection=0.; failure=0.
        for xa in range(ma+1):
            for ya in range(ma-xa+1):
                wa=multinomial.pmf([xa,ya,ma-xa-ya],ma,pa)
                for xb in range(mb+1):
                    for yb in range(mb-xb+1):
                        wb=multinomial.pmf([xb,yb,mb-xb-yb],mb,pb)
                        # IID category draws can repeat exact receptors; compute the
                        # CP-category oracle directly, unlike unique-index za in {0,1}.
                        eps=PAIRED_BOUND_DELTA/4
                        def u(k,m): return 1. if k==m else beta.ppf(1-eps,k+1,m-k)
                        upper=min(1,u(xa,ma)*u(yb,mb)+u(ya,ma)*u(xb,mb)+u(ya,ma)*u(yb,mb))
                        if upper<true: failure+=wa*wb
                        tails=PAIRED_BOUND_DELTA+binom.sf(range(-1,n),n,upper)
                        rejection+=wa*wb*sum(binom.pmf(k,n,true) for k in range(n+1) if tails[k]<alpha)
        assert failure<=PAIRED_BOUND_DELTA+1e-12
        assert rejection<=alpha+1e-12


def test_exact_identity_joint_filter_and_observation_retention():
    ann=Annotator.from_frame(reference())
    q=pl.DataFrame({'query_id':['q'], 'cdr3_alpha_aa':['CAVVF'], 'cdr3_beta_aa':['CASSF']})
    ctrl={'TRA':Index.build(['CAVVF','CAVAF'],'aa'),'TRB':Index.build(['CASSF','CASRF'],'aa')}
    # Joint identities must retain one-exact / other-near reference pairs.
    sp=search_params('1',matrix='',gap_open=1,gap_extend=1)
    detail,c=ann.paired_candidates(q,scope=sp,control=ctrl,exclude_exact=True,return_hits=True)
    assert set(detail['complex_id'])=={'aexact','bexact','bothnear','replicate'}
    r=c.row(0,named=True)
    assert r['n_hits']==3 and r['n_records']==4 and r['n_reference']==4
    assert r['n_control_joint']==3 and r['n_exact_control_alpha']==r['n_exact_control_beta']==1
    assert r['E_raw']==3
    assert r['calibration']=='paired_independent_joint_punctured_binomial_bound'
    assert r['finite_control_delta']==1e-6 and not r['rule_of_three']
    assert c.equals(ann.paired_candidates(q,scope=sp,control=ctrl,exclude_exact=True))
    _,empty=ann.paired_candidates(q,scope='0',control=ctrl,exclude_exact=True,return_hits=True)
    assert empty.height==0 and empty.schema==c.schema
    out=ann.annotate_paired(q,scope='0',control=ctrl,exclude_exact=True)
    assert out['vdjmatch_status'].to_list()==['no_joint_hit']


@pytest.mark.parametrize('args', [(1,2,1,3,1,2,1,1),(1,1,1,3,2,2,2,1),
                                 (1,0,1,3,2,2,1,0),(1,1,1,3,0,2,1,1),
                                 (1,1,1,3,2,2,1.0,1),(4,1,1,3,2,2,1,1)])
def test_invalid_counts(args):
    with pytest.raises(ValueError): result(*args)


def test_default_controls_do_not_extract_exact_memberships():
    class Control:
        def __init__(self, sequences): self.index=Index.build(sequences,'aa')
        def __len__(self): return len(self.index)
        def search_batch(self,*args): return self.index.search_batch(*args)
        def ref_seq(self,*args): raise AssertionError('default calibration changed its work')
    ann=Annotator.from_frame(reference())
    q=pl.DataFrame({'cdr3_alpha_aa':['CAVVF'],'cdr3_beta_aa':['CASSF']})
    ctrl={'TRA':Control(['CAVVF','CAVAF']),'TRB':Control(['CASSF','CASRF'])}
    c=ann.paired_candidates(q,scope='1',control=ctrl)
    assert c['calibration'].to_list()==['paired_independent_fixed_ball']
    assert 'E_raw' not in c.columns


def test_seqtrie_joint_puncture_does_not_use_zero_edit_metadata():
    ann=Annotator.from_frame(reference())
    q=pl.DataFrame({'cdr3_alpha_aa':['CAVVF'],'cdr3_beta_aa':['CASSF']})
    ctrl={'TRA':Index.build(['CAVVF','CAVAF'],'aa'),'TRB':Index.build(['CASSF','CASRF'],'aa')}
    h,c=ann.paired_candidates(q,scope=search_params('1',engine='seqtrie'),control=ctrl,
                              exclude_exact=True,return_hits=True)
    assert set(h['complex_id'])=={'aexact','bexact','bothnear','replicate'}
    assert c['n_control_joint'].to_list()==[3]


def test_published_cli_joint_puncture_and_manifest(tmp_path):
    import json
    import subprocess
    import sys
    reference().write_csv(tmp_path/'reference.tsv',separator='\t')
    pl.DataFrame({'junction_aa':['CAVVF','CASSF'],'locus':['TRA','TRB'],
                  'cell_id':['cell','cell']}).write_csv(tmp_path/'query.tsv',separator='\t')
    for locus,seqs in [('TRA',['CAVVF','CAVAF']),('TRB',['CASSF','CASRF'])]:
        pl.DataFrame({'junction_aa':seqs}).write_csv(tmp_path/f'{locus}.tsv',separator='\t')
    command=[sys.executable,'-m','vdjmatch.cli','match',str(tmp_path/'query.tsv'),
             '--vdjdb',str(tmp_path/'reference.tsv'),'--species','any','--control-species','human',
             '--paired','--link','cell_id','--exclude-exact','--fresh-control',
             '--control','TRA='+str(tmp_path/'TRA.tsv'),'--control','TRB='+str(tmp_path/'TRB.tsv'),
             '--scope','1','--threads','1','--no-align','--output-prefix',str(tmp_path/'out')]
    run=subprocess.run(command,capture_output=True,text=True,timeout=30)
    assert run.returncode==0,run.stderr+run.stdout
    candidate=pl.read_csv(tmp_path/'out.query.candidates.txt',separator='\t')
    assert candidate['n_hits'].to_list()==[3]
    assert candidate['n_control_joint'].to_list()==[3]
    manifest=json.loads((tmp_path/'out.query.manifest.json').read_text())
    assert manifest['paired_calibration_model']['delta']==1e-6
    assert not manifest['paired_calibration_model']['deduplicated_controls_establish_iid']
