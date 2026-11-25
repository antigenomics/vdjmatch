"""Exercise raw declared-target CLI inputs, calibrated ranks and linked exact-hit semantics."""
import json
import math

import polars as pl
import pytest

from vdjmatch.cli.__main__ import main


def _files(tmp_path):
    aa = ['CAVAAAAAAGKLIF', 'CAVGGGGGAGKLIF', 'CAVPPPPPAGKLIF', 'CAVWWWWWAGKLIF']
    bb = ['CASSAAAAAYEQYF', 'CASSGGGGGYEQYF', 'CASSPPPPPYEQYF', 'CASSWWWWWYEQYF']
    reference = []
    for i, (a,b) in enumerate(zip(aa,bb)):
        for locus,v,j,seq in [('TRA','TRAV12-2','TRAJ42',a),('TRB','TRBV19','TRBJ2-7',b)]:
            reference.append(dict(complex_id=i+1, cdr3=seq, v=v, j=j, gene=locus,
                epitope='A' if i<2 else 'B', species='HomoSapiens',
                mhc_a='HLA-A*02:01' if i<3 else 'HLA-A*03:01', mhc_b='B2M', mhc_class='MHCI'))
    r=tmp_path/'reference.tsv'
    pl.DataFrame(reference).write_csv(r,separator='\t')
    rows=[]
    for k,(a,b,ep,restriction) in enumerate([(aa[0],bb[0],'A','HLA-A*02'),
                                               (aa[0],bb[1],'A','HLA-A*02'),
                                               (aa[0],bb[0],'B','HLA-A*03'),
                                               (aa[0],bb[0],'NO','HLA-A*02')]):
        for locus,v,j,seq in [('TRA','TRAV12-2','TRAJ42',a),('TRB','TRBV19','TRBJ2-7',b)]:
            rows.append(dict(query_id=f'{k}-{locus}',pair_id=f'p{k}',junction_aa=seq,v_call=v,j_call=j,
                             locus=locus,epitope=ep,mhc_a=restriction,mhc_b='B2M',mhc_class='MHCI'))
    sample=tmp_path/'sample.tsv'
    pl.DataFrame(rows).write_csv(sample,separator='\t')
    return sample,r


@pytest.mark.parametrize('locus', ['TRA','TRB','paired'])
def test_raw_equal_tail_query_and_reference_memberships(tmp_path,locus):
    sample,r=_files(tmp_path)
    def run(label,*options):
        prefix=tmp_path/label
        main(['equal-tail',str(sample),'--vdjdb',str(r),'--locus',locus,'--output-prefix',str(prefix),*options])
        assert json.loads(prefix.with_suffix('.manifest.json').read_text())['queries']==4
        return pl.read_csv(prefix.with_suffix('.scores.tsv'),separator='\t')
    keep=run('keep')
    exc=run('exc','--exclude-exact')
    compat=run('compatible','--background-scope','compatible-hla','--exclude-exact')
    assert keep['query_id'].to_list()==([f'p{i}' for i in range(4)] if locus=='paired'
                                        else [f'{i}-{locus}' for i in range(4)])
    assert keep['N'].to_list()==exc['N'].to_list()==[2,2,1,0]
    assert keep['M'].to_list()==exc['M'].to_list()==[2,2,2,4]
    assert compat['M'].to_list()==[1,1,0,3]
    assert keep['d'][0]==0 and exc['d'][0]>0
    if locus=='paired':
        assert keep['d'][1]>0
    assert compat['status'].to_list()==['available','available','empty_background','no_target_reference']
    assert compat['equal_tail_score'][2] is None and keep['equal_tail_score'][3] is None
    for g in (keep,exc,compat):
        for z in g.filter(pl.col('status')=='available').iter_rows(named=True):
            assert z['product_numerator']==z['D_tail_n']*z['K_tail_n']
            assert z['tail_denominator']==1+z['M']
            assert z['equal_tail_score']==pytest.approx(-math.log(z['product_numerator']/z['tail_denominator']**2))
    assert exc['d'].equals(compat['d']) and exc['k'].equals(compat['k'])


def test_equal_tail_rejects_missing_candidate_and_duplicate_paired_roles(tmp_path):
    s,r=_files(tmp_path);base=pl.read_csv(s,separator='\t')
    cmd=['equal-tail',str(s),'--vdjdb',str(r),'--locus','paired','--output-prefix',str(tmp_path/'out')]
    for frame in [base.drop('epitope'),pl.concat([base,base.head(1).with_columns(pl.lit('duplicate').alias('query_id'))])]:
        frame.write_csv(s,separator='\t')
        with pytest.raises(SystemExit):
            main(cmd)


def test_equal_tail_class_alias_and_unsupported_query(tmp_path):
    sample, ref = _files(tmp_path)
    base = pl.read_csv(sample, separator='\t').filter(pl.col('locus') == 'TRB').head(1)
    alias = base.with_columns(pl.lit('alias').alias('query_id'), pl.lit(' MHC-I ').alias('mhc_class'),
                              pl.lit(' A ').alias('epitope'))
    unsupported = base.with_columns(pl.lit('unsupported').alias('query_id'),
                                    pl.lit(None, dtype=pl.String).alias('junction_aa'))
    pl.concat([base, alias, unsupported]).write_csv(sample, separator='\t')
    prefix = tmp_path / 'checked'
    main(['equal-tail', str(sample), '--vdjdb', str(ref), '--locus', 'TRB', '--output-prefix', str(prefix)])
    scores = pl.read_csv(prefix.with_suffix('.scores.tsv'), separator='\t')
    assert scores['N'].to_list() == [2, 2, 2]
    assert scores['equal_tail_score'][0] == scores['equal_tail_score'][1]
    assert scores['calibration_available'].all()
    assert scores['has_evidence'].to_list() == [True, True, False]
    item = scores.row(2, named=True)
    assert item['status'] == 'no_evidence' and item['equal_tail_score'] == 0
    assert item['D_tail_n'] == item['K_tail_n'] == item['M'] + 1


@pytest.mark.parametrize('scope', ['all-hla', 'compatible-hla'])
def test_double_cdr12_ranks_queries_and_background_with_same_geometry(tmp_path, scope):
    import numpy as np
    from vdjmatch.match.geometry import distance_matrix

    sample, ref = _files(tmp_path)
    records = pl.read_csv(ref, separator='\t').filter(pl.col('gene') == 'TRB').with_columns(
        pl.Series('v', ['TRBV19*01', 'TRBV7-9*01', 'TRBV19*01', 'TRBV7-9*01']))
    records.write_csv(ref, separator='\t')
    q = pl.read_csv(sample, separator='\t').filter(pl.col('query_id') == '0-TRB').with_columns(
        pl.lit('TRBV7-9*01').alias('v_call'))
    q = pl.concat([q, q.with_columns(pl.lit('short').alias('query_id'),
                                      pl.lit('CASSF').alias('junction_aa')),
                   q.with_columns(pl.lit('unknown').alias('query_id'), pl.lit('unknown').alias('v_call'))])
    q.write_csv(sample, separator='\t')
    cmd = ['equal-tail', str(sample), '--vdjdb', str(ref), '--locus', 'TRB', '--exclude-exact',
           '--background-scope', scope, '--output-prefix', str(tmp_path/'rank')]
    main(cmd)
    baseline = pl.read_csv(tmp_path/'rank.scores.tsv', separator='\t')
    main(cmd + ['--distance-geometry', 'blosum62-double-cdr12'])
    new = pl.read_csv(tmp_path/'rank.scores.tsv', separator='\t')
    assert new.select('N', 'M', 'K_tail_n', 'k', 'kernel_missing').equals(
        baseline.select('N', 'M', 'K_tail_n', 'k', 'kernel_missing'))
    assert new['d'][0] != baseline['d'][0]
    assert new['distance_missing'].to_list() == [False, False, True]
    target = records.head(2)
    background = records.tail(2) if scope == 'all-hla' else records.slice(2, 1)
    def minimum(junction, v):
        hits = target.filter(pl.col('cdr3') != junction)
        return distance_matrix([junction], hits['cdr3'].to_list(), [v], hits['v'].to_list(),
                               matrix='blosum62', v_regions='double-cdr12', allow_short=True).min()
    others = np.array([minimum(seq, v) for seq, v in background.select('cdr3', 'v').iter_rows()])
    for i in range(2):
        distance = minimum(q['junction_aa'][i], q['v_call'][i])
        assert new['d'][i] == distance
        assert new['D_tail_n'][i] == 1 + (others <= distance).sum()
    assert new['D_tail_n'][2] == len(background) + 1
