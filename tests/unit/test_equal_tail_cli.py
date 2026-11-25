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


@pytest.mark.parametrize('locus', ['paired','TRB'])
def test_database_any_call_fit_apply_and_binding(tmp_path, locus):
    sample, ref = _files(tmp_path)
    raw = pl.read_csv(sample, separator='\t').drop('epitope','mhc_a','mhc_b','mhc_class')
    # Drop one linked partner, preserving its full marginal scoring population.
    raw = raw.filter(~((pl.col('pair_id') == 'p1') & (pl.col('locus') == 'TRA')))
    raw.write_csv(sample, separator='\t')
    prefix = tmp_path / 'cal'
    cmd = ['equal-tail', str(sample), '--vdjdb', str(ref), '--locus', locus, '--threads','1',
           '--all-targets','--exclude-exact','--distance-geometry','blosum62-double-cdr12']
    assert main(cmd + ['--fit-any-call','.25','--min-null-target-records','1','--output-prefix',str(prefix)]) == 0
    recipe = prefix.with_suffix('.calibration.json')
    from vdjmatch.db.cache import sha256
    from vdjmatch.cli import calibration
    from vdjmatch.evalue import directional
    from pathlib import Path
    pins = json.loads(recipe.read_text())['equal_tail_binding']['source_sha256']
    for owner in (calibration, directional):
        assert pins[owner.__name__] == sha256(Path(owner.__file__))
    cases = pl.read_csv(prefix.with_suffix('.null-cases.tsv'),separator='\t')
    assert cases.group_by('null_source').len().sort('null_source')['len'].to_list() == [4,4]
    run = cmd + ['--apply-any-call',str(recipe),'--output-prefix',str(tmp_path/'result')]
    assert main(run) == 0
    selected = pl.read_csv(tmp_path/'result.selection.tsv',separator='\t')
    assert len(selected) == 4 and selected['score'].null_count() == 0
    for row in selected.iter_rows(named=True):
        tails = [g['score'].ge(row['score']).sum()/g.height for g in cases.partition_by('null_source')]
        assert row['E_null'] == max(tails)
    with pytest.raises(SystemExit):
        main(run + ['--background-scope','compatible-hla'])
    with pytest.raises(SystemExit):
        main(run + ['--distance-geometry','native'])
    prefix.with_suffix('.null-cases.tsv').write_text('changed')
    with pytest.raises(SystemExit):
        main(run)


def test_complete_search_ties_no_hits_and_withholding(tmp_path):
    from fractions import Fraction
    from vdjmatch import db
    from vdjmatch.cli.equal_tail import _queries
    from vdjmatch.evalue.equal_tail import equal_tail_scores
    from vdjmatch.evalue.equal_tail_selection import equal_tail_selection

    sample, ref = _files(tmp_path)
    # B and C have the same target model; C remains eligible when B is withheld.
    r = db.load(ref, species='HomoSapiens').filter('reference_valid')
    r = pl.concat([r, r.filter(pl.col('epitope') == 'B').with_columns(pl.lit('C').alias('epitope'))])
    q, _ = _queries(sample, 'TRB', 'human', declared=False)
    q = q.head(1).with_columns(pl.lit(None, dtype=pl.String).alias('cdr3b'))
    q = pl.concat([q, q.with_columns(pl.lit('seen').alias('query_id'), pl.lit('CASSWWWWWYEQYF').alias('cdr3b'))])
    opts = dict(locus='TRB', species='human', threads=1, distance_geometry='blosum62-double-cdr12')
    selected = equal_tail_selection(q, r, **opts)
    assert selected['selection_status'].to_list() == ['no_support','ambiguous']
    assert selected['score'][0] is None and selected['selected_epitope'][1] is None
    held = equal_tail_selection(q.with_columns(pl.lit('B').alias('held')), r, withheld_column='held', **opts)
    assert held['selection_status'].to_list() == ['no_support','unique']
    assert held['selected_epitope'][1] == 'C'
    truth = []
    for task in r.select('epitope','mhc_a','mhc_b','mhc_class').unique().iter_rows(named=True):
        declared = q.with_columns(*[pl.lit(v).alias(k) for k,v in task.items()])
        z = equal_tail_scores(declared,r,**opts).row(1,named=True)
        if z['calibration_available'] and z['has_evidence']:
            truth.append(Fraction(z['product_numerator'], z['tail_denominator']**2))
    assert min(truth) == Fraction(selected['product_numerator'][1], selected['tail_denominator'][1]**2)
    assert selected['score'][1] == held['score'][1]



def test_paired_all_target_search_can_use_unlinked_beta_alone(tmp_path):
    from vdjmatch import db
    from vdjmatch.cli.equal_tail import _queries
    from vdjmatch.evalue.equal_tail_selection import equal_tail_selection

    s,r = _files(tmp_path)
    ref = db.load(r,species='HomoSapiens').filter(pl.col('reference_valid') & (pl.col('gene')=='TRB'))
    q,_ = _queries(s,'paired','human',declared=False)
    options=dict(species='human',threads=1,distance_geometry='blosum62-double-cdr12')
    pairs=equal_tail_selection(q,ref,locus='paired',**options)
    beta=equal_tail_selection(q,ref,locus='TRB',**options)
    assert pairs.equals(beta) and pairs['score'].null_count()==0


def test_kernel_score_slab_budget_preserves_declared_and_null_ranks(tmp_path):
    sample, ref = _files(tmp_path)
    # Repeat A as C to retain exact multi-target ties. Include a no-hit query.
    raw = pl.read_csv(sample, separator='\t')
    raw = pl.concat([raw,raw.tail(2).with_columns(pl.lit('nohit').alias('pair_id'),
                         pl.concat_str(pl.lit('x'), 'query_id').alias('query_id'),
                         pl.lit(None,dtype=pl.String).alias('junction_aa'))])
    raw.write_csv(sample,separator='\t')
    records=pl.read_csv(ref,separator='\t')
    extra=records.filter(pl.col('epitope')=='A').with_columns(pl.lit('C').alias('epitope'),
                                                               (pl.col('complex_id')+4).alias('complex_id'))
    pl.concat([records,extra]).write_csv(ref,separator='\t')
    cmd=['equal-tail',str(sample),'--vdjdb',str(ref),'--locus','paired',
         '--distance-geometry','blosum62-double-cdr12','--exclude-exact']
    for name, budget in [('split',1),('single',2**30)]:
        args=cmd+['--kernel-score-budget-bytes',str(budget)]
        def run(suffix, *options):
            prefix=tmp_path/(name+suffix)
            main(args+list(options)+['--output-prefix',str(prefix)])
            return prefix
        run('declared')
        run('fit','--all-targets','--fit-any-call','.25','--min-null-target-records','1')
        run('apply','--all-targets','--apply-any-call',str(tmp_path/(name+'fit.calibration.json')))
    for suffix in ['declared.scores.tsv','fit.null-cases.tsv','fit.operating_points.tsv','apply.selection.tsv','apply.null-tails.tsv']:
        left=pl.read_csv(tmp_path/('split'+suffix),separator='\t')
        right=pl.read_csv(tmp_path/('single'+suffix),separator='\t')
        assert left.equals(right)
    for pref in ['split','single']:
        result=pl.read_csv(tmp_path/(pref+'apply.selection.tsv'),separator='\t')
        assert result.filter(pl.col('query_id')=='nohit')['E_null'].item()==0
        assert 'ambiguous' in result['selection_status'].to_list()
    r0=json.loads((tmp_path/'splitfit.calibration.json').read_text())
    r1=json.loads((tmp_path/'singlefit.calibration.json').read_text())
    assert r0==r1
    with pytest.raises(SystemExit):
        main(cmd+['--kernel-score-budget-bytes','0','--output-prefix',str(tmp_path/'invalid')])


def test_previous_fixed_budget_binding_requires_the_verified_sources_and_event():
    from vdjmatch.cli.equal_tail_compatibility import _BUDGET_REVISION, compatible_binding
    from vdjmatch.db.cache import sha256
    from vdjmatch.cli import equal_tail as cli, equal_tail_calibration
    from vdjmatch.evalue import equal_tail, equal_tail_selection
    from pathlib import Path
    files={m.__name__:Path(m.__file__) for m in [cli,equal_tail,equal_tail_selection]}
    files['calibration']=Path(equal_tail_calibration.__file__)
    old={name:ids[0] for name,ids in _BUDGET_REVISION.items()}
    live={name:sha256(path) for name,path in files.items()}
    frozen=dict(source_sha256=old|{'geometry':'same'},locus='TRB')
    current=dict(source_sha256=live|{'geometry':'same'},locus='TRB')
    assert compatible_binding(frozen,current)
    assert not compatible_binding(frozen,current|{'locus':'TRA'})
    assert not compatible_binding(frozen,current|{'source_sha256':live|{'geometry':'different'}})
    assert not compatible_binding(frozen,current|{'source_sha256':current['source_sha256']|{'calibration':'unknown'}})


@pytest.mark.parametrize('exclude_exact', [False, True])
def test_nearest_batch_handles_shared_empty_and_unsupported_targets(exclude_exact):
    import numpy as np
    from vdjmatch.evalue.equal_tail import _nearest, _NO_DISTANCE
    from vdjmatch.match.geometry import paired_distance_matrix, precompute_v_loop_distances

    ctx = precompute_v_loop_distances(matrix='blosum62', v_regions='double-cdr12')
    q = pl.DataFrame(dict(cdr3a=['CAVAF', 'CAVGF', 'CAVAF', None],
                          cdr3b=['CASSF', 'CASSF', 'CASRF', 'CASSF'],
                          nativea=['TRAV12-2*01'] * 4, nativeb=['TRBV19*01'] * 4,
                          _junction=[0, 1, 2, 3]))
    targets = [q[:2], q[1:3], q[3:], q[:0]]
    batch = _nearest(q, targets, 'ab', ctx, 'human', 1, exclude_exact)
    expected = np.full((4, len(q)), _NO_DISTANCE, dtype=np.int32)
    for i, refs in enumerate(targets[:2]):
        d = paired_distance_matrix(*[f['cdr3' + c] for f in [q[:3], refs] for c in 'ab'],
                                   *[f['native' + c] for f in [q[:3], refs] for c in 'ab'],
                                   matrix='blosum62', v_loop_distances=ctx, allow_short=True)
        if exclude_exact:
            d[q['_junction'][:3].to_numpy()[:, None] == refs['_junction'].to_numpy()] = _NO_DISTANCE
        expected[i, :3] = d.min(axis=1)
    np.testing.assert_array_equal(batch, expected)
