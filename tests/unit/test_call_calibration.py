"""Empirical any-call thresholds: ties, separate null sources and no-hit denominators."""
import argparse
import json
from pathlib import Path

import polars as pl
import pytest

from vdjmatch.cli import calibration


def test_cli_separate_nulls_ties_and_no_hits(tmp_path):
    rows = pl.DataFrame({
        'query_id': ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'],
        'score': [10., 8., None, 8., 8., None, 7., None],
        'outcome': ['correct', 'incorrect', 'incorrect', 'null', 'null', 'null', 'null', 'null'],
        'null_source': [None, None, None, 'withheld', 'withheld', 'withheld', 'generated', 'generated'],
    })
    path = tmp_path / 'cases.tsv'
    rows.write_csv(path, separator='\t')
    parser = argparse.ArgumentParser()
    calibration.register(parser.add_subparsers())
    args = parser.parse_args(['calibrate-calls', str(path), '--target-fpr', '.34',
                              '--output-prefix', str(tmp_path / 'out')])
    assert calibration.main(args) == 0
    result = json.loads((tmp_path / 'out.calibration.json').read_text())
    assert result['cutoff'] == 10.
    assert result['vdjdb_queries'] == 3
    assert result['correct_calls'] == 1
    assert result['vdjdb_correct_call_recall'] == pytest.approx(1/3)
    assert result['vdjdb_call_precision'] == 1.
    assert result['empirical_E_false_calls_per_query'] == 0.
    assert result['nulls'] == [
        {'source': 'generated', 'queries': 2, 'calls': 0, 'call_rate': 0.},
        {'source': 'withheld', 'queries': 3, 'calls': 0, 'call_rate': 0.},
    ]
    curve = pl.read_csv(tmp_path / 'out.operating_points.tsv', separator='\t')
    assert curve.filter(pl.col('cutoff') == 8.)['eligible'][0] is False
    assert curve.filter(pl.col('cutoff') == 8.)['empirical_E_false_calls_per_query'][0] == pytest.approx(2/3)
    # Inverting score direction preserves every operating point and selected call.
    lower = rows.with_columns((-pl.col('score')).alias('score'))
    same, _ = calibration.operating_points(lower, .34, higher_is_better=False)
    assert same['cutoff'] == -10.
    assert same['nulls'] == result['nulls']
    # A null tied with the highest correct hit forbids splitting that tie.
    tied = rows.with_columns(pl.when(pl.col('query_id') == 'd').then(10.)
                             .otherwise(pl.col('score')).alias('score'))
    abstain, _ = calibration.operating_points(tied, 0.)
    assert abstain['cutoff'] is None and abstain['action'] == 'abstain_all'
    with pytest.raises(ValueError, match='unique'):
        calibration.operating_points(pl.concat([rows, rows.head(1)]), .01)
    with pytest.raises(ValueError, match='finite'):
        calibration.operating_points(rows.with_columns(pl.lit(float('nan')).alias('score')), .01)
    with pytest.raises(ValueError, match='null_source'):
        calibration.operating_points(rows.drop('null_source'), .01)


@pytest.mark.parametrize('directory',[False,True])
def test_reference_calls_retains_truth_and_withholds_all_restrictions(tmp_path,directory):
    from vdjmatch.cli import reference_calls
    ref = tmp_path / 'ref.tsv'
    if directory:
        release=tmp_path/'release';release.mkdir();ref=release/'vdjdb.txt'
    control = tmp_path / 'control.tsv'
    pl.DataFrame({'gene': ['TRB']*3, 'cdr3': ['CASSF', 'CASSF', 'CATSF'],
        'v': ['TRBV1']*3, 'j': ['TRBJ1']*3, 'species': ['HomoSapiens']*3,
        'epitope': ['ONE', 'TWO', 'ONE'], 'mhc_a': ['A', 'B', 'A']}).write_csv(ref, separator='\t')
    pl.DataFrame({'junction_aa': ['CASRF', 'CWWWF']}).write_csv(control, separator='\t')
    parser = argparse.ArgumentParser()
    reference_calls.register(parser.add_subparsers())
    for mode in ('labelled', 'withheld'):
        prefix = tmp_path / mode
        args = parser.parse_args(['reference-calls', '--vdjdb', str(ref.parent if directory else ref), '--locus', 'TRB',
            '--min-target-records', '1',
            '--control', str(control), '--output-prefix', str(prefix),
            *(['--leave-epitope-out'] if mode == 'withheld' else [])])
        assert reference_calls.main(args) == 0
        calls = pl.read_csv(str(prefix)+'.calls.tsv', separator='\t')
        cases = pl.read_csv(str(prefix)+'.calibration-cases.tsv', separator='\t')
        assert cases['score'].to_list() == calls.select(
            pl.when(pl.col('vdjmatch_status') == 'matched').then(pl.col('vdjmatch_score'))
            .otherwise(None))['vdjmatch_score'].to_list()
        assert calls.height == (3 if mode == 'withheld' else 2)
        truth = pl.read_csv(str(prefix)+'.truth.tsv', separator='\t')
        assert truth.height == 3
        if mode == 'withheld':
            assert not calls.filter(pl.col('epitope') == pl.col('vdjmatch_epitope')).height
        else:
            assert truth.filter(pl.col('query_id') == 0)['epitope'].n_unique() == 2
    # Eligibility counts source records before self-exclusion, not distinct junctions.
    pl.DataFrame({'gene': ['TRB']*59, 'cdr3': ['CASSF']*30+['CATSF']*29,
        'v': ['TRBV1']*59, 'j': ['TRBJ1']*59, 'species': ['HomoSapiens']*59,
        'epitope': ['ONE']*30+['TWO']*29, 'reference_id': [f'study{i}' for i in range(59)]}).write_csv(ref, separator='\t')
    prefix = tmp_path / 'min30'
    args = parser.parse_args(['reference-calls', '--vdjdb', str(ref.parent if directory else ref), '--locus', 'TRB',
        '--control', str(control), '--leave-epitope-out', '--output-prefix', str(prefix)])
    assert reference_calls.main(args) == 0
    calls = pl.read_csv(str(prefix)+'.calls.tsv', separator='\t')
    assert calls.height == 1 and calls['epitope'].to_list() == ['ONE']
    assert calls['vdjmatch_epitope'].to_list() == ['TWO']  # complete searchable reference
    manifest = json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert manifest['eligible_benchmark_records'] == 30
    assert manifest['ineligible_benchmark_records'] == 29


def test_apply_frozen_calibration_rejects_without_using_query_outcomes(tmp_path):
    frozen = {'action': 'apply_cutoff', 'cutoff': 8., 'comparison': '>='}
    queries = pl.DataFrame({'query_id': ['a', 'b', 'c', 'd'],
        'score': [9., 8., 7., None], 'outcome': ['incorrect'] * 4})
    out = calibration.apply_calibration(queries, frozen)
    assert out['call_status'].to_list() == ['accepted', 'accepted', 'no_supported_match', 'no_supported_match']
    assert out['call_score'].to_list() == [9., 8., None, None]
    assert out.select(queries.columns).equals(queries)
    assert out.select('call_status', 'call_score').equals(
        calibration.apply_calibration(queries.drop('outcome'), frozen).select('call_status', 'call_score'))
    lower = queries.with_columns((-pl.col('score')).alias('score'))
    assert calibration.apply_calibration(lower, frozen | {'cutoff': -8., 'comparison': '<='})['call_status'].equals(out['call_status'])
    assert calibration.apply_calibration(queries, frozen | {'action': 'abstain_all', 'cutoff': None})['call_score'].null_count() == 4
    assert calibration.apply_calibration(queries.head(0), frozen)['call_status'].dtype == pl.String
    with pytest.raises(ValueError, match='finite'):
        calibration.apply_calibration(queries.with_columns(pl.lit(float('nan')).alias('score')), frozen)
    with pytest.raises(ValueError, match='unique'):
        calibration.apply_calibration(pl.concat([queries, queries.head(1)]), frozen)
    with pytest.raises(ValueError, match='cutoff'):
        calibration.apply_calibration(queries, frozen | {'cutoff': None})
    with pytest.raises(ValueError, match='comparison'):
        calibration.apply_calibration(queries, frozen | {'comparison': '>'})
    with pytest.raises(ValueError, match='output'):
        calibration.apply_calibration(out, frozen)
    path = tmp_path / 'queries.tsv'
    queries.write_csv(path, separator='\t')
    recipe = tmp_path / 'frozen.json'
    recipe.write_text(json.dumps(frozen))
    parser = argparse.ArgumentParser()
    calibration.register(parser.add_subparsers())
    target = tmp_path / 'applied.tsv'
    args = parser.parse_args(['apply-call-calibration', str(path), '--calibration', str(recipe), '--output', str(target)])
    assert calibration.apply_main(args) == 0
    assert pl.read_csv(target, separator='\t')['call_status'].to_list() == out['call_status'].to_list()
    manifest = json.loads(Path(str(target) + '.manifest.json').read_text())
    assert manifest['queries'] == 4 and manifest['accepted_calls'] == 2
    assert manifest['rejected_queries'] == 2 and manifest['no_score_queries'] == 1
    args.output = str(tmp_path / 'frozen')
    args.calibration = args.output + '.manifest.json'
    Path(args.calibration).write_text(json.dumps(frozen))
    with pytest.raises(ValueError, match='overwrite'):
        calibration.apply_main(args)


def test_frozen_null_tails_keep_ties_nohits_and_sources_separate():
    queries = pl.DataFrame({'query_id':['x','y','z'], 'score':[8.,11.,None],
                            'outcome':['incorrect']*3})
    nulls = pl.DataFrame({'query_id':['a','b','c','d','e'],
        'score':[8.,8.,None,10.,None], 'outcome':['null']*5,
        'null_source':['withheld']*3+['generated']*2})
    scored, tails = calibration.null_tail_probabilities(queries, nulls)
    assert scored.select(queries.columns).equals(queries)
    assert scored['p_null'].to_list() == pytest.approx([.75,1/3,1.])
    assert scored['E_null'].to_list() == pytest.approx([2/3,0.,0.])
    withheld = tails.filter(pl.col('null_source')=='withheld')
    assert withheld['null_queries'].to_list() == [3]*3
    assert withheld['null_tail_count'].to_list() == [2,0,0]
    lower_queries = queries.with_columns((-pl.col('score')).alias('score'))
    lower_nulls = nulls.with_columns((-pl.col('score')).alias('score'))
    lower, _ = calibration.null_tail_probabilities(lower_queries,lower_nulls,higher_is_better=False)
    assert lower.select('p_null','E_null').equals(scored.select('p_null','E_null'))
    with pytest.raises(ValueError,match='null'):
        calibration.null_tail_probabilities(queries,nulls.with_columns(pl.lit('correct').alias('outcome')))
    with pytest.raises(ValueError,match='finite'):
        calibration.null_tail_probabilities(queries.with_columns(pl.lit(float('inf')).alias('score')),nulls)


@pytest.mark.parametrize('score_column',['score','likelihood_score'])
def test_cli_null_tails_are_bound_to_frozen_calibration_input(tmp_path,score_column):
    parser=argparse.ArgumentParser();calibration.register(parser.add_subparsers())
    cases=tmp_path/'calibration.tsv'
    pl.DataFrame({'query_id':['label','n1','n2','n3'],score_column:[9.,8.,8.,None],
        'outcome':['correct','null','null','null'],
        'null_source':[None,'withheld','withheld','withheld']}).write_csv(cases,separator='\t')
    prefix=tmp_path/'recipe'
    assert calibration.main(parser.parse_args(['calibrate-calls',str(cases),
        '--score-column',score_column,'--target-fpr','.7','--output-prefix',str(prefix)]))==0
    query=tmp_path/'query.tsv';out=tmp_path/'out.tsv'
    query_frame=pl.DataFrame({'query_id':['x','y'],score_column:[8.,None]})
    if score_column!='score':query_frame=query_frame.with_columns(pl.lit(100.).alias('score'))
    query_frame.write_csv(query,separator='\t')
    args=parser.parse_args(['apply-call-calibration',str(query),'--calibration',
        str(prefix)+'.calibration.json','--null-cases',str(cases),'--output',str(out)])
    assert calibration.apply_main(args)==0
    result=pl.read_csv(out,separator='\t')
    assert result.select(query_frame.columns).equals(query_frame)
    assert result['p_null'].to_list()==[.75,1.]
    assert result['E_null'].to_list()==pytest.approx([2/3,0.])
    assert pl.read_csv(str(out)+'.null-tails.tsv',separator='\t').height==2
    cases.write_text(cases.read_text()+'\n')
    with pytest.raises(ValueError,match='hash'):
        calibration.apply_main(args)


def test_frozen_recipe_selects_likelihood_column_and_preserves_historical_score():
    query=pl.DataFrame({'query_id':['a','b'],'score':[100.,100.],'likelihood_score':[9.,7.]})
    recipe={'action':'apply_cutoff','cutoff':8.,'comparison':'>=','score_column':'likelihood_score'}
    scored=calibration.apply_calibration(query,recipe)
    assert scored['call_score'].to_list()==[9.,None]
    assert scored.select(query.columns).equals(query)
    nulls=pl.DataFrame({'query_id':['n1','n2'],'score':[100.,100.],
        'likelihood_score':[8.,None],'outcome':['null']*2,'null_source':['withheld']*2})
    out,_=calibration.null_tail_probabilities(scored,nulls,score_column='likelihood_score')
    assert out['p_null'].to_list()==pytest.approx([1/3,2/3])
    with pytest.raises(ValueError,match='likelihood_score'):
        calibration.apply_calibration(query.drop('likelihood_score'),recipe)
