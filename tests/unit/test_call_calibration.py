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


def test_reference_calls_retains_truth_and_withholds_all_restrictions(tmp_path):
    from vdjmatch.cli import reference_calls
    ref = tmp_path / 'ref.tsv'
    control = tmp_path / 'control.tsv'
    pl.DataFrame({'gene': ['TRB']*3, 'cdr3': ['CASSF', 'CASSF', 'CATSF'],
        'v': ['TRBV1']*3, 'j': ['TRBJ1']*3, 'species': ['HomoSapiens']*3,
        'epitope': ['ONE', 'TWO', 'ONE'], 'mhc_a': ['A', 'B', 'A']}).write_csv(ref, separator='\t')
    pl.DataFrame({'junction_aa': ['CASRF', 'CWWWF']}).write_csv(control, separator='\t')
    parser = argparse.ArgumentParser()
    reference_calls.register(parser.add_subparsers())
    for mode in ('labelled', 'withheld'):
        prefix = tmp_path / mode
        args = parser.parse_args(['reference-calls', '--vdjdb', str(ref), '--locus', 'TRB',
            '--min-target-records', '1',
            '--control', str(control), '--output-prefix', str(prefix),
            *(['--leave-epitope-out'] if mode == 'withheld' else [])])
        assert reference_calls.main(args) == 0
        calls = pl.read_csv(str(prefix)+'.calls.tsv', separator='\t')
        cases = pl.read_csv(str(prefix)+'.calibration-cases.tsv', separator='\t')
        assert cases['score'].to_list() == calls.select(
            pl.when(pl.col('vdjmatch_status') == 'matched').then(pl.col('vdjmatch_E'))
            .otherwise(None))['vdjmatch_E'].to_list()
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
    args = parser.parse_args(['reference-calls', '--vdjdb', str(ref), '--locus', 'TRB',
        '--control', str(control), '--leave-epitope-out', '--output-prefix', str(prefix)])
    assert reference_calls.main(args) == 0
    calls = pl.read_csv(str(prefix)+'.calls.tsv', separator='\t')
    assert calls.height == 1 and calls['epitope'].to_list() == ['ONE']
    assert calls['vdjmatch_epitope'].to_list() == ['TWO']  # complete searchable reference
    manifest = json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert manifest['eligible_benchmark_records'] == 30
    assert manifest['ineligible_benchmark_records'] == 29
