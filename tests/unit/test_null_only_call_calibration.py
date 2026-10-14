"""Null-only fitting preserves every source denominator without labelled outcomes."""
from types import SimpleNamespace
import json
import polars as pl
import pytest
from vdjmatch.cli import calibration


def test_null_only_fit_and_cli_roundtrip_keep_recall_unavailable(tmp_path):
    frame=pl.DataFrame({'query_id':[str(i) for i in range(150)],
        'score':[.8]+[None]*149,'outcome':['null']*150,
        'null_source':['withheld']*100+['generated']*50})
    result,curve=calibration.operating_points(frame,.01)
    assert result['cutoff']==.8
    assert result['vdjdb_queries']==0
    assert result['correct_calls']==result['incorrect_calls']==0
    assert result['vdjdb_correct_call_recall'] is None
    assert result['vdjdb_call_precision'] is None
    assert result['nulls']==[{'source':'generated','queries':50,'calls':0,'call_rate':0.},
        {'source':'withheld','queries':100,'calls':1,'call_rate':.01}]
    assert curve.schema['vdjdb_correct_call_recall']==pl.Float64
    assert curve['vdjdb_correct_call_recall'].null_count()==curve.height
    cases=tmp_path/'cases.tsv';frame.write_csv(cases,separator='\t')
    calibration.main(SimpleNamespace(cases=str(cases),target_fpr=.01,lower_is_better=False,
        output_prefix=str(tmp_path/'fit'),score_column='score',input_manifest=None))
    fitted=json.loads((tmp_path/'fit.calibration.json').read_text())
    assert fitted['vdjdb_correct_call_recall'] is None
    applied=calibration.apply_calibration(frame,fitted)
    assert applied.height==150
    assert applied.filter(pl.col('call_status')=='accepted').height==1
    with pytest.raises(ValueError,match='null queries'):
        calibration.operating_points(frame.with_columns(pl.lit('incorrect').alias('outcome')),.01)
