"""Recorded mock evidence remains unchanged by publication rendering."""
import pytest
import polars as pl


def test_report_preserves_recorded_p_q_and_all_mocks(tmp_path):
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    from vdjmatch.cli.repertoire_report import render
    counts=pl.DataFrame({'sample_id':['S']*201,'group_id':['T']*201,'replicate':range(-1,200),
                         'matched_curve_mean':[10.]+[1.]*200})
    summary=pl.DataFrame({'sample_id':['S'],'group_id':['T'],'original_p_curve_mean':[1/201],
                          'original_q_curve_mean':[.03]})
    result=render(counts,summary,'curve_mean',tmp_path/'report').row(0,named=True)
    assert result['recorded_p']==1/201 and result['recorded_q']==.03
    assert result['mock_panels']==200 and result['real_mock_ratio']==10
    assert (tmp_path/'report.pdf').exists()
    with pytest.raises(ValueError,match='200 distinct mocks'):
        render(counts.head(200),summary,'curve_mean',tmp_path/'invalid')
