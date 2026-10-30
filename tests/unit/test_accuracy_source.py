"""The manuscript sorted contrast preserves original gene calls and source labels."""
import importlib.util
import gzip
from pathlib import Path
import polars as pl
import pytest


def test_original_sorted_source_is_not_replaced_by_mirror(tmp_path):
    spec=importlib.util.spec_from_file_location('accuracy_example',Path(__file__).resolve().parents[2]/'examples/accuracy.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    airr=tmp_path/'vdjmatch/airr';airr.mkdir(parents=True)
    for number,labels in [(1,['NLVPMVATV']*2),(2,['LLWNGPMAV','LLLGIGILV'])]:
        frame=pl.DataFrame({'junction_aa':['CASSIRSSYEQYF','CASSIRSSYDQYF'],
            'v_call':['TRBV19']*2,'j_call':['TRBJ2-7']*2,'binder':[1,0] if number==1 else [1,1],
            'epitope':labels,'species':['human']*2,
            'mhc_a':['HLA-A*02']*2,'mhc_b':['B2M']*2,'mhc_class':['MHCI']*2})
        with gzip.open(airr/f'dataset_{number}_TRB.tsv.gz','wt') as stream:
            stream.write(frame.write_csv(separator='\t'))
    source=tmp_path/'sorted.tsv'
    pl.DataFrame({'cdr3':['CASSIRSSYEQYF','CASSIRSSYDQYF'],'v.segm':['TRBV20-1','TRBV7-9'],
        'j.segm':['TRBJ2-7']*2,'antigen.epitope':['LLWNGPMAV','LLLGIGILV']}).write_csv(source,separator='\t')
    raw=module.prepare(tmp_path,source)
    contrast=raw.filter(pl.col('input_origin')=='original-sorted-source')
    assert 'cdr3' not in raw.columns and 'source_cdr3' in raw.columns
    assert raw.height==6 and raw['query_id'].n_unique()==6
    assert raw['locus'].to_list()==['TRB']*6
    from vdjmatch.io import read_rearrangement
    assembled=tmp_path/'assembled.tsv';raw.write_csv(assembled,separator='\t')
    normalized=read_rearrangement(assembled,source='airr')
    assert normalized.height==6 and normalized['locus'].to_list()==['TRB']*6
    assert contrast['v_call'].to_list()==['TRBV20-1','TRBV7-9']
    assert contrast['source_epitope'].to_list()==['LLWNGPMAV','LLLGIGILV']
    assert contrast['epitope'].to_list()==['LLWNGPMAV']*2
    assert contrast['positive'].to_list()==[1,0] and contrast['binder'].to_list()==[1,1]
    assert contrast['task'].str.contains('original V calls').all()
    pl.read_csv(source,separator='\t').with_columns(pl.lit('OTHER').alias('antigen.epitope')).write_csv(source,separator='\t')
    with pytest.raises(ValueError,match='no rows are silently removed'): module.prepare(tmp_path,source)
