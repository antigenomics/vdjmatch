"""Fresh old/new detailed CLI equality on synthetic full pMHC/source identities."""
import importlib.util
from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal


def inputs(tmp_path,paired):
    path=Path(__file__).with_name('test_grouped_kernel_selection.py')
    spec=importlib.util.spec_from_file_location('detail_fixture',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    sample,reference=module.fixture(tmp_path,paired)
    raw=pl.read_csv(sample,separator='\t')
    first=raw.head(2 if paired else 1)
    orphan=first.with_columns((pl.lit('orphan-')+pl.col('sequence_id')).alias('sequence_id'),
        pl.lit('orphan').alias('clone_id'),
        pl.when(pl.col('locus')=='TRB').then(pl.lit('C'+'W'*90+'F')).otherwise(pl.col('junction_aa')).alias('junction_aa'))
    invalid=first.with_columns((pl.lit('invalid-')+pl.col('sequence_id')).alias('sequence_id'),
        pl.lit('invalid').alias('clone_id'),pl.lit('UNKNOWN').alias('v_call'))
    pl.concat([raw,orphan,invalid]).write_csv(sample,separator='\t')
    return sample,reference


@pytest.mark.parametrize('paired',[False,True])
@pytest.mark.parametrize('assayed',[False,True])
@pytest.mark.parametrize('exclude',[False,True])
def test_detailed_old_new_all_fields_and_native_batch_budget(tmp_path,monkeypatch,paired,assayed,exclude):
    from vdjmatch.cli import tcrdist as cli
    from vdjmatch.match import tcrdist as td
    from vdjmatch.cli.__main__ import main
    sample,reference=inputs(tmp_path,paired)
    args=['annotate','--details',str(sample),'--vdjdb',str(reference),'--locus','paired' if paired else 'TRB']
    if assayed:args+=['--targets-from-sample']
    if exclude:args+=['--exclude-exact']
    eligibility=cli._grouped_detailed_eligible
    native=td.grouped_total_distance_count_batch;calls=[]
    def counted(*args,**kwargs):
        calls.append(kwargs['threads'])
        return native(*args,**kwargs)
    monkeypatch.setattr(td,'grouped_total_distance_count_batch',counted)
    for threads in (1,4):
        old,new=tmp_path/('old'+str(threads)),tmp_path/('new'+str(threads))
        monkeypatch.setattr(cli,'_grouped_detailed_eligible',lambda *a:False)
        assert main(args+['--threads',str(threads),'--output-prefix',str(old)])==0
        calls.clear()
        monkeypatch.setattr(cli,'_grouped_detailed_eligible',eligibility)
        assert main(args+['--threads',str(threads),'--output-prefix',str(new)])==0
        assert calls==[threads]*(3 if paired else 1)
        old_files=sorted(tmp_path.glob(old.name+'*.tsv'))
        assert old_files
        for path in old_files:
            other=path.with_name(path.name.replace(old.name,new.name,1))
            expected=pl.read_csv(path,separator='\t',infer_schema_length=0)
            actual=pl.read_csv(other,separator='\t',infer_schema_length=0)
            assert_frame_equal(actual,expected,check_exact=True)
            if threads==4:
                serial=other.with_name(other.name.replace('new4','new1',1))
                assert_frame_equal(actual,pl.read_csv(serial,separator='\t',infer_schema_length=0),check_exact=True)
