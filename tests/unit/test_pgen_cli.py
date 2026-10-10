"""Marginal full-junction CLI preserves raw rows and delegates to native Pgen."""
import json
from pathlib import Path
from importlib.resources import files

import numpy as np
import polars as pl
import pytest

from vdjmatch.cli.__main__ import main


@pytest.fixture(scope='module')
def model_path():
    return str(files('vdjtools.model').joinpath('_bundled/olga/TRA'))


def run(tmp_path, model_path, frame, *options):
    sample=tmp_path/'input.tsv'
    frame.write_csv(sample,separator='\t')
    prefix=tmp_path/'result'
    result=main(['pgen',str(sample),'--model-path',model_path,'--species','human',
                 '--locus','TRA','--output-prefix',str(prefix),*options])
    return result,pl.read_csv(str(prefix)+'.scores.tsv',separator='\t',infer_schema_length=0),json.loads(
        (tmp_path/'result.manifest.json').read_text())


def test_marginal_oracle_and_raw_identity(tmp_path,model_path,monkeypatch):
    from vdjtools.model import load_model
    from vdjtools.model import native
    seq='CAVRDSNYQLIW'
    expected=native.pgen_aa_batch(load_model(model_path,validate=True),[seq,seq],threads=1)
    original=native.pgen_aa_batch
    calls=[]
    def batch(model,seqs,**kwargs):
        calls.append((seqs,kwargs))
        return original(model,seqs,**kwargs)
    monkeypatch.setattr(native,'pgen_aa_batch',batch)
    frame=pl.DataFrame({'query_id':['001','002','003','004','005','006'],
        'sequence_id':['a','b','c','d','e','f'],'clone_id':['pair']*6,
        'junction_aa':[seq,seq,None,'CAXF',seq,seq],
        'v_call':[None]*6,'locus':['TRA']*5+['TRB'],
        'species':['human']*4+['mouse','human']})
    result,out,manifest=run(tmp_path,model_path,frame,'--threads','4')
    assert result==0 and out['query_id'].to_list()==frame['query_id'].to_list()
    assert out['sequence_id'].to_list()==frame['sequence_id'].to_list()
    assert out['clone_id'].to_list()==frame['clone_id'].to_list()
    assert out['status'].to_list()==['scored','scored','missing_junction','invalid_junction',
                                     'unselected_species','unselected_locus']
    assert np.array_equal(out['pgen'].head(2).cast(pl.Float64).to_numpy(),expected)
    assert out['pgen'].tail(4).null_count()==4
    assert calls==[([seq,seq],{'mismatches':0,'threads':4})]
    assert manifest['model']['collapse_alleles'] is False
    assert len(manifest['model']['sha256'])==64 and manifest['counts']['input_rows']==6


def test_missing_v_locus_and_zero_mass(tmp_path,model_path,monkeypatch):
    from vdjtools.model import native
    monkeypatch.setattr(native,'pgen_aa_batch',lambda *args,**kwargs:[0.0])
    _,out,_=run(tmp_path,model_path,pl.DataFrame({'junction_aa':['CAVRDSNYQLIW',None]}))
    assert out['query_id'].to_list()==['0','1']
    assert out['status'].to_list()==['zero_mass','missing_junction']
    assert out['log10_pgen'].null_count()==2


@pytest.mark.parametrize('kwargs,match',[
    ({'species':'mouse'},'model species'),({'locus':'TRB'},'model locus'),
    ({'threads':'0'},'positive threads')])
def test_model_and_resource_validation(tmp_path,model_path,kwargs,match,capsys):
    path=tmp_path/'input.tsv';pl.DataFrame({'junction_aa':['CAVRDSNYQLIW']}).write_csv(path,separator='\t')
    options={'species':'human','locus':'TRA','threads':'1',**kwargs}
    argv=['pgen',str(path),'--model-path',model_path,'--output-prefix',str(tmp_path/'result')]
    for key,value in options.items():argv.extend(['--'+key,value])
    with pytest.raises(SystemExit):main(argv)
    assert match in capsys.readouterr().err


@pytest.mark.parametrize('probability',[float('nan'),float('inf'),-0.1,1.1])
def test_invalid_backend_probability_fails(tmp_path,model_path,monkeypatch,probability,capsys):
    from vdjtools.model import native
    monkeypatch.setattr(native,'pgen_aa_batch',lambda *args,**kwargs:[probability])
    with pytest.raises(SystemExit):
        run(tmp_path,model_path,pl.DataFrame({'junction_aa':['CAVRDSNYQLIW']}))
    assert 'invalid probability' in capsys.readouterr().err


def test_bare_cdr3_not_fabricated(tmp_path,model_path,capsys):
    with pytest.raises(SystemExit):run(tmp_path,model_path,pl.DataFrame({'cdr3_aa':['AVRDSNYQLI']}))
    assert 'junction_aa' in capsys.readouterr().err


def test_duplicate_query_ids_rejected(tmp_path,model_path,capsys):
    with pytest.raises(SystemExit):
        run(tmp_path,model_path,pl.DataFrame({'query_id':['01','01'],'junction_aa':['CAVRDSNYQLIW']*2}))
    assert 'query_id' in capsys.readouterr().err


def test_empty_batch_keeps_schema_and_never_calls_native(tmp_path,model_path,monkeypatch):
    from vdjtools.model import native
    def forbidden(*args,**kwargs):
        raise AssertionError('empty batch must not invoke native')
    monkeypatch.setattr(native,'pgen_aa_batch',forbidden)
    _,out,manifest=run(tmp_path,model_path,pl.DataFrame({'junction_aa':pl.Series([],dtype=pl.String)}))
    assert out.height==0 and 'pgen' in out.columns
    assert manifest['counts']['native_rows']==0


def test_native_serial_parallel_and_normalization(tmp_path,model_path):
    frame=pl.DataFrame({'junction_aa':[' cavrdsnyqliw ']*64})
    _,serial,_=run(tmp_path,model_path,frame,'--threads','1')
    _,parallel,_=run(tmp_path,model_path,frame,'--threads','4')
    assert serial.equals(parallel)
    assert parallel['junction_aa'].to_list()==frame['junction_aa'].to_list()
    assert parallel['junction_aa_normalized'].unique().to_list()==['CAVRDSNYQLIW']


def test_explicit_locus_survives_missing_sequence(tmp_path,model_path):
    _,out,_=run(tmp_path,model_path,pl.DataFrame({'junction_aa':[None,None], 'locus':['TRB','BAD']}))
    assert out['status'].to_list()==['unselected_locus','unselected_locus']


def test_unknown_v_does_not_condition_marginal_pgen(tmp_path,model_path):
    _,out,_=run(tmp_path,model_path,pl.DataFrame({'junction_aa':['CAVRDSNYQLIW'],'v_call':['unknown']}))
    assert out['status'].to_list()==['scored']


def test_reserved_internal_column_rejected(tmp_path,model_path,capsys):
    with pytest.raises(SystemExit):
        run(tmp_path,model_path,pl.DataFrame({'junction_aa':['CAVRDSNYQLIW'],'_pgen_locus':['TRA']}))
    assert 'reserved' in capsys.readouterr().err


@pytest.mark.parametrize('rich',[False,True])
def test_vdjdb_zip_preserves_reference_observations(tmp_path,model_path,rich):
    import zipfile
    from vdjmatch import db
    seq='CAVRDSNYQLIW'
    if rich:
        # Reuse the existing tiny release producer, with two alpha observations.
        from test_releases import rich_tables
        directory=tmp_path/'tables';rich_tables(directory)
        chains=pl.read_csv(directory/'chains.tsv',separator='\t').with_columns(
            pl.when(pl.col('record_id')=='R2').then(pl.lit('TRA')).otherwise(pl.col('gene')).alias('gene'),
            pl.when(pl.col('gene')!='TRB').then(pl.lit(seq)).otherwise(pl.col('cdr3')).alias('cdr3'))
        # The second expression sees original TRB; set its newly selected alpha explicitly.
        chains=chains.with_columns(pl.when(pl.col('gene')=='TRA').then(pl.lit(seq)).otherwise(pl.col('cdr3')).alias('cdr3'))
        chains.write_csv(directory/'chains.tsv',separator='\t')
        pl.DataFrame({'record_id':['R1'],'evidence_id':['E1'],'evidence_type':['independent_study']}).write_csv(
            directory/'evidence.tsv',separator='\t')
        archive=tmp_path/'rich.zip'
        with zipfile.ZipFile(archive,'w') as zf:
            for path in directory.iterdir():zf.write(path,path.name)
    else:
        archive=tmp_path/'legacy.zip'
        with zipfile.ZipFile(archive,'w') as zf:
            zf.writestr('release/vdjdb.txt','gene\tcdr3\tv.segm\tj.segm\tspecies\tantigen.epitope\treference.id\tmethod\n'+
                f'TRA\t{seq}\tTRAV1\tTRAJ1\tHomoSapiens\tPEP\tPMID:1\tfirst\n'+
                f'TRA\t{seq}\tTRAV2\tTRAJ2\tHomoSapiens\tPEP\tPMID:2\tsecond\n'+
                'TRB\tCASSF\tTRBV1\tTRBJ1\tHomoSapiens\tPEP\tPMID:3\tother\n')
    expected=db.load(archive,species='HomoSapiens',gene='TRA')
    prefix=tmp_path/'reference-pgen'
    assert main(['pgen',str(archive),'--vdjdb','--model-path',model_path,'--species','human',
                 '--locus','TRA','--output-prefix',str(prefix)])==0
    out=pl.read_csv(str(prefix)+'.scores.tsv',separator='\t',infer_schema_length=0)
    assert out.height==expected.height==2 and out['query_id'].to_list()==['0','1']
    assert out['junction_aa'].to_list()==[seq,seq] and out['pgen'].n_unique()==1
    assert out['reference_id'].to_list()==expected['reference_id'].to_list()
    assert out['v_call'].to_list()==expected['v'].to_list()
    for column in ('record_id','clonotype_id','method_verification') if rich else ('method',):
        assert out[column].to_list()==expected[column].cast(pl.String).to_list()
    if rich:
        assert json.loads(out['evidence'][0])['evidence'][0]['evidence_id']=='E1'
    manifest=json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert manifest['input_format']=='vdjdb' and manifest['reference']['input_sha256']==db.provenance(archive)['input_sha256']
    assert manifest['counts']['input_rows']==2
