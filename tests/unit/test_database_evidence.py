"""Whole-database presence and matched nearest-distance null tails, no label choice."""
import json
from pathlib import Path
import polars as pl
import pytest
from seqtree import Index
from vdjmatch.cli.__main__ import main


def test_presence_and_null_tail_use_one_batch_and_retain_every_draw(tmp_path,monkeypatch):
    reference=tmp_path/'reference.tsv'
    pl.DataFrame({'cdr3':['CASSF','CATSF','CASSF'],'gene':['TRB']*3,
        'species':['HomoSapiens']*3,'v':['TRBV1']*3,'j':['TRBJ1-1']*3,
        'epitope':['A','B','C']}).write_csv(reference,separator='\t')
    sample=tmp_path/'query.tsv'
    pl.DataFrame({'query_id':['exact','near','far','bad','missing','alpha'],
        'junction_aa':['CASSF','CARSF','WWWWW','CA*F',None,'CASSF'],
        'v_call':['TRBV1']*6,'j_call':['TRBJ1-1']*6,
        'locus':['TRB']*5+['TRA']}).write_csv(sample,separator='\t')
    null=tmp_path/'null.tsv'
    pl.DataFrame({'junction_aa':['CASSF','CARSF','CARSF','WWWWW'],'locus':['TRB']*4}).write_csv(null,separator='\t')
    calls=[];native=Index.search_batch
    def spy(self,queries,params,threads=0):
        calls.append((len(queries),threads));return native(self,queries,params,threads)
    monkeypatch.setattr(Index,'search_batch',spy)
    prefix=tmp_path/'result'
    args=['database-evidence',str(sample),'--vdjdb',str(reference),'--null-sample',str(null),
        '--locus','TRB','--max-edits','1','--threads','1','--output-prefix',str(prefix)]
    assert main(args)==0
    q=pl.read_csv(str(prefix)+'.database.tsv',separator='\t')
    assert calls==[(3,1)]
    assert q['query_id'].to_list()==['exact','near','far','bad','missing','alpha']
    assert q['database_exact_junction'].to_list()==[True,False,False,None,None,None]
    assert q['database_exact_gene_clonotype'].to_list()==[True,False,False,None,None,None]
    assert q['database_nearest_edits'].to_list()==[0,1,None,None,None,None]
    assert q['database_p'].to_list()==[2/5,4/5,1.,None,None,None]
    assert q['database_status'].to_list()==['matched','matched','no_hit','invalid_sequence','missing_sequence','locus_mismatch']
    draws=pl.read_csv(str(prefix)+'.null.tsv',separator='\t')
    assert draws.height==4 and draws['database_nearest_edits'].null_count()==1
    tails=pl.read_csv(str(prefix)+'.null-tails.tsv',separator='\t')
    assert tails['null_tail_count'].to_list()==[1,3,None]
    assert 'E_null' not in tails.columns
    manifest=json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert manifest['unique_reference_junctions']==2 and manifest['null_sources'][0]['draws']==4
    # Replacing all specificity labels changes neither membership nor the P-value.
    pl.read_csv(reference,separator='\t').with_columns(pl.lit('OTHER').alias('epitope')).write_csv(reference,separator='\t')
    assert main(args)==0
    assert pl.read_csv(str(prefix)+'.database.tsv',separator='\t')['database_p'].equals(q['database_p'])
    other=tmp_path/'other-null.tsv'
    pl.DataFrame({'junction_aa':['CASSF','CASSF']}).write_csv(other,separator='\t')
    assert main(args+['--null-sample',str(other)])==0
    multi=pl.read_csv(str(prefix)+'.database.tsv',separator='\t')
    assert multi['database_p'].head(3).to_list()==[1.,1.,1.]
    pl.read_csv(reference,separator='\t').with_columns(pl.when(pl.col('cdr3')=='CASSF')
        .then(pl.lit('TRBV1*01,TRBV2*01')).otherwise(pl.col('v')).alias('v')).write_csv(reference,separator='\t')
    assert main(args)==0
    assert pl.read_csv(str(prefix)+'.database.tsv',separator='\t')['database_exact_gene_clonotype'][0] is False
    pl.read_csv(sample,separator='\t').with_columns(pl.when(pl.col('query_id')=='exact')
        .then(pl.lit('TRBV1*01,TRBV2*01')).otherwise(pl.col('v_call')).alias('v_call')).write_csv(sample,separator='\t')
    assert main(args)==0
    assert pl.read_csv(str(prefix)+'.database.tsv',separator='\t')['database_exact_gene_clonotype'][0] is None
    pl.read_csv(sample,separator='\t').with_columns(pl.when(pl.col('query_id')=='exact')
        .then(pl.lit('TRBV1/TRBV2')).otherwise(pl.col('v_call')).alias('v_call')).write_csv(sample,separator='\t')
    assert main(args)==0
    assert pl.read_csv(str(prefix)+'.database.tsv',separator='\t')['database_exact_gene_clonotype'][0] is None
    # Invalid null draws must not silently reduce the denominator.
    pl.DataFrame({'junction_aa':['CA*F']}).write_csv(null,separator='\t')
    with pytest.raises(SystemExit): main(args)
    from vdjmatch.cli.database import prepare
    from argparse import Namespace
    pl.read_csv(sample,separator='\t').with_columns(pl.lit('source metadata').alias('_database_v')).write_csv(sample,separator='\t')
    with pytest.raises(ValueError,match='reserved _database_'):
        prepare(sample,Namespace(input_format='auto',sequence_convention=None,species='human',locus='TRB'))
    pl.DataFrame({'junction_aa':['CASSF']*2,'locus':['BETA']*2,
        'species':['human']*2,'organism':['HomoSapiens','MusMusculus']}).write_csv(sample,separator='\t')
    prepared,_=prepare(sample,Namespace(input_format='auto',sequence_convention=None,species='human',locus='TRB'))
    assert prepared['database_status'].to_list()==['ready','species_mismatch']
