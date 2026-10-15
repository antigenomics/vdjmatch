"""All local donors/replicates and rare eligible receptors reach the CLI batch."""
import gzip
import importlib.util
from pathlib import Path

import polars as pl
import pytest


def test_complete_repertoires_and_private_metadata(tmp_path):
    spec=importlib.util.spec_from_file_location('vaccination_example',Path(__file__).resolve().parents[2]/'examples/vaccination.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    metadata=pl.DataFrame({'file_name':['first.airr.tsv.gz','second.airr.tsv.gz','absent.airr.tsv.gz'],
        'donor':['private-a','private-b','private-a'],'day':[0,0,7],'replica':['r1','r1','r1']})
    metadata.write_csv(tmp_path/'metadata.txt',separator='\t')
    raw=pl.DataFrame({'locus':['TRB','TRB','TRB','TRA','TRB','TRB','TRB'],
        'junction_aa':['CASSF','CASSF','CASAF','CAVF','CA*F','CASKF','QASSIRSSYEQYF'],
        'v_gene':['TRBV19','TRBV19','TRBV19','TRAV1','TRBV19',None,'TRBV19'],
        'j_gene':['TRBJ2-7','TRBJ2-7',None,'TRAJ1','TRBJ2-7','TRBJ2-7','TRBJ2-7'],
        'duplicate_count':[100,20,1,10,2,3,2]})
    for name in metadata['file_name'][:2]:
        with gzip.open(tmp_path/name,'wt') as stream:stream.write(raw.write_csv(separator='\t'))
    observations,queries,coverage,selected,donors=module.prepare(tmp_path)
    assert observations.height==6 and queries.height==3
    assert selected.height==2 and len(donors)==2
    assert queries['mhc_a'].to_list()==['HLA-A*02']*3
    assert queries['locus'].to_list()==['TRB']*3
    assert observations['duplicate_count'].sort().to_list()==[1,1,2,2,120,120]
    assert observations['query_id'].n_unique()==3 and observations['query_id'].null_count()==0
    assert queries.filter(pl.col('junction_aa')=='QASSIRSSYEQYF')['canonical_anchor_residues'].to_list()==[False]
    assert set(observations['donor'])=={'Donor 1','Donor 2'}
    assert coverage.select('raw_rows','eligible_rows','non_trb_rows','invalid_junction_rows','missing_v_rows').rows()==[(7,4,1,1,1)]*2
    assert coverage['noncanonical_anchor_rows'].to_list()==[1,1]
    assert coverage['eligible_read_counts'].to_list()==[123,123]
    metadata.with_columns(pl.lit('private-a').alias('donor'),
        pl.when(pl.col('file_name')=='first.airr.tsv.gz').then(10).otherwise(2).alias('day')).write_csv(
            tmp_path/'metadata.txt',separator='\t')
    assert module.prepare(tmp_path)[2]['day'].to_list()==[2,10]
    with gzip.open(tmp_path/'first.airr.tsv.gz','wt') as stream:
        stream.write(raw.with_columns(pl.lit(-1).alias('duplicate_count')).write_csv(separator='\t'))
    with pytest.raises(ValueError,match='non-negative integer counts'):module.prepare(tmp_path)
