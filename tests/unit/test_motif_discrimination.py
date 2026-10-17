"""Epitope document frequency is conditional and ignores record multiplicity."""
import math

import polars as pl
import pytest

from vdjmatch.aggregate.motif_discrimination import epitope_motif_documents


@pytest.mark.parametrize('mode,token', [('exact3', 'AAA'), ('wild4', 'AXAA')])
def test_epitope_idf_uses_eligible_epitopes_within_each_stratum(mode, token):
    source = pl.DataFrame({'epitope': ['A', 'B', 'C', 'tail'],
        'cdr3': ['CAAAAAAF', 'CAAAAAAF', 'CGGGGGGF', 'CAAAAAAF'], 'v': ['TRBV19']*4})
    source = pl.concat([source.filter(pl.col('epitope') != 'tail')]*30+[source.filter(pl.col('epitope') == 'tail')])
    docs = epitope_motif_documents(source, locus='TRB', mode=mode)
    for scope in ['all_lengths', 'exact_length']:
        row = docs.filter((pl.col('length_scope') == scope) & (pl.col('motif') == token)).row(0, named=True)
        assert row['epitope_documents'] == 3 and row['epitope_document_frequency'] == 2
        assert row['epitope_idf'] == pytest.approx(1+math.log(4/3))
    assert docs.equals(epitope_motif_documents(pl.concat([source, source.filter(pl.col('epitope') == 'A')]), locus='TRB', mode=mode))
    assert epitope_motif_documents(source.head(1), locus='TRB', mode=mode).is_empty()


def test_cli_full20_idfs_use_complete_corpus_despite_target_filter(tmp_path):
    import json
    from vdjmatch.cli.__main__ import main
    reference, control, output = tmp_path/'ref.tsv', tmp_path/'control.tsv', tmp_path/'out'
    source = pl.DataFrame({'epitope': ['A', 'B', 'C'],
        'cdr3': ['CAAAAAAF', 'CAAAAAAF', 'CGGGGGGF'], 'v': ['TRBV19']*3,
        'j': ['TRBJ1-1']*3, 'gene': ['TRB']*3, 'species': ['HomoSapiens']*3})
    pl.concat([source]*30).write_csv(reference, separator='\t')
    source.select(pl.col('cdr3').alias('junction_aa'), pl.col('v').alias('v_call'),
        pl.col('j').alias('j_call')).write_csv(control, separator='\t')
    assert main(['reference-structure', '--vdjdb', str(reference), '--control', str(control),
        '--locus', 'TRB', '--motif-mode', 'full20', '--epitope-idf', '--epitope', 'A',
        '--output-prefix', str(output)]) == 0
    motifs = pl.read_csv(str(output)+'.motifs.tsv', separator='\t')
    assert set(motifs['motif_mode']) == {'exact3', 'wild4'}
    assert set(motifs['motif_alphabet']) == {'full20'}
    row = motifs.filter((pl.col('motif') == 'AAA') & (pl.col('length_scope') == 'exact_length')).row(0, named=True)
    assert row['epitope_documents'] == 3 and row['epitope_document_frequency'] == 2
    assert row['epitope_idf'] == pytest.approx(1+math.log(4/3))
    assert row['conditional_epitope_tfidf'] == row['epitope_idf']
    assert row['control_idf'] == row['idf']
    assert row['conditional_control_idf_bits'] == 1.
    manifest = json.loads((tmp_path/'out.manifest.json').read_text())
    assert manifest['motif_information']['epitope_idf'] is not None
    assert manifest['motif_information']['k'] is None
