"""Streamed source rows and controls retain eager semantics across boundaries."""
import gzip
import polars as pl
import pytest
from vdjmatch import io


@pytest.mark.parametrize('header,rows', [
    ('junction_aa\tv_call\tj_call\tduplicate_count\tpair_id',
     ['CASSF\tTRBV1\tTRBJ1-1\t2\t001', 'CATSF\tTRBV2\tTRBJ1-1\t3\t002', 'CASSF\tTRBV1\tTRBJ1-1\t2\t001']),
    ('count\tfreq\tcdr3nt\tcdr3aa\tv\td\tj',
     ['2\t.5\tTGT\tCASSF\tTRBV1\tTRBD1\tTRBJ1-1', '2\t.5\tTGT\tCATSF\tTRBV2\tTRBD1\tTRBJ1-1']),
    ('Read.count\tRead.proportion\tCDR3.nucleotide.sequence\tCDR3.amino.acid.sequence\tbestVGene\tbestJGene',
     ['2\t.5\tTGT\tCASSF\tTRBV1\tTRBJ1-1', '2\t.5\tTGT\tCATSF\tTRBV2\tTRBJ1-1']),
])
def test_stream_matches_eager(tmp_path, header, rows):
    path = tmp_path/'sample.tsv.gz'
    with gzip.open(path, 'wt') as stream:
        stream.write(header+'\n'+'\n'.join(rows)+'\n')
    expected = io.read_rearrangement(path)
    for size in (1, 2, 50):
        result = pl.concat(io.iter_rearrangement(path, batch_size=size))
        assert result.equals(expected)


def test_occurrence_control_projection_preserves_keys_and_weights(tmp_path):
    from vdjmatch.cli.repertoire_reference import _read_control
    from vdjmatch.evalue.mock_reference import control_events
    raw = pl.DataFrame(dict(cdr3aa=['CASSF','CASSF','CATSF'], v=['TRBV1']*3,
        j=['TRBJ1-1']*3, cdr3nt=['TGT',None,'TGA'], occurrences=['2','3','0'],
        sample_id=['a','b','a'], incidence=['1','1','1'], convergence=['1','1','1']))
    path=tmp_path/'control.tsv'
    raw.write_csv(path,separator='\t')
    for batch_size in (1,2,65536):
        assert _read_control(path,species='human',locus='TRB',batch_size=batch_size).equals(control_events(raw,species='human',locus='TRB'))
    raw=raw.with_columns(pl.Series('v',[None,'TRBV2','TRBV2']))
    raw.write_csv(path,separator='\t')
    assert _read_control(path,species='human',locus='TRB',batch_size=1).equals(control_events(raw,species='human',locus='TRB'))


def test_stream_preserves_canonical_junction_over_legacy(tmp_path):
    path=tmp_path/'mixed.tsv'
    path.write_text('junction_aa\tcdr3aa\tcount\tv\tj\nCASSF\tCWWWF\t2\tTRBV1\tTRBJ1-1\n')
    assert pl.concat(io.iter_rearrangement(path,batch_size=1)).equals(io.read_rearrangement(path))


def test_stream_control_skips_other_populations(tmp_path):
    from vdjmatch.cli.repertoire_reference import _read_control
    from vdjmatch.evalue.mock_reference import control_events
    raw=pl.DataFrame(dict(junction_aa=['CASSF','CATSF'],v=['TRBV1','TRBV2'],
        j=['TRBJ1-1']*2,species=['mouse','human'],occurrences=['1','2']))
    path=tmp_path/'mixed-control.tsv'
    raw.write_csv(path,separator='\t')
    assert _read_control(path,species='human',locus='TRB',batch_size=1).equals(control_events(raw,species='human',locus='TRB'))
