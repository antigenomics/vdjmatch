"""Small user-owned labeled references reuse ordinary matching and pairing."""
import polars as pl
import pytest
from polars.testing import assert_frame_equal

from vdjmatch import Annotator


@pytest.mark.parametrize('kind', ['rows', 'frame', 'epitope', 'csv'])
def test_single_labels_keep_exact_hits_and_arbitrary_labels(tmp_path, kind):
    rows = [('CASSQYF', 'KLGGALQAK'), ('CATSRYF', 'control-label')]
    data = rows
    if kind != 'rows':
        data = pl.DataFrame(rows, schema=['junction_aa', 'label'], orient='row')
        if kind == 'epitope':
            data = data.rename({'label': 'epitope'})
        if kind == 'csv':
            path = tmp_path / 'reference.csv'
            data.write_csv(path)
            data = path
    ann = Annotator.from_labels(data, locus='TRB')
    assert ann.reference_report == dict(reference_rows=2, searchable_rows=2, unsearchable_rows=0)
    got = ann.annotate(['CASSQYF', 'CWWWWWWWF'], scope='0', threads=1,
                       sequence_convention='junction')
    assert got['vdjmatch_epitope'].to_list() == ['KLGGALQAK', None]
    assert got['vdjmatch_nearest_edits'].to_list()[0] == 0
    excluded = ann.annotate(['CASSQYF'], scope='0', threads=1,
                            sequence_convention='junction', exclude_exact=True)
    assert excluded['vdjmatch_epitope'].to_list() == [None]


def test_label_pairs_keep_real_links_and_conflicting_labels():
    rows = [('CAVQYF', 'CASSQYF', 'first'), ('CAARYF', 'CATSRYF', 'second')]
    ann = Annotator.from_labels(rows)
    assert set(ann.loci) == {'TRA', 'TRB'}
    reference = ann._index.reference
    assert reference['complex_id'].n_unique() == 2
    assert reference.group_by('complex_id').len()['len'].to_list() == [2, 2]
    queries = pl.DataFrame(rows, schema=['junction_aa_alpha', 'junction_aa_beta', 'label'], orient='row')
    calls = ann.annotate_paired(queries, scope='0', threads=1)
    assert calls['vdjmatch_epitope'].to_list() == ['first', 'second']
    expected = Annotator.from_frame(pl.DataFrame(rows,
        schema=['junction_aa_alpha', 'junction_aa_beta', 'epitope'], orient='row'))
    assert_frame_equal(reference.drop('label').sort('complex_id', 'gene'),
                       expected._index.reference.sort('complex_id', 'gene'), check_column_order=False)


def test_single_label_locus_is_explicit_and_ambiguous_headers_fail():
    with pytest.raises(ValueError, match='locus'):
        Annotator.from_labels([('CASSQYF', 'label')])
    with pytest.raises(ValueError, match='locus'):
        Annotator.from_labels([('CASSQYF', 'label')], locus='IGH')
    with pytest.raises(ValueError, match='label.*epitope'):
        Annotator.from_labels(pl.DataFrame({'junction_aa': ['CASSQYF'],
            'label': ['a'], 'epitope': ['b']}), locus='TRB')


def test_single_alpha_reference_supplies_locus_but_query_metadata_wins():
    ann = Annotator.from_labels([('CAVQYF', 'alpha')], locus='TRA')
    assert ann.annotate(['CAVQYF'], scope='0', threads=1)['vdjmatch_epitope'].to_list() == ['alpha']
    query = pl.DataFrame({'junction_aa': ['CAVQYF'], 'gene': ['TRB']})
    assert ann.annotate(query, scope='0', threads=1)['vdjmatch_epitope'].to_list() == [None]


@pytest.mark.parametrize('locus', ['TRA', 'TRB'])
def test_cli_label_reference_needs_no_vdjdb_metadata_or_preparation(tmp_path, monkeypatch, locus):
    import json
    from vdjmatch import db
    from vdjmatch.db import index
    from vdjmatch.cli.__main__ import main

    junction = 'CAVQYF' if locus == 'TRA' else 'CASSQYF'
    reference, query = tmp_path / 'labels.csv', tmp_path / 'queries.tsv'
    pl.DataFrame({'junction_aa': [junction], 'label': ['user-class']}).write_csv(reference)
    pl.DataFrame({'sequence_id': ['known', 'unknown'], 'junction_aa': [junction, 'CWWWWWWWF']}).write_csv(query, separator='\t')
    monkeypatch.setattr(db, 'fetch_latest', lambda **kw: pytest.fail('custom reference fetched VDJdb'))
    monkeypatch.setattr(index, 'load_or_build_index', lambda *a, **kw: pytest.fail('label query built diagnostics'))
    prefix = tmp_path / 'result'
    args = ['match', str(query), '--reference', str(reference), '--reference-format', 'labels',
            '--locus', locus, '--no-align', '-o', str(prefix)]
    assert main(args) == 0
    calls = pl.read_csv(tmp_path / 'result.queries.calls.txt', separator='\t')
    assert calls['vdjmatch_epitope'].to_list() == ['user-class', None]
    manifest = json.loads((tmp_path / 'result.queries.manifest.json').read_text())
    assert manifest['exclude_exact'] is False and manifest['reference_format'] == 'labels'
    assert manifest['parameters']['evalue'] is False and manifest['parameters']['locus'] == locus
    assert main(args + ['--exclude-exact', '-o', str(tmp_path / 'excluded')]) == 0
    excluded = pl.read_csv(tmp_path / 'excluded.queries.calls.txt', separator='\t')
    assert excluded['vdjmatch_epitope'].to_list() == [None, None]


def test_cli_paired_labels_use_the_same_wide_columns_for_reference_and_queries(tmp_path):
    from vdjmatch.cli.__main__ import main
    frame = pl.DataFrame({'junction_aa_alpha': ['CAVQYF', 'CAARYF'],
        'junction_aa_beta': ['CASSQYF', 'CATSRYF'], 'label': ['first', 'second']})
    reference, query = tmp_path / 'labels.tsv', tmp_path / 'queries.tsv'
    frame.write_csv(reference, separator='\t')
    frame.drop('label').write_csv(query, separator='\t')
    assert main(['match', str(query), '--reference', str(reference), '--reference-format', 'labels',
                 '--paired', '--no-align', '-o', str(tmp_path / 'paired')]) == 0
    calls = pl.read_csv(tmp_path / 'paired.queries.calls.txt', separator='\t')
    assert calls['vdjmatch_epitope'].to_list() == ['first', 'second']


@pytest.mark.parametrize('official', [False, True])
def test_ordinary_vdjdb_routes_retain_exact_hits(tmp_path, monkeypatch, official):
    from vdjmatch import db
    from vdjmatch.cli.__main__ import main
    reference, query = tmp_path / 'custom-vdjdb.csv', tmp_path / 'query.tsv'
    pl.DataFrame({'cdr3': ['CASSQYF'], 'gene': ['TRB'], 'species': ['HomoSapiens'],
        'antigen.epitope': ['CUSTOM']}).write_csv(reference)
    pl.DataFrame({'junction_aa': ['CASSQYF'], 'locus': ['TRB']}).write_csv(query, separator='\t')
    fetched = []
    def fetch(**kwargs):
        fetched.append(kwargs)
        return reference
    monkeypatch.setattr(db, 'fetch_latest', fetch)
    args = ['match', str(query), '--fresh-reference', '--no-evalue', '--no-align',
            '-o', str(tmp_path / 'out')]
    if not official:
        args += ['--reference', str(reference)]
    assert main(args) == 0
    assert len(fetched) == int(official)
    calls = pl.read_csv(tmp_path / 'out.query.calls.txt', separator='\t')
    assert calls['vdjmatch_epitope'].to_list() == ['CUSTOM']


def test_conflicting_exact_labels_stay_ambiguous_and_bare_cdr3_is_not_reinterpreted():
    ann = Annotator.from_labels([('CASSQYF', 'A'), ('CASSQYF', 'B')], locus='TRB')
    candidates = ann.candidates(['CASSQYF'], scope='0', threads=1)
    assert set(candidates['epitope']) == {'A', 'B'}
    assert candidates['n_exact'].to_list() == [1, 1]
    assert ann.annotate(['CASSQYF'], scope='0', threads=1)['vdjmatch_epitope'].to_list() == [None]
    with pytest.raises(ValueError, match='anchors'):
        Annotator.from_labels(pl.DataFrame({'cdr3_aa': ['ASSQY'], 'label': ['A']}), locus='TRB')


def test_label_reference_mhc_selector_retains_distinct_two_field_alleles():
    reference = pl.DataFrame({'junction_aa': ['CASSQYF', 'CASSQYF'],
        'label': ['first', 'second'], 'mhc_a': ['B*27:02', 'B*27:05']})
    ann = Annotator.from_labels(reference, locus='TRB', mhc_a='B*27:05')
    assert ann.reference_report['reference_rows'] == 1
    got = ann.annotate(['CASSQYF'], scope='0', threads=1)
    assert got['vdjmatch_epitope'].to_list() == ['second']
    assert got['vdjmatch_mhc_a'].to_list() == ['B*27:05']


def test_label_cli_rejects_missing_reference_and_conflicting_pair_locus_before_io(monkeypatch, capsys):
    from vdjmatch import db
    from vdjmatch.cli.__main__ import main
    monkeypatch.setattr(db, 'fetch_latest', lambda **kw: pytest.fail('invalid label request fetched VDJdb'))
    with pytest.raises(SystemExit):
        main(['match', 'missing.tsv', '--reference-format', 'labels'])
    assert 'requires --reference' in capsys.readouterr().err
    with pytest.raises(SystemExit):
        main(['match', 'missing.tsv', '--paired', '--locus', 'TRB'])
    assert 'omit it with --paired' in capsys.readouterr().err
