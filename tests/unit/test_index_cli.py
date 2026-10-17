"""Public reference preparation: custom, fetched and reusable native assets."""
import polars as pl


def test_custom_and_fetched_index_and_api_share_the_saved_bundle(tmp_path, monkeypatch, capsys):
    from vdjmatch import Annotator, db
    from vdjmatch.cli.__main__ import main
    from vdjmatch.match.engine import VdjdbIndex

    source = tmp_path / 'reference.tsv'
    pl.DataFrame({'junction_aa': ['CASSF', 'CASAF'], 'locus': ['TRB'] * 2,
        'epitope': ['DEMO'] * 2, 'v_call': ['TRBV19'] * 2}).write_csv(source, separator='\t')
    assert main(['index', '--vdjdb', str(source), '--threads', '1']) == 0
    bundle = next(source.with_name(source.name + '.vdjmatch').glob('*/manifest.json')).parent
    assert pl.read_csv(bundle / 'profiles.tsv', separator='\t').height == 2
    joint = pl.read_csv(bundle / 'joint-profiles.tsv', separator='\t')
    assert joint.height == 1 and joint['both_chains_observed'].to_list() == [False]
    assert 'profiles.tsv' in capsys.readouterr().out

    def forbidden(*args, **kwargs):
        raise AssertionError('reopening a saved reference must not rebuild it')
    monkeypatch.setattr(VdjdbIndex, 'build', forbidden)
    annotation = Annotator.from_path(source, threads=1)
    assert annotation.reference_profiles.height == 2
    assert annotation.reference_index_manifest['artifact_path'] == str(bundle)
    assert annotation.hits(['CASSF'], scope='0', threads=1).height == 1
    monkeypatch.setattr(db, 'fetch_latest', lambda **kwargs: source)
    assert main(['update', '--threads', '1']) == 0
    assert main(['index', '--threads', '1']) == 0
