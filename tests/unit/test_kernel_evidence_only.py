"""Evidence-only exports retain detailed scores, ties and dispositions without dense candidates."""
import importlib.util
import json
from pathlib import Path

import polars as pl
import pytest


def fixture(tmp_path, paired):
    path = Path(__file__).with_name('test_grouped_kernel_details.py')
    spec = importlib.util.spec_from_file_location('evidence_fixture', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sample, reference = module.inputs(tmp_path, paired)
    raw = pl.read_csv(sample, separator='\t')
    invalid_target = raw.head(2 if paired else 1).with_columns(
        (pl.col('sequence_id')+'-target').alias('sequence_id'),
        pl.lit('invalid-target').alias('clone_id'), pl.lit(None, dtype=pl.String).alias('mhc_a'))
    pl.concat([raw, invalid_target]).write_csv(sample, separator='\t')
    refs = pl.read_csv(reference, separator='\t')
    twin = refs.filter(pl.col('epitope') == 'C').with_columns(
        pl.lit('TIE').alias('epitope'), (pl.col('complex_id').cast(pl.String)+'-tie').alias('complex_id'))
    refs = refs.with_columns(pl.col('complex_id').cast(pl.String))
    if paired:
        unlinked = refs.filter(pl.col('gene') == 'TRA').head(1).with_columns(
            pl.lit('D-unlinked').alias('epitope'), pl.lit('0').alias('complex_id'))
        refs = pl.concat([refs, unlinked])
    pl.concat([refs, twin]).write_csv(reference, separator='\t')
    return sample, reference


@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('assayed', [False, True])
@pytest.mark.parametrize('exclude', [False, True])
def test_evidence_only_exact_exports_without_distance_matrices(tmp_path, monkeypatch, paired, assayed, exclude):
    from vdjmatch.cli.__main__ import main
    from vdjmatch.match import geometry
    sample, reference = fixture(tmp_path, paired)
    args = ['evidence', str(sample), '--vdjdb', str(reference),
            '--locus', 'paired' if paired else 'TRB']
    args += ['--targets-from-sample'] if assayed else ['--all-targets']
    if exclude:
        args += ['--exclude-exact']
    full, compact = tmp_path/'full', tmp_path/'evidence'
    assert main(args + ['--output-prefix', str(full)]) == 0
    def forbidden(*args, **kwargs):
        raise AssertionError('evidence-only must not build candidate distance matrices')
    monkeypatch.setattr(geometry, 'distance_matrix', forbidden)
    monkeypatch.setattr(geometry, 'paired_distance_matrix', forbidden)
    assert main(args + ['--evidence-only', '--output-prefix', str(compact)]) == 0
    for lane in ['', '.TRA', '.TRB'] if paired else ['']:
        for kind in ['queries', 'enrichment', 'selection']:
            assert Path(str(full)+lane+'.'+kind+'.tsv').read_bytes() == Path(str(compact)+lane+'.'+kind+'.tsv').read_bytes()
        assert not Path(str(compact)+lane+'.candidates.tsv').exists()
        manifest = json.loads(Path(str(compact)+lane+'.manifest.json').read_text())
        assert manifest['evidence_only'] is True
        assert 'candidates' not in manifest['output_tables']
    if paired:
        assert Path(str(full)+'.paired-evidence.tsv').read_bytes() == Path(str(compact)+'.paired-evidence.tsv').read_bytes()
    status = pl.read_csv(str(compact)+'.queries.tsv', separator='\t')
    assert 'no_neighbours' in status['status']
    assert any('invalid' in s or 'unknown' in s for s in status['status'])
    if assayed:
        assert 'invalid_target' in status['status'] and 'no_reference' in status['status']
    if not assayed:
        scores = pl.read_csv(str(compact)+'.enrichment.tsv', separator='\t')
        for label in ['C', 'TIE']:
            assert scores.filter(pl.col('epitope') == label).height
        ties = pl.read_csv(str(compact)+'.selection.tsv', separator='\t')
        assert ties['tied_peptides'].max() > 1
        if paired:
            union = pl.read_csv(str(compact)+'.paired-evidence.tsv', separator='\t')
            marginal = union.filter(pl.col('epitope') == 'D-unlinked')
            assert marginal.height and marginal['kernel_linked_score'].null_count() == marginal.height
            assert marginal['kernel_alpha_score'].null_count() < marginal.height


@pytest.mark.parametrize('option', [['--selection-only'], ['--unbounded-nearest'], ['--shortlist-audit'],
                                  ['--background', 'real'], ['--position-weighting', 'ndn']])
def test_evidence_only_rejects_incompatible_diagnostics(tmp_path, option):
    from vdjmatch.cli.__main__ import main
    sample, reference = fixture(tmp_path, False)
    with pytest.raises((ValueError, SystemExit)):
        main(['evidence', str(sample), '--vdjdb', str(reference), '--locus', 'TRB',
              '--evidence-only', *option, '--output-prefix', str(tmp_path/'invalid')])


def test_evidence_only_rejects_null_assayed_peptide_instead_of_dense_fallback(tmp_path):
    from vdjmatch.cli.__main__ import main
    sample, reference = fixture(tmp_path, False)
    pl.read_csv(sample, separator='\t').with_columns(
        pl.lit(None, dtype=pl.String).alias('epitope')).write_csv(sample, separator='\t')
    with pytest.raises(SystemExit):
        main(['evidence', str(sample), '--vdjdb', str(reference), '--locus', 'TRB',
              '--targets-from-sample', '--evidence-only', '--output-prefix', str(tmp_path/'invalid')])


def test_gzip_evidence_preserves_plain_paired_exports(tmp_path):
    import gzip
    from vdjmatch.cli.__main__ import main
    sample, reference = fixture(tmp_path, True)
    args = ['evidence', str(sample), '--vdjdb', str(reference), '--locus', 'paired',
            '--all-targets', '--exclude-exact', '--evidence-only']
    plain, compressed = tmp_path/'plain', tmp_path/'gzip'
    assert main(args+['--output-prefix', str(plain)]) == 0
    assert main(args+['--gzip-evidence', '--output-prefix', str(compressed)]) == 0
    from vdjmatch.cli import neighbours
    from vdjmatch.runner.annotation import _core_defaults
    served = tmp_path/'served'
    options = _core_defaults(reference, served)
    options.sample, options.locus, options.reference_population = str(sample), 'paired', 'curated-benchmark'
    options.all_targets = options.exclude_exact = options.evidence_only = options.gzip_evidence = True
    marginals = {locus: {kind: pl.read_csv(str(plain)+'.'+locus+'.'+kind+'.tsv',
        separator='\t', schema_overrides={'query_id': pl.UInt32})
        for kind in ['queries', 'selection', 'enrichment']} for locus in ['TRA', 'TRB']}
    assert neighbours.main(options, marginal_tables=marginals) == 0
    for exported in [compressed, served]:
        for lane in ['', '.TRA', '.TRB']:
            for kind in ['queries', 'selection']:
                assert Path(str(plain)+lane+'.'+kind+'.tsv').read_bytes() == Path(str(exported)+lane+'.'+kind+'.tsv').read_bytes()
            for kind in ['enrichment'] + (['paired-evidence'] if not lane else []):
                expected = Path(str(plain)+lane+'.'+kind+'.tsv').read_bytes()
                actual = Path(str(exported)+lane+'.'+kind+'.tsv.gz')
                assert gzip.decompress(actual.read_bytes()) == expected
                assert not actual.with_suffix('').exists()
            if exported == compressed or not lane:
                manifest = json.loads(Path(str(exported)+lane+'.manifest.json').read_text())
                assert manifest['evidence_compression'] == dict(codec='gzip', level=1,
                    tables=['enrichment']+(['paired-evidence'] if not lane else []))
    ties = pl.read_csv(str(compressed)+'.selection.tsv', separator='\t')
    assert ties['tied_peptides'].max() > 1
    with pytest.raises(SystemExit):
        main(['evidence', str(sample), '--vdjdb', str(reference), '--locus', 'paired',
              '--gzip-evidence', '--output-prefix', str(tmp_path/'invalid')])
