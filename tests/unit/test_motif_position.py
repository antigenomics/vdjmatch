"""Hand-enumerated V/J/length-conditioned junction contrasts."""
import polars as pl
import pytest


def test_position_contrast_pools_alleles_deduplicates_and_matches_both_genes_length():
    from vdjmatch.aggregate.motif_position import position_contrast
    reference = pl.DataFrame({'cdr3': ['CAGF', 'CAGF', 'CATF', 'CCCC'],
        'v': ['TRBV1*01', 'TRBV1*02', 'TRBV1', 'TRBV2'], 'j': ['TRBJ1-1']*4})
    control = pl.DataFrame({'junction_aa': ['CAGF', 'CGGF', 'CTTF', 'CAAF', 'CAAAF'],
        'v': ['TRBV1']*5, 'j': ['TRBJ1-1', 'TRBJ1-1', 'TRBJ1-1', 'TRBJ2-1', 'TRBJ1-1'],
        'event_weight': [1000, 1, 1, 100, 100]})
    positions, motifs = position_contrast(reference, control, locus='TRB',
        v_gene='trbv1*03', j_gene='TRBJ1-1', junction_length=4, motif='AG')
    assert positions.height == 80
    row = positions.filter((pl.col('position') == 2) & (pl.col('amino_acid') == 'A')).row(0, named=True)
    assert row['reference_keys'] == 2 and row['control_keys'] == 3
    assert row['reference_count'] == 2 and row['control_count'] == 1
    assert row['reference_fraction'] == 1 and row['control_fraction'] == pytest.approx(1/3)
    assert row['difference'] == pytest.approx(2/3)
    assert positions.filter(pl.col('position') == 1)['difference'].sum() == 0
    start = motifs.filter(pl.col('position') == 2).row(0, named=True)
    assert start['reference_count'] == 1 and start['control_count'] == 1
    assert motifs.filter(pl.col('position') == 0)['reference_count'].item() == 1
    unsupported, _ = position_contrast(reference, control.head(0), locus='TRB',
        v_gene='TRBV1', j_gene='TRBJ1-1', junction_length=4)
    assert unsupported['control_status'].unique().to_list() == ['unsupported_control_stratum']
    assert unsupported['difference'].null_count() == 80
    with pytest.raises(ValueError, match='no reference'):
        position_contrast(reference, control, locus='TRB', v_gene='TRBV1',
            j_gene='TRBJ2-1', junction_length=4)


@pytest.mark.parametrize('generative', [False, True])
def test_cli_prepared_events_only_records_provenance_and_curates_reference(tmp_path, generative):
    import json
    from vdjmatch.cli.__main__ import main
    from vdjmatch.db.cache import sha256
    raw, reference, model = tmp_path/'raw.tsv', tmp_path/'reference.tsv', tmp_path/'model'
    pl.DataFrame({'junction_aa': ['CAGF', 'CGGF'], 'v_call': ['TRBV1']*2,
        'j_call': ['TRBJ1-1']*2, 'occurrences': [100, 1]}).write_csv(raw, separator='\t')
    assert main(['prepare-control', '--control', str(raw), '--species', 'human',
        '--locus', 'TRB', '--motif-mode', 'exact3', '--output-dir', str(model)]) == 0
    for path in model.iterdir():
        if path.name not in {'events.parquet', 'manifest.json'}:
            path.unlink()
    pl.DataFrame({'cdr3': ['CAGF', 'CATF', 'CCCC'], 'v': ['TRBV1']*3,
        'j': ['TRBJ1-1']*3, 'epitope': ['DEMO']*3, 'species': ['HomoSapiens']*3,
        'gene': ['TRB']*3, 'mhc_a': ['HLA-A*02:01']*3, 'mhc_b': ['B2M']*3,
        'mhc_class': ['MHCI']*3, 'reference_id': ['PMID:1', 'PMID:2', 'PMID:40498839']}
        ).write_csv(reference, separator='\t')
    prefix = tmp_path/'contrast'
    assert main(['motif-position', '--vdjdb', str(reference), '--control', str(raw),
        '--control-model', str(model), '--species', 'human', '--locus', 'TRB',
        '--epitope', 'DEMO', '--v-gene', 'TRBV1', '--j-gene', 'TRBJ1-1',
        '--junction-length', '4', '--motif', 'AG', '--output-prefix', str(prefix)] +
        (['--generative-background', 'olga'] if generative else [])) == 0
    positions = pl.read_csv(str(prefix)+'.positions.tsv', separator='\t')
    assert positions['reference_keys'].unique().to_list() == [2]
    assert positions['control_keys'].unique().to_list() == [2]
    manifest = json.loads((tmp_path/'contrast.manifest.json').read_text())
    assert manifest['control']['control_model']['manifest_sha256'] == sha256(model/'manifest.json')
    assert manifest['native_batches'] == int(generative)
    if generative:
        assert manifest['generative_background']['native_sha256']
        assert positions['generative_status'].unique().to_list() == ['unsupported_generative_stratum']
    assert manifest['tables']['positions']['sha256'] == sha256(tmp_path/'contrast.positions.tsv')


def test_exact_generative_masks_partition_mass_and_match_singleton_pgen(monkeypatch):
    from vdjmatch.aggregate.motif_position import generative_position_contrast, position_contrast
    from vdjtools.model.bundled import load_bundled
    from vdjtools.model.native import pack, pgen_aa
    from importlib import import_module
    import vdjtools.model.native as native
    calls = []
    original = native.pgen_aa_degenerate_batch
    def batched(*args, **kwargs):
        calls.append(len(args[1]))
        return original(*args, **kwargs)
    monkeypatch.setattr(native, 'pgen_aa_degenerate_batch', batched)
    alleles_for = import_module('vdjtools.overlap.alice')._alleles_for
    sequence = 'CAAAAAAF'
    reference = pl.DataFrame({'cdr3': [sequence], 'v': ['TRBV2'], 'j': ['TRBJ1-1']})
    control = reference.rename({'cdr3': 'junction_aa'})
    positions, motifs = position_contrast(reference, control, locus='TRB',
        v_gene='TRBV2', j_gene='TRBJ1-1', junction_length=8, motif=sequence)
    positions, motifs, identity = generative_position_contrast(positions, motifs,
        species='human', locus='TRB', threads=1)
    assert positions['generative_status'].unique().to_list() == ['available']
    for row in positions.group_by('position').agg(pl.col('generative_fraction').sum()).iter_rows():
        assert row[1] == pytest.approx(1, abs=1e-10)
    model = load_bundled('TRB', 'olga', organism='human', collapse=False)
    _, vi, ji = pack(model)
    exact = sum(pgen_aa(model, sequence, va, ja) for va in alleles_for(vi, 'TRBV2')
        for ja in alleles_for(ji, 'TRBJ1-1'))
    denominator = positions['generative_stratum_pgen'][0]
    assert denominator > 0
    assert exact > 0
    assert motifs.filter(pl.col('position') == 1)['generative_fraction'].item() == pytest.approx(exact/denominator)
    assert motifs.filter(pl.col('position') == 0)['generative_status'].item() == 'unsupported_anywhere_union'
    assert motifs.filter(pl.col('position') == 0)['generative_fraction'].item() is None
    assert identity['native_batches'] == 1 and identity['source'] == 'olga'
    assert identity['native_sha256'] and identity['files']['manifest.json']
    assert len(identity['v_alleles']) == 2
    assert calls == [2*(1+20*8+1)]
    unsupported, _, _ = generative_position_contrast(positions.with_columns(
        pl.lit('TRBV1').alias('v_gene')), None, species='human', locus='TRB', threads=1)
    assert unsupported['generative_status'].unique().to_list() == ['unsupported_generative_stratum']
    assert unsupported['generative_fraction'].null_count() == unsupported.height
