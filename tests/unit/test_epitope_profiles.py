"""The primary cohort pools MHC, counts observations and leaves the small tail out."""
import json
import polars as pl


def test_combined_profiles_pool_mhc_and_use_shared_modal_pair(tmp_path, monkeypatch):
    from vdjmatch.cli.__main__ import main
    import vdjmatch.evalue.mock_reference as mocks
    sequences = ['CASSLAAAF', 'CASSLAAAAF', 'CAGGGGGGF']
    rows = []
    for locus, count in [('TRB', 30), ('TRA', 40)]:
        for i in range(count):
            rows.append(dict(cdr3=sequences[i % 3], v='TRBV1' if locus == 'TRB' else 'TRAV1-1',
                j='TRBJ1-1' if locus == 'TRB' else 'TRAJ1', species='HomoSapiens', gene=locus,
                epitope='MAIN', mhc_a='HLA-B*07:02' if locus == 'TRB' and i < 20 else 'HLA-A*02:01',
                mhc_b='B2M', mhc_class='MHCI'))
    rows.extend([rows[0] | dict(epitope='TAIL') for _ in range(29)])
    rows.extend([rows[0] | dict(cdr3='INVALIDX') for _ in range(2)])
    reference, control, prefix = [tmp_path/name for name in ['reference.tsv', 'control.tsv', 'profile']]
    pl.DataFrame(rows).write_csv(reference, separator='\t')
    pl.DataFrame(dict(junction_aa=sequences, v_call=['TRBV1']*3, j_call=['TRBJ1-1']*3,
        occurrences=[100]*3)).write_csv(control, separator='\t')
    monkeypatch.setattr(mocks, 'compute_pgen', lambda ref, pool, **kwargs:
        pool.select('junction_aa').with_columns(pl.lit(.25).alias('pgen')))
    prepare = mocks._prepare_mock_control
    prepared_records = []
    def checked_prepare(ref, *args, **kwargs):
        prepared_records.append(ref.height)
        return prepare(ref, *args, **kwargs)
    monkeypatch.setattr(mocks, '_prepare_mock_control', checked_prepare)
    assert main(['reference-diagnostics', '--vdjdb', str(reference), '--locus', 'TRB',
        '--control', str(control), '--geometry', 'levenshtein', '--mock-replicates', '2',
        '--nboot', '3', '--output-prefix', str(prefix)]) == 0
    assert prepared_records == [30]
    structure = pl.read_csv(str(prefix)+'.structures.tsv', separator='\t')
    assert structure.height == 1 and structure['epitope'].item() == 'MAIN'
    assert structure['mhc_a'].item() == 'HLA-A*02:01'
    manifest = json.loads((tmp_path/'profile.manifest.json').read_text())
    assert manifest['epitope_profiles'] is True and manifest['min_records'] == 30
    inventory = pl.DataFrame(manifest['epitope_inventory'])
    counts = inventory.filter(pl.col('gene') == 'TRB').sort('epitope')
    assert counts['source_records'].to_list() == [32, 29]
    assert counts['valid_source_records'].to_list() == [30, 29]
    assert counts['classification_eligible'].to_list() == [True, False]
    coverage = pl.read_csv(str(prefix)+'.levenshtein.coverage.tsv', separator='\t').filter(
        pl.col('population') == 'reference')
    assert coverage['source_records'].to_list() == [32, 32]
    assert coverage['retained_records'].to_list() == [30, 30]
    owner = json.loads((tmp_path/'profile.levenshtein.manifest.json').read_text())
    assert owner['epitope_labels'] == manifest['epitope_labels']
