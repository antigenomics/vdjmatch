"""Discrete spectratype effects preserve joint-key mass and unsupported lengths."""
import json

import polars as pl
import pytest


def test_length_information_identical_disjoint_and_partial_support():
    from vdjmatch.aggregate.length_profile import length_control, length_summary

    control = length_control(pl.DataFrame({'junction_length': [8, 10], 'keys': [1, 1]}))
    def summary(lengths, counts):
        return length_summary(pl.DataFrame({'junction_length': lengths, 'keys': counts}), control)
    same = summary([8, 10], [2, 2])
    assert same['length_kl_bits'] == same['length_js_bits'] == 0
    selected = summary([8], [4])
    assert selected['length_kl_bits'] == 1
    assert selected['length_js_bits'] == pytest.approx(.3112781244591328)
    absent = summary([9], [4])
    assert absent['length_kl_bits'] == float('inf')
    assert absent['length_js_bits'] == 1
    partial = summary([8, 9], [1, 1])
    assert partial['length_kl_bits'] == float('inf')
    assert partial['length_js_bits'] == .5
    empty = length_summary(pl.DataFrame(schema={'junction_length': pl.Int64, 'keys': pl.Int64}), control)
    assert empty['length_kl_bits'] is empty['length_js_bits'] is None


def test_cdf_effects_keep_unequal_grids_and_discrete_tail_boundaries():
    from vdjmatch.aggregate.length_profile import length_control, length_summary

    control = length_control(pl.DataFrame({'junction_length': [8, 10, 12], 'keys': [1, 2, 1]}))
    row = length_summary(pl.DataFrame({'junction_length': [7, 9, 13], 'keys': [1, 1, 2]}), control)
    assert row['cdf_max_distance'] == pytest.approx(.5)
    assert row['wasserstein1_residues'] == pytest.approx(1.5)
    assert row['mean_length_shift_residues'] == pytest.approx(.5)
    assert row['median_length_shift_residues'] == pytest.approx(1.)
    assert row['reference_keys'] == row['control_keys'] == 4
    assert row['reference_tail_keys'] == 3
    assert row['reference_tail_fraction'] == .75
    assert row['control_tail_keys'] == 0
    tied = length_control(pl.DataFrame({'junction_length': [8, 10, 12], 'keys': [1, 38, 1]}))
    row = length_summary(pl.DataFrame({'junction_length': [7, 8, 10, 12], 'keys': [1, 1, 1, 1]}), tied)
    assert (row['control_central95_lower_length'], row['control_central95_upper_length']) == (8, 10)
    assert row['reference_tail_keys'] == 2  # endpoints stay inside the interval
    assert row['control_tail_keys'] == 1
    assert row['control_tail_fraction'] == .025


def test_joint_key_histogram_pools_alleles_and_keeps_missing_gene_keys():
    from vdjmatch.aggregate.length_profile import length_histogram

    rows = pl.DataFrame({'cdr3': ['CASSF'] * 5 + ['CAXF', None],
        'v': ['TRBV19*01', 'TRBV19*02', None, 'TRBV19,TRBV20', 'TRBV19', 'TRBV19', 'TRBV19'],
        'j': ['TRBJ2-7'] * 4 + ['TRBJ1-1', 'TRBJ2-7', 'TRBJ2-7']})
    histogram, counts = length_histogram(rows, locus='TRB')
    assert histogram.to_dicts() == [{'junction_length': 5, 'keys': 3, 'unresolved_gene_keys': 1}]
    assert counts == {'source_rows': 7, 'selected_rows': 5, 'dropped_rows': 2,
        'selected_keys': 3, 'unresolved_gene_keys': 1}


def test_reference_structure_exports_spectratype_without_changing_native_batches(tmp_path):
    from vdjmatch.cli.__main__ import main

    reference, control = tmp_path/'reference.tsv', tmp_path/'control.tsv'
    pl.DataFrame({'cdr3': ['CASSAAAF', 'CASSAAAAYEQYF', 'CASSAAAF'],
        'v': ['TRBV19'] * 3, 'j': ['TRBJ2-7'] * 3, 'epitope': ['A'] * 3,
        'species': ['HomoSapiens'] * 3, 'gene': ['TRB'] * 3,
        'mhc_a': ['HLA-A*02:01'] * 3, 'mhc_b': ['B2M'] * 3,
        'mhc_class': ['MHCI'] * 3}).write_csv(reference, separator='\t')
    pl.DataFrame({'junction_aa': ['CASSAAAAAF', 'CASSAAAAAYQF'],
        'v_call': ['TRBV19'] * 2, 'j_call': ['TRBJ2-7'] * 2,
        'occurrences': [99, 1]}).write_csv(control, separator='\t')
    prefix = tmp_path/'profile'
    assert main(['reference-structure', '--vdjdb', str(reference), '--locus', 'TRB',
        '--control', str(control), '--publicity', '--threads', '1', '--output-prefix', str(prefix)]) == 0
    row = pl.read_csv(str(prefix)+'.spectratype.tsv', separator='\t').row(0, named=True)
    assert row['epitope'] == 'A' and row['source_observations'] == 3
    assert row['reference_keys'] == row['control_keys'] == 2  # event weights are not key weights
    assert row['cdf_max_distance'] == .5
    assert row['wasserstein1_residues'] == 1.5
    assert row['length_kl_bits'] == float('inf')
    assert row['length_js_bits'] == 1
    degrees = pl.read_csv(str(prefix)+'.degrees.tsv', separator='\t')
    assert degrees['inclusive_degree'].to_list() == (degrees['degree']+1).to_list()
    acquisition = pl.read_csv(str(prefix)+'.graph-acquisition.tsv', separator='\t')
    assert acquisition['sample_junctions'].to_list() == [1, 2]
    assert acquisition['expected_connected_fraction'].to_list() == [0., 0.]
    structure = pl.read_csv(str(prefix)+'.structures.tsv', separator='\t').row(0, named=True)
    assert structure['publicity_valid_source_records'] == 3
    assert structure['publicity_control_junctions'] == 2
    assert structure['publicity_matched_source_records'] == 0
    assert structure['publicity_record_fraction'] == 0
    bins = pl.read_csv(str(prefix)+'.spectratype-lengths.tsv', separator='\t')
    assert bins.filter(pl.col('population') == 'reference')['junction_length'].to_list() == [8, 13]
    assert bins.filter(pl.col('population') == 'reference')['keys'].to_list() == [1, 1]
    manifest = json.loads((tmp_path/'profile.manifest.json').read_text())
    assert manifest['native_batches'] == 1
    assert 'profile.spectratype.tsv' in manifest['tables']
    assert 'profile.spectratype-lengths.tsv' in manifest['tables']
    assert manifest['spectratype']['control_cdf_builds'] == 1
