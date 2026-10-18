"""Finite-size, selection-aware reference-only gene/motif evidence."""
import polars as pl


AA = 'ACDEFGHIKLMNPQRSTVWY'


def sequences(core, n=30, prefix='CAS'):
    return [prefix+core+AA[i//20]+AA[i%20]+'QYF' for i in range(n)]


def rows(groups):
    return pl.DataFrame([(target, seq, gene, str(target))
        for target, seqs, gene in groups for seq in seqs],
        schema=['_target', 'cdr3', 'v', 'epitope'], orient='row')


def test_planted_interior_motif_survives_selection_aware_null_and_duplicates():
    from vdjmatch.aggregate.profile_noise import profile_noise
    source = rows([(0, sequences('WGG'), 'TRBV19'), (1, sequences('LLL', 90), 'TRBV19')])
    result = profile_noise(source, locus='TRB', eligible_targets=[0])
    assert result == profile_noise(pl.concat([source, source]), locus='TRB', eligible_targets=[0])
    value = result[0]
    legacy = profile_noise(source, locus='TRB', eligible_targets=[0], shuffles=199)[0]
    assert legacy['noise_shuffles'] == 199 and legacy['motif_p'] == .005
    assert legacy['motif'] == value['motif']
    assert value['noise_shuffles'] == 999
    assert value['gene_kl_bits'] == value['gene_excess_bits'] == 0.
    assert value['gene_p'] == 1. and not value['gene_supported']
    assert value['motif_p'] == .001 and value['motif_supported']
    assert value['motif_excess_bits'] > 0
    assert value['motif_reference_with_motif'] == value['motif_reference_v_keys'] == 30
    assert value['motif_background_v_keys'] == 120
    assert value['motif'] == 'WGG' and value['motif_v_gene'] == 'TRBV19'
    assert value['motif_background_with_motif'] == 30
    assert value['motif_kl_bits'] > value['motif_null_max_bits']


def test_v_only_contrast_does_not_manufacture_conditional_motif_support():
    from vdjmatch.aggregate.profile_noise import profile_noise
    source = rows([(0, sequences('GGG'), 'TRBV19'), (1, sequences('GGG'), 'TRBV20-1')])
    result = profile_noise(source, locus='TRB', eligible_targets=[0])[0]
    assert result['gene_kl_bits'] == 1.
    assert result['gene_p'] == .001 and result['gene_supported']
    assert result['motif_kl_bits'] == 0.
    assert result['motif_p'] == 1. and not result['motif_supported']
    assert result['motif'] is None
    assert result['noise_status'] == 'observed'


def test_shared_membership_vectors_and_germline_prefix_are_not_motif_evidence():
    from vdjmatch.aggregate.profile_noise import profile_noise
    common = sequences('GGG')
    shared = profile_noise(rows([(0, common, 'TRBV19'), (1, common, 'TRBV19')]),
        locus='TRB', eligible_targets=[0, 1])
    assert shared[0] == shared[1]
    assert shared[0]['gene_p'] == shared[0]['motif_p'] == 1.
    prefix = profile_noise(rows([(0, common, 'TRBV19'),
        (1, sequences('GGG', prefix='CAV'), 'TRBV19')]), locus='TRB', eligible_targets=[0])[0]
    assert prefix['motif_kl_bits'] == 0. and prefix['motif_p'] == 1.


def test_unknown_genes_and_degenerate_populations_remain_explicit():
    from vdjmatch.aggregate.profile_noise import profile_noise
    source = rows([(0, sequences('WGG'), None), (1, sequences('LLL'), 'TRBV19')])
    value = profile_noise(source, locus='TRB', eligible_targets=[0])[0]
    assert value['noise_status'] == 'observed'
    assert value['gene_supported'] is False and value['motif_supported'] is True
    background = profile_noise(source, locus='TRB', eligible_targets=[1])[1]
    assert background['noise_status'] == 'observed'
    assert isinstance(background['gene_supported'],bool) and isinstance(background['motif_supported'],bool)
    assert source['v'].null_count()==30
    assert value==profile_noise(pl.concat([source,source]).reverse(),locus='TRB',eligible_targets=[0])[0]
    missing=profile_noise(source.with_columns(pl.lit(None,dtype=pl.String).alias('v')),locus='TRB',eligible_targets=[0])[0]
    assert missing['noise_status']=='observed' and missing['gene_supported'] is False
    assert missing['gene_p']==1. and missing['motif_supported'] is True
    singleton = profile_noise(rows([(0, ['CASSF'], 'TRBV19')]), locus='TRB', eligible_targets=[0])[0]
    assert singleton['noise_status'] == 'insufficient_unique_junctions'
    assert singleton['gene_supported'] is singleton['motif_supported'] is None
    assert profile_noise(source.head(0), locus='TRB', eligible_targets=[]) == {}
    empty = profile_noise(source.head(0), locus='TRB', eligible_targets=[0])[0]
    assert empty['noise_status'] == 'empty_reference'
    assert empty['gene_supported'] is empty['motif_supported'] is None


def test_short_junctions_have_no_core_motif_and_absent_targets_are_retained():
    from vdjmatch.aggregate.profile_noise import profile_noise
    source = rows([(0, ['CASSF', 'CASAF'], 'TRBV19'), (1, ['CASTF', 'CASGF'], 'TRBV19')])
    result = profile_noise(source, locus='TRB', eligible_targets=[0, 2])
    assert result[0]['motif'] is None and result[0]['motif_p'] == 1.
    assert result[0]['motif_kl_bits'] == result[0]['motif_null_max_bits'] == 0.
    assert result[2]['noise_status'] == 'empty_target'
    assert result[2]['gene_supported'] is result[2]['motif_supported'] is None


def test_length_divergence_is_descriptive_and_preserves_gene_motif_nulls():
    import pytest
    from vdjmatch.aggregate.profile_noise import profile_noise
    source = rows([(0, sequences('WGG'), 'TRBV19'), (1, sequences('LLLL', 90), 'TRBV19')])
    row = profile_noise(source, locus='TRB', eligible_targets=[0, 2])
    assert row[0]['length_database_kl_bits'] == 2.
    assert row[0]['length_database_js_bits'] == pytest.approx(.5487949406953986)
    assert row[2]['length_database_kl_bits'] is row[2]['length_database_js_bits'] is None


def test_full20_union_selects_gapped_fourmer_inside_every_null():
    from vdjmatch.aggregate.profile_noise import profile_noise
    target = ['CASW'+AA[i%20]+'GG'+AA[i//20]+AA[i%20]+'QYF' for i in range(30)]
    source = rows([(0, target, 'TRBV19'), (1, sequences('LVCC', 90), 'TRBV19')])
    exact = profile_noise(source, locus='TRB', eligible_targets=[0])[0]
    union = profile_noise(source, locus='TRB', eligible_targets=[0], motif_modes=('exact3', 'wild4'))[0]
    assert union['motif'] == 'WXGG'
    assert union['motif_p'] == .001 and union['motif_supported']
    assert union['motif_kl_bits'] > exact['motif_kl_bits']
    assert union['motif_reference_with_motif'] == union['motif_reference_v_keys'] == 30
    assert union['noise_shuffles'] == 999
    shared = source.with_columns(pl.lit(0).alias('_target'))
    assert not profile_noise(shared, locus='TRB', eligible_targets=[0],
        motif_modes=('exact3', 'wild4'))[0]['motif_supported']


def test_imputed_v_blocks_preserve_conditional_shuffle_counts():
    import math
    from vdjmatch.aggregate.profile_noise import profile_noise
    source=rows([(0,sequences('WGG'),'TRBV19'),(1,sequences('LLL'),'TRBV20-1'),(2,sequences('AAA'),None)])
    before=source.clone()
    result=profile_noise(source,locus='TRB',eligible_targets=[0,1,2],shuffles=19)
    assert source.equals(before)
    assert all(math.isfinite(row[field]) for row in result.values() for field in ['gene_kl_bits','motif_kl_bits','motif_null_mean_bits','motif_null_max_bits'])
    assert result==profile_noise(source.reverse(),locus='TRB',eligible_targets=[0,1,2],shuffles=19)
