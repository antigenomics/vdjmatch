"""Tiny oracles for fixed-reference morphology, never outcome-fitted classes."""
import polars as pl
import pytest

from vdjmatch.db.schema import normalize


def reference(rows):
    return normalize(pl.DataFrame(rows, schema=['species', 'gene', 'epitope', 'cdr3', 'v', 'mhc_a'], orient='row'))


def test_joint_marginal_support_keeps_missing_chains_unknown():
    from vdjmatch.aggregate.fast_profiles import joint_reference_profiles
    rows = []
    for epitope, a, b in [('A', True, False), ('B', False, False), ('C', None, False), ('D', True, None)]:
        for locus, supported in [('TRA', a), ('TRB', b)]:
            if epitope == 'D' and locus == 'TRB':
                continue
            rows.append(dict(species='HomoSapiens', epitope=epitope, level='epitope', locus=locus,
                valid_source_records=30, classification_eligible=True,
                motif_supported=supported, profile_class='example'))
    result = joint_reference_profiles(pl.DataFrame(rows))
    assert result['any_chain_motif_supported'].to_list() == [True, False, None, True]
    assert result['both_chains_observed'].to_list() == [True, True, True, False]
    assert result['both_chains_eligible'].to_list() == [True, True, True, False]
    assert result['joint_class'][-1] == 'alpha: example; beta: missing_chain'


def test_identical_label_sets_have_enrichment_one_and_duplicate_records_do_not_add_edges():
    from vdjmatch.aggregate.fast_profiles import reference_profiles
    raw = reference([('HomoSapiens', 'TRB', ep, seq, 'TRBV19', 'HLA-A*02:01')
        for ep in ['A', 'B'] for seq in ['CASSF', 'CASAF'] for _ in range(15)])
    profile = reference_profiles(raw, threads=1)
    assert profile.height == 4
    for row in profile.to_dicts():
        assert row['source_records'] == row['valid_source_records'] == 30
        assert row['junctions'] == 2
        assert row['within_neighbor_pairs'] == row['within_opportunities'] == 1
        assert row['between_neighbor_pairs'] == row['between_opportunities'] == 2
        assert row['exact_shared_junctions'] == 2
        assert row['neighbor_enrichment'] == 1.
        assert row['neighbor_coverage'] == 1.
        assert row['max_v_gene_share'] == 1.


def test_indel_neighbors_use_one_batch_per_species_locus_without_species_leak(monkeypatch):
    from vdjmatch.aggregate.fast_profiles import reference_profiles
    import vdjmatch.cluster as cluster
    native, calls = cluster.overlap, []
    def record(*args, **kwargs):
        calls.append((len(args[0]), kwargs['scope'], kwargs['threads']))
        return native(*args, **kwargs)
    monkeypatch.setattr(cluster, 'overlap', record)
    raw = reference([('HomoSapiens', 'TRB', 'A', 'CASSF', 'TRBV19', 'HLA-A*02:01'),
        ('HomoSapiens', 'TRB', 'A', 'CASSSF', 'TRBV19', 'HLA-A*02:01'),
        ('HomoSapiens', 'TRB', 'B', 'CALAF', 'TRBV19', 'HLA-A*02:01'),
        ('MusMusculus', 'TRB', 'B', 'CASAF', 'TRBV19', 'H-2Kb')])
    result = reference_profiles(raw, threads=1)
    assert calls == [(3, '1,1,1,1', 1), (1, '1,1,1,1', 1)]
    human = result.filter((pl.col('species') == 'HomoSapiens') & (pl.col('epitope') == 'A'))
    assert human['within_neighbor_pairs'].to_list() == [1, 1]
    assert human['neighbor_enrichment'].to_list() == [float('inf')] * 2
    assert human['neighbor_coverage'].to_list() == [1., 1.]
    mouse = result.filter(pl.col('species') == 'MusMusculus')
    assert mouse['background_junctions'].to_list() == [0, 0]
    assert mouse['neighbor_enrichment'].to_list() == [None, None]


def test_missing_v_keeps_mass_and_full_pmhc_background_excludes_the_whole_peptide():
    from vdjmatch.aggregate.fast_profiles import reference_profiles
    raw = reference([('HomoSapiens', 'TRB', 'A', 'CASSF', 'TRBV19*01', 'HLA-A*02:01'),
        ('HomoSapiens', 'TRB', 'A', 'CASAF', None, 'HLA-A*02:01'),
        ('HomoSapiens', 'TRB', 'A', 'CASGF', None, 'HLA-A*02:01'),
        ('HomoSapiens', 'TRB', 'A', 'CASTF', None, 'HLA-A*02:01'),
        ('HomoSapiens', 'TRB', 'A', 'CAAAF', 'TRBV19*02', 'HLA-A*03:01'),
        ('HomoSapiens', 'TRB', 'B', 'CATTF', 'TRBV19', None)])
    result = reference_profiles(raw, threads=1)
    target = result.filter((pl.col('level') == 'pmhc') & (pl.col('mhc_a') == 'HLA-A*02:01')).row(0, named=True)
    assert target['reference_v_keys'] == 4 and target['unresolved_v_keys'] == 3
    assert target['max_v_gene_share'] == .25
    assert target['background_junctions'] == 1  # other restriction of A is not background
    assert target['support_status'] == 'insufficient_source_records'
    assert result.filter(pl.col('epitope') == 'B')['mhc_context_status'].to_list() == ['modal_annotation', 'incomplete_mhc']


def test_no_neighbours_and_all_invalid_targets_remain_explicit():
    from vdjmatch.aggregate.fast_profiles import reference_profiles
    raw = reference([('HomoSapiens', 'TRA', 'A', 'CASSF', 'TRAV1-1', None),
        ('HomoSapiens', 'TRA', 'B', 'CWWWYF', 'TRAV1-1', None),
        ('HomoSapiens', 'TRA', 'C', 'CAXF', 'TRAV1-1', None)])
    result = reference_profiles(raw, threads=1)
    assert result.height == 6
    assert result.filter(pl.col('epitope') != 'C')['neighbor_enrichment'].null_count() == 4
    bad = result.filter(pl.col('epitope') == 'C')
    assert bad['excluded_source_records'].to_list() == [1, 1]
    assert bad['junctions'].to_list() == [0, 0]
    assert bad['background_junctions'].to_list() == [2, 2]
    assert bad['max_v_gene_share'].null_count() == 2
    assert reference_profiles(raw.head(0), threads=1).height == 0
    assert reference_profiles(raw.filter(pl.col('epitope') == 'C'), threads=1)['junctions'].to_list() == [0, 0]
    single = reference_profiles(raw.filter(pl.col('epitope') == 'A'), threads=1)
    assert single['background_junctions'].to_list() == [0, 0]
    with pytest.raises(ValueError):
        reference_profiles(raw, threads=0)


def test_hand_counted_shared_background_and_adjacent_length_opportunities():
    from vdjmatch.aggregate.fast_profiles import reference_profiles
    raw = reference([('HomoSapiens', 'TRB', ep, seq, 'TRBV19', mhc)
        for ep, seq, mhc in [('A', 'CASSF', 'HLA-A*02:01'), ('A', 'CASAF', 'HLA-A*02:01'),
            ('A', 'CASSSF', 'HLA-A*03:01'), ('B', 'CASSF', 'HLA-A*02:01'),
            ('B', 'CALAF', 'HLA-A*02:01'), ('C', 'CASTF', 'HLA-A*02:01')]])
    result = reference_profiles(raw, threads=1)
    pooled = result.filter((pl.col('level') == 'epitope') & (pl.col('epitope') == 'A')).row(0, named=True)
    assert pooled['within_opportunities'] == 3 and pooled['within_neighbor_pairs'] == 2
    assert pooled['background_junctions'] == 3 and pooled['exact_shared_junctions'] == 1
    assert pooled['between_opportunities'] == 8 and pooled['between_neighbor_pairs'] == 5
    assert pooled['neighbor_enrichment'] == pytest.approx(16/15)
    assert pooled['neighbor_coverage'] == 1.
    context = result.filter((pl.col('level') == 'pmhc') & (pl.col('epitope') == 'A') &
        (pl.col('mhc_a') == 'HLA-A*02:01')).row(0, named=True)
    assert context['within_opportunities'] == context['within_neighbor_pairs'] == 1
    assert context['between_opportunities'] == 5 and context['between_neighbor_pairs'] == 4
    assert context['neighbor_enrichment'] == pytest.approx(1.25)


def test_a_twofold_ratio_does_not_supply_a_noise_validated_class():
    from vdjmatch.aggregate.fast_profiles import reference_profiles
    raw = reference([('HomoSapiens', 'TRB', 'A', s, 'TRBV19', 'HLA-A*02:01')
        for s in ['CASSF', 'CASAF', 'CASSSF'] for _ in range(10)] +
        [('HomoSapiens', 'TRB', 'B', 'CALAF', 'TRBV19', 'HLA-A*02:01')])
    result = reference_profiles(raw, threads=1).filter(pl.col('epitope') == 'A')
    assert result['neighbor_enrichment'].to_list() == [2., 2.]
    assert result['profile_class'].to_list() == ['no_supported_axis'] * 2
    assert result['gene_p'].to_list() == [1., 1.]
    assert result['motif_p'].to_list() == [1., 1.]
