"""Actual retained in-memory component adapters with tiny raw model populations."""
import numpy as np
import polars as pl
import pytest
from seqtree import Index

from vdjmatch.runner.combined import score_components
from vdjmatch.aggregate.candidates import PMHC


def _single():
    seqs = ['CASSLGQETQYF', 'CASSIGQETQYF', 'CATSLGQETQYF']
    reference = pl.DataFrame({'cdr3': seqs, 'v': ['TRBV1'] * 3, 'j': ['TRBJ1-1'] * 3,
        'gene': ['TRB'] * 3, 'species': ['HomoSapiens'] * 3, 'epitope': ['A', 'A', 'B'],
        'mhc_a': ['HLA-A*02:01'] * 3, 'mhc_b': ['B2M'] * 3, 'mhc_class': ['MHCI'] * 3})
    queries = pl.DataFrame({'query_id': ['q', 'missing', 'far'],
        'cdr3': ['CASSMGQETQYF', 'CASSLGQETQYF', 'W' * 200],
        'v': ['TRBV1', None, 'TRBV1'], 'j': ['TRBJ1-1'] * 3, 'locus': ['TRB'] * 3})
    context = {'density_index': Index.build(seqs + ['CASSAGQETQYF'], 'aa'),
               'prior_background': reference.select('cdr3', 'v', 'j')}
    return queries, reference, context


def test_single_epitope_custom_reference_keeps_bounded_kernel_evidence():
    queries, reference, context = _single()
    # The bounded kernel divides by 1+m, so an empty other-epitope population
    # has m=0 and a defined score. Final control calibration supplies its null.
    reference = reference.with_columns(pl.lit('A').alias('epitope'))
    result = score_components(queries.head(1), reference, context,
        species='human', locus='TRB', exclude_exact=True)
    assert result.height == 1
    assert result['gapblock_score'].item() > 0
    assert result['stem4_score'].item() > 0


def test_actual_single_components_preserve_original_oracle_and_missing_rows():
    from vdjmatch.match.historical import density, germline_prior, unpaired_score
    queries, reference, context = _single()
    scored = score_components(queries, reference, context, species='human', locus='TRB', threads=1, exclude_exact=True)
    assert scored.height == 6
    q = scored.filter((pl.col('query_id') == 'q') & (pl.col('epitope') == 'A'))
    raw = density(queries.head(1), reference.filter(pl.col('epitope') == 'A'), context['density_index'],
        threads=1, pool_reference=True, puncture_populations=True)
    prior = germline_prior(queries.head(1), reference.filter(pl.col('epitope') == 'A'), context['prior_background'],
        exclude_exact=True, include_length=True)
    expected = unpaired_score(raw['score'], prior['germline_lr'], raw['n_reference'], log_scale=True)
    assert q['original_score'].item() == expected.item()
    assert q['gapblock_score'].item() > 0 and q['stem4_score'].item() > 0
    assert scored.filter(pl.col('query_id') == 'missing')['original_score'].null_count() == 2
    assert scored.filter(pl.col('query_id') == 'missing')['gapblock_score'].null_count() == 2
    far = scored.filter(pl.col('query_id') == 'far')
    assert far['original_score'].is_finite().all()
    assert far['original_supported'].to_list() == [False, False]
    assert far['gapblock_score'].null_count() == far.height
    mask = scored.filter(pl.col('query_id') == 'q').select('query_id', *PMHC).head(1)
    restricted = score_components(queries, reference, context, species='human', locus='TRB', threads=1,
        exclude_exact=True, target_masks=mask)
    assert restricted.select(scored.columns).equals(scored.filter(pl.col('query_id') == 'q').head(1))


def test_paired_marginals_survive_absent_linked_model_and_original_is_log_sum():
    reference = pl.DataFrame({'cdr3': ['CAVRDSNYQLIW', 'CASSLGQETQYF', 'CAVRDTNYQLIW', 'CASSIGQETQYF'],
        'v': ['TRAV1-2', 'TRBV1', 'TRAV1-2', 'TRBV1'], 'j': ['TRAJ33', 'TRBJ1-1', 'TRAJ33', 'TRBJ1-1'],
        'gene': ['TRA', 'TRB', 'TRA', 'TRB'], 'species': ['HomoSapiens'] * 4,
        'complex_id': ['p1', 'p1', '0', '0'], 'epitope': ['A', 'A', 'B', 'B'],
        'mhc_a': ['HLA-A*02:01'] * 4, 'mhc_b': ['B2M'] * 4, 'mhc_class': ['MHCI'] * 4})
    query = pl.DataFrame({'query_id': ['pair'], 'cdr3a': ['CAVRDSNYQLIW'], 'cdr3b': ['CASSLGQETQYF'],
        'va': ['TRAV1-2'], 'vb': ['TRBV1'], 'ja': ['TRAJ33'], 'jb': ['TRBJ1-1']})
    context = {'density_index': {chain: Index.build(reference.filter(pl.col('gene') == chain)['cdr3'].to_list(), 'aa')
        for chain in ('TRA', 'TRB')}, 'prior_background': reference.filter(pl.col('gene') == 'TRB').select('cdr3', 'v', 'j')}
    scored = score_components(query, reference, context, species='human', locus='paired', threads=1)
    marginal = scored.filter(pl.col('epitope') == 'B')
    assert marginal['gapblock_linked_score'].item() is None
    assert marginal['gapblock_alpha_score'].item() > 0 and marginal['gapblock_beta_score'].item() > 0
    expected = (np.log(marginal['original_density_alpha'].item() + 1e-6) +
        np.log(marginal['original_density_beta'].item() + 1e-6) +
        marginal['germline_lr_alpha'].item() + marginal['germline_lr_beta'].item())
    assert marginal['original_score'].item() == pytest.approx(expected)
    with pytest.raises(ValueError, match='accepted internal'):
        score_components(query, reference, pl.DataFrame(), species='human', locus='paired')


def test_actual_component_bank_combines_in_memory_with_graph_weights_and_nohit():
    from vdjmatch.runner.combined import combine_components
    query, reference, context = _single()
    query = query.with_columns(pl.lit(-20).alias('pgen_bin'))
    bank = pl.DataFrame({'query_id': ['normalizer1', 'normalizer2', 'normalizer3'],
        'cdr3': ['CASSAGQETQYF', 'CASSVGQETQYF', 'W' * 200], 'v': ['TRBV1'] * 3,
        'j': ['TRBJ1-1'] * 3, 'locus': ['TRB'] * 3, 'pgen_bin': [-20] * 3})
    components = score_components(query, reference, context, species='human', locus='TRB', exclude_exact=True)
    normalizer = score_components(bank, reference, context, species='human', locus='TRB', exclude_exact=True)
    combined = combine_components(components, normalizer, reference, query, species='human', locus='TRB', exclude_exact=True)
    assert combined.height == 6
    assert combined.filter(pl.col('query_id') == 'q')['score'].is_finite().all()
    assert combined.filter(pl.col('query_id') != 'q')['score'].null_count() == 4
    assert combined.filter(pl.col('query_id') == 'q')['axis'].to_list() == ['single', 'single']
    assert np.allclose(combined.select(pl.sum_horizontal('original_weight', 'gapblock_weight', 'stem4_weight')).to_series(), 1)
    with pytest.raises(ValueError, match='Pgen stratum'):
        combine_components(components, normalizer.with_columns(pl.lit(-10).alias('pgen_bin')),
            reference, query, species='human', locus='TRB')


def test_ordinary_exact_singleton_closed_form_retains_historical_density():
    from vdjmatch.runner.combined import combine_components
    query, reference, context = _single()
    reference = reference.with_columns(pl.col('cdr3').alias('epitope'))
    query = query.head(1).with_columns(pl.lit(reference['cdr3'][0]).alias('cdr3'))
    ordinary = score_components(query, reference, context, species='human', locus='TRB', exclude_exact=False)
    exact = ordinary.filter(pl.col('epitope') == reference['cdr3'][0])
    assert exact['historical_original_density'].item() == 0
    # One unique reference AA, four control AAs, one exact control match.
    expected_exact = 1 / max(1 / 4 * 1, 0.01)
    assert exact['original_exact_contribution'].item() == expected_exact
    assert exact['original_density'].item() == expected_exact
    assert exact['original_score'].item() == pytest.approx(np.log(expected_exact + 1e-6) + exact['germline_lr'].item())
    assert exact['original_supported'].item()
    benchmark = score_components(query, reference, context, species='human', locus='TRB', exclude_exact=True)
    exact = benchmark.filter(pl.col('epitope') == reference['cdr3'][0])
    assert not exact['original_supported'].item()
    assert exact['original_score'].item() is None


def test_actual_paired_complete_calibration_retains_linked_roles_and_nohits(monkeypatch):
    from vdjmatch.evalue import combined as calibration
    reference = pl.DataFrame({'cdr3': ['CAVRDSNYQLIW'] * 2 + ['CASSLGQETQYF'] * 2 + ['CAVRDTNYQLIW', 'CASSIGQETQYF'],
        'v': ['TRAV1-2'] * 2 + ['TRBV1'] * 2 + ['TRAV1-2', 'TRBV1'],
        'j': ['TRAJ33'] * 2 + ['TRBJ1-1'] * 2 + ['TRAJ33', 'TRBJ1-1'],
        'gene': ['TRA'] * 2 + ['TRB'] * 2 + ['TRA', 'TRB'], 'species': ['HomoSapiens'] * 6,
        'complex_id': ['p1'] * 4 + ['0', '0'], 'epitope': ['A'] * 4 + ['B', 'B'],
        'mhc_a': ['HLA-A*02:01'] * 6, 'mhc_b': ['B2M'] * 6, 'mhc_class': ['MHCI'] * 6})
    query = pl.DataFrame({'query_id': ['pair', 'far'], 'cdr3a': ['CAVRDSNYQLIW', 'W' * 200],
        'cdr3b': ['CASSLGQETQYF', 'W' * 200], 'va': ['TRAV1-2'] * 2, 'vb': ['TRBV1'] * 2,
        'ja': ['TRAJ33'] * 2, 'jb': ['TRBJ1-1'] * 2})
    control = reference.unique('cdr3').select(pl.col('cdr3').alias('junction_aa'),
        pl.col('v').alias('v_call'), pl.col('j').alias('j_call'), pl.col('gene').alias('locus'), 'species').with_columns(
            pl.lit(1).alias('occurrences'))
    def pgen(reference, control, **options):
        return pl.concat([reference.select('junction_aa'), control.select('junction_aa')]).unique().with_columns(
            pl.lit(0.01).alias('pgen'))
    monkeypatch.setattr(calibration, 'compute_pgen', pgen)
    detail, tails = calibration.calibrate_combined(query, reference, control,
        species='human', locus='paired', replicates=2, seed=17, threads=1)
    assert detail.height == 4
    assert detail.filter(pl.col('query_id') == 'pair')['score'].is_finite().all()
    assert detail.filter(pl.col('query_id') == 'far')['score'].null_count() == 2
    assert tails.height == 6
    assert tails.filter(pl.col('query_id') == 'far')['p_both'].to_list() == [1.0] * 3
    assert tails['reference_draws'].to_list() == [2] * 6
    assert tails['query_draws'].to_list() == [2] * 6
    assert set(tails['test_scope']) == {'epitope', 'database'}
    assert tails['p_both'].to_list() == tails.select(pl.max_horizontal('p_reference', 'p_query')).to_series().to_list()
