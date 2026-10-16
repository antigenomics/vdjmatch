"""Focused synthetic checks for event-weighted, Pgen-matched mock databases."""
import numpy as np
import polars as pl
import pytest

from vdjmatch.evalue.mock_reference import normalize_control, compute_pgen, sample_mock_references, iter_mock_references, pgen_strata


def test_mock_event_annotations_stay_joint_with_rearrangement_weights():
    from vdjmatch.evalue.mock_reference import control_events
    raw = pl.DataFrame({'junction_aa': ['CAAF', 'CAAF', 'CACF'],
        'junction': ['AAA', 'CCC', 'GGG'], 'v_call': ['TRBV1', 'TRBV2', 'TRBV3'],
        'j_call': ['TRBJ1', 'TRBJ2', 'TRBJ1'], 'incidence': [2, 3, 1]})
    events = control_events(raw)
    assert events.select('junction_aa', 'v', 'j', 'event_weight').rows() == [
        ('CAAF', 'TRBV1', 'TRBJ1', 2), ('CAAF', 'TRBV2', 'TRBJ2', 3),
        ('CACF', 'TRBV3', 'TRBJ1', 1)]
    reference = pl.DataFrame({'group_id': [1] * 6, 'junction_aa': ['CREF'] * 6})
    pgen = pl.DataFrame({'junction_aa': ['CAAF', 'CACF', 'CREF'], 'pgen': [0.01] * 3})
    mock = sample_mock_references(reference, events, replicates=1, pgen=pgen)
    assert mock.group_by('junction_aa', 'v', 'j').len().sort('junction_aa', 'v')['len'].to_list() == [2, 3, 1]
    with pytest.raises(ValueError, match='convergence totals'):
        control_events(raw.drop('junction', 'incidence').with_columns(pl.lit(3).alias('convergence')))
    with pytest.raises(ValueError, match='requires supplied V/J'):
        control_events(raw.select('junction_aa'))


def test_event_weights_ignore_reads_and_deduplicate_nt_events():
    raw = pl.DataFrame({'junction_aa': ['CASSF'] * 4 + ['CATF'],
        'junction': ['AAA', 'AAA', 'CCC', 'CCC', 'GGG'],
        'sample_id': ['a', 'a', 'a', 'b', 'a'], 'v_call': ['TRBV1'] * 5,
        'j_call': ['TRBJ1'] * 5, 'duplicate_count': [100, 200, 1, 2, 1000]})
    assert normalize_control(raw).select('junction_aa', 'event_weight').rows() == [('CASSF', 3), ('CATF', 1)]
    assert normalize_control(raw.with_columns(pl.Series('occurrences', [2, 3, 4, 5, 6])))['event_weight'].to_list() == [14, 6]
    assert normalize_control(raw.select('junction_aa', 'duplicate_count'))['event_weight'].to_list() == [1, 1]
    assert normalize_control(raw.select('junction_aa', 'sample_id').with_columns(pl.lit(2).alias('convergence')))['event_weight'].to_list() == [4, 2]


def test_sequence_and_population_contract():
    assert normalize_control(pl.DataFrame({'cdr3aa': [' cassf ']}))['junction_aa'].to_list() == ['CASSF']
    with pytest.raises(ValueError, match='AIRR'):
        normalize_control(pl.DataFrame({'cdr3_aa': ['ASS']}))
    for fields in ({'species': ['mouse']}, {'locus': ['TRA']}):
        with pytest.raises(ValueError, match='match'):
            normalize_control(pl.DataFrame({'junction_aa': ['CASSF'], **fields}))


def test_fixed_bins_source_row_counts_zero_and_seed():
    reference = pl.DataFrame({'group_id': ['x', 'x', '__database__', 'zero'],
                              'junction_aa': ['CREF', 'CREF', 'CREF', 'CZERF']})
    control = normalize_control(pl.DataFrame({'junction_aa': ['CAAF', 'CACF', 'CADF', 'CAEF', 'CAFF', 'CAGF'],
                                             'occurrences': [1, 2, 3, 4, 5, 6]}))
    probabilities = pl.DataFrame({'junction_aa': ['CREF', 'CZERF', *control['junction_aa']],
                                 'pgen': [2**-4, 0, 2**-4.49, 2**-3.51, 2**-4, 0, 0, 2**-6]})
    actual = sample_mock_references(reference, control, replicates=3, seed=3, pgen=probabilities)
    assert actual.group_by('replicate', 'group_id').len().sort('replicate', 'group_id')['len'].to_list() == [1, 2, 1] * 3
    assert set(actual.filter(pl.col('group_id') == 'zero')['junction_aa']) <= {'CAEF', 'CAFF'}
    assert set(actual.filter(pl.col('group_id') != 'zero')['junction_aa']) <= {'CAAF', 'CACF', 'CADF'}
    assert actual.equals(sample_mock_references(reference.reverse(), control.reverse(), replicates=3, seed=3, pgen=probabilities.reverse()))
    capacities = dict(control.select('junction_aa', 'event_weight').iter_rows())
    for _, _, sequence, count in actual.group_by('replicate', 'group_id', 'junction_aa').len().iter_rows():
        assert count <= capacities[sequence]
    with pytest.raises(ValueError, match='insufficient'):
        sample_mock_references(reference, control.head(1), replicates=1, pgen=probabilities)


def test_native_pgen_one_union_batch_and_ignores_supplied_scores(monkeypatch):
    from vdjtools import model
    calls = []
    monkeypatch.setattr(model, 'load_bundled', lambda *a, **kw: calls.append((a, kw)) or 'model')
    def batch(m, sequences, **kw):
        calls.append((m, sequences, kw))
        return np.full(len(sequences), 0.01)
    monkeypatch.setattr(model.native, 'pgen_aa_batch', batch)
    reference = pl.DataFrame({'junction_aa': ['CATF', 'CATF'], 'pgen': [0, 1]})
    control = pl.DataFrame({'junction_aa': ['CAAF', 'CATF'], 'pgen': [0, 1]})
    result = compute_pgen(reference, control, species='mouse', locus='TRA', threads=2)
    assert result['pgen'].to_list() == [0.01, 0.01]
    assert calls == [(('TRA', 'arda'), {'organism': 'mouse', 'collapse': False}),
                     ('model', ['CAAF', 'CATF'], {'mismatches': 0, 'threads': 2})]


def test_mixed_population_filter_and_weighted_stream():
    raw = pl.DataFrame({'junction_aa': ['CAAF', 'CACF', 'CADF', 'CAEF'],
        'species': ['human', 'human', 'mouse', None], 'locus': ['TRB'] * 4,
        'occurrences': [1000000, 1, 4, 5]})
    control = normalize_control(raw)
    assert control['junction_aa'].to_list() == ['CAAF', 'CACF']
    reference = pl.DataFrame({'group_id': [1], 'junction_aa': ['CREF']}).with_columns(pl.col('group_id').cast(pl.UInt32))
    probabilities = pl.DataFrame({'junction_aa': ['CREF', 'CAAF', 'CACF'], 'pgen': [0.1] * 3})
    stream = iter_mock_references(reference, control, replicates=3, seed=0, pgen=probabilities)
    assert iter(stream) is stream
    samples = list(stream)
    assert [part['replicate'].item() for part in samples] == [0, 1, 2]
    assert all(part['junction_aa'].item() == 'CAAF' for part in samples)
    assert all(part.schema['group_id'] == pl.UInt32 for part in samples)


def test_shared_pgen_bins_and_incidental_exact_mock_matches():
    probabilities = pl.DataFrame({'junction_aa': ['CAAF', 'CACF', 'CADF', 'CAEF'],
        'pgen': [0.0, 2**-4.49, 2**-3.51, 2**-6]})
    assert pgen_strata(probabilities)['pgen_bin'].to_list() == [-2147483648, -4, -4, -6]
    reference = pl.DataFrame({'group_id': [0], 'junction_aa': ['CAAF']})
    control = normalize_control(pl.DataFrame({'junction_aa': ['CAAF']}))
    actual = sample_mock_references(reference, control, replicates=2, pgen=probabilities)
    assert actual['junction_aa'].to_list() == ['CAAF', 'CAAF']


def test_pooled_nt_incidence_counts_events_without_read_weights():
    raw = pl.DataFrame({'junction_aa': ['CASSF'] * 3, 'junction': ['AAA', 'CCC', 'AAA'],
        'v_call': ['TRBV1'] * 3, 'j_call': ['TRBJ1'] * 3,
        'incidence': [3, 4, 3], 'convergence': [2, 2, 2], 'reads': [100, 20, 100]})
    assert normalize_control(raw)['event_weight'].to_list() == [7]
    # Sample identities establish distinct observed events; incidence does not multiply them.
    assert normalize_control(raw.with_columns(pl.Series('sample_id', ['a', 'a', 'a'])))['event_weight'].to_list() == [2]


def test_event_slots_preserve_aa_rearrangement_multiplicity():
    reference = pl.DataFrame({'group_id': [0] * 7, 'junction_aa': ['CREF'] * 7})
    control = normalize_control(pl.DataFrame({'junction_aa': ['CAAF', 'CACF'], 'occurrences': [3, 4]}))
    pgen = pl.DataFrame({'junction_aa': ['CREF', 'CAAF', 'CACF'], 'pgen': [0.01] * 3})
    actual = sample_mock_references(reference, control, replicates=2, seed=1, pgen=pgen)
    assert actual.height == 14
    assert actual.group_by('replicate', 'junction_aa').len().sort('replicate', 'junction_aa')['len'].to_list() == [3, 4, 3, 4]
    with pytest.raises(ValueError, match='insufficient control events'):
        sample_mock_references(pl.concat([reference, reference.head(1)]), control, replicates=1, pgen=pgen)
    with pytest.raises(ValueError, match='positive integers'):
        sample_mock_references(reference, control.with_columns(pl.lit(1.5).alias('event_weight')), replicates=1, pgen=pgen)


def test_independent_paired_sampler_source_counts_and_chain_bins(monkeypatch):
    from vdjmatch.evalue import mock_reference as module
    reference = pl.DataFrame({'group_id': [0] * 8 + [1] * 4,
        'junction_aa_alpha': ['CREF'] * 6 + ['CRAF'] * 2 + ['CREF'] * 3 + ['CRAF'],
        'junction_aa_beta': ['CREF'] * 6 + ['CRGF'] * 2 + ['CREF'] * 3 + ['CRGF']}).with_columns(
            pl.col('group_id').cast(pl.UInt32))
    control = pl.DataFrame({'junction_aa': ['CAAF', 'CACF', 'CADF', 'CAEF'] * 2 + ['CMOUSEF'],
        'locus': ['TRA'] * 4 + ['TRB'] * 5, 'species': ['human'] * 8 + ['mouse'],
        'occurrences': [4, 4, 4, 4] * 2 + [100]})
    calls = []
    def compute(raw, pool, **options):
        calls.append((options['locus'], set(pool['junction_aa'])))
        sequences = pl.concat([raw.select('junction_aa'), pool.select('junction_aa')]).unique()
        return sequences.with_columns(pl.when(pl.col('junction_aa').is_in(['CRAF', 'CRGF', 'CADF', 'CAEF']))
            .then(0.0).otherwise(0.01).alias('pgen'))
    monkeypatch.setattr(module, 'compute_pgen', compute)
    samples = list(module.iter_mock_paired_references(reference, control,
        species='human', replicates=2, seed=13, threads=1))
    assert calls == [('TRA', {'CAAF', 'CACF', 'CADF', 'CAEF'}), ('TRB', {'CAAF', 'CACF', 'CADF', 'CAEF'})]
    for sample in samples:
        assert sample.schema['group_id'] == pl.UInt32
        assert sample.group_by('group_id').len().sort('group_id')['len'].to_list() == [8, 4]
        for lane in ('alpha', 'beta'):
            assert sample.with_columns(pl.col(f'junction_aa_{lane}').is_in(['CADF', 'CAEF']).alias('zero')).group_by('group_id').agg(
                pl.col('zero').sum()).sort('group_id')['zero'].to_list() == [2, 1]
    assert any(part['junction_aa_alpha'].to_list() != sorted(part['junction_aa_alpha'].to_list())
               for sample in samples for part in sample.partition_by('group_id'))
    assert any(sample['junction_aa_alpha'].to_list() != sample['junction_aa_beta'].to_list() for sample in samples)
    reversed_samples = list(module.iter_mock_paired_references(reference.reverse(), control.reverse(),
        species='human', replicates=2, seed=13, threads=1))
    assert all(a.equals(b) for a, b in zip(samples, reversed_samples))


def test_paired_sampler_requires_locus_and_complete_genuine_pairs():
    from vdjmatch.evalue.mock_reference import iter_mock_paired_references
    reference = pl.DataFrame({'group_id': [0], 'junction_aa_alpha': ['CREF'], 'junction_aa_beta': [None]})
    control = pl.DataFrame({'junction_aa': ['CAAF', 'CACF'], 'locus': ['TRA', 'TRB'], 'species': ['mouse', 'mouse']})
    with pytest.raises(ValueError, match='both valid'):
        list(iter_mock_paired_references(reference, control, replicates=1))
    reference = reference.with_columns(pl.lit('CREF').alias('junction_aa_beta'))
    with pytest.raises(ValueError, match='locus metadata'):
        list(iter_mock_paired_references(reference, control.drop('locus'), replicates=1))
    with pytest.raises(ValueError, match='matching requested'):
        list(iter_mock_paired_references(reference, control, species='human', replicates=1))


def test_control_locus_inference_and_conflicting_gene_metadata():
    raw = pl.DataFrame({'junction_aa': ['CAAF', 'CACF'], 'v_call': [' trav1 ', 'trbv1'],
        'j_call': ['traj1', 'trbj1']})
    assert normalize_control(raw, locus='TRB')['junction_aa'].to_list() == ['CACF']
    assert normalize_control(raw, locus='TRA')['junction_aa'].to_list() == ['CAAF']
    with pytest.raises(ValueError, match='conflicts with explicit'):
        normalize_control(raw.with_columns(pl.lit('TRB').alias('locus')))
    with pytest.raises(ValueError, match='conflicting V/J'):
        normalize_control(raw.with_columns(pl.lit('TRBJ1').alias('j_call')))


def test_control_nt_keys_canonicalized_and_incomplete_events_fail():
    raw = pl.DataFrame({'junction_aa': ['CASSF'] * 3, 'junction': ['AAA', ' aaa ', 'CCC'],
        'v_call': ['TRBV1', ' trbv1 ', 'TRBV1'], 'j_call': ['TRBJ1', ' trbj1 ', 'TRBJ1'],
        'sample_id': ['SampleA'] * 3})
    assert normalize_control(raw)['event_weight'].to_list() == [2]
    assert normalize_control(raw.select('junction_aa', 'junction', 'sample_id'))['event_weight'].to_list() == [2]
    incomplete = raw.with_columns(pl.Series('junction', ['AAA', 'CCC', None]))
    with pytest.raises(ValueError, match='incomplete NT event metadata'):
        normalize_control(incomplete)
    with pytest.raises(ValueError, match='incomplete NT event metadata'):
        normalize_control(raw.with_columns(pl.Series('v_call', ['TRBV1', None, 'TRBV1']),
            pl.Series('j_call', ['TRBJ1', None, 'TRBJ1'])))
    assert normalize_control(incomplete.with_columns(pl.lit(3).alias('convergence')))['event_weight'].to_list() == [3]
    # IDs are opaque: case differences remain distinct observed samples.
    assert normalize_control(raw.with_columns(pl.Series('sample_id', ['SampleA', 'samplea', 'SampleA'])))['event_weight'].to_list() == [3]


def test_mock_models_restore_identical_linked_roles_without_inventing_pairs():
    from vdjmatch.evalue.mock_reference import iter_mock_models
    from vdjmatch.evalue.paired import build_paired_ref
    reference = pl.DataFrame({'gene': ['TRA', 'TRA', 'TRB', 'TRB', 'TRA', 'TRB'],
        'cdr3': ['CAAF', 'CAAF', 'CACF', 'CACF', 'CADF', 'CAEF'],
        'v': ['TRAV1', 'TRAV1', 'TRBV1', 'TRBV1', 'TRAV1', 'TRBV1'],
        'j': ['TRAJ1', 'TRAJ1', 'TRBJ1', 'TRBJ1', 'TRAJ1', 'TRBJ1'],
        'complex_id': ['p1', 'p1', 'p1', 'p1', '0', '0'],
        'epitope': ['A', 'A', 'A', 'A', 'B', 'B']})
    controls = {chain: pl.DataFrame({'junction_aa': ['CAAF', 'CADF'] if chain == 'TRA' else ['CACF', 'CAEF'],
        'v': [chain + 'V1'] * 2, 'j': [chain + 'J1'] * 2, 'event_weight': [1, 1]}) for chain in ('TRA', 'TRB')}
    pgen = {chain: pool.select('junction_aa').with_columns(pl.lit(0.01).alias('pgen')) for chain, pool in controls.items()}
    assert build_paired_ref(reference).height == 1
    for mock in iter_mock_models(reference, controls, pgen=pgen, replicates=2, seed=3):
        assert mock.height == reference.height
        assert build_paired_ref(mock).height == 1
        assert mock.group_by('complex_id', 'gene').agg(pl.col('cdr3').n_unique()).filter(
            pl.col('complex_id') == 'p1')['cdr3'].to_list() == [1, 1]
        assert mock.group_by('epitope').len().sort('epitope')['len'].to_list() == [4, 2]


def test_event_slot_sampler_matches_exhaustive_hypergeometric_law():
    """Every equally likely event subset maps to the exact multivariate law."""
    from collections import Counter
    from itertools import combinations
    from math import comb
    from vdjmatch.evalue.mock_reference import _sample_event_indices

    weight = np.array([1, 3, 2], dtype=np.int64)
    cumulative = weight.cumsum()
    observed = Counter()

    class Draw:
        def choice(self, total, count, *, replace, shuffle):
            assert (total, count, replace, shuffle) == (6, 3, False, False)
            return np.array(self.slots)

    draw = Draw()
    for draw.slots in combinations(range(6), 3):
        indices = _sample_event_indices(cumulative, 3, draw)
        observed[tuple(np.bincount(indices, minlength=3))] += 1
    assert sum(observed.values()) == comb(6, 3)
    for counts, ways in observed.items():
        assert ways == comb(1, counts[0]) * comb(3, counts[1]) * comb(2, counts[2])
    for count in (0, 1, 6):
        a = _sample_event_indices(cumulative, count, np.random.default_rng(7))
        b = _sample_event_indices(cumulative, count, np.random.default_rng(7))
        assert np.array_equal(a, b)
        assert len(a) == count
        assert np.all(np.bincount(a, minlength=3) <= weight)


def test_prepared_control_preserves_seeded_draws_and_rejects_wrong_inputs():
    from vdjmatch.evalue.mock_reference import _prepare_mock_control
    reference = pl.DataFrame({'group_id': [2, 0, 2, 0],
        'junction_aa': ['CAAF', 'CACF', 'CAAF', 'CACF']})
    control = pl.DataFrame({'junction_aa': ['CACF', 'CAAF', 'CAAF', 'CADF'],
        'v': ['TRBV2', 'TRBV1', 'TRBV2', 'TRBV1'], 'j': ['TRBJ1'] * 4,
        'event_weight': [7, 5, 3, 9], 'pgen': [0.] * 4})
    pgen = pl.DataFrame({'junction_aa': ['CAAF', 'CACF', 'CADF'],
        'pgen': [.25, .25, .125]})
    prepared = _prepare_mock_control(reference, control, pgen=pgen)
    options = dict(replicates=4, seed=19, pgen=pgen)
    original = list(iter_mock_references(reference, control, **options))
    reused = list(iter_mock_references(reference, control, prepared_control=prepared, **options))
    assert all(a.equals(b) for a, b in zip(original, reused))
    assert all(not cumulative.flags.writeable for _, cumulative in prepared[-1].values())
    for other_control, other_pgen, species, locus in (
            (control.clone(), pgen, 'human', 'TRB'), (control, pgen.clone(), 'human', 'TRB'),
            (control, pgen, 'mouse', 'TRB'), (control, pgen, 'human', 'TRA')):
        with pytest.raises(ValueError, match='same control/Pgen'):
            list(iter_mock_references(reference, other_control, replicates=2,
                pgen=other_pgen, species=species, locus=locus, prepared_control=prepared))
