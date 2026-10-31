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
    mocked=control_events(raw.select('junction_aa'))
    assert mocked.height==2 and mocked['event_weight'].sum()==2
    assert mocked['v_imputed'].all() and mocked['j_imputed'].all()


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
    calls.clear()
    compute_pgen(reference,control,species='human',locus='TRA',threads=2,source='isalgo')
    assert calls[0] == (('TRA','isalgo'), {'organism':'human','collapse':False})


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


def test_joint_control_imputes_conflicting_v_for_the_requested_locus():
    from vdjmatch.evalue.mock_reference import control_events
    raw = pl.DataFrame({'junction_aa': ['CASSF']*2, 'junction': ['AAA']*2,
        'v_call': ['TRBV19', 'TRAV9D-1'], 'j_call': ['TRBJ1-1']*2})
    original = raw.clone()
    events = control_events(raw, species='mouse', locus='TRB')
    assert raw.equals(original)
    assert events['event_weight'].sum() == 2
    assert events['v'].to_list() == ['TRBV19']
    assert events['v_imputed'].all() and not events['j_imputed'].any()
    assert events.equals(control_events(raw.with_columns(pl.lit('TRB').alias('locus')),
        species='mouse', locus='TRB'))


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


def test_zero_reference_mass_uses_lowest_occupied_control_bin():
    reference=pl.DataFrame({'group_id':[0,0,0], 'junction_aa':['CREF','CREF','CACF']})
    control=pl.DataFrame({'junction_aa':['CAAF','CACF','CADF'],'event_weight':[1,2,10]})
    pgen=pl.DataFrame({'junction_aa':['CREF','CAAF','CACF','CADF'],
        'pgen':[0.,2**-8,2**-6,2**-4]})
    sampled=sample_mock_references(reference,control,replicates=2,pgen=pgen)
    assert sampled.height==6
    assert sampled.group_by('replicate','junction_aa').len().sort('replicate','junction_aa')['len'].to_list()==[1,2,1,2]
    assert set(sampled['junction_aa'])=={'CACF','CADF'}
    assert pgen['pgen'].to_list()==[0.,2**-8,2**-6,2**-4]


def test_absent_positive_and_zero_bins_share_one_capacity_aware_floor():
    from vdjmatch.evalue.mock_reference import _prepare_mock_control

    reference = pl.DataFrame({'group_id': [0]*7,
        'junction_aa': ['CZERO']*2 + ['CPOS']*2 + ['CFLOOR'] + ['CMATCH']*2})
    control = pl.DataFrame({'junction_aa': ['CLOW', 'CFLOOR', 'CMATCH', 'CHIGH'],
        'event_weight': [1, 4, 3, 10]})
    pgen = pl.DataFrame({'junction_aa': ['CZERO', 'CPOS', 'CLOW', 'CFLOOR', 'CMATCH', 'CHIGH'],
        'pgen': [0., 2.**-108, 2.**-8, 2.**-6, 2.**-4, 2.**-2]})
    original = pgen.clone()
    for floor_capacity, expected_floor in [(4, -2), (6, -6)]:
        events = control.with_columns(pl.when(pl.col('junction_aa') == 'CFLOOR')
            .then(floor_capacity).otherwise(pl.col('event_weight')).alias('event_weight'))
        prepared = _prepare_mock_control(reference, events, pgen=pgen)
        bins = dict(prepared[4].iter_rows())
        assert bins['CZERO'] == bins['CPOS'] == expected_floor
        assert bins['CFLOOR'] == -6 and bins['CMATCH'] == -4
        draws = list(iter_mock_references(reference, events, pgen=pgen, prepared_control=prepared,
            replicates=2, seed=7))
        for draw in draws:
            counts = dict(draw.group_by('junction_aa').len().iter_rows())
            assert counts == ({'CHIGH': 4, 'CFLOOR': 1, 'CMATCH': 2} if floor_capacity == 4 else
                {'CFLOOR': 5, 'CMATCH': 2})
        assert pgen.equals(original)
    exhausted = control.with_columns(pl.lit(3).alias('event_weight'))
    with pytest.raises(ValueError, match='insufficient control events'):
        _prepare_mock_control(reference, exhausted, pgen=pgen)
    # An occupied bin keeps its ordinary matching rule; insufficient slots still fail.
    thin_match = control.with_columns(pl.when(pl.col('junction_aa') == 'CMATCH')
        .then(1).otherwise(pl.col('event_weight')).alias('event_weight'))
    with pytest.raises(ValueError, match='Pgen bin -4: need 2, have 1'):
        list(iter_mock_references(reference, thin_match, pgen=pgen, replicates=1))


def test_missing_control_genes_are_imputed_without_losing_source_events(tmp_path):
    from vdjmatch.evalue.mock_reference import control_events
    raw=pl.DataFrame({'junction_aa':['CASSF']*3,'junction':['AAA','AAA','CCC'],'v_call':['TRBV19','TRBV19',None],'j_call':['TRBJ2-7',None,'TRBJ2-7'],'sample_id':['one']*3})
    original=raw.clone();events=control_events(raw)
    from vdjmatch.cli.combined import _raw_control
    path=tmp_path/'partial-control.tsv';raw.write_csv(path,separator='\t')
    retained,report=_raw_control(path,retain_missing_genes=True)
    assert report['control_valid_rows']==3 and control_events(retained).equals(events)
    assert raw.equals(original)
    assert events['event_weight'].sum()==3
    assert events['v'].to_list()==['TRBV19'] and events['j'].to_list()==['TRBJ2-7']
    assert events['v_imputed'].item() and events['j_imputed'].item()
    assert events.equals(control_events(raw.reverse()))
    assert control_events(raw.with_columns(pl.col('j_call').fill_null('')))['event_weight'].sum()==3
    mocked=control_events(raw.with_columns(pl.lit(None,dtype=pl.String).alias('v_call'),pl.lit(None,dtype=pl.String).alias('j_call')))
    assert mocked['event_weight'].sum()==2 and mocked['v_imputed'].all() and mocked['j_imputed'].all()


def test_imputation_retains_original_aa_event_keys_and_infers_missing_declared_locus():
    from vdjmatch.evalue.mock_reference import control_events
    raw=pl.DataFrame({'junction_aa':['CASSF','CASSF'],'v_call':['TRBV19',None],'j_call':['TRBJ2-7']*2})
    assert control_events(raw)['event_weight'].sum()==2
    counted=raw.with_columns(pl.Series('locus',['TRB',None]),pl.Series('occurrences',[2,3]))
    assert control_events(counted)['event_weight'].sum()==5
    assert control_events(counted.with_columns(pl.col('locus').fill_null('')))['event_weight'].sum()==5


def test_epitope_mock_policies_preserve_bins_and_exclude_reference_keys():
    from vdjmatch.evalue import mock_reference as m
    assert hasattr(m, 'iter_epitope_mocks'), 'M3 reference mock sampler is missing'
    reference = pl.DataFrame(dict(group_id=['X','X'], junction_aa=['CASSF','CATSF'],
                                 v=['TRBV2','TRBV3'], j=['TRBJ1-1']*2, pgen_bin=[-10]*2))
    control = pl.DataFrame(dict(junction_aa=['CASSF','CASRF','CAVRF','CARRF','CAARF','CAGRF'],
                               v=['TRBV2','TRBV2','TRBV3','TRBV3','TRBV4','TRBV4'],
                               j=['TRBJ1-1']*6,pgen_bin=[-10]*6))
    for policy in ['pgen','pgen-vj','pgen-different-v']:
        draws=list(m.iter_epitope_mocks(reference,control,replicates=3,seed=7,policy=policy))
        assert all(x.height==2 and x['pgen_bin'].to_list()==[-10,-10] for x in draws)
        assert all('CASSF' not in x['junction_aa'].to_list() for x in draws)
        assert all(x.select('junction_aa','v','j').n_unique()==2 for x in draws)
        if policy=='pgen':
            import hashlib,numpy as np
            keys=['junction_aa','v','j']
            eligible=control.sort(keys).join(reference.select(keys),on=keys,how='anti')
            identity=int.from_bytes(hashlib.sha256(b'X').digest()[:8],'little')
            for b,draw in enumerate(draws):
                rng=np.random.default_rng(np.random.SeedSequence([7,identity,b]))
                expected=eligible[rng.choice(eligible.height,2,replace=False)]
                assert draw.select(keys).equals(expected.select(keys))
        if policy=='pgen-vj':
            assert all(x['v'].to_list()==reference['v'].to_list() for x in draws)
        if policy=='pgen-different-v':
            assert all((x['v']!=reference['v']).all() for x in draws)
        again=list(m.iter_epitope_mocks(reference,control,replicates=3,seed=7,policy=policy))
        assert all(a.equals(b) for a,b in zip(draws,again))
    with pytest.raises(ValueError,match='insufficient'):
        list(m.iter_epitope_mocks(reference,control.head(1),replicates=1))
    with pytest.raises(ValueError,match='bin'):
        list(m.iter_epitope_mocks(reference,control.with_columns(pl.lit(-11).alias('pgen_bin')),replicates=1))


def test_mock_sampler_keeps_pool_index_dtype(monkeypatch):
    """Mixed signed/unsigned search keys copy the entire large control stratum."""
    import numpy as np
    from vdjmatch.evalue.mock_reference import iter_epitope_mocks
    searchsorted = np.searchsorted
    def checked(a, v, *args, **kwargs):
        if a.dtype == np.dtype('uint32'):
            assert v.dtype == a.dtype
        return searchsorted(a, v, *args, **kwargs)
    monkeypatch.setattr(np, 'searchsorted', checked)
    control = pl.DataFrame(dict(junction_aa=['CAAF','CABF','CACF','CADF'],
        v=['V1','V2','V1','V2'], j=['J1']*4, pgen_bin=[-10]*4))
    reference = control.head(2).with_columns(pl.lit('X').alias('group_id'))
    for policy in ['pgen', 'pgen-vj', 'pgen-different-v']:
        assert len(list(iter_epitope_mocks(reference, control, replicates=2, policy=policy))) == 2


def test_compiled_reference_strata_preserves_compiled_r_grouping():
    from vdjmatch.evalue.mock_reference import compiled_reference_strata
    probes=pl.DataFrame(dict(junction_aa=['CAAF','CAAF','CACF'],
        v=['V1','V2','V1'],j=['J1','J1','J2'],pgen=[2.**-20,2.**-18,0.]))
    sample=pl.DataFrame(dict(v=['V1','V1','V2'],j=['J1','J2','J1'],records=[3,1,100]))
    model=pl.DataFrame(dict(v=['V1','V1','V2','V2'],j=['J1','J2','J1','J2'],records=[1,3,1,1]))
    result,audit=compiled_reference_strata(probes,sample,model)
    # Original AA-only join/max; V/J correction is normalized within V because
    # dplyr summarise retains the V grouping. Zero uses the compiled 2^-80 floor.
    assert result['pcdr'].to_list()==[2.**-18,2.**-18,0.]
    assert result['pgen_bin'].to_list()==[-17,-17,-80]
    selected,_=compiled_reference_strata(probes,sample,model,selection_q=4)
    assert selected['pgen_bin'].to_list()==[-15,-15,-80]
    with pytest.raises(ValueError,match='selection_q'):
        compiled_reference_strata(probes,sample,model,selection_q=0)
    assert result['vj_factor'].to_list()==pytest.approx([2.,202/102,.5])
    assert audit.filter(pl.col('v')=='V1')['sample_probability'].sum()==pytest.approx(1.)
    with pytest.raises(ValueError,match='unique'):
        compiled_reference_strata(pl.concat([probes,probes.head(1)]),sample,model)


def test_reference_panel_r_sampler_matches_base_r_seed42():
    import shutil
    if not shutil.which('Rscript'):pytest.skip('optional R reproduction sampler')
    from vdjmatch.evalue.mock_reference import sample_stratified_panels
    probes=pl.DataFrame(dict(junction_aa=['CAAF','CACF','CADF','CAEF','CAFF','CAGF'],
        v=['V1']*6,j=['J1']*6,pgen_bin=[-10]*6))
    reference=probes.head(2).with_columns(pl.lit('X').alias('group_id'))
    result=sample_stratified_panels(reference,probes,replicates=3,seed=42,sampler='r')
    # R: set.seed(42); replicate(3, sample.int(4,2)) => (1,4),(1,4),(2,4).
    assert result['junction_aa'].to_list()==['CADF','CAGF','CADF','CAGF','CAEF','CAGF']
    assert result['replicate'].to_list()==[0,0,1,1,2,2]


def test_profile_strata_corrects_against_selected_vj_probabilities():
    from vdjmatch.evalue.mock_reference import compiled_reference_strata
    probes=pl.DataFrame({'junction_aa':['CASSF','CATSF'], 'v':['V','V'], 'j':['J1','J2'], 'pgen':[.01,.02]})
    sample=pl.DataFrame({'v':['V','V'], 'j':['J1','J2'], 'records':[3,1]})
    selected=pl.DataFrame({'v':['V','V'], 'j':['J1','J2'], 'p':[.8,.2]})
    scored,audit=compiled_reference_strata(probes,sample,model_vj_probabilities=selected,selection_q=2.)
    assert audit['model_probability'].to_list()==[.8,.2]
    assert scored['vj_factor'].to_list()==pytest.approx([(4/6)/.8,(2/6)/.2])
    assert scored['adjusted_pgen'].to_list()==pytest.approx([.02*(4/6)/.8,.04*(2/6)/.2])
    with pytest.raises(ValueError,match='positive'):
        compiled_reference_strata(probes,sample,model_vj_probabilities=selected.with_columns(pl.lit(0.).alias('p')))
