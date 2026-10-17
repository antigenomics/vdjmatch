"""Complete-score calibration uses fresh models and independent internal banks."""
import polars as pl


def test_complete_score_two_directions_keep_source_roles_and_nohits(monkeypatch):
    import weakref
    from vdjmatch.evalue import combined as module
    from vdjmatch.runner import combined as runner
    from vdjmatch.match import historical
    from vdjmatch.evalue.mock_reference import control_events, iter_mock_models
    query = pl.DataFrame({'query_id': ['q', 'far'],
        'cdr3': ['CASSMGQETQYF', 'W' * 200], 'v': ['TRBV1'] * 2,
        'j': ['TRBJ1-1'] * 2, 'locus': ['TRB'] * 2})
    reference = pl.DataFrame({'cdr3': ['CASSLGQETQYF', 'CASSIGQETQYF', 'CATSLGQETQYF'],
        'v': ['TRBV1'] * 3, 'j': ['TRBJ1-1'] * 3, 'gene': ['TRB'] * 3,
        'species': ['HomoSapiens'] * 3, 'epitope': ['A', 'A', 'B'],
        'mhc_a': ['HLA-A*02:01'] * 3, 'mhc_b': ['B2M'] * 3, 'mhc_class': ['MHCI'] * 3})
    raw = pl.DataFrame({'junction_aa': ['CASSAGQETQYF', 'CASSVGQETQYF', 'CASSTGQETQYF'],
        'v_call': ['TRBV1'] * 3, 'j_call': ['TRBJ1-1'] * 3, 'occurrences': [3, 4, 5]})
    def pgen(ref, control, **kwargs):
        return pl.concat([ref.select('junction_aa'), control.select('junction_aa')]).unique().with_columns(
            pl.lit(0.01).alias('pgen'))
    monkeypatch.setattr(module, 'compute_pgen', pgen)
    transient = []
    released = []
    def track(function):
        def call(*args, **kwargs):
            result = function(*args, **kwargs)
            if result is not None:
                transient.append(weakref.ref(result))
            return result
        return call
    def models(*args, **kwargs):
        for model in iter_mock_models(*args, **kwargs):
            assert all(ref() is None for ref in transient), 'wide scores retained across calibration models'
            transient.clear()
            released.append(True)
            yield model
        assert all(ref() is None for ref in transient), 'last mock scores retained during tail export'
    monkeypatch.setattr(runner, 'score_components', track(runner.score_components))
    monkeypatch.setattr(runner, 'combine_components', track(runner.combine_components))
    monkeypatch.setattr(module, 'iter_mock_models', models)
    prepared = []
    prepare = historical._prepare_prior_background
    def prepare_once(background):
        prepared.append(background.height)
        return prepare(background)
    monkeypatch.setattr(historical, '_prepare_prior_background', prepare_once)
    detail, tails = module.calibrate_combined(query, reference, raw,
        species='human', locus='TRB', replicates=2, seed=7, exclude_exact=True)
    assert prepared == [3]  # Fixed control shared by real and both mock models.
    assert len(released) == 2
    assert detail.height == 4 and tails.height == 6
    assert tails.filter(pl.col('query_id') == 'far')['p_both'].to_list() == [1.] * 3
    assert tails['reference_draws'].to_list() == [2] * 6
    assert tails['query_draws'].to_list() == [2] * 6
    assert tails.select((pl.col('p_both') == pl.max_horizontal('p_reference', 'p_query')).all()).item()
    assert tails.select((pl.col('null_query_match_rate') == pl.col('query_hit_draws') / 2).all()).item()
    assert tails.select((pl.col('null_reference_match_rate') == pl.col('reference_hit_draws') / 2).all()).item()
    events = control_events(raw)
    probabilities = {'TRB': pgen(reference.select(pl.col('cdr3').alias('junction_aa')), events)}
    models = list(iter_mock_models(reference, {'TRB': events}, pgen=probabilities, replicates=2))
    for model in models:
        assert model.select('gene', 'species', 'epitope', 'mhc_a', 'mhc_b', 'mhc_class').equals(
            reference.select('gene', 'species', 'epitope', 'mhc_a', 'mhc_b', 'mhc_class'))
        assert set(model['cdr3']) <= set(raw['junction_aa'])


def test_query_banks_are_separate_and_pairs_are_random_joint_events():
    from vdjmatch.evalue.combined import _query_banks
    from vdjmatch.evalue.mock_reference import control_events
    query = pl.DataFrame({'query_id': ['q'], 'cdr3a': ['CAVRDSNYQLIW'],
        'cdr3b': ['CASSLGQETQYF'], 'va': ['TRAV1-2'], 'vb': ['TRBV1'], 'ja': ['TRAJ33'], 'jb': ['TRBJ1-1']})
    controls = {}
    probabilities = {}
    for chain, sequences, genes, j in [('TRA', ['CAVRDTNYQLIW', 'CAVRDMNYQLIW'], ['TRAV1-2', 'TRAV2'], 'TRAJ33'),
            ('TRB', ['CASSIGQETQYF', 'CASSTGQETQYF'], ['TRBV1', 'TRBV2'], 'TRBJ1-1')]:
        controls[chain] = control_events(pl.DataFrame({'junction_aa': sequences, 'v_call': genes, 'j_call': [j] * 2}), locus=chain)
        probabilities[chain] = pl.DataFrame({'junction_aa': sequences + [query['cdr3a' if chain == 'TRA' else 'cdr3b'][0]],
                                            'pgen': [0.01] * 3})
    tagged, normalizer, test = _query_banks(query, controls, probabilities, locus='paired', replicates=40, seed=3)
    assert tagged['pgen_bin'].item() == '-7:-7'
    assert not set(normalizer['query_id']) & set(test['query_id'])
    for bank in (normalizer, test):
        assert bank.select('cdr3a', 'cdr3b').n_unique() == 4
        assert bank.filter((pl.col('cdr3a') == 'CAVRDTNYQLIW') & (pl.col('va') != 'TRAV1-2')).height == 0


def test_hand_enumerated_complete_score_tails_repeat_label_selection(monkeypatch):
    from vdjmatch.evalue import combined as module
    from vdjmatch.runner import combined as runner
    query = pl.DataFrame({'query_id': ['q', 'none'], 'cdr3': ['CAAF', 'CACF'],
                         'v': ['TRBV1'] * 2, 'j': ['TRBJ1-1'] * 2, 'locus': ['TRB'] * 2})
    reference = pl.DataFrame({'cdr3': ['CADF', 'CAEF'], 'v': ['TRBV1'] * 2,
        'j': ['TRBJ1-1'] * 2, 'gene': ['TRB'] * 2, 'epitope': ['A', 'B'],
        'mhc_a': ['HLA-A*02:01'] * 2, 'mhc_b': ['B2M'] * 2, 'mhc_class': ['MHCI'] * 2,
        'model': [0, 0]})
    raw_control = pl.DataFrame({'junction_aa': ['CAAF', 'CACF', 'CADF', 'CAEF'],
        'v_call': ['TRBV1'] * 4, 'j_call': ['TRBJ1-1'] * 4})
    monkeypatch.setattr(module, 'compute_pgen', lambda ref, control, **kw: pl.concat([
        ref.select('junction_aa'), control.select('junction_aa')]).unique().with_columns(pl.lit(.01).alias('pgen')))
    def banks(q, *args, **kw):
        q = q.with_columns(pl.lit('bin').alias('pgen_bin'))
        def bank(name):
            return q.head(1).drop('query_id').join(pl.DataFrame({
                'query_id': [f'__calibration_{name}_{i}' for i in range(3)], 'replicate': range(3)}), how='cross')
        return q, bank('normalize'), bank('test')
    monkeypatch.setattr(module, '_query_banks', banks)
    monkeypatch.setattr(module, 'iter_mock_models', lambda model, *a, **kw:
        (model.with_columns(pl.lit(i).alias('model')) for i in (1, 2, 3)))
    def components(q, model, context, **kw):
        return q.select('query_id', 'pgen_bin').join(model.select(
            'epitope', 'mhc_a', 'mhc_b', 'mhc_class', 'model'), how='cross')
    monkeypatch.setattr(runner, 'score_components', components)
    seen_models = []
    def combine(scored, norm, model, q, **kw):
        marker = model['model'][0]
        seen_models.append(marker)
        assert norm['query_id'].str.starts_with('__calibration_normalize_').all()
        assert norm['model'].unique().to_list() == [marker]
        values = []
        for row in scored.iter_rows(named=True):
            name, epitope = row['query_id'], row['epitope']
            if name == 'q':
                value = {0: (5., 7.), 1: (8., 0.), 2: (5., 9.), 3: (None, None)}[marker][epitope == 'B']
            elif name.startswith('__calibration_test_'):
                value = [(6., 1.), (4., 8.), (None, None)][int(name[-1])][epitope == 'B']
            else:
                value = None
            values.append(value)
        return scored.with_columns(pl.Series('score', values, dtype=pl.Float64),
            pl.Series('call_score', values, dtype=pl.Float64))
    def streaming_combine(*args, **kwargs):
        consume = kwargs.pop('_consume', None)
        frame = combine(*args, **kwargs)
        if consume is None:
            return frame
        for part in frame.iter_slices(1):
            consume(part)
        return None
    monkeypatch.setattr(runner, 'combine_components', streaming_combine)
    context = {'density_index': None, 'prior_background': raw_control.rename(
        {'junction_aa': 'cdr3', 'v_call': 'v', 'j_call': 'j'})}
    _, tails = module.calibrate_combined(query, reference, raw_control, species='human', locus='TRB',
        replicates=3, model_context=context)
    q = tails.filter(pl.col('query_id') == 'q')
    for epitope, p_reference, p_query in [('A', .75, .5), ('B', .5, .5), (None, .75, .5)]:
        row = q.filter(pl.col('epitope').eq_missing(epitope))
        assert row['p_reference'].item() == p_reference
        assert row['p_query'].item() == p_query
        assert row['p_both'].item() == max(p_reference, p_query)
    assert tails.filter(pl.col('query_id') == 'none')['p_both'].to_list() == [1.] * 3
    assert seen_models == [0, 1, 2, 3]
    # A whole missing null group must fail, rather than look like zero exceedances.
    import pytest
    def missing_group(*args, **kwargs):
        frame = combine(*args, **kwargs)
        return frame.filter(~(pl.col('query_id').str.starts_with('__calibration_test_') &
                              (pl.col('epitope') == 'A')))
    monkeypatch.setattr(runner, 'combine_components', missing_group)
    with pytest.raises(RuntimeError, match='every testing group'):
        module.calibrate_combined(query, reference, raw_control, species='human', locus='TRB',
            replicates=3, model_context=context)
    # Streaming mocks must reject missing cases and repeats across batches.
    for duplicate in (False, True):
        def broken_mock(*args, **kwargs):
            consume = kwargs.pop('_consume', None)
            frame = combine(*args, **kwargs)
            if consume is None:
                return frame
            consume(frame if duplicate else frame.filter(pl.col('epitope') != 'B'))
            if duplicate:
                consume(frame.filter(pl.col('query_id') == 'q').head(1))
        monkeypatch.setattr(runner, 'combine_components', broken_mock)
        with pytest.raises(RuntimeError, match='complete observed case grid'):
            module.calibrate_combined(query, reference, raw_control, species='human', locus='TRB',
                replicates=3, model_context=context)


def test_current_model_graph_matches_union_induced_weights_and_tails(monkeypatch):
    """Unused control vertices cannot join components or enlarge graph searches."""
    from vdjmatch.evalue import combined as module
    from vdjmatch.runner import combined as runner
    monkeypatch.setattr(runner, '_COMBINATION_BATCH_ROWS', 1)
    from vdjmatch import cluster
    from vdjmatch.aggregate.candidates import PMHC
    from vdjmatch.aggregate.reference_combination import graph_mixture_weights
    query = pl.DataFrame({'query_id': ['q'], 'cdr3': ['CASSMGQETQYF'],
        'v': ['TRBV1'], 'j': ['TRBJ1-1'], 'locus': ['TRB']})
    reference = pl.DataFrame({'cdr3': ['CASSLGQETQYF', 'CASSIGQDTQYF', 'CATSLGQETQYVF'],
        'v': ['TRBV1'] * 3, 'j': ['TRBJ1-1'] * 3, 'gene': ['TRB'] * 3,
        'species': ['HomoSapiens'] * 3, 'epitope': ['A', 'A', 'B'],
        'mhc_a': [None] * 3, 'mhc_b': [None] * 3, 'mhc_class': [None] * 3}).with_columns(
            pl.col('mhc_a', 'mhc_b', 'mhc_class').cast(pl.String))
    control = pl.DataFrame({'junction_aa': ['CASSIGQETQYF', 'CASSVGQETQYF',
        'CASSAGQETQYVF', 'CATSVWGQETQYF', 'CASSLGQETQYF', 'WWWWWWWW'],
        'v_call': ['TRBV1'] * 6, 'j_call': ['TRBJ1-1'] * 6, 'occurrences': [3, 4, 5, 6, 3, 2]})
    monkeypatch.setattr(module, 'compute_pgen', lambda ref, events, **kw: pl.concat([
        ref.select('junction_aa'), events.select('junction_aa')]).unique().with_columns(pl.lit(.01).alias('pgen')))
    nodes = pl.concat([reference.select('cdr3'), control.select(pl.col('junction_aa').alias('cdr3'))]).unique().sort('cdr3').with_row_index('_node')
    overlap = cluster.overlap
    pairs = overlap(nodes['cdr3'].to_list(), scope='1,0,0,1', threads=1).filter(pl.col('a_idx') < pl.col('b_idx'))
    union = pairs.join(nodes.select(pl.col('_node').cast(pl.Int64).alias('a_idx'),
        pl.col('cdr3').alias('cdr3_a')), on='a_idx').join(nodes.select(
        pl.col('_node').cast(pl.Int64).alias('b_idx'), pl.col('cdr3').alias('cdr3_b')), on='b_idx').select('cdr3_a', 'cdr3_b')
    combine = runner.combine_components
    histories = []
    for use_union in (True, False):
        weights = []; searched = []
        def combine_and_record(parts, normalizer, model, queries, **kw):
            groups = model.select(PMHC).unique().sort(PMHC).with_row_index('group_id')
            records = model.join(groups, on=PMHC, nulls_equal=True).select('group_id', 'cdr3', 'v')
            edges = union if use_union else None
            weights.append(graph_mixture_weights(records, threads=1, hamming_edges=edges)[0].select(
                'group_id', 'original_weight', 'gapblock_weight', 'stem4_weight'))
            kw['hamming_edges'] = edges
            return combine(parts, normalizer, model, queries, **kw)
        def record_overlap(sequences, **kw):
            searched.append(len(sequences))
            return overlap(sequences, **kw)
        monkeypatch.setattr(runner, 'combine_components', combine_and_record)
        monkeypatch.setattr(cluster, 'overlap', record_overlap)
        detail, tails = module.calibrate_combined(query, reference, control,
            species='human', locus='TRB', replicates=3, seed=7, threads=1)
        histories.append((detail.sort(['query_id', *PMHC]), tails.sort(['query_id', 'test_scope', *PMHC]), weights))
        if not use_union:
            assert searched and max(searched) <= reference.height
    a, b = histories
    assert a[0].equals(b[0]) and a[1].equals(b[1])
    assert len(a[2]) == len(b[2]) == 4
    assert all(x.equals(y) for x, y in zip(a[2], b[2]))
    assert any(not a[2][0].equals(mock) for mock in a[2][1:])


def test_case_budget_rejects_before_pgen_and_before_density_index(monkeypatch):
    import pytest
    import seqtree
    from types import SimpleNamespace
    from vdjmatch.evalue import combined as module
    from vdjmatch.runner import combined as runner
    queries = pl.DataFrame({'query_id': ['q'], 'cdr3': ['CASSMGQETQYF'],
        'v': ['TRBV1'], 'j': ['TRBJ1-1'], 'locus': ['TRB']})
    reference = pl.DataFrame({'cdr3': ['CASSLGQETQYF', 'CASSIGQETQYF'],
        'v': ['TRBV1'] * 2, 'j': ['TRBJ1-1'] * 2, 'gene': ['TRB'] * 2, 'epitope': ['A', 'B']})
    control = pl.DataFrame({'junction_aa': ['CASSAGQETQYF', 'CASSVGQETQYF'],
        'v_call': ['TRBV1'] * 2, 'j_call': ['TRBJ1-1'] * 2, 'occurrences': [4, 4]})
    pgen_calls = []
    def pgen(ref, events, **kwargs):
        pgen_calls.append(1)
        return pl.concat([ref.select('junction_aa'), events.select('junction_aa')]).unique().with_columns(pl.lit(.01).alias('pgen'))
    monkeypatch.setattr(module, 'compute_pgen', pgen)
    monkeypatch.setattr(seqtree, 'Index', SimpleNamespace(build=lambda *a, **kw: pytest.fail('density index built before case preflight')))
    # Minimum observed+unconditional-normalizer grid: (1+2)*2=6 cases.
    monkeypatch.setattr(runner, '_CASE_LIMIT', 5)
    with pytest.raises(ValueError, match='case budget exceeded') as error:
        module.calibrate_combined(queries, reference, control, species='human', locus='TRB', replicates=2)
    assert not pgen_calls and 'masks' not in str(error.value)
    # The first check now passes, but the actual test bank expands the grid to10.
    monkeypatch.setattr(runner, '_CASE_LIMIT', 6)
    with pytest.raises(ValueError, match='case budget exceeded'):
        module.calibrate_combined(queries, reference, control, species='human', locus='TRB', replicates=2)
    assert len(pgen_calls) == 1
