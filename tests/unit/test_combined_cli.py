"""One-control annotation preserves raw rows and runs the complete score on demand."""
import polars as pl
import pytest


def test_control_eligibility_reports_cross_chain_rows_and_keeps_shared_genes(tmp_path):
    from vdjmatch.cli.combined import _raw_control
    path = tmp_path / 'control.tsv'
    pl.DataFrame({'junction_aa': ['CAVRDSNYQLIW'] * 4,
        'v_call': ['TRAV14/DV4', 'TRAV9-2', 'TRBV1,TRBV2', 'TRBV1'],
        'j_call': ['TRAJ33', 'TRBJ2-7', 'TRBJ1-1', 'TRBJ1-1'],
        'locus': ['TRA', 'TRA', 'TRB', 'TRA']}).write_csv(path, separator='\t')
    eligible, counts = _raw_control(path)
    assert eligible['v_call'].to_list() == ['TRAV14/DV4']
    assert counts['control_source_rows'] == 4
    assert counts['control_excluded_invalid_gene_rows'] == 3
    assert counts['control_valid_rows'] == 1


def inputs(tmp_path):
    reference = tmp_path / 'reference.tsv'
    pl.DataFrame({'cdr3': ['CASSLGQETQYF', 'CASSIGQETQYF', 'CATSLGQETQYF'],
        'v': ['TRBV1'] * 3, 'j': ['TRBJ1-1'] * 3, 'gene': ['TRB'] * 3,
        'epitope': ['A', 'A', 'KLGGALQAK'], 'evidence_type': ['binding'] * 3,
        'mhc_a': ['HLA-A*02:01'] * 3, 'mhc_b': ['B2M'] * 3, 'mhc_class': ['MHCI'] * 3}).write_csv(reference, separator='\t')
    sample = tmp_path / 'sample.tsv'
    pl.DataFrame({'query_id': ['near', 'far', 'invalid', 'other-species'],
        'junction_aa': ['CASSMGQETQYF', 'W' * 80, 'CA*F', 'CASSMGQETQYF'],
        'v_call': ['TRBV1'] * 3 + ['TRAV1-2'], 'j_call': ['TRBJ1-1'] * 3 + ['TRAJ33'],
        'locus': ['TRB'] * 3 + ['TRA'],
        'species': ['human', 'human', 'human', 'mouse']}).write_csv(sample, separator='\t')
    control = tmp_path / 'control.tsv'
    pl.DataFrame({'junction_aa': ['CASSAGQETQYF', 'CASSVGQETQYF', 'CASSTGQETQYF', 'CA*F'],
        'v_call': ['TRBV1'] * 3 + [None], 'j_call': ['TRBJ1-1'] * 3 + [None],
        'occurrences': [3, 4, 5, 1], 'productive': ['T', 'T', 'T', 'F']}).write_csv(control, separator='\t')
    return sample, reference, control


@pytest.mark.parametrize('missing_reference_chain', [False, True])
@pytest.mark.parametrize('details', [False, True])
def test_one_raw_control_runs_actual_native_sum_without_json(tmp_path, monkeypatch, missing_reference_chain, details):
    from vdjmatch.evalue import combined
    from vdjmatch.cli.__main__ import main
    sample, reference, control = inputs(tmp_path)
    pl.read_csv(sample, separator='\t').with_columns(pl.lit('A').alias('epitope'),
        pl.lit('HLA-A*02:01').alias('mhc_a'),pl.lit('B2M').alias('mhc_b'),
        pl.lit('MHCI').alias('mhc_class')).write_csv(sample, separator='\t')
    if missing_reference_chain:
        pl.read_csv(sample, separator='\t').with_columns(pl.lit('human').alias('species')).write_csv(sample, separator='\t')
    def constant_pgen(reference, events, **kwargs):
        return pl.concat([reference.select('junction_aa'), events.select('junction_aa')]).unique().with_columns(
            pl.lit(0.01).alias('pgen'))
    monkeypatch.setattr(combined, 'compute_pgen', constant_pgen)
    prefix = tmp_path / 'result'
    assert main(['annotate', str(sample), '--vdjdb', str(reference), '--control', str(control),
        '--locus', 'auto' if missing_reference_chain else 'TRB', '--mock-replicates', '2', '--seed', '7',
        '-o', str(prefix)] + (['--details'] if details else [])) == 0
    queries = pl.read_csv(str(prefix) + '.queries.tsv', separator='\t')
    assert queries.height == 4
    assert queries['query_id'].to_list() == ['near', 'far', 'invalid', 'other-species']
    assert queries['calibration_status'].to_list()[2:] == ['invalid_sequence', 'no_reference' if missing_reference_chain else 'species_mismatch']
    assert queries['calibration_query_id'].n_unique() == 4
    database = pl.read_csv(str(prefix) + '.database.tsv', separator='\t')
    epitopes = pl.read_csv(str(prefix) + '.epitopes.tsv', separator='\t')
    calibration = pl.read_csv(str(prefix) + '.calibration.tsv', separator='\t')
    assert database.height == 4 and epitopes.height == (4 if details else 2)
    assert set(epitopes['epitope']) == ({'A', 'KLGGALQAK'} if details else {'A'})
    if not details:
        assert epitopes['declared_target'].all() and not epitopes['reported_specificity'].any()
        assert epitopes['p_any'].min() >= 1/3
    assert database['p_both'].null_count() == 2
    assert database.filter(pl.col('calibration_query_id') == queries['calibration_query_id'][1])['p_both'].item() == 1
    assert epitopes.select((pl.col('p_both') == pl.max_horizontal('p_reference', 'p_query')).all()).item()
    assert calibration['mock_replicates'].unique().to_list() == [2]
    assert calibration['reference_event_sampler'].unique().to_list() == ['uniform-event-slots-v1']
    assert calibration['control_productive_rows'].unique().to_list() == [3]
    assert calibration['reference_species_policy'].unique().to_list() == ['caller_declared']
    assert calibration['seqtree_native_sha256'].str.len_chars().unique().to_list() == [64]
    assert not list(tmp_path.glob('result*.json'))


def test_declared_targets_expand_compatible_groups_without_pooling():
    from vdjmatch.cli.combined import _detail_targets
    raw = pl.DataFrame({'calibration_query_id':['q','q','other'], 'epitope':['E','F','E'],
        'mhc_a':['HLA-A*02:01']*3,'mhc_b':['B2M']*3,'mhc_class':['MHCI']*3})
    reference = pl.DataFrame({'epitope':['E']*4,'mhc_a':['HLA-A*02','HLA-A*02:01','HLA-A*02:02',None],
        'mhc_b':['B2M']*4,'mhc_class':['MHCI']*4})
    actual = _detail_targets(raw,pl.DataFrame({'query_id':['q']}),reference)
    assert actual.height == 3 and actual['query_id'].unique().to_list() == ['q']
    assert set(actual.filter(pl.col('epitope')=='E')['mhc_a']) == {'HLA-A*02','HLA-A*02:01'}
    assert actual.filter(pl.col('epitope')=='F')['mhc_a'].item() == 'HLA-A*02:01'
    assert reference.height == 4


def test_preflight_rejects_same_control_and_unsupported_geometry(tmp_path):
    from vdjmatch.cli.__main__ import main
    sample, reference, control = inputs(tmp_path)
    base = ['annotate', str(sample), '--vdjdb', str(reference), '--locus', 'TRB', '-o', str(tmp_path / 'out')]
    for extra in (['--control', str(sample)], ['--control', str(control), '--stem-depth', '0'],
                  ['--control', str(control), '--targets-from-sample']):
        with pytest.raises(SystemExit):
            main(base + extra)
    assert not list(tmp_path.glob('out.*'))


@pytest.mark.parametrize('fpr', [None, .025])
def test_mixed_wide_pairs_keep_missing_partner_and_global_identity(tmp_path, monkeypatch, fpr):
    from vdjmatch.evalue import combined
    from vdjmatch.cli.__main__ import main
    inputs_dir = tmp_path / 'inputs'; inputs_dir.mkdir()
    pl.DataFrame({'query_id': ['complete', 'partial'],
        'junction_aa_alpha': ['CAVRDSNYQLIW', None], 'junction_aa_beta': ['CASSLGQETQYF', 'CASSMGQETQYF'],
        'v_alpha': ['TRAV1-2', None], 'v_beta': ['TRBV1'] * 2,
        'j_alpha': ['TRAJ33', None], 'j_beta': ['TRBJ1-1'] * 2}).write_csv(inputs_dir / 'wide.tsv', separator='\t')
    pl.DataFrame({'query_id': ['complete'], 'junction_aa': ['CASSIGQETQYF'],
        'v_call': ['TRBV1'], 'j_call': ['TRBJ1-1'], 'locus': ['TRB']}).write_csv(inputs_dir / 'single.tsv', separator='\t')
    reference = tmp_path / 'reference.tsv'
    pl.DataFrame({'cdr3': ['CAVRDSNYQLIW', 'CASSLGQETQYF'], 'gene': ['TRA', 'TRB'],
        'v': ['TRAV1-2', 'TRBV1'], 'j': ['TRAJ33', 'TRBJ1-1'], 'epitope': ['E'] * 2}).write_csv(reference, separator='\t')
    control = tmp_path / 'control.tsv'
    pl.DataFrame({'junction_aa': ['CAVRDSNYQLIW', 'CASSLGQETQYF'], 'locus': ['TRA', 'TRB'],
        'v_call': ['TRAV1-2', 'TRBV1'], 'j_call': ['TRAJ33', 'TRBJ1-1'], 'occurrences': [4, 4]}).write_csv(control, separator='\t')
    calls = []
    def calibrate(queries, reference, raw, **kw):
        assert kw['report_p'] == (.01 if fpr is None else fpr)
        calls.append((kw['locus'], queries))
        detail = queries.select('query_id').with_columns(pl.lit('E').alias('epitope'),
            *[pl.lit(None, dtype=pl.String).alias(c) for c in ('mhc_a', 'mhc_b', 'mhc_class')],
            pl.lit(1.).alias('score'))
        tail = detail.with_columns(pl.lit('epitope').alias('test_scope'))
        database = detail.with_columns(*[pl.lit(None, dtype=pl.String).alias(c)
            for c in ('epitope', 'mhc_a', 'mhc_b', 'mhc_class')], pl.lit('database').alias('test_scope'))
        tails = pl.concat([tail, database]).with_columns(*[pl.lit(.5).alias(c) for c in ('p_reference', 'p_query', 'p_both')])
        return detail, tails
    monkeypatch.setattr(combined, 'calibrate_combined', calibrate)
    prefix = tmp_path / 'mixed'
    args = ['annotate', str(inputs_dir), '--vdjdb', str(reference), '--control', str(control), '-o', str(prefix)]
    assert main(args + ([] if fpr is None else ['--fpr', str(fpr)])) == 0
    assert pl.read_csv(str(prefix) + '.calibration.tsv', separator='\t')['report_p'].unique().to_list() == [.01 if fpr is None else fpr]
    assert [(locus, q.height) for locus, q in calls] == [('TRB', 2), ('paired', 1)]
    ids = [i for _, q in calls for i in q['query_id']]
    assert len(ids) == len(set(ids)) == 3
    query = pl.read_csv(str(prefix) + '.queries.tsv', separator='\t')
    assert query.height == 5 and query['calibration_query_id'].n_unique() == 4
    assert query['calibration_pair_status'].to_list().count('complete_pair') == 2
    assert query['calibration_pair_status'].to_list().count('incomplete_pair_marginal') == 1
    assert query['calibration_status'].to_list().count('missing_sequence') == 1
    assert pl.read_csv(str(prefix) + '.database.tsv', separator='\t').height == 4
    assert not list(tmp_path.glob('mixed*.json'))


@pytest.mark.parametrize('fpr', ['0', '1', '-.01', 'nan', 'inf'])
def test_invalid_reporting_threshold_fails_before_inputs(fpr):
    from vdjmatch.cli.__main__ import main
    with pytest.raises(SystemExit):
        main(['annotate', 'missing.tsv', '--control', 'missing-control.tsv', '--fpr=' + fpr])


def test_reporting_threshold_is_not_silently_ignored_without_calibration(monkeypatch):
    from vdjmatch.cli.__main__ import main
    from vdjmatch.runner import annotation
    calls = []
    monkeypatch.setattr(annotation, 'annotate_repertoires', lambda *args, **kw: calls.append(kw))
    assert main(['annotate', 'sample.tsv']) == 0
    assert 'fpr' not in calls[0]
    with pytest.raises(SystemExit):
        main(['annotate', 'sample.tsv', '--fpr', '.025'])
    assert len(calls) == 1


def test_control_filter_counts_and_reference_species_aliases(tmp_path):
    from vdjmatch.cli.combined import _raw_control
    sample, reference, control = inputs(tmp_path)
    raw = pl.read_csv(control, separator='\t').drop('productive')
    raw.with_columns(pl.col('junction_aa').str.to_lowercase()).write_csv(control, separator='\t')
    filtered, counts = _raw_control(control)
    assert filtered['junction_aa'].to_list() == ['CASSAGQETQYF', 'CASSVGQETQYF', 'CASSTGQETQYF']
    assert counts == dict(control_source_rows=4, control_productive_rows=4,
        control_excluded_nonproductive_rows=0, control_excluded_invalid_sequence_rows=1,
        control_excluded_invalid_gene_rows=0, control_valid_rows=3)
    # Use real Pgen and the actual native score: control covers every tested Pgen stratum.
    ref = pl.read_csv(reference, separator='\t').with_columns(pl.lit('human').alias('species'))
    ref.write_csv(reference, separator='\t')
    seqs = pl.concat([ref.select(pl.col('cdr3').alias('junction_aa')),
        pl.DataFrame({'junction_aa': ['CASSMGQETQYF', 'W' * 80]})])
    seqs.with_columns(pl.lit('TRBV1').alias('v_call'), pl.lit('TRBJ1-1').alias('j_call'),
        pl.lit(10).alias('occurrences')).write_csv(control, separator='\t')
    from vdjmatch.cli.__main__ import main
    prefix = tmp_path / 'real-pgen'
    assert main(['annotate', str(sample), '--vdjdb', str(reference), '--control', str(control),
        '--locus', 'TRB', '--mock-replicates', '2', '-o', str(prefix)]) == 0
    metadata = pl.read_csv(str(prefix) + '.calibration.tsv', separator='\t')
    assert metadata['pgen_model_source'].item() == 'olga'
    assert metadata['reference_species_policy'].item() == 'source_annotated'
    assert metadata['reference_excluded_species_rows'].item() == 0
    assert not list(tmp_path.glob('real-pgen*.json'))


def test_invalid_queries_survive_and_valid_control_requires_genes(tmp_path):
    from vdjmatch.cli.__main__ import main
    sample, reference, control = inputs(tmp_path)
    pl.read_csv(sample, separator='\t').with_columns(pl.lit('CA*F').alias('junction_aa')).write_csv(sample, separator='\t')
    prefix = tmp_path / 'invalid'
    args = ['annotate', str(sample), '--vdjdb', str(reference), '--control', str(control),
        '--locus', 'TRB', '--mock-replicates', '2', '-o', str(prefix)]
    assert main(args) == 0
    query = pl.read_csv(str(prefix) + '.queries.tsv', separator='\t')
    assert query.height == 4 and query['calibration_status'].to_list() == ['invalid_sequence'] * 4
    assert pl.read_csv(str(prefix) + '.database.tsv', separator='\t')['p_both'].null_count() == 4
    assert pl.read_csv(str(prefix) + '.epitopes.tsv', separator='\t').height == 0
    pl.read_csv(control, separator='\t').drop('j_call').write_csv(control, separator='\t')
    with pytest.raises(SystemExit):
        main(args[:-1] + [str(tmp_path / 'missing-gene')])
    assert not list(tmp_path.glob('missing-gene.*'))
