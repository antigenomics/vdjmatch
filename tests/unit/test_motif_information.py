"""Scalar presence oracle and exact binary-information identities for motifs."""
import math

import polars as pl
import pytest


@pytest.mark.parametrize('mode,size', [('exact3', 8000), ('imgt3', 1331), ('wild4', 32000)])
def test_motif_views_reuse_tokenizer_and_preserve_key_denominators(mode, size):
    from vdjmatch.match.motifs import masked_motif_presence_counts
    from vdjmatch.evalue.epitope_rules import _tokens
    from vdjmatch.aggregate.motif_information import motif_population, motif_selection_information
    frame = pl.DataFrame({'_id': [0, 1, 2, 3], 'v0': [0, 0, 1, 1],
        'seq0': ['CAAAAAAF', 'CAAAAACF', 'CAF', 'CAAAAAAF']})
    counts, vocabulary = masked_motif_presence_counts(frame, mode=mode)
    expected = _tokens(frame, [0], mode, full_junction=True,
        wildcard_positions=(0, 1, 2, 3)).group_by('v0', 't0').len()
    lookup = {token: i for i, token in enumerate(vocabulary['t0'])}
    assert counts.shape == (2, size)
    assert counts[:, lookup['XAAA' if 'wild' in mode else 'AAA']].tolist() == [2, 1]
    assert counts.sum() == expected['len'].sum()
    for group, token, n in expected.iter_rows():
        assert counts[group, lookup[token]] == n
    records = frame.select(pl.col('seq0').alias('cdr3'), pl.lit('TRBV1').alias('v'))
    control = motif_population(records, locus='TRB', mode=mode)
    result = motif_selection_information(records, control, locus='TRB', mode=mode)
    assert set(result['motif_mode']) == {mode}
    assert result['total_reference_keys'].unique().to_list() == [3]
    assert result['motif_information_bits'].to_list() == pytest.approx([0.] * result.height)
    with pytest.raises(ValueError, match='motif mode'):
        motif_selection_information(records, control, locus='TRB')


@pytest.mark.parametrize('mode', ['exact3', 'imgt3', 'wild4', 'imgt_wild4'])
def test_reference_structure_cli_labels_motif_views(tmp_path, monkeypatch, mode):
    import json
    import vdjmatch.cluster as cluster
    from vdjmatch.cli.__main__ import main
    reference, control, prefix = (tmp_path / name for name in ('reference.tsv', 'control.tsv', 'out'))
    junctions = ['CAAAAAAAF', 'CAAAAAACF', 'CAGGGGGGF']
    pl.DataFrame({'cdr3': junctions, 'v': ['TRBV1'] * 3, 'j': ['TRBJ1-1'] * 3,
        'species': ['HomoSapiens'] * 3, 'gene': ['TRB'] * 3, 'epitope': ['DEMO'] * 3}).write_csv(reference, separator='\t')
    pl.DataFrame({'junction_aa': junctions, 'v_call': ['TRBV1'] * 3,
        'j_call': ['TRBJ1-1'] * 3}).write_csv(control, separator='\t')
    calls, native = [], cluster.overlap
    def batched(*args, **kwargs):
        calls.append(len(args[0]))
        return native(*args, **kwargs)
    monkeypatch.setattr(cluster, 'overlap', batched)
    assert main(['reference-structure', '--vdjdb', str(reference), '--control', str(control),
        '--locus', 'TRB', '--motif-mode', mode, '--output-prefix', str(prefix)]) == 0
    rows = pl.read_csv(str(prefix) + '.motifs.tsv', separator='\t')
    assert rows.height and set(rows['motif_mode']) == {mode}
    assert set(rows['motif_alphabet']) == {'imgt11' if 'imgt' in mode else 'full20'}
    assert rows['motif_information_bits'].to_list() == pytest.approx([0.] * rows.height)
    manifest = json.loads((tmp_path / 'out.manifest.json').read_text())
    assert calls == [3] and manifest['native_batches'] == 1
    assert manifest['motif_information']['mode'] == mode
    assert manifest['motif_information']['k'] == (4 if 'wild' in mode else 3)
    assert manifest['motif_information']['masks'] == ([0, 1, 2, 3] if 'wild' in mode else [])


def test_motif_modes_validate_before_count_allocation_or_reference_loading(tmp_path):
    from vdjmatch.match.motifs import masked_motif_presence_counts
    from vdjmatch.cli.__main__ import main
    frame = pl.DataFrame({'_id': [0], 'v0': [2000], 'seq0': ['CAAAAAAAF']})
    with pytest.raises(ValueError, match='unknown motif mode'):
        masked_motif_presence_counts(frame, mode='bad')
    with pytest.raises(ValueError, match='256 MiB'):
        masked_motif_presence_counts(frame, mode='wild4')
    with pytest.raises(SystemExit) as exc:
        main(['reference-structure', '--vdjdb', str(tmp_path / 'absent.tsv'), '--locus', 'TRB',
            '--motif-mode', 'exact3', '--output-prefix', str(tmp_path / 'out')])
    assert exc.value.code == 2


def test_bounded_motif_counts_match_existing_tokenizer_and_count_presence():
    from vdjmatch.match.motifs import masked_motif_presence_counts
    from vdjmatch.evalue.epitope_rules import _tokens
    frame = pl.DataFrame({'_id': list(range(7)),
        'seq0': ['CAAAAAAF', 'CAAAAAF', 'CASSF', 'CWF', 'CWRHKF', 'CASSAAAAF', 'CAAAAAAF'],
        'v0': [0, 0, 1, 1, 2, 2, 1]})
    counts, vocabulary = masked_motif_presence_counts(frame)
    expected = _tokens(frame, [0], 'imgt_wild4', full_junction=True,
        wildcard_positions=(0, 1, 2, 3)).group_by('v0', 't0').len()
    lookup = {token: i for i, token in enumerate(vocabulary['t0'])}
    assert vocabulary.height == 5324 and counts.shape == (3, 5324)
    assert counts.sum() == expected['len'].sum()
    for group, token, n in expected.iter_rows():
        assert counts[group, lookup[token]] == n
    assert counts[:, lookup['XAAA']].tolist() == [2, 1, 1]


def test_motif_selection_uses_matching_key_units_and_retains_missing_support():
    from vdjmatch.aggregate.motif_information import motif_population, motif_selection_information
    reference = pl.DataFrame({'cdr3': ['CAAAAAAF', 'CAAAAAAF', 'CAAAAACF', 'CAGGGGGF', 'CAAAAAAAF'],
        'v': ['TRBV1*01', 'TRBV1*02', 'TRBV1', 'TRBV1', 'TRBV2']})
    control = pl.DataFrame({'cdr3': ['CAAAAAAF', 'CAAAAAAF', 'CAGGGGGF', 'CATTTTTF'],
        'v': ['TRBV1*01', 'TRBV1*02', 'TRBV1', 'TRBV1']})
    background = motif_population(control, locus='TRB')
    result = motif_selection_information(reference, background, locus='TRB')
    row = result.filter((pl.col('v_gene') == 'TRBV1') & (pl.col('motif') == 'XAAA') &
        (pl.col('length_scope') == 'exact_length')).row(0, named=True)
    assert row['reference_keys'] == row['control_keys'] == 3
    assert row['reference_with_motif'] == 2 and row['control_with_motif'] == 1
    assert row['presence_bits'] == pytest.approx(1.)
    assert row['total_reference_keys'] == 4 and row['total_control_keys'] == 3
    assert row['stratum_bits'] == pytest.approx(math.log2(3/4))
    assert row['joint_presence_bits'] == pytest.approx(math.log2(1.5))
    supported = result.filter(pl.col('joint_presence_bits').is_not_null())
    assert supported['joint_presence_bits'].to_list() == pytest.approx(
        supported.select(pl.col('stratum_bits') + pl.col('presence_bits')).to_series().to_list())
    expected = 2/3 * math.log2(2) + 1/3 * math.log2(.5)
    assert row['motif_information_bits'] == pytest.approx(expected)
    assert row['information_status'] == 'available'
    # Missing control family retains the reference mass and produces no infinite information.
    missing = result.filter(pl.col('v_gene') == 'TRBV2')
    assert missing.height and missing['reference_keys'].unique().to_list() == [1]
    assert missing['motif_information_bits'].null_count() == missing.height
    assert missing['stratum_bits'].null_count() == missing.height
    assert missing['joint_presence_bits'].null_count() == missing.height
    assert not result.select(pl.col('presence_bits').is_infinite().any()).item()
    identical = motif_selection_information(control, background, locus='TRB')
    assert identical['motif_information_bits'].to_list() == pytest.approx([0.] * identical.height)
    assert identical['joint_presence_bits'].to_list() == pytest.approx([0.] * identical.height)


def test_binary_motif_depletion_and_unobserved_absence_are_distinct():
    from vdjmatch.aggregate.motif_information import motif_population, motif_selection_information
    target = pl.DataFrame({'cdr3': ['CAAAAAAF', 'CAAAAACF', 'CAGGGGGF'], 'v': ['TRBV1'] * 3})
    background = pl.DataFrame({'cdr3': ['CAAAAAAF'], 'v': ['TRBV1']})
    result = motif_selection_information(target, motif_population(background, locus='TRB'), locus='TRB')
    row = result.filter((pl.col('motif') == 'XAAA') & (pl.col('length_scope') == 'all_lengths')).row(0, named=True)
    assert row['presence_bits'] == pytest.approx(math.log2(2/3))
    assert row['motif_information_bits'] is None
    assert row['information_status'] == 'absence_not_observed_in_control'
    short = pl.DataFrame({'cdr3': ['CAF'], 'v': [None]}, schema={'cdr3': pl.String, 'v': pl.String})
    population = motif_population(short, locus='TRB')
    assert population['exact_length'][0]['keys'].item() == 1
    assert population['exact_length'][1].sum() == 0
    empty = motif_selection_information(short, population, locus='TRB')
    assert empty.is_empty()
    assert empty.write_csv().startswith('v_gene,junction_length,length_scope,motif,')
