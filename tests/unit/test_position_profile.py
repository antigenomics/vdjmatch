"""Opt-in per-position substitution-cost multipliers, synthetic only."""
import numpy as np
import polars as pl
import pytest
from vdjmatch.match import geometry as td
from vdjmatch.match.regions import load_position_profile, profile_weights, position_profile_summary


def write_profile(path, multipliers, loci=('TRA', 'TRB'), lengths=range(8, 20)):
    """``multipliers(locus, length, position) -> float`` over complete junction frames."""
    with open(path, 'w') as fh:
        fh.write('locus\tlength\tposition\tmultiplier\n')
        for locus in loci:
            for length in lengths:
                for position in range(length):
                    fh.write(f'{locus}\t{length}\t{position}\t{multipliers(locus, length, position)}\n')
    return path


def frames(paired):
    q = pl.DataFrame({'seq0': ['CASSLGQAYEQYF', 'CASSLGRAYEQYF'],
                      'v0': ['TRAV12-1', 'TRAV12-2'] if paired else ['TRBV19', 'TRBV7-9']})
    r = pl.DataFrame({'seq0': ['CASSLGQAYEQYF', 'WASSLGQAYEQYW', 'CASSLWGQAAYEQYF'],
                      'v0': ['TRAV12-2', 'TRAV12-1', 'TRAV12-1'] if paired else ['TRBV7-9', 'TRBV19', 'TRBV19']})
    if paired:
        q = q.with_columns(pl.Series('seq1', ['CASSIVGGNEQFF', 'CASSMRSTGELFF']),
                           pl.Series('v1', ['TRBV19', 'TRBV7-9']))
        r = r.with_columns(pl.Series('seq1', ['CASSIVGGNEQFF', 'WASSIVGGNEQFW', 'CASSIVGAGNEQFF']),
                           pl.Series('v1', ['TRBV7-9', 'TRBV19', 'TRBV19']))
    return q, r


MEMBERSHIPS = dict(reference_target_offsets=[0, 1, 2, 3], reference_target_ids=[0, 0, 0],
                   target_peptide_ids=[0], reference_exclusive_peptide_ids=[0, 0, 0])


def counts(q, r, cuts, paired, **kwargs):
    iterator = td.grouped_total_distance_count_batch(q, r, cuts, paired=paired, ctrim=3,
        stem_depth=4, **MEMBERSHIPS, **kwargs)
    return np.concatenate([np.asarray(chunk.counts)[:, 0, :] for chunk in iterator])


@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('threads', [1, 4])
def test_all_ones_profile_reproduces_the_unweighted_counts_bitwise(tmp_path, paired, threads):
    profile = load_position_profile(write_profile(tmp_path/'ones.tsv', lambda *_: 1))
    q, r = frames(paired)
    cuts = list(range(0, 12001, 37))
    baseline = counts(q, r, cuts, paired, threads=threads)
    weighted = counts(q, r, cuts, paired, threads=threads, position_profile=profile)
    np.testing.assert_array_equal(weighted, baseline)
    assert baseline.sum()  # the ladder actually admits pairs, so equality is not vacuous


def test_profile_weights_match_the_uniform_table_and_zero_the_trimmed_flanks(tmp_path):
    profile = load_position_profile(write_profile(tmp_path/'ones.tsv', lambda *_: 1, lengths=[13]))
    assert profile_weights(profile, 'TRB', [13], 3) == {13: [0]*3 + [100]*7 + [0]*3}
    edged = load_position_profile(write_profile(tmp_path/'edge.tsv',
        lambda locus, length, position: 9 if position < 3 else 1, lengths=[13]))
    assert profile_weights(edged, 'TRB', [13], 3) == {13: [0]*3 + [100]*7 + [0]*3}


def test_a_heavier_column_moves_one_known_pair_away_and_a_lighter_one_moves_it_closer(tmp_path):
    # One Q->R substitution at junction position6: clipped BLOSUM62 cost3, native distance 3*cost.
    q = pl.DataFrame({'seq0': ['CASSLGQAYEQYF'], 'v0': ['TRBV19']})
    r = pl.DataFrame({'seq0': ['CASSLGRAYEQYF'], 'v0': ['TRBV19']})
    memberships = dict(reference_target_offsets=[0, 1], reference_target_ids=[0],
                       target_peptide_ids=[0], reference_exclusive_peptide_ids=[0])
    def distance(**kwargs):
        ladder = list(range(0, 20001, 1))
        iterator = td.grouped_total_distance_count_batch(q, r, ladder, ctrim=3, stem_depth=4,
            **memberships, **kwargs)
        found = np.concatenate([np.asarray(chunk.counts)[:, 0, :] for chunk in iterator])[0]
        return int(np.argmax(found > 0))
    baseline = distance()
    assert baseline == 1800  # 3 (CDR3 weight) * 3 (cost) * 200 native units
    for multiplier in (0, 2, 4):
        profile = load_position_profile(write_profile(tmp_path/f'm{multiplier}.tsv',
            lambda locus, length, position, m=multiplier: m if position == 6 else 1, lengths=[13]))
        assert distance(position_profile=profile) == multiplier*baseline


def test_profile_rejects_gaps_conflicting_weights_and_malformed_rows(tmp_path):
    q, r = frames(False)
    short = load_position_profile(write_profile(tmp_path/'short.tsv', lambda *_: 1, lengths=[13]))
    with pytest.raises(ValueError, match=r'does not cover TRB junction lengths \[15\]'):
        counts(q, r, [0, 5000], False, position_profile=short)
    full = load_position_profile(write_profile(tmp_path/'full.tsv', lambda *_: 1))
    with pytest.raises(ValueError, match='uniform'):
        counts(q, r, [0], False, position_profile=full, position_weighting='significance-trimmed')
    (tmp_path/'bad-header.tsv').write_text('locus\tlength\tmultiplier\n')
    with pytest.raises(ValueError, match='header'):
        load_position_profile(tmp_path/'bad-header.tsv')
    (tmp_path/'negative.tsv').write_text('locus\tlength\tposition\tmultiplier\nTRB\t9\t0\t-1\n')
    with pytest.raises(ValueError, match='nonnegative'):
        load_position_profile(tmp_path/'negative.tsv')
    (tmp_path/'partial.tsv').write_text('locus\tlength\tposition\tmultiplier\nTRB\t9\t0\t1\n')
    with pytest.raises(ValueError, match=r'covers 1/9 positions'):
        load_position_profile(tmp_path/'partial.tsv')
    (tmp_path/'dup.tsv').write_text('locus\tlength\tposition\tmultiplier\nTRB\t9\t0\t1\nTRB\t9\t0\t2\n')
    with pytest.raises(ValueError, match='duplicate'):
        load_position_profile(tmp_path/'dup.tsv')


def test_summary_reports_the_declared_multiplier_range_per_locus(tmp_path):
    profile = load_position_profile(write_profile(tmp_path/'mixed.tsv',
        lambda locus, length, position: 4 if position == 5 else 1, loci=['TRB'], lengths=[13]))
    summary = position_profile_summary(profile)
    assert summary['TRB']['lengths'] == [13] and summary['TRB']['positions'] == 13
    assert (summary['TRB']['min'], summary['TRB']['max']) == (1.0, 4.0)
    assert summary['TRB']['mean'] == pytest.approx((12 + 4)/13)


def test_kernel_selection_forwards_the_profile_into_the_native_batch(monkeypatch, tmp_path):
    from vdjmatch.evalue.kernel import grouped_kernel_selection
    profile = load_position_profile(write_profile(tmp_path/'ones.tsv', lambda *_: 1))
    q, r = frames(False)
    q = q.rename({'seq0': 'cdr3', 'v0': '_allele'}).with_row_index('query_id')
    r = r.rename({'seq0': 'cdr3', 'v0': '_allele'}).with_columns(pl.Series('epitope', ['A', 'A', 'B']))
    seen = []
    native = td.grouped_total_distance_count_batch
    def capture(*args, **kwargs):
        seen.append(kwargs.get('position_profile'))
        return native(*args, **kwargs)
    monkeypatch.setattr(td, 'grouped_total_distance_count_batch', capture)
    groups = [({'epitope': epitope}, q, r.filter(pl.col('epitope') == epitope)) for epitope in ('A', 'B')]
    state, _ = grouped_kernel_selection(q, r, groups, query_keys=['cdr3', '_allele'],
        search_keys=['cdr3', '_allele'], paired=False, cutoffs=[0, 10, 30, 100, 300],
        radius=100, kernel_scale=90, scorer='v-bounded', ctrim=3, stem_depth=4,
        position_profile=profile)
    assert seen == [profile]
    assert state['kernel_score'].is_finite().all()


def test_default_runs_carry_no_profile():
    import inspect
    from vdjmatch.evalue.reference_rule_bank import stream_reference_kernels
    assert inspect.signature(stream_reference_kernels).parameters['position_profile'].default is None
    assert inspect.signature(td.grouped_total_distance_count_batch
                             ).parameters['position_profile'].default is None
