"""Opt-in end-anchored per-position substitution-cost multipliers, synthetic only."""
import numpy as np
import polars as pl
import pytest
from vdjmatch.match import geometry as td
from vdjmatch.match.regions import (load_position_profile, profile_multipliers, profile_weights,
                                    position_profile_summary)


def write_profile(path, cys, fw, loci=('TRA', 'TRB')):
    """``cys``/``fw`` are multiplier lists indexed by offset from their anchor."""
    with open(path, 'w') as fh:
        fh.write('locus\tside\toffset\tmultiplier\n')
        for locus in loci:
            for side, row in (('cys', cys), ('fw', fw)):
                for offset, multiplier in enumerate(row):
                    fh.write(f'{locus}\t{side}\t{offset}\t{multiplier}\n')
    return str(path)


def ones(path, loci=('TRA', 'TRB')):
    return load_position_profile(write_profile(path, [1]*9, [1]*9, loci))


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
    q, r = frames(paired)
    cuts = list(range(0, 12001, 37))
    baseline = counts(q, r, cuts, paired, threads=threads)
    weighted = counts(q, r, cuts, paired, threads=threads, position_profile=ones(tmp_path/'ones.tsv'))
    np.testing.assert_array_equal(weighted, baseline)
    assert baseline.sum()  # the ladder actually admits pairs, so equality is not vacuous


def test_anchor_counts_meet_in_the_middle_and_converge_on_a_pooled_core(tmp_path):
    profile = load_position_profile(write_profile(tmp_path/'ramp.tsv', [.2, .5, .9, 1.4], [.1, .4, 1.1]))
    # Short: both ramps still climbing when they meet; the equidistant centre averages the sides.
    assert profile_multipliers(profile, 'TRB', 9) == [.2, .5, .9, 1.4, 1.25, 1.1, 1.1, .4, .1]
    # Long: every interior position takes its nearer side's pooled core, so no length is special.
    assert profile_multipliers(profile, 'TRB', 21) == (
        [.2, .5, .9] + [1.4]*7 + [1.25] + [1.1]*8 + [.4, .1])
    # Lengths no cohort can estimate are covered by construction, not enumeration. An even
    # length has no equidistant centre, so the two pooled cores abut without averaging.
    longest = profile_multipliers(profile, 'TRB', 126)
    assert len(longest) == 126 and set(longest[4:-4]) == {1.4, 1.1}
    assert longest[62] == 1.4 and longest[63] == 1.1


def test_profile_weights_match_the_uniform_table_and_zero_the_trimmed_flanks(tmp_path):
    assert profile_weights(ones(tmp_path/'ones.tsv'), 'TRB', [13], 3) == {13: [0]*3 + [100]*7 + [0]*3}
    # Multipliers inside the trimmed flanks never reach the native table.
    edged = load_position_profile(write_profile(tmp_path/'edge.tsv', [9, 9, 9, 1], [9, 9, 9, 1]))
    assert profile_weights(edged, 'TRB', [13], 3) == {13: [0]*3 + [100]*7 + [0]*3}


def test_a_heavier_column_moves_one_known_pair_away_and_a_lighter_one_moves_it_closer(tmp_path):
    # One Q->R substitution at junction position6 of13: clipped BLOSUM62 cost3, distance 3*cost.
    # Position6 is the equidistant centre, so both sides carry the multiplier.
    q = pl.DataFrame({'seq0': ['CASSLGQAYEQYF'], 'v0': ['TRBV19']})
    r = pl.DataFrame({'seq0': ['CASSLGRAYEQYF'], 'v0': ['TRBV19']})
    memberships = dict(reference_target_offsets=[0, 1], reference_target_ids=[0],
                       target_peptide_ids=[0], reference_exclusive_peptide_ids=[0])
    def distance(**kwargs):
        iterator = td.grouped_total_distance_count_batch(q, r, list(range(0, 20001)), ctrim=3,
            stem_depth=4, **memberships, **kwargs)
        found = np.concatenate([np.asarray(chunk.counts)[:, 0, :] for chunk in iterator])[0]
        return int(np.argmax(found > 0))
    baseline = distance()
    assert baseline == 1800  # 3 (CDR3 weight) * 3 (cost) * 200 native units
    for multiplier in (0, 2, 4):
        core = [1, 1, 1, 1, 1, 1, multiplier]
        profile = load_position_profile(write_profile(tmp_path/f'm{multiplier}.tsv', core, core))
        assert distance(position_profile=profile) == multiplier*baseline


def test_profile_rejects_one_sided_noncontiguous_and_malformed_files(tmp_path):
    q, r = frames(False)
    with pytest.raises(ValueError, match='does not declare TRB'):
        counts(q, r, [0, 5000], False, position_profile=ones(tmp_path/'tra.tsv', loci=['TRA']))
    with pytest.raises(ValueError, match='uniform'):
        counts(q, r, [0], False, position_profile=ones(tmp_path/'full.tsv'),
               position_weighting='significance-trimmed')
    def bad(name, body):
        (tmp_path/name).write_text('locus\tside\toffset\tmultiplier\n'+body)
        return tmp_path/name
    (tmp_path/'head.tsv').write_text('locus\tlength\tposition\tmultiplier\n')
    with pytest.raises(ValueError, match='header'):
        load_position_profile(tmp_path/'head.tsv')
    with pytest.raises(ValueError, match='side must be cys or fw'):
        load_position_profile(bad('side.tsv', 'TRB\tV\t0\t1\n'))
    with pytest.raises(ValueError, match='nonnegative'):
        load_position_profile(bad('negative.tsv', 'TRB\tcys\t0\t-1\n'))
    with pytest.raises(ValueError, match='contiguously'):
        load_position_profile(bad('gap.tsv', 'TRB\tcys\t0\t1\nTRB\tcys\t2\t1\nTRB\tfw\t0\t1\n'))
    with pytest.raises(ValueError, match='both the cys and fw sides'):
        load_position_profile(bad('half.tsv', 'TRB\tcys\t0\t1\n'))
    with pytest.raises(ValueError, match='duplicate'):
        load_position_profile(bad('dup.tsv', 'TRB\tcys\t0\t1\nTRB\tcys\t0\t2\n'))


def test_summary_reports_the_anchor_depths_and_multiplier_range_per_locus(tmp_path):
    profile = load_position_profile(write_profile(tmp_path/'mixed.tsv', [1, 1, 4], [1, 2], loci=['TRB']))
    summary = position_profile_summary(profile)
    assert summary['TRB']['cys_offsets'] == 3 and summary['TRB']['fw_offsets'] == 2
    assert summary['TRB']['core'] == dict(cys=4.0, fw=2.0)
    assert (summary['TRB']['min'], summary['TRB']['max']) == (1.0, 4.0)
    assert summary['TRB']['mean'] == pytest.approx((1 + 1 + 4 + 1 + 2)/5)


def test_kernel_selection_forwards_the_profile_into_the_native_batch(monkeypatch, tmp_path):
    from vdjmatch.evalue.kernel import grouped_kernel_selection
    profile = ones(tmp_path/'ones.tsv')
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
