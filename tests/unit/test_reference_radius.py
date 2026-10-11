"""Reference-derived boundaries preserve units, full identities and true pairs."""
import numpy as np
import polars as pl
import pytest

from vdjmatch.match import tcrdist


def test_kneedle_scale_ties_and_no_knee():
    distances = [*range(11), *range(30, 211, 20)]
    assert tcrdist.kneedle_radius(distances) == 10
    assert tcrdist.kneedle_radius(np.array(distances)[::-1] * 200) == 2000
    assert tcrdist.kneedle_radius([0] * 10 + [10, 20]) == 0
    for values in ([], [1, 2], [3] * 10, range(20)):
        assert tcrdist.kneedle_radius(values) is None
    for values in ([0, float('nan'), 2], [-1, 2, 3], [[1, 2], [3, 4]]):
        with pytest.raises(ValueError):
            tcrdist.kneedle_radius(values)


def test_degenerate_knee_floor_is_scale_equivariant_and_geometry_local():
    raw = np.array([0, .01, 10, 20, np.nan])
    got, floor, source = tcrdist.regularize_reference_radii(raw, [1, 2, 3])
    np.testing.assert_array_equal(got, [5, 5, 10, 20, np.nan])
    assert floor == 5 and source == 'positive_knee_median'
    scaled, _, _ = tcrdist.regularize_reference_radii(raw*200, [200, 400, 600])
    np.testing.assert_array_equal(scaled, got*200)
    got, floor, source = tcrdist.regularize_reference_radii([0, np.nan], [0, 2, 4, 6])
    np.testing.assert_array_equal(got, [2, np.nan])
    assert floor == 2 and source == 'positive_nearest_median'
    got, floor, source = tcrdist.regularize_reference_radii([0, np.nan], [0, 0])
    assert np.isnan(got).all() and floor is None and source == 'no_positive_scale'


@pytest.mark.parametrize('paired', [False, True])
def test_reference_nearest_matches_exhaustive_full_identity_oracle(paired):
    seq = ['CASSLGQAYEQYF', 'CASSLGRAYEQYF', 'CASSLGQAPAYEQYF',
           'AASSLGQAYEQYF', 'CASSLGQAYEQYF']
    frame = pl.DataFrame({'seq0': seq, 'v0': ['TRAV1-1'] * 5})
    scores = tcrdist.distance_matrix(seq, seq, frame['v0'], frame['v0'], ctrim=3)
    same = np.asarray(seq)[:, None] == np.asarray(seq)[None, :]
    if paired:
        beta = [seq[0], seq[2], seq[1], seq[3], seq[2]]
        frame = frame.with_columns(pl.Series('seq1', beta), pl.lit('TRBV19').alias('v1'))
        scores += tcrdist.distance_matrix(beta, beta, frame['v1'], frame['v1'], ctrim=3)
        same &= np.asarray(beta)[:, None] == np.asarray(beta)[None, :]
    want = np.where(same, np.inf, scores).min(axis=1)
    for threads in (1, 4):
        got = tcrdist.reference_nearest_distances(frame, paired=paired, ctrim=3, threads=threads)
        np.testing.assert_array_equal(got, want)
    lone = frame.head(1)
    assert np.isnan(tcrdist.reference_nearest_distances(lone, paired=paired)).all()


def test_reference_radius_cli_keeps_unresolved_and_no_knee_dispositions(tmp_path):
    from vdjmatch.cli.__main__ import main
    source = tmp_path/'reference.tsv'
    prefix = tmp_path/'radius'
    pl.DataFrame({'cdr3': ['CASSLGQAYEQYF', 'CASSLGRAYEQYF', 'CASSLGQAPAYEQYF'],
                  'v': ['TRBV19', 'TRBV19', 'unknown'], 'gene': ['TRB']*3,
                  'epitope': ['PEPTIDE']*3, 'species': ['HomoSapiens']*3}).write_csv(source, separator='\t')
    assert main(['reference-radius', '--vdjdb', str(source), '--locus', 'TRB', '--species', 'human',
                 '--output-prefix', str(prefix)]) == 0
    out = pl.read_csv(str(prefix)+'.radii.tsv', separator='\t')
    assert out['source_observations'].to_list() == [3]
    assert out['unavailable_observations'].to_list() == [1]
    assert out['reference_keys'].to_list() == [2]
    assert out['nearest_keys'].to_list() == [2]
    assert out['radius'].null_count() == 1
    assert out['radius_status'].to_list() == ['no_knee']
