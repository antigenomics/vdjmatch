"""Reference-derived boundaries preserve units, full identities and true pairs."""
import numpy as np
import polars as pl
import pytest

from vdjmatch.match import geometry


def test_kneedle_scale_ties_and_no_knee():
    distances = [*range(11), *range(30, 211, 20)]
    assert geometry.kneedle_radius(distances) == 10
    assert geometry.kneedle_radius(np.array(distances)[::-1] * 200) == 2000
    assert geometry.kneedle_radius([0] * 10 + [10, 20]) == 0
    for values in ([], [1, 2], [3] * 10, range(20)):
        assert geometry.kneedle_radius(values) is None
    for values in ([0, float('nan'), 2], [-1, 2, 3], [[1, 2], [3, 4]]):
        with pytest.raises(ValueError):
            geometry.kneedle_radius(values)


def test_reference_structure_uses_matrix_units_and_all_keys():
    distances=np.array([0.,2.,4.,20.,100.,np.nan])
    coverage,shape=geometry.reference_structure(distances,radius_factor=5)
    assert coverage==4/6 and 0<=shape<=1
    costs=geometry._substitution_costs()
    scaled=geometry.reference_structure(distances*7,radius_factor=5,costs=costs*7)
    np.testing.assert_allclose(scaled,[coverage,shape],rtol=0,atol=1e-15)
    assert geometry.reference_structure([np.nan])==(0.,0.)


def test_kneedle_uses_strongest_confirmed_bend_not_minor_leading_bump():
    curve=[0]*5+[15,15]+list(range(16,36))+[70,100,140]
    assert geometry.kneedle_radius(curve)==35
    assert geometry.kneedle_radius(np.array(curve)*200)==7000


def test_degenerate_knee_floor_is_scale_equivariant_and_geometry_local():
    raw = np.array([0, .01, 10, 20, np.nan])
    got, floor, source = geometry.regularize_reference_radii(raw, [1, 2, 3])
    np.testing.assert_array_equal(got, [5, 5, 10, 20, np.nan])
    assert floor == 5 and source == 'positive_knee_median'
    scaled, _, _ = geometry.regularize_reference_radii(raw*200, [200, 400, 600])
    np.testing.assert_array_equal(scaled, got*200)
    got, floor, source = geometry.regularize_reference_radii([0, np.nan], [0, 2, 4, 6])
    np.testing.assert_array_equal(got, [2, np.nan])
    assert floor == 2 and source == 'positive_nearest_median'
    got, floor, source = geometry.regularize_reference_radii([0, np.nan], [0, 0])
    assert np.isnan(got).all() and floor is None and source == 'no_positive_scale'


@pytest.mark.parametrize('paired', [False, True])
def test_reference_nearest_matches_exhaustive_full_identity_oracle(paired):
    seq = ['CASSLGQAYEQYF', 'CASSLGRAYEQYF', 'CASSLGQAPAYEQYF',
           'AASSLGQAYEQYF', 'CASSLGQAYEQYF']
    frame = pl.DataFrame({'seq0': seq, 'v0': ['TRAV1-1'] * 5})
    scores = geometry.distance_matrix(seq, seq, frame['v0'], frame['v0'], ctrim=3)
    same = np.asarray(seq)[:, None] == np.asarray(seq)[None, :]
    if paired:
        beta = [seq[0], seq[2], seq[1], seq[3], seq[2]]
        frame = frame.with_columns(pl.Series('seq1', beta), pl.lit('TRBV19').alias('v1'))
        scores += geometry.distance_matrix(beta, beta, frame['v1'], frame['v1'], ctrim=3)
        same &= np.asarray(beta)[:, None] == np.asarray(beta)[None, :]
    want = np.where(same, np.inf, scores).min(axis=1)
    for threads in (1, 4):
        got = geometry.reference_nearest_distances(frame, paired=paired, ctrim=3, threads=threads)
        np.testing.assert_array_equal(got, want)
    lone = frame.head(1)
    assert np.isnan(geometry.reference_nearest_distances(lone, paired=paired)).all()


@pytest.mark.parametrize('matrix', ['clipped', 'blosum62'])
def test_reference_radius_cli_keeps_unresolved_and_no_knee_dispositions(tmp_path, matrix):
    from vdjmatch.cli.__main__ import main
    source = tmp_path/'reference.tsv'
    prefix = tmp_path/'radius'
    pl.DataFrame({'cdr3': ['CASSLGQAYEQYF', 'CASSLGRAYEQYF', 'CASSLGQAPAYEQYF'],
                  'v': ['TRBV19', 'TRBV19', 'unknown'], 'gene': ['TRB']*3,
                  'epitope': ['PEPTIDE']*3, 'species': ['HomoSapiens']*3}).write_csv(source, separator='\t')
    assert main(['reference-radius', '--vdjdb', str(source), '--locus', 'TRB', '--species', 'human',
                 '--matrix', matrix, '--output-prefix', str(prefix)]) == 0
    out = pl.read_csv(str(prefix)+'.radii.tsv', separator='\t')
    assert out['source_observations'].to_list() == [3]
    assert out['unavailable_observations'].to_list() == [1]
    assert out['reference_keys'].to_list() == [2]
    assert out['nearest_keys'].to_list() == [2]
    assert out['radius'].null_count() == 1
    assert out['radius_status'].to_list() == ['no_knee']


def test_reference_structure_cli_profiles_compatible_union(tmp_path):
    from vdjmatch.cli.__main__ import main
    source=tmp_path/'reference.tsv';sample=tmp_path/'queries.tsv';prefix=tmp_path/'profile'
    pl.DataFrame({'cdr3':['CASSLGQAYEQYF','CASSLGRAYEQYF'], 'v':['TRBV19']*2,
        'gene':['TRB']*2, 'epitope':['PEPTIDE']*2, 'species':['HomoSapiens']*2,
        'mhc_a':['HLA-A*02:01','HLA-A*02:02'],'mhc_b':['B2M']*2,'mhc_class':['MHCI']*2
    }).write_csv(source,separator='\t')
    pl.DataFrame({'query_id':['q1','q2'],'species':['human']*2,
        'epitope':['PEPTIDE','ABSENT'],'mhc_a':['HLA-A*02']*2,
        'mhc_b':['B2M']*2,'mhc_class':['MHCI']*2}).write_csv(sample,separator='\t')
    assert main(['reference-radius','--vdjdb',str(source),'--species','human','--locus','TRB',
        '--targets-from-sample',str(sample),'--output-prefix',str(prefix)])==0
    rows=pl.read_csv(str(prefix)+'.radii.tsv',separator='\t').sort('epitope')
    assert rows['reference_keys'].to_list()==[0,2]
    assert rows['source_observations'].to_list()==[0,2]
    assert rows['reference_coverage'].to_list()==[0.,1.]


def test_reference_radius_cli_uses_selected_matrix_units(tmp_path, monkeypatch):
    import json
    from vdjmatch.cli.__main__ import main
    source, prefix = tmp_path/'reference.tsv', tmp_path/'radius'
    pl.DataFrame({'cdr3': ['CASSLGQAYEQYF', 'CASSLGRAYEQYF', 'CASSIGQAYEQYF'],
        'v': ['TRBV19']*3, 'gene': ['TRB']*3, 'epitope': ['PEPTIDE']*3,
        'species': ['HomoSapiens']*3}).write_csv(source, separator='\t')
    scale = geometry.distance_units(geometry._substitution_costs(matrix='blosum62'))[0]
    def nearest(frame, **options):
        assert options['matrix'] == 'blosum62'
        assert frame.height == 3
        return np.array([10., 20., 30.])*scale
    monkeypatch.setattr(geometry, 'reference_nearest_distances', nearest)
    monkeypatch.setattr(geometry, 'kneedle_radius', lambda _: 10.*scale)
    assert main(['reference-radius', '--vdjdb', str(source), '--locus', 'TRB',
        '--matrix', 'blosum62', '--output-prefix', str(prefix)]) == 0
    row = pl.read_csv(str(prefix)+'.radii.tsv', separator='\t').row(0, named=True)
    assert row['reference_coverage'] == pytest.approx(2/3)
    assert row['knee_radius_ratio'] == pytest.approx(10/32.5)
    manifest = json.loads(prefix.with_suffix('.manifest.json').read_text())
    assert manifest['structure_features']['matrix_scale'] == scale
    assert manifest['options']['matrix'] == 'blosum62'
    assert 'blosum62' in manifest['geometry']
