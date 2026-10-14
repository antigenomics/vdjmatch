"""Ordered class geometry preserves matrix units and raw receptor identity."""
import numpy as np
import polars as pl
import pytest

from vdjmatch.evalue.epitope_rules import IMGT_CLASSES
from vdjmatch.match.tcrdist import (
    _cdr3_options, _substitution_costs, amino_acids, distance_matrix,
    distance_units, load_v_loops, paired_total_distance_mass_batch,
    reference_nearest_distances, total_distance_count_batch,
)


def test_class_matrix_and_v_loops_share_geometry():
    alphabet = amino_acids()
    full = _substitution_costs()
    coarse = _substitution_costs(alphabet='imgt')
    for i, a in enumerate(alphabet):
        ca = IMGT_CLASSES.get(a, a)
        ai = [k for k, x in enumerate(alphabet) if IMGT_CLASSES.get(x, x) == ca]
        for j, b in enumerate(alphabet):
            cb = IMGT_CLASSES.get(b, b)
            bi = [k for k, x in enumerate(alphabet) if IMGT_CLASSES.get(x, x) == cb]
            assert coarse[i, j] == (0 if ca == cb else full[np.ix_(ai, bi)].min())
    np.testing.assert_array_equal(coarse, coarse.T)
    assert np.all(coarse <= full)
    scale, _ = distance_units(coarse)
    options = _cdr3_options(costs=coarse)
    assert options['gap_open'] == options['gap_extend'] == 2 * scale

    model = load_v_loops('human')
    a, b = 'TRBV19*01', 'TRBV20-1*01'
    loops = lambda v: ''.join(model[v][c] for c in ('cdr1', 'cdr2', 'cdr25'))
    expected = sum(coarse[alphabet.index(x), alphabet.index(y)]
                   for x, y in zip(loops(a), loops(b)) if x in IMGT_CLASSES and y in IMGT_CLASSES)
    junction = 'CASSLTQETQYF'
    observed = distance_matrix([junction], [junction], [a], [b], alphabet='imgt')
    assert observed[0, 0] == expected
    np.testing.assert_array_equal(_substitution_costs(alphabet='full'), full)
    with pytest.raises(ValueError):
        _substitution_costs(alphabet='unknown')
    frame = pl.DataFrame({'seq0': ['CASSLTQETQYF', 'CATTLTQETQYF'], 'v0': [a, a]})
    np.testing.assert_array_equal(reference_nearest_distances(frame, alphabet='imgt'), [0, 0])
    assert np.all(reference_nearest_distances(frame) > 0)


@pytest.mark.parametrize('paired', [False, True])
def test_class_zero_distance_excludes_only_raw_exact_receptor(paired):
    q, alias = 'CASSLTQETQYF', 'CATTLTQETQYF'
    assert q != alias
    if paired:
        count, _ = paired_total_distance_mass_batch(
            [q], [q], [q, alias], [q, q], ['TRAV1-1*01'], ['TRBV19*01'],
            ['TRAV1-1*01'] * 2, ['TRBV19*01'] * 2, [[0]],
            alphabet='imgt', exclude_exact=True)
        # Marginal zero-distance matches from two different references do not
        # constitute a linked match when each reference's other chain differs.
        wrong = 'CASSFFQETQYF'
        unlinked, _ = paired_total_distance_mass_batch(
            [q], [q], [alias, wrong], [wrong, alias], ['TRAV1-1*01'], ['TRBV19*01'],
            ['TRAV1-1*01'] * 2, ['TRBV19*01'] * 2, [[0]], alphabet='imgt')
        assert unlinked == [[0]]
    else:
        count = total_distance_count_batch([q], [q, alias], ['TRBV19*01'],
            ['TRBV19*01'] * 2, [[0]], alphabet='imgt', exclude_exact=True)
    assert count == [[1]]


def test_class_cli_retains_raw_exclusion_and_matrix_recipe(tmp_path):
    import json
    from vdjmatch.cli.__main__ import main
    q, alias = 'CASSLTQETQYF', 'CATTLTQETQYF'
    meta = dict(species='human', epitope='E', mhc_a='HLA-A*02', mhc_b='B2M', mhc_class='MHCI')
    sample, ref = tmp_path/'q.tsv', tmp_path/'r.tsv'
    pl.DataFrame([dict(junction_aa=q, v_call='TRBV19*01', locus='TRB', **meta)]).write_csv(sample, separator='\t')
    pl.DataFrame([dict(cdr3=s, v='TRBV19*01', gene='TRB', **(meta | dict(species='HomoSapiens', epitope=e)))
                  for s, e in [(q, 'E'), (alias, 'E'), ('CASSFTQETQYF', 'OTHER')]]).write_csv(ref, separator='\t')
    args = ['tcrdist-neighbours', str(sample), '--vdjdb', str(ref), '--locus', 'TRB',
            '--targets-from-sample', '--exclude-exact', '--unbounded-nearest', '--nearest-rank',
            '--background', 'vdjdb-other', '--neighbour-weighting', 'local-rank']
    for name, extra in [('full', []), ('imgt', ['--alphabet', 'imgt']),
                        ('parallel', ['--alphabet', 'imgt', '--threads', '2'])]:
        assert main(args + extra + ['--output-prefix', str(tmp_path/name)]) == 0
    for suffix in ('queries', 'candidates', 'enrichment'):
        a = pl.read_csv(tmp_path/f'imgt.{suffix}.tsv', separator='\t')
        b = pl.read_csv(tmp_path/f'parallel.{suffix}.tsv', separator='\t')
        assert a.equals(b)
    full = pl.read_csv(tmp_path/'full.candidates.tsv', separator='\t')
    coarse = pl.read_csv(tmp_path/'imgt.candidates.tsv', separator='\t')
    assert full['distance'][0] > 0 and coarse['distance'][0] == 0
    e = pl.read_csv(tmp_path/'imgt.enrichment.tsv', separator='\t')
    assert e['reference_population'][0] == 1
    manifest = json.loads((tmp_path/'imgt.manifest.json').read_text())
    scale = manifest['distance_scale']
    assert scale['alphabet'] == 'imgt'
    assert scale['radius'] == scale['radius_factor'] * scale['matrix_unit']
    assert scale['gap_cost'] == 3 * scale['matrix_unit']
