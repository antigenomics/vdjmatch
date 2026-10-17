"""Independent expression fixtures for retained advanced geometry, no cohort data."""
import numpy as np
import pytest
from vdjmatch.match import geometry as td


@pytest.mark.parametrize('species', ['human', 'mouse'])
@pytest.mark.parametrize('matrix', ['clipped', 'blosum62'])
def test_cdr12_context_matches_independent_aligned_expression(species, matrix):
    model = td.load_v_loops(species)
    costs = td._substitution_costs(matrix=matrix)
    context = td.precompute_v_loop_distances(species=species, matrix=matrix, v_regions='cdr12')
    alphabet = td.amino_acids()
    for locus in ('TRA', 'TRB'):
        alleles, actual = context['tables'][locus]
        # Independent aligned columns: omit CDR2.5, unknown symbols cost zero.
        loops = [model[v]['cdr1'] + model[v]['cdr2'] for v in alleles]
        lookup = {a: i for i, a in enumerate(alphabet)}
        expected = np.zeros(actual.shape, dtype=np.int32)
        for i, left in enumerate(loops):
            for j, right in enumerate(loops):
                expected[i, j] = sum(int(costs[lookup[a], lookup[b]])
                    for a, b in zip(left, right) if a in 'ACDEFGHIKLMNPQRSTVWY' and b in 'ACDEFGHIKLMNPQRSTVWY')
        np.testing.assert_array_equal(actual, expected)


def test_cdr12_dense_native_counts_and_default_preserved():
    seq = ['CASSIVGGNEQFF', 'CASSMRSTGELFF', 'CASSIVGGNEQFF']
    v = ['TRBV19', 'TRBV6-5', 'TRBV19']
    default = td.distance_matrix(seq, seq, v, v, ctrim=3, stem_depth=4)
    np.testing.assert_array_equal(default,
        td.distance_matrix(seq, seq, v, v, ctrim=3, stem_depth=4, v_regions='cdr12plus'))
    projected = {k: dict(meta, cdr25='') for k, meta in td.load_v_loops().items()}
    costs = td._substitution_costs()
    loop, ids = td._v_loop_groups(td.resolve_v_alleles(v, projected).to_list(), projected, costs)
    _, cdr3 = td.distance_matrix(seq, seq, v, v, ctrim=3, return_cdr3=True)
    expected = cdr3 + loop[ids[:, None], ids[None, :]]
    actual = td.distance_matrix(seq, seq, v, v, ctrim=3, v_regions='cdr12')
    np.testing.assert_array_equal(actual, expected)
    cutoffs = list(range(0, 24001, 200))
    for threads in (1, 4):
        count = td.total_distance_count_batch(seq, seq, v, v, cutoffs,
            ctrim=3, v_regions='cdr12', exclude_exact=True, shared_thresholds=True, threads=threads)
        oracle = ((200 * expected[:, :, None] <= np.asarray(cutoffs)) &
            (np.asarray(seq)[:, None, None] != np.asarray(seq)[None, :, None])).sum(axis=1)
        np.testing.assert_array_equal(count, oracle)


def test_v_region_validation_and_context_coherence():
    for bad in ('cdr25', '', 1):
        with pytest.raises(ValueError, match='v_regions'):
            td.precompute_v_loop_distances(v_regions=bad)
    context = td.precompute_v_loop_distances(v_regions='cdr12')
    with pytest.raises(ValueError, match='matching V regions'):
        td.distance_matrix([], [], [], [], v_regions='cdr12plus', v_loop_distances=context)


def test_recorded_trimmed_profile_needs_rounding_before_v_addition():
    """An exact floor boundary cannot be replaced by fractional weighted counts."""
    from seqtree import SubstitutionMatrix, gapblock
    from vdjmatch.match.regions import significance_weights
    length = 13
    profile = significance_weights(length)[3:-3]
    mean = sum(profile) / len(profile)
    weights = [max(1, round(100 * value / mean)) for value in profile]
    costs = td._substitution_costs()
    symbols = td.amino_acids()
    # Find a genuine quantized profile residual, independent of gene/reference data.
    site, residue, residual = next((i, aa, 2 * int(costs[symbols.index('A'), symbols.index(aa)]) * w)
        for i, w in enumerate(weights) for aa in 'RNDCEQGHILKMFPSTWYV'
        if (2 * int(costs[symbols.index('A'), symbols.index(aa)]) * w) % 200)
    q = 'A' * len(weights)
    r = q[:site] + residue + q[site + 1:]
    raw = gapblock.score_matrix([q], [r],
        matrix=SubstitutionMatrix.from_similarity((-costs).tolist()),
        position_weights_by_length={len(q): weights}, threads=1)[0, 0]
    assert raw == residual
    cutoff = 3 * (raw // 200)
    assert cutoff < 3 * raw / 200
    # Existing integer native weighting retains the residual, whereas the recorded
    # expression floors the CDR3 component before the threefold composition.
    fractional = gapblock.count_batch([q], [r], [[200 * cutoff]],
        matrix=SubstitutionMatrix.from_similarity((-3 * costs).tolist()),
        position_weights_by_length={len(q): weights}, threads=1)
    assert fractional == [[0]]
    assert int(3 * (raw // 200) <= cutoff) == 1


@pytest.mark.parametrize('paired', [False, True])
def test_cdr12_grouped_counts_match_dense_geometry(paired):
    import polars as pl
    q = pl.DataFrame({'seq0': ['CASSIVGGNEQFF'], 'v0': ['TRBV19']})
    r = pl.DataFrame({'seq0': ['CASSMRSTGELFF', 'CASSIVGGNEQFF'],
                      'v0': ['TRBV6-5', 'TRBV19']})
    if paired:
        q = q.rename({'seq0': 'seq1', 'v0': 'v1'}).with_columns(
            pl.lit('CAVRDSNYQLIW').alias('seq0'), pl.lit('TRAV12-2').alias('v0'))
        r = r.rename({'seq0': 'seq1', 'v0': 'v1'}).with_columns(
            pl.Series('seq0', ['CAVRDSNYQLVW', 'CAVRDSNYQLIW']),
            pl.Series('v0', ['TRAV14/DV4', 'TRAV12-2']))
    context = td.precompute_v_loop_distances(v_regions='cdr12')
    if paired:
        dense = td.paired_distance_matrix(q['seq0'], q['seq1'], r['seq0'], r['seq1'],
            q['v0'], q['v1'], r['v0'], r['v1'], ctrim=3, v_regions='cdr12')
    else:
        dense = td.distance_matrix(q['seq0'], r['seq0'], q['v0'], r['v0'],
            ctrim=3, v_regions='cdr12')
    cuts = [0, 6000, 12000, 24000, 48000]
    expected = (dense[:, :, None] * 200 <= np.asarray(cuts)).sum(axis=1)
    for threads in (1, 4):
        iterator = td.grouped_total_distance_count_batch(q, r, cuts, paired=paired,
            v_regions='cdr12', v_loop_distances=context, ctrim=3, threads=threads,
            reference_target_offsets=[0, 1, 2], reference_target_ids=[0, 0],
            target_peptide_ids=[0], reference_exclusive_peptide_ids=[0, 0])
        chunks = list(iterator)
        np.testing.assert_array_equal(np.asarray(chunks[0].counts)[:, 0, :], expected)
