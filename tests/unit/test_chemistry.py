"""Independent receptor-mixture oracle, including linked-pair hit exclusion."""
import itertools
import math

import numpy as np
import polars as pl
import pytest
from vdjtools.features.physchem import load_property_table

from vdjmatch.match.chemistry import local_chemistry_evidence


@pytest.mark.parametrize('paired', [False, True])
def test_local_chemistry_matches_scalar_mixture(paired):
    table = load_property_table()
    values = dict(zip(table['amino_acid'], table['kf4']))
    sd = np.std(list(values.values()))
    def junction(core):
        return 'CAA' + core + 'AAF'
    def frame(cores, genes, query=False):
        rows = []
        for i, (core, gene) in enumerate(zip(cores, genes)):
            row = {'seq0': junction(core[0]), 'v0': gene}
            if paired:
                row.update(seq1=junction(core[1]), v1='TRBV1*01')
            if query:
                row['query_id'] = str(i)
            rows.append(row)
        return pl.DataFrame(rows)
    query = [('AAILV', 'DDDDE'), ('AAILV', 'DDDDE'), ('AAA', 'DDD')]
    ref = [('AAILV', 'DDDDE'), ('AAILV', 'DDDDE'), ('AAILV', 'DDDDE'),
           ('AVVIL', 'DDDEE'), ('DDDEE', 'AAILV')]
    control = [('AAILV', 'DDDDE'), ('AAILV', 'DDDDE'), ('DDDDE', 'AAILV'),
               ('DDDDD', 'AAAIL')]
    q = frame(query, ['TRAV1*02', 'TRAV1*01', 'TRAV1'], True)
    r = frame(ref, ['TRAV1*01', 'TRAV1*02', 'TRAV2', 'TRAV1', 'TRAV1'])
    b = frame(control, ['TRAV1', 'TRAV2', 'TRAV1', 'TRAV1'])
    axes = range(2 if paired else 1)
    def windows(cores):
        parts = [[tuple(values[a] / sd for a in core[i:i+4]) for i in range(len(core)-3)]
                 for core in cores[:len(axes)]]
        return [tuple(itertools.chain.from_iterable(t)) for t in itertools.product(*parts)]
    def density(pop, exclude):
        keys = {(tuple(x[f'seq{i}'] for i in axes), x['v0'].split('*')[0]) for x in pop.to_dicts()}
        eligible = [s for s, v in keys if v == 'TRAV1' and
                    (not exclude or s != tuple(q[f'seq{i}'][0] for i in axes))]
        Q = windows(query[0])
        total = 0
        for seqs in eligible:
            R = windows(tuple(s[3:-3] for s in seqs))
            total += sum(math.exp(-sum((a-b)**2 for a, b in zip(x, y))/2)
                         for x in Q for y in R) / (len(Q)*len(R))
        n = sum(v == 'TRAV1' for _, v in keys)
        return total / n, len(keys), n
    prior_values = []
    for exclude in (False, True):
        out = local_chemistry_evidence(q, r, b, paired=paired, exclude_exact=exclude)
        expected_r, N, n = density(r, exclude)
        expected_b, M, m = density(b, exclude)
        x = out.row(0, named=True)
        assert x['chemistry_reference_density'] == pytest.approx(expected_r, rel=1e-12)
        assert x['chemistry_control_density'] == pytest.approx(expected_b, rel=1e-12)
        assert (x['chemistry_reference_population'], x['chemistry_control_population']) == (N, M)
        assert (x['chemistry_reference_v_keys'], x['chemistry_control_v_keys']) == (n, m)
        prior_values.append(x['chemistry_log_v_ratio'])
        assert x['chemistry_log_conditional_ratio'] == pytest.approx(math.log(expected_r / expected_b))
        assert out.row(1, named=True)['chemistry_log_score'] == pytest.approx(x['chemistry_log_score'])
        assert out['chemistry_status'][2] == 'no_core_window'
    assert prior_values[0] == prior_values[1]
    with pytest.raises(ValueError, match='unique'):
        local_chemistry_evidence(pl.concat([q, q]), r, b)
    empty = local_chemistry_evidence(q.head(0), r.head(0), b.head(0), paired=paired)
    assert empty.height == 0 and 'chemistry_log_score' in empty.columns
    missing = local_chemistry_evidence(q, r, b.head(0), paired=paired)
    assert missing['chemistry_log_score'].null_count() == q.height
    assert missing['chemistry_status'][0] == 'no_control_v'
    with pytest.raises(ValueError, match='bandwidth'):
        local_chemistry_evidence(q, r, b, bandwidth_factor=0)
