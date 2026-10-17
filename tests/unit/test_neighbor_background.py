"""Uniform same-N graph draws retain exact-junction links between clones."""
import numpy as np
import polars as pl
import pytest


def test_uniform_draws_match_direct_hamming_oracle_and_one_native_batch(monkeypatch):
    from vdjmatch.aggregate.neighbor_background import uniform_neighbor_draws
    from vdjmatch import cluster
    pool = ['CASSF', 'CASSF', 'CASAF', 'CWWWF', 'CWWAF']
    calls, native = [], cluster.overlap
    def counted(*args, **kwargs):
        calls.append(len(args[0]))
        return native(*args, **kwargs)
    monkeypatch.setattr(cluster, 'overlap', counted)
    result = uniform_neighbor_draws(pool, [1, 3, 5], replicates=3, seed=7, threads=1)
    assert len(calls) == 1
    rng = np.random.default_rng(7)
    for replicate in range(3):
        draw = rng.choice(len(pool), size=5, replace=False)
        for n in (1, 3, 5):
            sample = [pool[i] for i in draw[:n]]
            degree = [sum(len(a) == len(b) and sum(x != y for x, y in zip(a, b)) <= 1
                for j, b in enumerate(sample) if i != j) for i, a in enumerate(sample)]
            row = result.filter((pl.col('replicate') == replicate) & (pl.col('clonotypes') == n)).row(0, named=True)
            assert row['graph_edges'] == sum(degree)//2
            assert row['connected_fraction'] == pytest.approx(np.mean(np.asarray(degree) > 0))
            assert row['pair_density'] == pytest.approx(sum(degree)/(n*(n-1))) if n > 1 else row['pair_density'] is None
    assert result.equals(uniform_neighbor_draws(pool, [1, 3, 5], replicates=3, seed=7, threads=2))
    for sizes in ([0], [6], [1.5]):
        with pytest.raises(ValueError):
            uniform_neighbor_draws(pool, sizes, replicates=1, seed=0, threads=1)


def test_structure_cli_exports_uniform_vdjdb_baseline(tmp_path):
    from vdjmatch.cli.__main__ import main
    raw = pl.DataFrame({'cdr3': ['CASSF', 'CASSF', 'CASAF', 'CWWWF'],
        'v': ['TRBV19', 'TRBV20-1', 'TRBV19', 'TRBV19'], 'j': ['TRBJ1-1']*4,
        'epitope': ['A']*4, 'species': ['HomoSapiens']*4, 'gene': ['TRB']*4})
    raw.write_csv(tmp_path/'reference.tsv', separator='\t')
    prefix = tmp_path/'fresh'
    assert main(['reference-structure', '--vdjdb', str(tmp_path/'reference.tsv'),
        '--locus', 'TRB', '--neighbor-draws', '3', '--threads', '1',
        '--output-prefix', str(prefix)]) == 0
    draws = pl.read_csv(str(prefix)+'.neighbor-draws.tsv', separator='\t')
    assert draws['graph_edges'].to_list() == [3]*3
    summary = pl.read_csv(str(prefix)+'.neighbor-background.tsv', separator='\t').row(0, named=True)
    assert summary['background'] == 'vdjdb' and summary['draws'] == 3
    assert summary['observed_clonotypes'] == 4
    assert summary['pair_density_ratio'] == 1 and summary['connected_fraction_ratio'] == 1
    assert summary['mean_degree_ratio'] == 1
    control = raw.select(pl.col('cdr3').alias('junction_aa'), pl.col('v').alias('v_call'),
        pl.col('j').alias('j_call')).with_columns(pl.lit(100).alias('occurrences'))
    control.write_csv(tmp_path/'control.tsv', separator='\t')
    controlled = tmp_path/'controlled'
    assert main(['reference-structure', '--vdjdb', str(tmp_path/'reference.tsv'),
        '--locus', 'TRB', '--control', str(tmp_path/'control.tsv'), '--neighbor-draws', '3',
        '--threads', '1', '--output-prefix', str(controlled)]) == 0
    compared = pl.read_csv(str(controlled)+'.neighbor-background.tsv', separator='\t')
    assert compared['background'].to_list() == ['vdjdb', 'control']
    assert compared['population_clonotypes'].to_list() == [4, 4]
    assert compared['pair_density_ratio'].to_list() == [1, 1]
