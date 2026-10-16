"""One hand-built reducer/CLI fixture; no biological calculations or native calls."""
import argparse
import json

import polars as pl
import pytest

from vdjmatch import db
from vdjmatch.db.cache import sha256
from vdjmatch.cli import reference_pairing as pairing


def test_linked_full_keys_mhc_separation_and_unlinked_marginals(tmp_path):
    aseq, bseq = 'CAAAAAAF', 'CAGGGGGF'
    rows = []
    def record(link, locus, seq, v, j, mhc='HLA-A*01:01'):
        rows.append(dict(complex_id=link, gene=locus, cdr3=seq, v=v, j=j,
            species='HomoSapiens', epitope='DEMO', mhc_a=mhc, mhc_b='B2M', mhc_class='MHCI'))
    record('p1', 'TRA', aseq, 'TRAV1*01', 'TRAJ1*01')
    record('p1', 'TRA', aseq, 'TRAV1*01', 'TRAJ1*01')  # duplicate chain observation
    record('p1', 'TRB', bseq, 'TRBV1*01', 'TRBJ1*01')
    record('p2', 'TRA', aseq, 'TRAV2*01', 'TRAJ2*01')  # same AA, distinct V/J clonotype
    record('p2', 'TRB', bseq, 'TRBV2*01', 'TRBJ2*01')
    record('p3', 'TRA', aseq, 'TRAV1*01', 'TRAJ1*01', 'HLA-A*02:01')
    record('p3', 'TRB', bseq, 'TRBV1*01', 'TRBJ1*01', 'HLA-A*02:01')
    record('0', 'TRA', 'CARRRRAF', 'TRAV3*01', 'TRAJ3*01')
    record('0', 'TRB', 'CAYYYYYF', 'TRBV3*01', 'TRBJ3*01')
    source = tmp_path/'reference.tsv'
    pl.DataFrame(rows).write_csv(source, separator='\t')
    raw = db.load(source, species='HomoSapiens')
    original = raw.rename({'gene': 'locus'}).select(pairing.KEY).unique()
    paths = []
    for locus in ['TRA', 'TRB']:
        frame = original.filter(pl.col('locus') == locus)
        # Deliberately reuse community IDs across loci; only full keys establish a link.
        frame = frame.with_columns(
            pl.when(pl.col('mhc_a') == 'HLA-A*02:01').then(pl.lit(30))
            .when(pl.col('v').str.contains('2')).then(pl.lit(11))
            .otherwise(pl.lit(10)).cast(pl.Int64).alias('community'),
            pl.lit(True).alias('geometry_available'),
            pl.lit('reference').alias('population'), pl.lit(-1).alias('replicate'))
        if locus == 'TRA':
            frame = frame.filter(pl.col('cdr3') != 'CARRRRAF')  # missing membership retained explicitly
        else:
            frame = frame.with_columns(
                pl.when(pl.col('v') == 'TRBV2*01').then(None).otherwise(pl.col('community')).alias('community'),
                (pl.col('v') != 'TRBV2*01').alias('geometry_available'))
        # A mock assignment must never leak into genuine observed-pair counts.
        mock = frame.head(1).with_columns(pl.lit('matched_control').alias('population'),
            pl.lit(0).alias('replicate'), pl.lit(999, dtype=pl.Int64).alias('community'))
        path = tmp_path/f'{locus}.communities.tsv'
        pl.concat([frame, mock]).write_csv(path, separator='\t')
        manifest = dict(reference=db.provenance(source), species='HomoSapiens', locus=locus,
            geometry='sw', reference_population='curated-benchmark', targets_sha256=None,
            tables={path.name: sha256(path)})
        path.with_name(f'{locus}.manifest.json').write_text(json.dumps(manifest))
        paths.append(path)
    parser = argparse.ArgumentParser()
    pairing.register(parser.add_subparsers())
    args = parser.parse_args(['reference-pairing', '--vdjdb', str(source),
        '--species', 'human', '--epitope', 'DEMO', '--alpha-communities', str(paths[0]),
        '--beta-communities', str(paths[1]), '--geometry', 'sw', '--output-prefix', str(tmp_path/'out')])
    assert args.func(args) == 0
    linked = pl.read_csv(tmp_path/'out.linked.tsv', separator='\t')
    assert linked.height == 3 and set(linked['complex_id']) == {'p1', 'p2', 'p3'}
    by_id = {row['complex_id']: row for row in linked.to_dicts()}
    assert by_id['p1']['community_alpha'] == 10
    assert by_id['p2']['community_alpha'] == 11
    assert by_id['p3']['community_alpha'] == 30
    assert by_id['p2']['community_beta'] is None
    assert by_id['p2']['community_status_beta'] == 'geometry_unavailable'
    pairs = pl.read_csv(tmp_path/'out.community-pairs.tsv', separator='\t')
    assert pairs['complete_complex_links'].sum() == 3
    assert pairs.filter(pl.col('mhc_a') == 'HLA-A*01:01')['distinct_paired_clonotypes'].sum() == 2
    assert 999 not in pairs['community_alpha']
    marginal = pl.read_csv(tmp_path/'out.marginals.tsv', separator='\t')
    assert marginal.height == 8
    unlinked = marginal.filter(pl.col('complete_complex_links') == 0)
    assert set(unlinked['cdr3']) == {'CARRRRAF', 'CAYYYYYF'}
    assert unlinked.filter(pl.col('locus') == 'TRA')['community_status'].item() == 'missing_community'
    coverage = pl.read_csv(tmp_path/'out.pairing-coverage.tsv', separator='\t')
    a = coverage.filter((pl.col('locus') == 'TRA') & (pl.col('mhc_a') == 'HLA-A*01:01')).row(0, named=True)
    assert a['chain_source_records'] == 4
    assert a['chain_full_clonotypes'] == 3 and a['paired_chain_full_clonotypes'] == 2
    assert a['complete_complex_links'] == 2 and a['distinct_paired_clonotypes'] == 2
    assert a['paired_chain_clonotype_fraction'] == pytest.approx(2/3)
    identity = db.provenance(source)
    manifest_path = paths[0].with_name('TRA.manifest.json')
    manifest = json.loads(manifest_path.read_text())
    for key, value in [('geometry', 'levenshtein'), ('locus', 'TRB'), ('targets_sha256', 'abc'),
            ('reference', {'input_sha256': '0'*64}), ('tables', {})]:
        manifest_path.write_text(json.dumps(manifest | {key: value}))
        with pytest.raises(ValueError):
            pairing._communities(paths[0], locus='TRA', species='HomoSapiens', geometry='sw',
                policy='curated-benchmark', reference=identity)
    manifest_path.write_text(json.dumps(manifest))
    # All-valid includes a canonical observation excluded only by shared curation.
    raw = raw.with_columns(pl.when((pl.col('cdr3') == 'CARRRRAF') | (pl.col('complex_id') == 'p2')).then(False)
        .otherwise(pl.col('reference_valid')).alias('reference_valid'))
    tables = pairing.pairing_tables(raw, *[pairing._communities(path, locus=locus,
        species='HomoSapiens', geometry='sw', policy='curated-benchmark', reference=identity)[0]
        for path, locus in zip(paths, ['TRA', 'TRB'])])
    assert tables['.marginals.tsv'].filter(pl.col('cdr3') == 'CARRRRAF').height == 1

    assert tables['.linked.tsv'].filter(pl.col('complex_id') == 'p2').height == 1
    empty = pairing.pairing_tables(raw.filter(pl.col('complex_id') == '0'),
        *[pairing._communities(path, locus=locus, species='HomoSapiens', geometry='sw',
            policy='curated-benchmark', reference=identity)[0]
            for path, locus in zip(paths, ['TRA', 'TRB'])])
    assert empty['.linked.tsv'].height == 0
    assert empty['.marginals.tsv'].height == 2
    assert empty['.pairing-coverage.tsv']['complete_complex_links'].sum() == 0
