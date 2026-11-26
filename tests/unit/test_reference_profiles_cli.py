"""Fresh classification CLI never reopens previously calculated diagnostics."""
import json

import polars as pl
import pytest


def test_reference_profile_cli_loads_raw_and_recalculates_each_invocation(tmp_path, monkeypatch):
    from vdjmatch.cli.__main__ import main
    from vdjmatch import db
    from vdjmatch.match.engine import VdjdbIndex
    from vdjmatch.db import index

    raw = pl.DataFrame({'cdr3': ['CASSLAAAYEQYF', 'CASSLACAYEQYF', 'CASSLADAYEQYF', 'CASSLAEAYEQYF'],
        'v': ['TRBV19']*4, 'j': ['TRBJ2-7']*4, 'gene': ['TRB']*4,
        'species': ['HomoSapiens']*4, 'epitope': ['A', 'A', 'B', 'B']})
    raw = pl.concat([raw]*15+[raw.head(2).with_columns(pl.lit('CAXF').alias('cdr3'))])
    source = tmp_path/'reference.tsv'
    raw.write_csv(source, separator='\t')
    calls, load = [], db.load
    def fresh(*args, **kwargs):
        calls.append(kwargs['species'])
        return load(*args, **kwargs)
    def forbidden(*args, **kwargs):
        raise AssertionError('fresh profiles must not reopen a saved reference diagnostic bundle')
    monkeypatch.setattr(db, 'load', fresh)
    monkeypatch.setattr(index, 'load_or_build_index', forbidden)
    monkeypatch.setattr(VdjdbIndex, 'build', forbidden)
    outputs = []
    for name in ('first', 'second'):
        prefix = tmp_path/name
        assert main(['reference-profiles', '--vdjdb', str(source), '--species', 'human',
            '--threads', '1', '--motif-depth', '10000', '--output-prefix', str(prefix)]) == 0
        outputs.append(pl.read_csv(str(prefix)+'.profiles.tsv', separator='\t'))
        manifest = json.loads((tmp_path/(name+'.manifest.json')).read_text())
        discovery = pl.read_csv(str(prefix)+'.motif-discovery.tsv', separator='\t')
        assert discovery.height == 8 and discovery['rediscovery_probability'].null_count() == 8
        assert manifest['source_frozen'] and manifest['reference_rows'] == 62
        assert manifest['profile_rows'] == 4 and manifest['joint_rows'] == 2
    assert calls == ['HomoSapiens', 'HomoSapiens']
    assert outputs[0].equals(outputs[1])
    assert outputs[0].filter(pl.col('level') == 'epitope')['valid_source_records'].to_list() == [30, 30]

def test_motif_discovery_samples_distinct_keys_without_replacement():
    from vdjmatch.aggregate.fast_profiles import motif_discovery

    profile = dict(level='epitope', species='HomoSapiens', locus='TRB', epitope='X',
        mhc_a='A', mhc_b='B2M', mhc_class='MHCI', reference_v_keys=4)
    for prefix in ('', 'full20_'):
        profile.update({prefix + name: value for name, value in dict(motif_v_gene='TRBV19',
            motif='AAA', motif_supported=True, motif_reference_with_motif=2).items()})
    rows = motif_discovery(pl.DataFrame([profile]), [1, 2, 4, 5])
    assert rows['rediscovery_probability'].to_list() == [0.5, 5/6, 1., None]*2


def test_motif_depth_cli_counts_carriers_and_retains_unresolved_keys(tmp_path):
    from vdjmatch.cli.__main__ import main

    source = tmp_path/'input.tsv'
    pl.DataFrame(dict(cdr3=['CASSAAAAYEQ', 'CASSAAAAYEQ', 'CASSAAAAFEQ', 'CASSAAAAYEQ'],
        v=['TRBV19*01', 'TRBV19*02', 'TRBV19', None], j=['TRBJ2']*4,
        gene=['TRB']*4, species=['HomoSapiens']*4, epitope=['X']*4)).write_csv(source, separator='\t')
    prefix = tmp_path/'results'
    assert main(['motif-depth', '--vdjdb', str(source), '--depth', '1', '--depth', '3',
        '--depth', '4', '--output-prefix', str(prefix)]) == 0
    result = pl.read_csv(str(prefix)+'.motif-depth.tsv', separator='\t')
    assert result['junction_v_keys'].to_list() == [3]*3
    assert result['unresolved_v_keys'].to_list() == [1]*3
    assert result['observed_motifs'].to_list() == [2]*3
    assert result['guaranteed_catalog_depth'].to_list() == [2]*3
    assert result['expected_90pct_depth'].to_list() == [2]*3
    assert result['expected_observed_motifs'][:2].to_list() == pytest.approx([4/3, 2])
    assert result['expected_observed_motifs'][2] is None
