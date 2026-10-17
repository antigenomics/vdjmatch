"""Fresh classification CLI never reopens previously calculated diagnostics."""
import json

import polars as pl


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
            '--threads', '1', '--output-prefix', str(prefix)]) == 0
        outputs.append(pl.read_csv(str(prefix)+'.profiles.tsv', separator='\t'))
        manifest = json.loads((tmp_path/(name+'.manifest.json')).read_text())
        assert manifest['source_frozen'] and manifest['reference_rows'] == 62
        assert manifest['profile_rows'] == 4 and manifest['joint_rows'] == 2
    assert calls == ['HomoSapiens', 'HomoSapiens']
    assert outputs[0].equals(outputs[1])
    assert outputs[0].filter(pl.col('level') == 'epitope')['valid_source_records'].to_list() == [30, 30]
