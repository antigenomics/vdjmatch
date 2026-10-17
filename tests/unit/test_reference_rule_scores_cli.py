"""The raw cohort and every continuous eligible case survive retrieval."""
import argparse
import json
from pathlib import Path

import polars as pl
import pytest


def _run(source, prefix, *options):
    from vdjmatch.cli.reference_rule_scores import register
    parser = argparse.ArgumentParser()
    register(parser.add_subparsers())
    args = parser.parse_args(['reference-rule-scores', '--vdjdb', str(source), '--locus', 'TRB',
        '--output-prefix', str(prefix), *options])
    return args.func(args)


def test_raw_cli_complete_cases_fixed_models_and_sparse_candidates(tmp_path, monkeypatch):
    from vdjmatch import db
    from vdjmatch.evalue import reference_rule_bank
    from vdjmatch.evalue.epitope_rules import reference_epitope_rules

    basic = pl.DataFrame({'cdr3': ['CASSLGQAYEQYF', 'CASSLGRAYEQYF', 'CASSWGRAYEQYF', 'CASSAAAAYEQYF'],
        'v': ['TRBV19', 'TRBV19', 'TRBV19', None], 'j': ['TRBJ2-7']*4,
        'gene': ['TRB']*4, 'species': ['HomoSapiens']*4, 'epitope': ['A', 'A', 'B', 'B'],
        'mhc_a': ['HLA-A*02:01']*4, 'mhc_b': ['B2M']*4, 'mhc_class': ['MHCI']*4})
    sparse = basic.head(1).with_columns(pl.lit('S').alias('epitope'))
    # Source observations, not deduplicated search keys, determine eligibility.
    raw = pl.concat([basic]*15+[sparse, basic.head(1).with_columns(pl.lit('Z').alias('epitope'),
        pl.lit('UNKNOWN').alias('v'))])
    source, prefix = tmp_path/'raw.tsv', tmp_path/'scores'
    raw.write_csv(source, separator='\t')
    loaded, load = [], db.load
    def fresh(*args, **kwargs):
        loaded.append(kwargs)
        return load(*args, **kwargs)
    monkeypatch.setattr(db, 'load', fresh)
    original = reference_rule_bank.stream_reference_rules
    comparisons, expected_winners = [], []
    def checked(queries, reference, **kwargs):
        assert set(reference['epitope']) == {'A', 'B', 'S'}
        assert reference.height == 46
        expected = []
        for epitope in sorted(reference['epitope'].unique()):
            expected.append(reference_epitope_rules(queries,
                reference.filter(pl.col('epitope') == epitope).drop('epitope'),
                reference.filter(pl.col('epitope') != epitope).drop('epitope'),
                exclude_exact=True).with_columns(pl.lit(epitope).alias('epitope')))
        expected = pl.concat(expected).sort('query_id', 'epitope')
        consume = kwargs.pop('consume')
        def compare(scorer, frame):
            assert frame.sort('query_id', 'epitope')['score'].to_list() == pytest.approx(
                expected[scorer+'_score'].to_list(), abs=1e-14)
            comparisons.append(scorer)
            expected_winners.append(frame.filter(pl.col('score') == pl.col('score').max().over('query_id'))
                .with_columns(pl.lit(scorer).alias('scorer')))
            consume(scorer, frame)
        return original(queries, reference, consume=compare, **kwargs)
    monkeypatch.setattr(reference_rule_bank, 'stream_reference_rules', checked)
    assert _run(source, prefix) == 0
    assert loaded == [{'species': 'HomoSapiens'}]
    assert len(comparisons) == 15
    queries = pl.read_csv(str(prefix)+'.queries.tsv', separator='\t')
    truth = pl.read_csv(str(prefix)+'.truth.tsv', separator='\t')
    scores = pl.read_csv(str(prefix)+'.scores.tsv', separator='\t')
    continuous = pl.read_parquet(str(prefix)+'.continuous/*.parquet')
    support = pl.read_csv(str(prefix)+'.target-support.tsv', separator='\t')
    manifest = json.loads((tmp_path/'scores.manifest.json').read_text())
    assert queries.height == 4 and queries['model_available'].sum() == 3
    assert truth.height == 5 and truth.filter(pl.col('epitope') == 'S').height == 1
    assert set(support['epitope']) == {'A', 'B', 'S', 'Z'}
    assert support.filter(pl.col('classification_eligible'))['valid_source_records'].to_list() == [30, 30]
    assert continuous.height == 4*2*15
    assert continuous.select('query_id', 'epitope', 'scorer').unique().height == continuous.height
    assert set(continuous['epitope']) == {'A', 'B'}
    assert continuous.group_by('scorer').len()['len'].to_list() == [8]*15
    missing = queries.filter(~pl.col('model_available'))['query_id'][0]
    assert continuous.filter(pl.col('query_id') == missing)['score'].null_count() == 30
    nohits = scores.filter(pl.col('query_id') == missing)
    assert nohits.height == nohits['epitope'].null_count() == 15
    assert not nohits['supported'].any()
    assert scores['query_id'].n_unique() == 4
    winner_columns = ['query_id', 'epitope', 'scorer', 'score', 'supported']
    assert scores.filter(pl.col('supported')).with_columns(pl.col('query_id').cast(pl.UInt32)).select(winner_columns).sort('scorer', 'query_id', 'epitope').equals(
        pl.concat(expected_winners).select(winner_columns).sort('scorer', 'query_id', 'epitope'),
        null_equal=True)
    assert manifest['source_frozen'] and manifest['reference_frozen']
    assert manifest['continuous_rows'] == manifest['expected_continuous_rows'] == 120
    assert manifest['candidate_targets'] == 4
    with pytest.raises(FileExistsError):
        _run(source, prefix)


def test_all_unavailable_reference_keeps_complete_null_outputs(tmp_path):
    raw = pl.DataFrame({'cdr3': ['CASSLGQAYEQYF', 'CASSLGRAYEQYF'], 'v': ['UNKNOWN', None],
        'j': ['TRBJ2-7']*2, 'gene': ['TRB']*2, 'species': ['HomoSapiens']*2,
        'epitope': ['A', 'B']})
    source, prefix = tmp_path/'raw.tsv', tmp_path/'unavailable'
    raw.write_csv(source, separator='\t')
    assert _run(source, prefix, '--min-records', '1') == 0
    continuous = pl.read_parquet(str(prefix)+'.continuous/*.parquet')
    assert continuous.height == 60 and continuous['score'].null_count() == 60
    scores = pl.read_csv(str(prefix)+'.scores.tsv', separator='\t')
    assert scores.height == scores['epitope'].null_count() == 30


def test_missing_target_short_query_and_curated_population(tmp_path):
    rows = pl.DataFrame({'cdr3': ['CASSLGQAYEQYF', 'CASSLGRAYEQYF', 'CASF', 'CASSLGQAYEQYF', 'CASSLGQAYEQYF', 'CASSZZQAYEQYF'],
        'v': ['TRBV19', None, 'TRBV19', 'TRBV19', 'TRBV19', 'TRBV19'],
        'j': [None, 'TRBJ2-7', 'TRBJ2-7', 'TRBJ1-1', 'TRBJ2-7', 'TRBJ2-7'],
        'gene': ['TRB']*6, 'species': ['HomoSapiens']*6,
        'epitope': ['A', 'B', 'C', 'A', 'KLGGALQAK', 'INVALID']})
    source, prefix = tmp_path/'raw.tsv', tmp_path/'curated'
    rows.write_csv(source, separator='\t')
    assert _run(source, prefix, '--min-records', '1') == 0
    queries = pl.read_csv(str(prefix)+'.queries.tsv', separator='\t')
    assert queries.height == 4 and queries['model_available'].sum() == 2
    assert queries.filter(pl.col('j').is_null())['model_available'].all()
    continuous = pl.read_parquet(str(prefix)+'.continuous/*.parquet')
    assert continuous.height == 4*3*15
    absent = continuous.filter(pl.col('epitope').is_in(['B', 'C']))
    assert absent.height == absent['score'].null_count() == 120
    override = tmp_path/'allvalid'
    assert _run(source, override, '--min-records', '1', '--reference-population', 'all-valid') == 0
    targets = pl.read_csv(str(override)+'.target-support.tsv', separator='\t')
    assert set(targets['epitope']) == {'A', 'B', 'C', 'KLGGALQAK'}


def test_continuous_export_parts_are_bounded_and_hashes_verified(tmp_path):
    from itertools import islice, product
    from vdjmatch.db.cache import sha256
    sequences = ['CASS'+''.join(chars)+'YEQYF'
        for chars in islice(product('ACDEFGHIKLMNPQRSTVWY', repeat=3), 4400)]
    source, prefix = tmp_path/'raw.tsv', tmp_path/'parts'
    pl.DataFrame({'cdr3': sequences, 'v': ['UNKNOWN']*len(sequences), 'j': ['TRBJ2-7']*len(sequences),
        'gene': ['TRB']*len(sequences), 'species': ['HomoSapiens']*len(sequences),
        'epitope': ['A']*len(sequences)}).write_csv(source, separator='\t')
    assert _run(source, prefix) == 0
    parts = sorted((tmp_path/'parts.continuous').glob('*.parquet'))
    assert [pl.read_parquet(p).height for p in parts] == [65536, 464]
    manifest = json.loads((tmp_path/'parts.manifest.json').read_text())
    assert manifest['min_records'] == 30
    for name, expected in manifest['tables'].items():
        assert sha256(tmp_path/name) == expected


@pytest.mark.parametrize('invalid', [('--components',), ('--control', 'unused'), ('--germline-background', 'unused')])
def test_component_options_require_explicit_prior(tmp_path, invalid):
    with pytest.raises(ValueError, match='requires?'):
        _run(tmp_path/'not-read.tsv', tmp_path/'out', *invalid)


@pytest.mark.parametrize('mutation', [None, 'control', 'prior', 'duplicate'])
def test_component_completion_and_raw_input_freeze(tmp_path, monkeypatch, mutation):
    from vdjmatch.evalue import control
    from vdjmatch.runner import retrieval_components
    from vdjmatch import io
    from vdjmatch.db.cache import sha256
    raw = pl.DataFrame({'cdr3': ['CASSLGQAYEQYF', 'CASSLGRAYEQYF', 'CASF', 'CASSLGQAYEQYF'],
        'v': ['TRBV19', None, 'TRBV19', 'TRBV19'], 'j': ['TRBJ2-7']*4,
        'gene': ['TRB']*4, 'species': ['HomoSapiens']*4, 'epitope': ['A', 'B', 'B', 'S']})
    source, prefix = tmp_path/'raw.tsv', tmp_path/'components'
    raw.write_csv(source, separator='\t')
    density, prior = tmp_path/'control.tsv', tmp_path/'prior.tsv'
    density.write_text('cdr3\nCASSLGQAYEQYF\n')
    prior.write_text('cdr3\tv\tj\nCASSLGQAYEQYF\tTRBV19\tTRBJ2-7\n')
    density_before, prior_before = sha256(density), sha256(prior)
    calls, density_index = [], object()
    def raw_control(locus, species, path):
        assert locus == 'TRB' and species == 'human' and path == str(density)
        calls.append('control')
        return density_index, {'sha256': sha256(density), 'path': str(density)}
    def read_prior(path, *, source):
        assert path == prior and source == 'legacy'
        calls.append('prior')
        return raw.head(1)
    component_names = ['accepted_original', 'current_original', 'current_original_log', 'gapblock', 'stem4', 'motif']
    def stream(queries, reference, context, *, consume, **kwargs):
        assert queries['query_id'].to_list() == [0, 1, 2]
        assert queries.height == 3 and set(queries.columns) == {'query_id', 'cdr3', 'v', 'j', 'locus'}
        assert reference.height == 5 and set(reference['epitope']) == {'A', 'B', 'S'}
        assert reference.filter(pl.col('v').is_null()).height == 1
        assert context['density_index'] is density_index and context['prior_background'].height == 1
        # Explicit finite zero survives the continuous export despite unsupported.
        frame = queries.head(1).select('query_id').with_columns(pl.lit('B').alias('epitope'),
            pl.lit(0.).alias('score'), pl.lit(False).alias('supported'))
        for scorer in component_names:
            consume(scorer, frame)
            # Sparse labels participate in winner selection even when not eligible.
            sparse = frame.with_columns(pl.lit('S').alias('epitope'), pl.lit(5.).alias('score'), pl.lit(True).alias('supported'))
            consume(scorer, sparse)
        if mutation in {'control', 'prior'}:
            (density if mutation == 'control' else prior).write_text('changed while scoring\n')
        if mutation == 'duplicate':
            consume(component_names[0], frame)
    monkeypatch.setattr(control, 'raw_background', raw_control)
    monkeypatch.setattr(io, 'read_rearrangement', read_prior)
    monkeypatch.setattr(retrieval_components, 'stream_components', stream)
    options = ['--components', '--control', str(density), '--germline-background', str(prior), '--min-records', '2']
    # B has two source observations; A's repeated record makes it eligible too.
    raw = pl.concat([raw, raw.head(1)])
    raw.write_csv(source, separator='\t')
    if mutation:
        with pytest.raises(ValueError, match='duplicate component|changed during fresh'):
            _run(source, prefix, *options)
        assert not (tmp_path/'components.manifest.json').exists()
        return
    assert _run(source, prefix, *options) == 0
    assert calls == ['control', 'prior']
    continuous = pl.read_parquet(str(prefix)+'.continuous/*.parquet')
    assert continuous.height == 3*2*21
    assert continuous.select('query_id', 'epitope', 'scorer').unique().height == continuous.height
    rows = continuous.filter(pl.col('scorer').is_in(component_names))
    assert rows.height == 36 and rows['score'].null_count() == 30
    assert rows.filter(pl.col('score') == 0.).height == 6 and not rows['supported'].any()
    winners = pl.read_csv(str(prefix)+'.scores.tsv', separator='\t').filter(pl.col('scorer').is_in(component_names))
    assert winners.filter(pl.col('supported'))['epitope'].to_list() == ['S']*6
    assert winners.filter(~pl.col('supported')).height == 12
    manifest = json.loads((tmp_path/'components.manifest.json').read_text())
    assert manifest['component_inputs_before'] == manifest['component_inputs_after']
    assert manifest['component_inputs_before']['control']['sha256'] == density_before
    assert manifest['component_inputs_before']['germline_background']['sha256'] == prior_before
    assert manifest['component_inputs_frozen'] and len(manifest['scorers']) == 21


def test_bundled_control_identity_and_all_valid_component_population(tmp_path, monkeypatch):
    from importlib.resources import files
    from vdjmatch.evalue import control
    from vdjmatch.runner import retrieval_components
    from vdjmatch import io
    from vdjmatch.db.cache import sha256
    raw = pl.DataFrame({'cdr3': ['CASSLGQAYEQYF', 'CASSLGRAYEQYF'], 'v': ['TRBV19']*2,
        'j': ['TRBJ2-7']*2, 'gene': ['TRB']*2, 'species': ['HomoSapiens']*2,
        'epitope': ['A', 'KLGGALQAK']})
    source, prior, prefix = tmp_path/'raw.tsv', tmp_path/'prior.tsv', tmp_path/'bundled'
    raw.write_csv(source, separator='\t')
    prior.write_text('raw prior identity\n')
    asset = Path(str(files('seqtree').joinpath('data/control_human_trb_aa.txt.gz')))
    asset_hash = sha256(asset)
    calls = []
    def bundled(locus, species, path):
        assert path is None and locus == 'TRB' and species == 'human'
        calls.append('bundled')
        return object(), {'kind': 'bundled_raw', 'asset': 'seqtree/data/control_human_trb_aa.txt.gz', 'sha256': asset_hash}
    def stream(queries, reference, context, *, consume, **kwargs):
        assert set(reference['epitope']) == {'A', 'KLGGALQAK'}
        assert reference['reference_valid'].all()
        calls.append('stream')
    monkeypatch.setattr(control, 'raw_background', bundled)
    monkeypatch.setattr(io, 'read_rearrangement', lambda *args, **kwargs: raw)
    monkeypatch.setattr(retrieval_components, 'stream_components', stream)
    assert _run(source, prefix, '--components', '--germline-background', str(prior),
        '--min-records', '1', '--reference-population', 'all-valid') == 0
    assert calls == ['bundled', 'stream']
    manifest = json.loads((tmp_path/'bundled.manifest.json').read_text())
    assert manifest['component_inputs_before']['control']['sha256'] == asset_hash
    assert manifest['control_info']['sha256'] == asset_hash
    continuous = pl.read_parquet(str(prefix)+'.continuous/*.parquet')
    assert continuous.height == 2*2*21
