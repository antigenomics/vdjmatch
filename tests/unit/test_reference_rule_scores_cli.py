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
        assert set(reference['epitope']) == {'A', 'B', 'S', 'Z'}
        assert reference.height == 62
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
    assert queries.height == 4 and queries['model_available'].all()
    assert queries.filter(pl.col('v').is_null())['v_geometry_imputed'].all()
    assert truth.height == 5 and truth.filter(pl.col('epitope') == 'S').height == 1
    assert set(support['epitope']) == {'A', 'B', 'S', 'Z'}
    assert support.filter(pl.col('classification_eligible'))['valid_source_records'].to_list() == [30, 30]
    assert continuous.height == 4*2*15
    assert continuous.select('query_id', 'epitope', 'scorer').unique().height == continuous.height
    assert set(continuous['epitope']) == {'A', 'B'}
    assert continuous.group_by('scorer').len()['len'].to_list() == [8]*15
    assert continuous['score'].null_count() == 0
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


def test_all_unresolved_reference_is_imputed_and_keeps_complete_outputs(tmp_path):
    raw = pl.DataFrame({'cdr3': ['CASSLGQAYEQYF', 'CASSLGRAYEQYF'], 'v': ['UNKNOWN', None],
        'j': ['TRBJ2-7']*2, 'gene': ['TRB']*2, 'species': ['HomoSapiens']*2,
        'epitope': ['A', 'B']})
    source, prefix = tmp_path/'raw.tsv', tmp_path/'unavailable'
    raw.write_csv(source, separator='\t')
    assert _run(source, prefix, '--min-records', '1') == 0
    continuous = pl.read_parquet(str(prefix)+'.continuous/*.parquet')
    assert continuous.height == 60 and continuous['score'].null_count() == 0
    queries = pl.read_csv(str(prefix)+'.queries.tsv', separator='\t')
    assert queries['v_geometry_imputed'].all() and queries['model_available'].all()
    assert queries['v'].to_list() == ['UNKNOWN', None]
    scores = pl.read_csv(str(prefix)+'.scores.tsv', separator='\t')
    assert scores.select('query_id', 'scorer').unique().height == 30
    assert scores['epitope'].null_count() == 0


@pytest.mark.parametrize('kernel', [False, True])
def test_imputed_reference_key_does_not_collapse_into_template_key(tmp_path, kernel):
    source, prefix = tmp_path/'raw.tsv', tmp_path/'scores'
    pl.DataFrame(dict(cdr3=['CASSLGQAYEQYF']*2+['CASSLGRAYEQYF'],
        v=['TRBV1', None, 'TRBV1'], j=['TRBJ2-7']*3, gene=['TRB']*3,
        species=['HomoSapiens']*3, epitope=['A', 'A', 'B'])).write_csv(source, separator='\t')
    assert _run(source, prefix, '--min-records', '1', *(['--kernel-only'] if kernel else [])) == 0
    queries = pl.read_csv(str(prefix)+'.queries.tsv', separator='\t')
    assert queries.height == 3 and queries['model_available'].all()
    assert queries['v_geometry_imputed'].sum() == 1
    support = pl.read_csv(str(prefix)+'.target-support.tsv', separator='\t')
    assert support.filter(pl.col('epitope') == 'A')['modeled_search_keys'].item() == 2
    manifest = json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert manifest['backend']['reference_keys'] == 3


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
    assert queries.height == 4 and queries['model_available'].sum() == 3
    assert queries.filter(pl.col('j').is_null())['model_available'].all()
    continuous = pl.read_parquet(str(prefix)+'.continuous/*.parquet')
    assert continuous.height == 4*3*15
    absent = continuous.filter(pl.col('epitope') == 'C')
    assert absent.height == absent['score'].null_count() == 60
    short = queries.filter(pl.col('cdr3') == 'CASF')['query_id'].item()
    assert continuous.filter(pl.col('query_id') == short)['score'].null_count() == 45
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
        assert queries.height == 3 and set(queries.columns) == {'query_id', 'cdr3', 'v', 'j', 'locus', 'original_v', 'v_geometry_imputed'}
        assert reference.height == 5 and set(reference['epitope']) == {'A', 'B', 'S'}
        assert reference['v'].null_count() == 0
        assert reference.filter(pl.col('original_v').is_null())['v_geometry_imputed'].to_list() == [True]
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


def test_kernel_only_exports_current_six_scores_without_rule_or_control_rebuild(tmp_path,monkeypatch):
    from vdjmatch.evalue import reference_rule_bank
    def forbidden(*args,**kwargs):raise AssertionError('completed rules must not rerun')
    monkeypatch.setattr(reference_rule_bank,'stream_reference_rules',forbidden)
    source=tmp_path/'raw.tsv';prefix=tmp_path/'kernels'
    pl.DataFrame(dict(cdr3=['CASSLGQAYEQYF','CASSLGRAYEQYF','CASSWGRAYEQYF','CASSAAAAYEQYF'],v=['TRBV19']*3+[None],j=['TRBJ2-7']*4,gene=['TRB']*4,species=['HomoSapiens']*4,epitope=['A','A','B','B'],mhc_a=['HLA-A*02:01']*4,mhc_b=['B2M']*4,mhc_class=['MHCI']*4)).write_csv(source,separator='\t')
    assert _run(source,prefix,'--kernel-only','--min-records','1')==0
    data=pl.read_parquet(str(prefix)+'.continuous/*.parquet')
    assert data.height==4*2*6
    assert set(data['scorer'])=={'kernel_'+x+'_'+y for x in ['v','all','vonly'] for y in ['bounded','rare']}
    assert data.select('query_id','epitope','scorer').unique().height==data.height
    queries=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t')
    unknown=queries.filter(pl.col('v_geometry_imputed'))['query_id'].item()
    assert queries['model_available'].all()
    assert data.filter(pl.col('query_id')==unknown)['score'].null_count()==0
    from vdjmatch.cli.neighbours import register
    parser=argparse.ArgumentParser();register(parser.add_subparsers())
    query_path=tmp_path/'queries.tsv'
    queries.rename({'cdr3':'junction_aa','v':'v_call','j':'j_call'}).write_csv(query_path,separator='\t')
    native_prefix=tmp_path/'native'
    args=parser.parse_args(['neighbours',str(query_path),'--vdjdb',str(source),'--locus','TRB','--species','human','--reference-population','curated-benchmark','--all-targets','--exclude-exact','--background','vdjdb-other','--neighbour-weighting','kernel','--distance-matrix','blosum62','--junction-ends','trim3','--stem-depth','4','--evidence-only','--stream-evidence','--output-prefix',str(native_prefix)])
    assert args.func(args)==0
    native=pl.read_csv(str(native_prefix)+'.enrichment.tsv',separator='\t')
    for scorer in data['scorer'].unique():
        observed=data.filter(pl.col('scorer')==scorer).sort('query_id','epitope')
        expected=native.select('query_id','epitope',pl.col(scorer).alias('score')).sort('query_id','epitope')
        assert observed.select('query_id','epitope','score').equals(expected)

    manifest=json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert manifest['source_frozen'] and manifest['reference_frozen'] and manifest['kernel_only']
    assert manifest['backend']['native_batches']==1 and manifest['expected_continuous_rows']==48


@pytest.mark.parametrize('kernel',[False,True])
def test_genuine_paired_reference_exports_use_max_chain_eligibility(tmp_path,kernel):
    rows=[]
    for complex_id,epitope,alpha,beta,vbeta in [(1,'A','CAVRLGQAGKSTF','CASSLGQAYEQYF','TRBV19'),(2,'A','CAVRLGQGGKSTF','CASSLGRAYEQYF',None),(3,'B','CAVRLGQSGKSTF','CASSWGRAYEQYF','TRBV19')]:
        for locus,seq,v,j in [('TRA',alpha,'TRAV1-1','TRAJ1'),('TRB',beta,vbeta,'TRBJ2-7')]:
            rows.append(dict(complex_id=complex_id,epitope=epitope,gene=locus,cdr3=seq,v=v,j=j,species='HomoSapiens',mhc_a='HLA-A*02:01',mhc_b='B2M',mhc_class='MHCI'))
    rows.extend([dict(rows[0],complex_id=0)]*30)
    rows.extend([dict(rows[0],complex_id=0,epitope='C')]*30)
    source=tmp_path/'raw-paired.tsv';prefix=tmp_path/'paired';pl.DataFrame(rows).write_csv(source,separator='\t')
    options=['--locus','paired']+(['--kernel-only'] if kernel else [])
    assert _run(source,prefix,*options)==0
    q=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t');support=pl.read_csv(str(prefix)+'.target-support.tsv',separator='\t')
    assert q.height==2 and q['model_available'].all()
    assert q['v_geometry_imputed1'].sum()==1
    assert support.filter(pl.col('classification_eligible'))['epitope'].to_list()==['A','C']
    assert support.filter(pl.col('epitope')=='A')['paired_source_records'].item()==2
    assert support.filter(pl.col('epitope')=='A')['valid_source_records'].item()==32
    data=pl.read_parquet(str(prefix)+'.continuous/*.parquet')
    assert support.filter(pl.col('epitope')=='C')['paired_source_records'].item()==0
    assert data.height==4*(6 if kernel else 15) and set(data['locus'])=={'PAIRED'}
    assert data.filter(pl.col('epitope')=='C')['score'].null_count()==2*(6 if kernel else 15)
    assert data.filter(pl.col('epitope')=='A')['score'].null_count()==0
    manifest=json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert manifest['source_frozen'] and manifest['reference_frozen'] and manifest['paired']
    assert manifest['paired_eligibility']=='maximum qualifying chain source-record count before exact-hit filtering'
    if kernel:
        from vdjmatch.cli.neighbours import register
        chains=[]
        for side,locus in [('alpha','TRA'),('beta','TRB')]:
            chains.append(q.select(pl.col('query_id').cast(pl.String).alias('clone_id'),pl.lit(locus).alias('locus'),pl.col(side).alias('junction_aa'),pl.col('v_'+side).alias('v_call'),pl.col('j_'+side).alias('j_call')))
        raw=pl.concat(chains).sort('clone_id','locus').with_row_index('sequence_id');query_path=tmp_path/'paired-queries.tsv';raw.write_csv(query_path,separator='\t')
        parser=argparse.ArgumentParser();register(parser.add_subparsers());native_prefix=tmp_path/'native-paired'
        args=parser.parse_args(['neighbours',str(query_path),'--vdjdb',str(source),'--locus','paired','--species','human','--reference-population','curated-benchmark','--all-targets','--exclude-exact','--background','vdjdb-other','--neighbour-weighting','kernel','--distance-matrix','blosum62','--junction-ends','trim3','--stem-depth','4','--evidence-only','--stream-evidence','--output-prefix',str(native_prefix)])
        assert args.func(args)==0
        native=pl.read_csv(str(native_prefix)+'.enrichment.tsv',separator='\t').filter(pl.col('epitope')=='A')
        for scorer in data['scorer'].unique():
            observed=data.filter((pl.col('scorer')==scorer)&(pl.col('epitope')=='A')&pl.col('score').is_not_null()).sort('query_id','epitope')
            expected=native.filter(pl.col(scorer).is_not_null()).select('query_id','epitope',pl.col(scorer).alias('score')).sort('query_id','epitope')
            assert observed.select('query_id','epitope','score').equals(expected)


def test_resolved_first_cli_retains_cases_and_records_policy(tmp_path):
    raw = pl.DataFrame(dict(cdr3=['CASSLGQAYEQYF', 'CASSLGRAYEQYF', 'CASSLGQAYEQWF'],
        v=['TRBV19', 'TRBV19', None], j=['TRBJ2-7']*3, gene=['TRB']*3,
        species=['HomoSapiens']*3, epitope=['A', 'B', 'C']))
    source, prefix = tmp_path/'raw.tsv', tmp_path/'scores'
    raw.write_csv(source, separator='\t')
    assert _run(source, prefix, '--min-records', '1', '--kernel-only',
        '--v-model-policy', 'resolved-first') == 0
    manifest = json.loads(Path(str(prefix)+'.manifest.json').read_text())
    assert manifest['v_model_policy'] == 'resolved-first'
    assert manifest['backend']['native_batches'] == 2
    assert manifest['source_frozen'] and manifest['reference_frozen']
    assert manifest['continuous_rows'] == manifest['expected_continuous_rows'] == 54
    queries = pl.read_csv(str(prefix)+'.queries.tsv', separator='\t')
    scores = pl.read_parquet(str(prefix)+'.continuous/*.parquet')
    known = queries.filter(pl.col('baseline_v0').is_not_null())['query_id']
    absent = scores.filter(pl.col('query_id').is_in(known.implode()) & (pl.col('epitope') == 'C'))
    assert absent.height == 12 and absent['score'].null_count() == 12
    assert not absent['supported'].any()
    with pytest.raises(ValueError, match='requires --kernel-only'):
        _run(source, tmp_path/'invalid', '--v-model-policy', 'resolved-first')


def test_position_profile_cli_is_opt_in_frozen_and_changes_only_the_weighted_arm(tmp_path):
    raw = pl.DataFrame(dict(cdr3=['CASSLGQAYEQYF', 'CASSLGRAYEQYF', 'CASSLGQAYEQWF'],
        v=['TRBV19', 'TRBV19', 'TRBV7-9'], j=['TRBJ2-7']*3, gene=['TRB']*3,
        species=['HomoSapiens']*3, epitope=['A', 'B', 'C']))
    source = tmp_path/'raw.tsv'
    raw.write_csv(source, separator='\t')
    profiles = tmp_path/'profiles'
    profiles.mkdir()
    def profile(name, core):
        path = profiles/name
        with open(path, 'w') as fh:
            fh.write('locus\tside\toffset\tmultiplier\n')
            for side in ('cys', 'fw'):
                for offset, multiplier in enumerate([1, 1, 1, 1, 1, 1, core]):
                    fh.write(f'TRB\t{side}\t{offset}\t{multiplier}\n')
        return str(path)
    def scores(prefix, *options):
        assert _run(source, prefix, '--min-records', '1', '--kernel-only', *options) == 0
        return pl.read_parquet(str(prefix)+'.continuous/*.parquet').sort('query_id', 'epitope', 'scorer')
    plain = scores(tmp_path/'plain')
    ones = scores(tmp_path/'ones', '--position-profile', profile('ones.tsv', 1))
    heavy = scores(tmp_path/'heavy', '--position-profile', profile('heavy.tsv', 8))
    assert plain.drop('score').equals(ones.drop('score'))
    assert plain['score'].to_list() == pytest.approx(ones['score'].to_list(), nan_ok=True)
    assert plain['score'].to_list() != pytest.approx(heavy['score'].to_list(), nan_ok=True)

    manifest = json.loads(Path(str(tmp_path/'heavy')+'.manifest.json').read_text())
    recorded = manifest['position_profile']
    assert recorded['path'] == str(profiles/'heavy.tsv') and recorded['default_changed'] is False
    assert recorded['multipliers']['TRB'] == dict(cys_offsets=7, fw_offsets=7,
        core=dict(cys=8.0, fw=8.0), min=1.0, max=8.0, mean=pytest.approx((6 + 8)*2/14))
    assert manifest['component_inputs_frozen'] and 'position_profile' in manifest['component_inputs_before']
    assert 'per-position cost multipliers' in manifest['backend']['geometry']
    assert json.loads(Path(str(tmp_path/'plain')+'.manifest.json').read_text())['position_profile'] is None

    with pytest.raises(ValueError, match='requires --kernel-only'):
        _run(source, tmp_path/'invalid', '--position-profile', profile('again.tsv', 1))
    (profiles/'tra.tsv').write_text('locus\tside\toffset\tmultiplier\n'
                                    'TRA\tcys\t0\t1\nTRA\tfw\t0\t1\n')
    with pytest.raises(ValueError, match='does not declare TRB'):
        scores(tmp_path/'wrong-locus', '--position-profile', str(profiles/'tra.tsv'))


def test_explicit_matrix_is_frozen_and_cannot_change_during_scoring(tmp_path, monkeypatch):
    from seqtree import amino_acids
    from vdjmatch.match.geometry import _substitution_costs
    from vdjmatch.evalue import reference_rule_bank
    from vdjmatch.db.cache import sha256

    source = tmp_path/'raw.tsv'
    pl.DataFrame(dict(cdr3=['CASSLGQAYEQYF', 'CASSLGRAYEQYF'],
        v=['TRBV19']*2, j=['TRBJ2-7']*2, gene=['TRB']*2,
        species=['HomoSapiens']*2, epitope=['A', 'B'])).write_csv(source, separator='\t')
    # Built-in names keep precedence over same-named files in the working directory.
    monkeypatch.chdir(tmp_path)
    (tmp_path/'blosum62').write_text('not a matrix')
    assert _run(source, tmp_path/'builtin', '--min-records', '1', '--kernel-only') == 0
    builtin = json.loads((tmp_path/'builtin.manifest.json').read_text())
    assert 'matrix' not in builtin['component_inputs_before']
    matrix = tmp_path/'costs.tsv'
    symbols, costs = amino_acids(), _substitution_costs('full', 'blosum62')
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    matrix.write_text('aa.1\taa.2\tcost\n'+''.join(
        f'{a}\t{b}\t{costs[symbols.index(a), symbols.index(b)]}\n'
        for i, a in enumerate(aa) for b in aa[i+1:]))
    options = ['--min-records', '1', '--kernel-only', '--matrix', str(matrix)]
    assert _run(source, tmp_path/'fixed', *options) == 0
    manifest = json.loads((tmp_path/'fixed.manifest.json').read_text())
    assert manifest['component_inputs_before']['matrix'] == dict(path=str(matrix), sha256=sha256(matrix))
    assert manifest['component_inputs_before'] == manifest['component_inputs_after']
    original = reference_rule_bank.stream_reference_kernels
    def changed(*args, **kwargs):
        result = original(*args, **kwargs)
        matrix.write_text(matrix.read_text()+'\n')
        return result
    monkeypatch.setattr(reference_rule_bank, 'stream_reference_kernels', changed)
    with pytest.raises(ValueError, match='changed during fresh scoring'):
        _run(source, tmp_path/'changed', *options)
    assert not (tmp_path/'changed.manifest.json').exists()
    with pytest.raises(ValueError, match='--matrix requires --kernel-only'):
        _run(source, tmp_path/'ignored', '--matrix', str(matrix))
