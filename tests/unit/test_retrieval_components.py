"""Streamed single-chain retrieval preserves the retained fixed-model laws."""
import polars as pl
import pytest
from polars.testing import assert_frame_equal
from seqtree import Index

from vdjmatch.aggregate.candidates import PMHC
from vdjmatch.runner.retrieval_components import stream_components


def _inputs(species='human'):
    reference = pl.DataFrame({
        'cdr3': ['CASSLGQETQYF', 'CASSIGQETQYF', 'CATSLGQETQYF'],
        'v': ['TRBV1'] * 3, 'j': ['TRBJ1-1'] * 3, 'gene': ['TRB'] * 3,
        'species': ['HomoSapiens' if species == 'human' else 'MusMusculus'] * 3,
        'epitope': ['A', 'A', 'B'], 'mhc_a': ['HLA-A*02:01'] * 3,
        'mhc_b': ['B2M'] * 3, 'mhc_class': ['MHCI'] * 3})
    query = pl.DataFrame({'query_id': [0, 1, 2, 3, 4],
        'cdr3': ['CASSLGQETQYF', 'CASSLGQETQYF', 'CASSMGQETQYF', 'W' * 200, '?'],
        'v': ['TRBV1', None, 'TRBV1', 'TRBV1', 'TRBV1'],
        'j': ['TRBJ1-1', 'TRBJ1-1', None, 'TRBJ1-1', 'TRBJ1-1'], 'locus': ['TRB'] * 5})
    control = {'density_index': Index.build(reference['cdr3'].to_list() + ['CASSAGQETQYF']),
               'prior_background': reference.select('cdr3', 'v', 'j')}
    return query, reference, control


@pytest.mark.parametrize('species', ['human', 'mouse'])
def test_stream_matches_scalar_oracle_fixed_exact_counts_and_batched_consumers(monkeypatch, species):
    from vdjmatch.runner import combined
    from vdjmatch.match.historical import density, germline_prior, unpaired_score

    query, reference, control = _inputs(species)
    expected = combined.score_components(query, reference, control, species=species,
                                         locus='TRB', exclude_exact=True)
    # Excluded exact hits stay in both fixed reference counts and germline priors.
    selected, ref = query.head(1), reference.filter(pl.col('epitope') == 'A')
    density_result = density(selected, ref, control['density_index'],
                             pool_reference=True, puncture_populations=False)
    assert density_result['n_reference'].item() == 2
    prior = germline_prior(selected, ref, control['prior_background'],
                           include_length=True, exclude_exact=False)
    oracle = unpaired_score(density_result['score'], prior['germline_lr'], 2)
    # Other species/locus must not contaminate either the target axis or counts.
    extra = reference.head(1).with_columns(pl.lit('OTHER').alias('epitope'),
        pl.lit('MusMusculus' if species == 'human' else 'HomoSapiens').alias('species'))
    extra_chain = reference.head(1).with_columns(pl.lit('TRA_ONLY').alias('epitope'),
        pl.lit('TRA').alias('gene'))
    reference = pl.concat([reference, extra, extra_chain])
    calls = []
    original, kernel = combined._original, combined._kernel

    def checked_original(*args, **kwargs):
        assert args[2] is None and kwargs['exclude_exact'] is True
        assert callable(kwargs['_consume'])
        calls.append('original')
        return original(*args, **kwargs)

    def checked_kernel(*args, **kwargs):
        assert args[2] is None and kwargs['exclude_exact'] is True
        assert callable(kwargs['_consume'])
        calls.append('kernel')
        return kernel(*args, **kwargs)

    monkeypatch.setattr(combined, '_original', checked_original)
    monkeypatch.setattr(combined, '_kernel', checked_kernel)
    monkeypatch.setattr(combined, '_cases', lambda *a, **k: pytest.fail('dense cases forbidden'))
    frames = {}

    def consume(name, frame):
        assert frame.columns == ['query_id', *PMHC, 'score', 'supported']
        assert frame.height <= query.height * 2
        frames.setdefault(name, []).append(frame)

    stream_components(query, reference, control, species=species, locus='TRB', consume=consume)
    assert calls == ['original', 'kernel']
    actual = {name: pl.concat(parts).sort('query_id', *PMHC) for name, parts in frames.items()}
    assert set(actual) == {'accepted_original', 'current_original', 'current_original_log', 'gapblock', 'stem4'}
    assert_frame_equal(actual['accepted_original'], actual['current_original'], check_exact=True)
    for name, column in [('current_original_log', 'original_score'), ('gapblock', 'gapblock_score'),
                         ('stem4', 'stem4_score')]:
        got = actual[name]
        want = expected.select('query_id', *PMHC, pl.col(column).alias('score'))
        want = got.select('query_id', *PMHC).join(want, on=['query_id', *PMHC],
            nulls_equal=True, validate='1:1').sort('query_id', *PMHC)
        assert_frame_equal(got.drop('supported'), want, check_exact=True)
    assert actual['accepted_original'].filter((pl.col('query_id') == 0) &
        (pl.col('epitope') == 'A'))['score'].item() == oracle.item()
    for name in ('accepted_original', 'current_original_log'):
        far = actual[name].filter(pl.col('query_id') == 3)
        assert far['score'].is_finite().all() and not far['supported'].any()
        missing = actual[name].filter(pl.col('query_id').is_in([1, 2]))
        assert missing['score'].null_count() == missing.height
        assert not missing['supported'].any()


def test_dense_raw_law_uses_fixed_unique_count_not_log_or_punctured_count(monkeypatch):
    from vdjmatch.runner import combined

    query, reference, control = _inputs()
    # Distinct canonical junctions create a 500-key fixed model including query 0.
    amino = 'ACDEFGHIKLMNPQRSTVWY'
    junctions = ['CASSLGQETQYF'] + ['C' + amino[i // 20] + amino[i % 20] + 'F'
        for i in range(399)] + ['CA' + amino[i // 20] + amino[i % 20] + 'F' for i in range(100)]
    reference = pl.concat([reference.head(1)] * 500).with_columns(pl.Series('cdr3', junctions))
    control['v_loop_distances'] = {}

    def original(queries, reference, cases, *args, **kwargs):
        frame = queries.select('query_id').join(reference.select(PMHC).unique(), how='cross').with_columns(
            pl.lit(3.).alias('original_density'), pl.lit(2.).alias('germline_lr'),
            pl.lit(1.0986126220013874).alias('original_score'), pl.lit(True).alias('original_supported'))
        kwargs['_consume'](frame)

    monkeypatch.setattr(combined, '_original', original)
    monkeypatch.setattr(combined, '_kernel', lambda *a, **k: None)
    frames = {}
    stream_components(query.head(1), reference, control, species='human', locus='TRB',
        consume=lambda name, frame: frames.setdefault(name, frame))
    assert frames['current_original']['score'].item() == 3.
    assert frames['current_original_log']['score'].item() == 1.0986126220013874


def test_validation_and_explicit_paired_rejection_emit_nothing():
    query, reference, control = _inputs()
    consume = lambda *a: pytest.fail('invalid input must fail before scoring')
    with pytest.raises(NotImplementedError, match='cohort-rank'):
        stream_components(query, reference, control, species='human', locus='paired', consume=consume)
    with pytest.raises(ValueError, match='contiguous integers'):
        stream_components(query.with_columns(pl.col('query_id') + 1), reference, control,
                          species='human', locus='TRB', consume=consume)
    with pytest.raises(ValueError, match='locus differs'):
        stream_components(query, reference, control, species='human', locus='TRA', consume=consume)
    stream_components(query.head(0), reference, control, species='human', locus='TRB', consume=consume)


def test_motif_adapter_calls_existing_all_label_consumer_once(monkeypatch):
    from vdjmatch.runner import combined
    from vdjmatch.match import motifs
    query,reference,control=_inputs()
    reference=pl.concat([reference.head(1)]*7).with_columns(
        pl.Series('cdr3',['CASSAAAA'+a+'F' for a in 'AILVDEN']))
    query=query.with_columns(pl.lit('CASSAAAAAF').alias('cdr3'))
    control['v_loop_distances']={}
    monkeypatch.setattr(combined,'_original',lambda *a,**k:None)
    monkeypatch.setattr(combined,'_kernel',lambda *a,**k:None)
    owner=motifs.masked_motif_scores
    expected=owner(query.filter(pl.col('v').is_not_null()),reference,exclude_exact=True)
    calls=[]
    def checked(queries,reference,**kwargs):
        assert kwargs['exclude_exact'] is True and callable(kwargs['_consume'])
        assert 'target_reference' not in kwargs
        calls.append(queries.height)
        return owner(queries,reference,**kwargs)
    monkeypatch.setattr(motifs,'masked_motif_scores',checked)
    parts=[]
    def consume(name,frame):
        assert name=='motif'
        assert frame['supported'].all()
        parts.append(frame)
    stream_components(query,reference,control,species='human',locus='TRB',consume=consume)
    assert calls==[4]
    assert_frame_equal(pl.concat(parts).drop('supported').sort('query_id',*PMHC),
        expected.select('query_id',*PMHC,pl.col('motif_score').alias('score')).sort('query_id',*PMHC),
        check_exact=True)


def test_missing_v_reference_is_excluded_only_from_fixed_motif_population(monkeypatch):
    from vdjmatch.runner import combined
    from vdjmatch.match import motifs
    query,reference,control=_inputs()
    valid=pl.concat([reference.head(1)]*7).with_columns(
        pl.Series('cdr3',['CASSAAAA'+a+'F' for a in 'AILVDEN']))
    missing=valid.head(2).with_columns(pl.Series('v',[None,'  ']),
        pl.Series('cdr3',['CASSAAAATF','CASSAAAAGF']))
    reference=pl.concat([valid,missing])
    query=query.head(1).with_columns(pl.lit('CASSAAAAAF').alias('cdr3'))
    control['density_index']=Index.build(reference['cdr3'].to_list())
    owner=motifs.masked_motif_scores
    expected=owner(query,valid,exclude_exact=True).select('query_id',*PMHC,
        pl.col('motif_score').alias('score'),pl.lit(True).alias('supported'))
    original,kernel=combined._original,combined._kernel
    populations=[]
    def checked_original(queries,ref,*args,**kwargs):
        populations.append(('original',ref.height))
        return original(queries,ref,*args,**kwargs)
    def checked_kernel(queries,ref,*args,**kwargs):
        populations.append(('kernel',ref.height))
        return kernel(queries,ref,*args,**kwargs)
    monkeypatch.setattr(combined,'_original',checked_original)
    monkeypatch.setattr(combined,'_kernel',checked_kernel)
    frames={}
    stream_components(query,reference,control,species='human',locus='TRB',
        consume=lambda name,frame:frames.setdefault(name,[]).append(frame))
    assert set(frames)=={'accepted_original','current_original','current_original_log','gapblock','stem4','motif'}
    assert populations==[('original',9),('kernel',9)]
    assert_frame_equal(pl.concat(frames['motif']).sort('query_id',*PMHC),
                       expected.sort('query_id',*PMHC),check_exact=True)
