"""Pinned linked-first label policy; kernel magnitudes never compete across axes."""
import gzip
import json

import polars as pl
import pytest


def inputs():
    selection = pl.DataFrame({'query_id': ['4294967297', '9223372036854775806', 'tie', 'agree', 'none'],
        'status': ['matched', 'no_neighbours', 'matched', 'matched', 'invalid_alpha_query'],
        'pair_id': ['01', '02', '03', '04', '05'],
        'kernel_score': [.01, None, .02, None, None],
        'kernel_alpha_score': [.9, .8, .8, .7, None],
        'kernel_beta_score': [.8, .3, .5, .4, None],
        'tied_peptides': [1, 0, 2, 0, 0], 'alpha_tied_peptides': [1, 1, 1, 1, 0],
        'beta_tied_peptides': [1, 1, 1, 1, 0],
        'selected_linked_epitope': ['A', None, 'A', None, None],
        'selected_alpha_epitope': ['B', 'A', 'B', 'A', None],
        'selected_beta_epitope': ['B', 'B', 'A', 'A', None],
        'kernel_any_axis_score': [.9, .8, .8, .7, None]})
    evidence = pl.DataFrame({'query_id': ['4294967297']*3+['9223372036854775806']*2+['tie']*2+['agree'],
        'epitope': ['A', 'B', 'A', 'A', 'B', 'A', 'B', 'A'],
        'kernel_linked_score': [.01, .001, .01, None, None, .02, .02, None],
        'kernel_alpha_score': [.001, .9, .001, .8, .1, .1, .8, .7],
        'kernel_beta_score': [.001, .8, .001, .1, .3, .5, .1, .4]})
    return evidence, selection


def test_linked_first_preserves_ties_fallback_and_metadata():
    from vdjmatch.aggregate.paired_labels import paired_labels
    evidence, selection = inputs()
    calls, scores = paired_labels(evidence, selection, species='HomoSapiens')
    assert calls.select(selection.columns).unique(maintain_order=True).equals(selection)
    labels = {q: part['paired_epitope'].to_list() for (q,), part in calls.group_by('query_id', maintain_order=True)}
    assert labels == {'4294967297': ['A'], '9223372036854775806': ['A', 'B'], 'tie': ['A', 'B'], 'agree': ['A'], 'none': [None]}
    assert calls.filter(pl.col('query_id')=='none')['paired_label_source'].item()=='none'
    assert scores.filter(pl.col('supported'))['score'].unique().to_list()==[1.0]
    assert scores.filter(~pl.col('supported'))['score'].to_list()==[None]
    for col in ['kernel_score', 'kernel_alpha_score', 'kernel_beta_score', 'kernel_any_axis_score']:
        assert calls.group_by('query_id', maintain_order=True).agg(pl.col(col).first())[col].equals(selection[col])


@pytest.mark.parametrize('change', ['lost_tie', 'scale', 'unknown'])
def test_native_winner_drift_is_rejected(change):
    from vdjmatch.aggregate.paired_labels import paired_labels
    e, s = inputs()
    if change=='lost_tie': e=e.filter(~((pl.col('query_id')=='tie')&(pl.col('epitope')=='B')))
    elif change=='scale': e=e.with_columns(pl.col('kernel_linked_score')*2)
    else: e=pl.concat([e, e.head(1).with_columns(pl.lit('unknown').alias('query_id'))])
    with pytest.raises(ValueError): paired_labels(e,s,species='HomoSapiens')


def test_cli_plain_and_gzip_equal(tmp_path):
    from vdjmatch.cli.__main__ import main
    e,s=inputs();e.write_csv(tmp_path/'e.tsv',separator='\t');s.write_csv(tmp_path/'s.tsv',separator='\t')
    with gzip.open(tmp_path/'e.tsv.gz','wt') as f: f.write((tmp_path/'e.tsv').read_text())
    (tmp_path/'m.json').write_text(json.dumps({'parameters': {'locus':'paired','species':'human','all_targets':True,'background':'vdjdb-other','neighbour_weighting':'kernel','stem_depth':4,'kernel_scorer':'v-bounded'}}))
    for name,path in [('plain','e.tsv'),('gzip','e.tsv.gz')]:
        assert main(['paired-labels',str(tmp_path/path),'--selection',str(tmp_path/'s.tsv'),'--input-manifest',str(tmp_path/'m.json'),'--output-prefix',str(tmp_path/name),'--threads','4'])==0
    for table in ['calls','scores']:
        assert (tmp_path/f'plain.{table}.tsv').read_bytes()==(tmp_path/f'gzip.{table}.tsv').read_bytes()
    native=pl.read_csv(tmp_path/'s.tsv',separator='\t',infer_schema=False)
    result=pl.read_csv(tmp_path/'plain.calls.tsv',separator='\t',infer_schema=False)
    assert result.select(native.columns).unique(maintain_order=True).equals(native)
    assert json.loads((tmp_path/'plain.manifest.json').read_text())['detection_changed'] is False


@pytest.mark.parametrize('change', ['lost_tie', 'scale', 'unknown'])
def test_cli_rejects_native_drift(tmp_path, change):
    from vdjmatch.cli.__main__ import main
    e,s=inputs()
    if change=='lost_tie': e=e.filter(~((pl.col('query_id')=='tie')&(pl.col('epitope')=='B')))
    elif change=='scale': e=e.with_columns(pl.col('kernel_linked_score')*2)
    else: e=pl.concat([e,e.head(1).with_columns(pl.lit('unknown').alias('query_id'))])
    e.write_csv(tmp_path/'e.tsv',separator='\t');s.write_csv(tmp_path/'s.tsv',separator='\t')
    (tmp_path/'m.json').write_text(json.dumps({'parameters': {'locus':'paired','species':'human','all_targets':True,'background':'vdjdb-other','neighbour_weighting':'kernel','stem_depth':4,'kernel_scorer':'v-bounded'}}))
    with pytest.raises(SystemExit):
        main(['paired-labels',str(tmp_path/'e.tsv'),'--selection',str(tmp_path/'s.tsv'),'--input-manifest',str(tmp_path/'m.json'),'--output-prefix',str(tmp_path/'bad')])
    assert not (tmp_path/'bad.calls.tsv').exists()


def test_empty_and_all_unsupported():
    from vdjmatch.aggregate.paired_labels import paired_labels
    e,s=inputs()
    for n in [0,1]:
        selection=s.tail(n) if n else s.head(0)
        calls,scores=paired_labels(e.head(0),selection,species='MusMusculus')
        assert calls.height==scores.height==n
        assert calls['paired_label_source'].to_list()==['none']*n
        assert scores['supported'].to_list()==[False]*n
