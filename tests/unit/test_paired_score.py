"""Continuous linked-first evidence keeps native components and marginal rescue."""
import json
import polars as pl
import pytest


def sample():
    return pl.DataFrame(dict(query_id=['01','02','03','04','05'],
        kernel_score=['0.2',None,'0',None,None],kernel_alpha_score=['0.9','0.7','0.4',None,None],
        kernel_beta_score=['0.3','0.5','0.8',None,None],kernel_any_axis_score=['0.9','0.7','0.8',None,None],
        pair_id=['001','002','003','004','005']))


def test_candidate_preserves_native_values_and_missing_linked_rescue():
    from vdjmatch.cli.paired_score import linked_first_scores
    native=sample();out=linked_first_scores(native)
    assert out.select(native.columns).equals(native)
    assert out['score'].to_list()==[.2,.7,.8,None,None]
    assert out['score_source'].to_list()==['linked','marginal_fallback','marginal_fallback','none','none']
    assert out['supported'].to_list()==[True,True,True,False,False]
    for change in [native.with_columns(pl.lit('0.4').alias('kernel_any_axis_score')),
                   native.with_columns(pl.lit('inf').alias('kernel_score')),
                   native.with_columns(pl.lit('bad').alias('kernel_score')),
                   pl.concat([native,native.head(1)]),native.with_columns(pl.lit('old').alias('score'))]:
        with pytest.raises(ValueError):linked_first_scores(change)


def test_cli_reducer_is_fresh_and_calibration_remains_unbound(tmp_path, capsys):
    from vdjmatch.cli.__main__ import main
    from vdjmatch.db.cache import sha256
    selection,manifest=tmp_path/'selection.tsv',tmp_path/'native.json'
    sample().write_csv(selection,separator='\t')
    options=dict(locus='paired',background='vdjdb-other',neighbour_weighting='kernel',stem_depth=4,kernel_scorer='v-bounded',all_targets=False,targets_from_sample=True)
    manifest.write_text(json.dumps(dict(parameters=options, query_pairs=5, selection_sha256=sha256(selection))))
    args=['paired-score',str(selection),'--input-manifest',str(manifest),'--output-prefix',str(tmp_path/'candidate')]
    with pytest.raises(SystemExit):main(args+['--threads','-1'])
    assert main(args)==0
    out=pl.read_csv(tmp_path/'candidate.scores.tsv',separator='\t',schema_overrides={'query_id':pl.String,'pair_id':pl.String})
    assert out['query_id'].to_list()==sample()['query_id'].to_list()
    assert out['pair_id'].to_list()==sample()['pair_id'].to_list()
    report=json.loads((tmp_path/'candidate.manifest.json').read_text())
    assert report['calibration_transfer']=='none' and report['default_changed'] is False and report['query_rows']==5
    with pytest.raises(FileExistsError):main(args)
    sample().with_columns(pl.col('query_id').reverse()).write_csv(selection, separator='\t')
    with pytest.raises(SystemExit):
        main(args[:-1] + [str(tmp_path/'permuted')])
    assert 'selection file differs from its producer manifest' in capsys.readouterr().err
    assert not list(tmp_path.glob('permuted.*'))
    sample().head(4).write_csv(selection, separator='\t')
    truncated = tmp_path / 'truncated'
    with pytest.raises(SystemExit) as exc:
        main(args[:-1] + [str(truncated)])
    assert exc.value.code == 2
    assert 'row count disagrees with its manifest' in capsys.readouterr().err
    assert not list(tmp_path.glob('truncated.*'))
    options['locus']='TRB';manifest.write_text(json.dumps(dict(parameters=options)))
    with pytest.raises(SystemExit) as exc:main(args[:-1]+[str(tmp_path/'invalid')])
    assert exc.value.code==2


def test_conditional_filter_requires_unique_winner_and_resolved_null_tail():
    from vdjmatch.cli.paired_score import conditional_linked_calls, linked_first_scores

    native = sample().with_columns(
        *[pl.lit('A').alias(c) for c in ('selected_linked_epitope', 'selected_alpha_epitope', 'selected_beta_epitope')],
        *[pl.lit('ranked_evidence').alias(c) for c in ('linked_selection_status', 'alpha_selection_status', 'beta_selection_status')],
        pl.lit(1).alias('tied_peptides'),
    )
    null = native.head(1).drop('query_id').join(pl.DataFrame({'query_id': range(99)}), how='cross')
    null = null.with_columns(pl.lit('0.1').alias('kernel_score'))
    result = conditional_linked_calls(linked_first_scores(native), null)
    assert result['call_status'].to_list() == ['accepted'] + ['no_supported_match'] * 4
    assert result['p_null'][0] == .01 and result['E_null'][0] == 0
    for more in [null.head(98), null.with_columns(pl.lit('0.2').alias('kernel_score'))]:
        assert conditional_linked_calls(linked_first_scores(native), more)['call_status'][0] == 'no_supported_match'
    conflicting = native.with_columns(pl.lit('ambiguous').alias('linked_selection_status'))
    unresolved = conditional_linked_calls(linked_first_scores(conflicting), null)
    assert unresolved.select('unique_winner', 'null_queries', 'p_null', 'E_null', 'call_status').row(0) == (
        False, None, None, None, 'no_supported_match')

    ambiguous = null.tail(1).with_columns(pl.lit('ambiguous').alias('linked_selection_status'))
    ranks = conditional_linked_calls(linked_first_scores(native), pl.concat([null.head(98), ambiguous]))
    assert ranks['null_queries'][0] == 98
    assert ranks['p_null'][0] == pytest.approx(1 / 99)
    assert ranks['call_status'][0] == 'no_supported_match'
    no_group = conditional_linked_calls(linked_first_scores(native), ambiguous)
    assert no_group.height == native.height
    assert no_group['p_null'].null_count() == native.height
    assert no_group['call_status'].to_list() == ['no_supported_match'] * native.height
