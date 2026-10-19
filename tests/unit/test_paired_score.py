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


def test_cli_reducer_is_fresh_and_calibration_remains_unbound(tmp_path):
    from vdjmatch.cli.__main__ import main
    selection,manifest=tmp_path/'selection.tsv',tmp_path/'native.json'
    sample().write_csv(selection,separator='\t')
    options=dict(locus='paired',background='vdjdb-other',neighbour_weighting='kernel',stem_depth=4,kernel_scorer='v-bounded',all_targets=False,targets_from_sample=True)
    manifest.write_text(json.dumps(dict(parameters=options)))
    args=['paired-score',str(selection),'--input-manifest',str(manifest),'--output-prefix',str(tmp_path/'candidate')]
    with pytest.raises(SystemExit):main(args+['--threads','0'])
    assert main(args)==0
    out=pl.read_csv(tmp_path/'candidate.scores.tsv',separator='\t',schema_overrides={'query_id':pl.String,'pair_id':pl.String})
    assert out['query_id'].to_list()==sample()['query_id'].to_list()
    assert out['pair_id'].to_list()==sample()['pair_id'].to_list()
    report=json.loads((tmp_path/'candidate.manifest.json').read_text())
    assert report['calibration_transfer']=='none' and report['default_changed'] is False and report['query_rows']==5
    with pytest.raises(FileExistsError):main(args)
    options['locus']='TRB';manifest.write_text(json.dumps(dict(parameters=options)))
    with pytest.raises(SystemExit) as exc:main(args[:-1]+[str(tmp_path/'invalid')])
    assert exc.value.code==2
