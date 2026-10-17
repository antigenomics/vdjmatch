"""One immutable V geometry preparation per mixed cohort invocation."""
from test_mixed_kernel_union import inputs
from vdjmatch.runner.annotation import annotate_repertoires
import json
import pytest


def test_mixed_cdr12_precomputes_once(tmp_path, monkeypatch):
    from vdjmatch.match import geometry
    paths, reference = inputs(tmp_path)
    prepare = geometry.precompute_v_loop_distances
    calls = []
    def recorded(**kwargs):
        calls.append(kwargs)
        return prepare(**kwargs)
    monkeypatch.setattr(geometry, 'precompute_v_loop_distances', recorded)
    result = annotate_repertoires(paths, reference=reference, v_regions='cdr12',
        output_prefix=tmp_path/'cdr12', threads=2)
    assert calls == [dict(species='human', alphabet='full', matrix='blosum62', v_regions='cdr12')]
    assert all(m['distance_scale']['v_regions'] == 'cdr12'
        for m in result['manifest']['core_manifests'])
    assert result['queries']['query_id'].n_unique() == result['queries'].height


@pytest.mark.parametrize('matrix,profile', [('pam250','uniform'),('mj','uniform'),
    ('vdjam','uniform'),('clipped','significance-trimmed'),('blosum62','significance-trimmed')])
def test_retained_geometry_cli_routes(tmp_path, matrix, profile):
    from vdjmatch.cli.__main__ import main
    paths, reference = inputs(tmp_path)
    prefix=tmp_path/'extra'
    metadata=tmp_path/'metadata.tsv'
    metadata.write_text('path\tsample_id\n'+''.join(str(p)+'\t'+p.stem+'\n' for p in paths))
    assert main(['annotate','--metadata',str(metadata),'--vdjdb',str(reference),
        '--distance-matrix',matrix,'--position-weighting',profile,'--stem-depth','0',
        '--junction-ends','trim3','--threads','2','-o',str(prefix)]) == 0
    manifest=json.loads(prefix.with_suffix('.manifest.json').read_text())
    assert manifest['mixed_chain_union'] and manifest['population_service_calls']==3
    assert all(m['parameters']['position_weighting']==profile and
        m['distance_scale']['matrix']==matrix for m in manifest['core_manifests'])


def test_trimmed_significance_details_reject_before_input_loading(tmp_path):
    with pytest.raises(ValueError, match='compact kernel'):
        annotate_repertoires(inputs(tmp_path)[0][0],reference='missing',selection_only=False,
            position_weighting='significance-trimmed',output_prefix=tmp_path/'out')
