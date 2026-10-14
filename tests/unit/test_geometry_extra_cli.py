"""Fixed advanced recipes are resolved before any raw input or V table is loaded."""
import argparse

import pytest

from vdjmatch.cli import neighbours as geometry


@pytest.mark.parametrize('extra', geometry.GEOMETRY_EXTRAS)
def test_geometry_extra_fixed_recipe(extra):
    parser = argparse.ArgumentParser()
    geometry.register(parser.add_subparsers(), name='annotate')
    args = parser.parse_args(['annotate', 'unused', '--vdjdb', 'unused',
        '--locus', 'TRB', '--output-prefix', 'unused', '--selection-only',
        '--geometry-extra', extra])
    geometry._apply_geometry_extra(args)
    assert args.distance_matrix == ('clipped' if extra in
        ('framev', 'framemid', 'framej', 'freegap') else 'blosum62')
    assert args.junction_ends == 'trim3'
    assert args.stem_depth == (4 if extra == 'stem4j' else 0)
    assert args.stem_weight == 1
    args.selection_only = False
    with pytest.raises(ValueError, match='compact uniform'):
        geometry._apply_geometry_extra(args)


@pytest.mark.parametrize('extra', geometry.GEOMETRY_EXTRAS)
def test_mixed_metadata_cli_fixed_geometry(tmp_path, extra):
    import json
    import polars as pl
    from test_mixed_kernel_union import inputs
    from vdjmatch.cli.__main__ import main
    paths,reference=inputs(tmp_path)
    refs=pl.read_csv(reference,separator='\t').with_columns(
        (pl.col('gene')+pl.lit('J1')).alias('j'))
    refs.write_csv(reference,separator='\t')
    metadata=tmp_path/'metadata.tsv'
    metadata.write_text('path\tsample_id\n'+''.join(str(p)+'\t'+p.stem+'\n' for p in paths))
    prefix=tmp_path/'extra'
    assert main(['annotate','--metadata',str(metadata),'--vdjdb',str(reference),
        '--geometry-extra',extra,'--threads','2','-o',str(prefix)])==0
    manifest=json.loads(prefix.with_suffix('.manifest.json').read_text())
    assert manifest['mixed_chain_union'] and manifest['population_service_calls']==3
    assert len(manifest['core_manifests'])==3
    assert all(m['geometry_extra']==extra and m['distance_scale']['matrix']==
        ('clipped' if extra in ('framev','framemid','framej','freegap') else 'blosum62')
        for m in manifest['core_manifests'])
    queries=pl.read_csv(prefix.with_suffix('.queries.tsv'),separator='\t')
    assert queries.height==217  # 73 standalone observations/chain +71 source pairs.
    assert queries['query_id'].n_unique()==217
    if extra in geometry.J_GEOMETRY_EXTRAS:
        assert all('literal J families' in m['background_enrichment']['background']['observation_unit']
            for m in manifest['core_manifests'])


@pytest.mark.parametrize('extra,locus,field', [('jw','TRB','additional_features'),
    ('framev','auto','gap_prior'),('sw','paired','alignment_mode')])
def test_native_capability_rejected_before_cohort_loading(monkeypatch,extra,locus,field):
    from types import SimpleNamespace
    from seqtree import gapblock
    from vdjmatch.runner import annotation
    capabilities={name:True for name in ('additional_features','gap_prior','score_quantum','alignment_mode')}
    capabilities.pop(field)
    monkeypatch.setattr(gapblock,'_GroupedGeometry',SimpleNamespace(**capabilities))
    def forbidden(*args,**kwargs): raise AssertionError('input loading preceded capability validation')
    monkeypatch.setattr(annotation,'_input_rows',forbidden)
    with pytest.raises(RuntimeError,match='Python/native'):
        annotation.annotate_repertoires('missing',reference='missing',threads=2,locus=locus,
            geometry_extra=extra)
