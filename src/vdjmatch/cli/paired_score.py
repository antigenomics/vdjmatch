"""Export an uncalibrated linked-first continuous candidate from native selections."""
from pathlib import Path


def register(subparsers):
    parser = subparsers.add_parser('paired-score', help=__doc__)
    parser.add_argument('selection', help='native paired selection TSV, optionally gzip compressed')
    parser.add_argument('--input-manifest', required=True)
    parser.add_argument('--output-prefix', required=True)
    parser.add_argument('--threads', type=int, default=1, help='positive table thread budget (default: 1)')
    parser.set_defaults(func=main)


def linked_first_scores(selection):
    """Retain positive linked evidence, with native marginal rescue if absent.

    Preserve every original row, identity and component. This is continuous
    ranking evidence, distinct from linked-first label decision indicators and
    from the calibrated native any-axis detection statistic.
    """
    import polars as pl
    columns = ['kernel_score', 'kernel_alpha_score', 'kernel_beta_score']
    required = {'query_id', 'kernel_any_axis_score', *columns}
    if not required <= set(selection.columns) or set(selection.columns) & {'score', 'scorer', 'score_source', 'supported'}:
        raise ValueError('native paired components and unreserved score columns required')
    if (selection['query_id'].null_count() or selection['query_id'].n_unique() != selection.height
            or selection.filter(pl.col('query_id').cast(pl.String).str.strip_chars() == '').height):
        raise ValueError('unique nonblank native query IDs required')
    numeric = selection.select(pl.col(c).cast(pl.Float64, strict=False) for c in [*columns, 'kernel_any_axis_score'])
    if any((selection[c].is_not_null() & numeric[c].is_null()).any() for c in numeric.columns):
        raise ValueError('numeric native scores required')
    if numeric.filter(pl.any_horizontal(pl.col(c).is_not_null() &
            (~pl.col(c).is_finite() | (pl.col(c) < 0)) for c in numeric.columns)).height:
        raise ValueError('finite nonnegative native scores required')
    if not numeric['kernel_any_axis_score'].eq_missing(numeric.select(pl.max_horizontal(columns)).to_series()).all():
        raise ValueError('native any-axis score disagrees with its components')
    linked = pl.col('kernel_score').cast(pl.Float64) > 0
    any_axis = pl.col('kernel_any_axis_score').cast(pl.Float64)
    supported = (any_axis > 0).fill_null(False)
    return selection.with_columns(
        pl.when(linked).then(pl.col('kernel_score').cast(pl.Float64))
            .when(supported).then(any_axis).otherwise(None).alias('score'),
        pl.when(linked).then(pl.lit('linked')).when(supported)
            .then(pl.lit('marginal_fallback')).otherwise(pl.lit('none')).alias('score_source'),
        pl.lit('stem4_linked_first_continuous').alias('scorer'), supported.alias('supported'))


def main(args):
    import json
    import polars as pl
    from ..db.cache import sha256
    from ..runner.annotation import _validate_output_paths

    paths = [Path(args.selection), Path(args.input_manifest)]
    prefix = _validate_output_paths(args.output_prefix, paths)
    outputs = [Path(str(prefix)+suffix) for suffix in ('.scores.tsv', '.manifest.json')]
    if any(p.exists() or p.is_symlink() for p in outputs):
        raise FileExistsError('fresh paired-score output paths required')
    before = {str(p): sha256(p) for p in paths}
    manifest = json.loads(paths[1].read_text())
    required = dict(locus='paired', background='vdjdb-other', neighbour_weighting='kernel',
                    stem_depth=4, kernel_scorer='v-bounded')
    if any(manifest.get('parameters', {}).get(k) != v for k,v in required.items()):
        raise ValueError('native paired stem4 v-bounded VDJdb-other selection required')
    native = pl.read_csv(paths[0], separator='\t', infer_schema=False)
    scores = linked_first_scores(native)
    if before != {str(p): sha256(p) for p in paths}:
        raise ValueError('paired-score inputs changed during reduction')
    prefix.parent.mkdir(parents=True, exist_ok=True)
    scores.write_csv(outputs[0], separator='\t')
    outputs[1].write_text(json.dumps(dict(inputs_sha256=before, source_sha256=sha256(Path(__file__)),
        outputs_sha256={outputs[0].name: sha256(outputs[0])}, query_rows=scores.height,
        definition='positive finite linked score; otherwise existing any-axis score; no supported axis: explicit no-score',
        score_kind='uncalibrated continuous ranking candidate; not a label indicator or binding probability',
        original_components_unchanged=True, labels_changed=False, default_changed=False,
        calibration_transfer='none', calibration_status='requires its own complete selected-event fit'), indent=2)+'\n')
    return 0
