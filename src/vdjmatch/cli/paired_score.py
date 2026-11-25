"""Export an uncalibrated linked-first continuous candidate from native selections."""
from pathlib import Path


from . import cpu_budget

def register(subparsers):
    parser = subparsers.add_parser('paired-score', help=__doc__)
    parser.add_argument('selection', help='native paired selection TSV, optionally gzip compressed')
    parser.add_argument('--input-manifest', required=True)
    parser.add_argument('--withheld-manifest', type=Path, help='experimental winner-conditional rank filter using sibling reference-calls .calls.tsv')
    parser.add_argument('--output-prefix', required=True)
    parser.add_argument('--threads', type=cpu_budget, default="0", help='table CPU budget (0: automatic, all allocated physical cores)')
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


def conditional_linked_calls(scores, withheld):
    """Rank a unique linked-first winner against null winners of its peptide and axis."""
    import polars as pl
    from ..evalue.directional import null_tail_probabilities

    def groups(frame):
        linked = pl.col('score_source') == 'linked'
        alpha = pl.col('kernel_alpha_score').cast(pl.Float64).fill_null(0)
        beta = pl.col('kernel_beta_score').cast(pl.Float64).fill_null(0)
        def select(ab, a, b):
            return pl.when(linked).then(ab).when(alpha >= beta).then(a).otherwise(b)
        ambiguous_fallback = (~linked & (alpha == beta) &
            (~pl.col('selected_alpha_epitope').eq_missing(pl.col('selected_beta_epitope')) |
             (pl.col('beta_tied_peptides').cast(pl.Int64) != 1)))
        return frame.with_columns(
            select(pl.lit('linked'), pl.lit('alpha'), pl.lit('beta')).alias('winning_axis'),
            select('selected_linked_epitope', 'selected_alpha_epitope', 'selected_beta_epitope').alias('winning_peptide'),
            (select('tied_peptides', 'alpha_tied_peptides', 'beta_tied_peptides').cast(pl.Int64).eq(1)
             & ~ambiguous_fallback).alias('unique_winner'))

    keys = ['winning_axis', 'winning_peptide']
    observed = groups(scores)
    nulls = groups(linked_first_scores(withheld)).select('query_id', 'score', *keys,
        pl.lit('null').alias('outcome'), pl.lit('withheld').alias('null_source'))
    depths = nulls.group_by(keys).len(name='null_queries')
    observed = observed.join(depths, on=keys, how='left', validate='m:1', maintain_order='left')
    eligible = observed.filter(pl.col('supported') & pl.col('unique_winner') & pl.col('null_queries').is_not_null())
    tails, _ = null_tail_probabilities(eligible, nulls, group_columns=keys)
    return observed.join(tails.select('query_id', 'p_null', 'E_null'), on='query_id',
        how='left', validate='1:1', maintain_order='left').with_columns(
        pl.when(pl.col('p_null') <= 0.01).then(pl.lit('accepted'))
          .otherwise(pl.lit('no_supported_match')).alias('call_status'))


def main(args):
    import json
    import polars as pl
    from ..db.cache import sha256
    from ..runner.annotation import _validate_output_paths

    paths = [Path(args.selection), Path(args.input_manifest)]
    withheld_manifest = getattr(args, 'withheld_manifest', None)
    if withheld_manifest is not None:
        if not withheld_manifest.name.endswith('.manifest.json'):
            raise ValueError('expected a reference-calls .manifest.json')
        withheld_calls = withheld_manifest.with_name(withheld_manifest.name.removesuffix('.manifest.json')+'.calls.tsv')
        paths.extend([withheld_manifest, withheld_calls])
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
    conditional = {}
    if withheld_manifest is not None:
        from .calibration import pipeline_binding
        null_manifest = json.loads(withheld_manifest.read_text())
        if (not null_manifest.get('options', {}).get('leave_epitope_out')
                or pipeline_binding(manifest) != pipeline_binding(null_manifest)):
            raise ValueError('a matching peptide-withheld reference-calls population is required')
        withheld = pl.read_csv(withheld_calls, separator='\t', infer_schema=False)
        if withheld.height != null_manifest['query_pairs'] or native.height != manifest['query_pairs']:
            raise ValueError('paired score input row count disagrees with its manifest')
        scores = conditional_linked_calls(scores, withheld)
        conditional = dict(calibration_status='experimental conditional ranks against peptide-withheld reference winners',
            conditioning=['winning_axis','winning_peptide'], cutoff='inclusive p_null <= 0.01; unique supported winners only',
            p_null='(1+k)/(1+n)', E_null='k/n; empirical tail conditional on the winning peptide and axis')
    if before != {str(p): sha256(p) for p in paths}:
        raise ValueError('paired-score inputs changed during reduction')
    prefix.parent.mkdir(parents=True, exist_ok=True)
    scores.write_csv(outputs[0], separator='\t')
    outputs[1].write_text(json.dumps(dict(inputs_sha256=before, source_sha256=sha256(Path(__file__)),
        outputs_sha256={outputs[0].name: sha256(outputs[0])}, query_rows=scores.height,
        definition='positive finite linked score; otherwise existing any-axis score; no supported axis: explicit no-score',
        score_kind='uncalibrated continuous ranking candidate; not a label indicator or binding probability',
        original_components_unchanged=True, labels_changed=False, default_changed=False,
        calibration_transfer='none', calibration_status='requires its own complete selected-event fit') | conditional, indent=2)+'\n')
    return 0
