"""Apply a frozen reference-only class policy to continuous component scores.

This experimental reducer does not fit a policy, score receptors, calibrate
calls or alter the annotation default. Its inputs are final component tables
and a separately frozen reference-model policy.
"""
from pathlib import Path

KEY = ['species', 'locus', 'epitope']
ROUTES = {'V_dominant': 'rule_v_score',
          'junction_projection_dominant': 'rule_rank_score',
          'general_fallback': 'original'}
OUTPUT = ['reference_class', 'selected_score_column', 'policy_score', 'policy_status']


def apply_reference_policy(scores, policy, *, original_column='original_score'):
    """Select existing components without using labels or dropping query rows."""
    import polars as pl
    required = {*KEY, original_column, 'rule_v_score', 'rule_rank_score'}
    if not required <= set(scores.columns) or not {*KEY, 'reference_class'} <= set(policy.columns):
        raise ValueError('component scores and policy require species/locus/epitope and declared score/class columns')
    if set(OUTPUT) & set(scores.columns) or any(c.startswith('_policy_') for c in scores.columns):
        raise ValueError('score input contains reserved policy output columns')
    if policy.select(KEY).is_duplicated().any() or policy.select(pl.any_horizontal(pl.col(KEY).is_null()).any()).item():
        raise ValueError('policy identities must be unique and non-null')
    if policy['reference_class'].null_count() or not set(policy['reference_class']) <= set(ROUTES):
        raise ValueError('unknown or missing reference policy class')
    join_keys = ['_policy_'+k for k in KEY]
    aliases = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}
    def keys(frame):
        return frame.with_columns(*[(pl.col(k).replace(aliases) if k=='species' else pl.col(k))
            .alias(j) for k,j in zip(KEY,join_keys)])
    selected = keys(policy.select(*KEY,'reference_class')).select(*join_keys,'reference_class')
    # Alias-normalised identities must remain injective as well.
    if selected.select(join_keys).is_duplicated().any():
        raise ValueError('policy species aliases collide')
    out = keys(scores).join(selected,on=join_keys,how='left',validate='m:1',maintain_order='left')
    desired = pl.col('reference_class').replace_strict(ROUTES,default='original')
    v = pl.col('rule_v_score').cast(pl.Float64)
    r = pl.col('rule_rank_score').cast(pl.Float64)
    original = pl.col(original_column).cast(pl.Float64)
    chosen = pl.when(desired=='rule_v_score').then(v).when(desired=='rule_rank_score').then(r).otherwise(original)
    usable = chosen.is_finite().fill_null(False)
    out = out.with_columns(
        pl.when(usable).then(chosen).alias('policy_score'),
        pl.when(desired!='original').then(desired).otherwise(pl.lit(original_column)).alias('selected_score_column'),
        pl.when(~usable & ~original.is_finite().fill_null(False)).then(pl.lit('no_finite_score'))
          .when(pl.col('reference_class').is_null()).then(pl.lit('no_reference_profile'))
          .when(desired=='original').then(pl.lit('reference_fallback'))
          .when(~usable).then(pl.lit('component_unavailable_no_score')).otherwise(pl.lit('selected_component')).alias('policy_status'),
        pl.col('reference_class').fill_null('general_fallback'))
    return out.drop(join_keys)


def register(subparsers):
    p = subparsers.add_parser('score-policy',help=__doc__.splitlines()[0])
    p.add_argument('scores',help='flat final component-score TSV/Parquet')
    p.add_argument('--policy',required=True,help='frozen reference-model reference-policy.tsv')
    p.add_argument('--original-column',default='original_score')
    p.add_argument('--output-prefix',required=True)
    p.set_defaults(func=main)


def main(a):
    import json
    from ..io.airr import _read_table
    from ..db.cache import sha256
    from ..runner.annotation import _validate_output_paths
    model = Path(a.policy).resolve()
    companion = model.parent/'manifest.json'
    info = json.loads(companion.read_text())
    if info.get('reference_model_asset') is not True or info.get('no_benchmark_inputs') is not True or info.get('outputs_sha256',{}).get(model.name)!=sha256(model):
        raise ValueError('policy must match a frozen reference-only model manifest')
    prefix = _validate_output_paths(a.output_prefix,[a.scores,model,companion])
    before = {str(Path(p).resolve()):sha256(Path(p)) for p in (a.scores,model,companion)}
    scores = _read_table(a.scores)
    result = apply_reference_policy(scores,_read_table(model),original_column=a.original_column)
    if any(sha256(Path(p))!=h for p,h in before.items()):
        raise ValueError('policy or score input changed during selection')
    prefix.parent.mkdir(parents=True,exist_ok=True)
    output = Path(str(prefix)+'.scores.tsv')
    result.write_csv(output,separator='\t')
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(dict(
        input_sha256=before,output_sha256={output.name:sha256(output)},rows=result.height,
        original_column=a.original_column,experimental=True,default_changed=False,
        labels_used=False,calibration_status='raw continuous selection only; component P/E/posteriors do not transfer',
        policy='frozen reference classes select V or complete conditional rank; otherwise original'),indent=2)+'\n')
    return 0
