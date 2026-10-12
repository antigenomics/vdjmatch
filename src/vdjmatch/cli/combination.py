"""Fit/apply a frozen global evidence formula with explicit training groups."""
from pathlib import Path


def register(subparsers):
    for command in ('fit-evidence', 'combine-evidence'):
        p = subparsers.add_parser(command, help=__doc__)
        p.add_argument('evidence', help='TSV of component scores and reference structure')
        p.add_argument('--model', required=command == 'fit-evidence')
        p.add_argument('--group-column', default='epitope', help='held-out group identity')
        if command == 'fit-evidence':
            p.add_argument('--label-column', default='label')
            p.add_argument('--weight-column')
        else:
            p.add_argument('--null-evidence',help='frozen full-pMHC/species/locus component null scores; normalize scales before selecting stronger evidence')
            p.add_argument('--fixed', action='store_true', help='apply the unfitted two-score sparsity formula')
            p.add_argument('--output', required=True)
            p.add_argument('--require-heldout', action='store_true',
                           help='benchmark: reject any overlap with model training groups')
        p.set_defaults(func=main, fit=command == 'fit-evidence')


def main(a):
    import json
    import numpy as np
    import polars as pl
    from ..aggregate.combination import fit_evidence_combination, combine_evidence, fixed_evidence_combination
    from ..db.cache import sha256

    null_evidence=getattr(a,'null_evidence',None)
    if not a.fit and sum(map(bool,[a.model,a.fixed,null_evidence]))!=1:
        raise ValueError('choose exactly one of --model, --fixed and --null-evidence')
    paths = [Path(a.evidence)] + ([Path(a.model)] if a.model else []) + ([Path(null_evidence)] if null_evidence else [])
    if not a.fit:
        paths += [Path(a.output), Path(str(a.output)+'.manifest.json')]
    for i, left in enumerate(paths):
        for right in paths[i+1:]:
            if left.resolve() == right.resolve() or (
                    left.exists() and right.exists() and left.samefile(right)):
                raise ValueError('evidence, model, output and manifest paths must be distinct')
    evidence_hash = sha256(Path(a.evidence))
    frame = pl.read_csv(a.evidence, separator='\t', infer_schema_length=None)
    if not frame.height:
        raise ValueError('nonempty evidence table required')
    if null_evidence:
        from ..aggregate.combination import null_normalized_evidence
        scored=null_normalized_evidence(frame,pl.read_csv(null_evidence,separator='\t',infer_schema_length=None))
        output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
        scored.with_columns(pl.col('normalized_score').alias('combined_score'),
            pl.when(pl.col('normalized_score').is_not_null()).then(pl.lit('scored')).otherwise(pl.lit('unsupported')).alias('combined_status')).write_csv(output,separator='\t')
        Path(str(output)+'.manifest.json').write_text(json.dumps({
            'evidence_sha256':evidence_hash,'null_evidence_sha256':sha256(Path(null_evidence)),
            'rows':frame.height,'unsupported_rows':scored['normalized_score'].null_count(),
            'formula':'-log(min(p_component_likelihood,p_component_rank))',
            'normalization':'max source (k+1)/(N+1), exact ties and no-hit null queries retained',
            'population_keys':['species','locus','epitope','mhc_a','mhc_b','mhc_class'],
            'query_labels_used':False,'interpretation':'normalized ranking statistic; final selected-score P/E requires complete-pipeline calibration'},indent=2)+'\n')
        return 0
    if a.group_column not in frame.columns or frame[a.group_column].null_count():
        raise ValueError('explicit nonmissing training/evaluation groups required')
    groups = sorted(set(frame[a.group_column].cast(pl.String).to_list()))
    if a.fit:
        model = fit_evidence_combination(frame, frame[a.label_column].to_numpy(),
            weights=frame[a.weight_column].to_numpy() if a.weight_column else None)
        model.update(training_groups=groups, group_column=a.group_column,
            training_sha256=evidence_hash, label_column=a.label_column,
            weight_column=a.weight_column)
        path = Path(a.model)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(model, indent=2)+'\n')
    else:
        model_hash = sha256(Path(a.model)) if a.model else None
        model = json.loads(Path(a.model).read_text()) if a.model else None
        overlap = []
        fixed_score = fixed_evidence_combination(frame)
        if model:
            if model.get('group_column') != a.group_column or 'training_groups' not in model:
                raise ValueError('model must record matching training group identities')
            overlap = sorted(set(groups).intersection(model['training_groups']))
            if a.require_heldout and overlap:
                raise ValueError('evaluation groups overlap model training groups')
        score = combine_evidence(frame, model) if model else fixed_score
        output = Path(a.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        frame.with_columns(pl.Series('combined_score', score, nan_to_null=True),
            pl.Series('fixed_score', fixed_score, nan_to_null=True),
            pl.Series('combined_status', np.where(np.isfinite(score), 'scored', 'unsupported'))
        ).write_csv(output, separator='\t')
        Path(str(output)+'.manifest.json').write_text(json.dumps({
            'evidence_sha256': evidence_hash, 'model_sha256': model_hash,
            'rows': frame.height, 'unsupported_rows': int((~np.isfinite(score)).sum()),
            'evaluation_groups': groups, 'group_column': a.group_column,
            'require_heldout': a.require_heldout, 'training_group_overlap': overlap,
            'formula': 'GLM' if model else 'fixed',
            'fixed_formula': 'L+s*log(B); retain the full original log term',
            'fixed_gate': 's=1-original_support*reference_coverage*(1-(1-knee_shape)*knee_radius_ratio)',
            'interpretation': model['interpretation'] if model else 'unfitted empirical log-score mixture; complete-pipeline P/E calibration required'}, indent=2)+'\n')
    return 0
