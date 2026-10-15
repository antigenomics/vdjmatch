"""Fit/apply a frozen global evidence formula with explicit training groups."""
from pathlib import Path


def register(subparsers):
    for command in ('fit-evidence', 'combine-evidence'):
        p = subparsers.add_parser(command, help=__doc__)
        p.add_argument('evidence', help='TSV of component scores and reference structure')
        p.add_argument('--model', required=command == 'fit-evidence')
        p.add_argument('--group-column', default='epitope', help='held-out group identity')
        p.add_argument('--threads',type=int,default=1,help='positive table thread budget (default: 1); BLAS stays serial')
        if command == 'fit-evidence':
            p.add_argument('--label-column', default='label')
            p.add_argument('--weight-column')
        else:
            p.add_argument('--null-evidence',help='frozen full-pMHC/species/locus component null scores; normalize scales before selecting stronger evidence')
            p.add_argument('--structure-conditioned',action='store_true',
                           help='with --null-evidence, discount rank selection by the reference sparsity gate')
            p.add_argument('--three-channel',action='store_true',
                           help='experimental: normalize original/gapblock/stem4 scores and sum their explicit reference-derived weights')
            p.add_argument('--fixed', action='store_true', help='apply the unfitted two-score sparsity formula')
            p.add_argument('--output', required=True)
            p.add_argument('--require-heldout', action='store_true',
                           help='benchmark: reject any overlap with model training groups')
        p.set_defaults(func=main, fit=command == 'fit-evidence')


def main(a):
    # Load both NumPy and SciPy BLAS before limiting an already initialized runtime.
    from ..aggregate import combination
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        return _main(a)


def _main(a):
    import json
    import numpy as np
    import polars as pl
    from ..aggregate.combination import fit_evidence_combination, combine_evidence, fixed_evidence_combination
    from ..db.cache import sha256

    null_evidence=getattr(a,'null_evidence',None)
    structure_conditioned=getattr(a,'structure_conditioned',False)
    three_channel=getattr(a,'three_channel',False)
    if structure_conditioned and not null_evidence:
        raise ValueError('--structure-conditioned requires --null-evidence')
    if three_channel and (not null_evidence or structure_conditioned):
        raise ValueError('--three-channel requires --null-evidence and explicit weights; no two-channel gate')
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
        from ..aggregate.combination import null_normalized_evidence, THREE_COMPONENTS, THREE_WEIGHTS, weighted_component_support
        scored=null_normalized_evidence(frame,pl.read_csv(null_evidence,separator='\t',infer_schema_length=None),
                                       structure_conditioned=structure_conditioned,three_channel=three_channel)
        output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
        supported = weighted_component_support(scored) if three_channel else pl.col('normalized_score').is_not_null()
        exported = scored.with_columns(pl.when(supported).then(pl.col('normalized_score')).otherwise(None).alias('combined_score'))
        exported = exported.with_columns(pl.when(pl.col('combined_score').is_not_null()).then(pl.lit('scored'))
            .otherwise(pl.lit('unsupported')).alias('combined_status'))
        exported.write_csv(output, separator='\t')
        Path(str(output)+'.manifest.json').write_text(json.dumps({
            'evidence_sha256':evidence_hash,'null_evidence_sha256':sha256(Path(null_evidence)),
            'rows':frame.height,'unsupported_rows':exported['combined_score'].null_count(),
            'support_policy':'original_supported plus positive weighted kernels; raw normalized_score retained' if three_channel and 'original_supported' in frame.columns else 'legacy component availability',
            'formula':'sum(w_k * -log(p_k)); missing components contribute zero without weight redistribution' if three_channel else 'max(-log(p_L),-log(p_R)+log(s)); missing L uses -log(p_R)' if structure_conditioned else '-log(min(p_component_likelihood,p_component_rank))',
            'three_channel':three_channel,
            'components':list(THREE_COMPONENTS) if three_channel else ['likelihood_score','conditional_rank_score'],
            'weight_columns':list(THREE_WEIGHTS) if three_channel else None,
            'structure_conditioned':structure_conditioned,
            'structure_gate':'s=1-reference_coverage*(1-(1-knee_shape)*knee_radius_ratio)' if structure_conditioned else None,
            'normalization':'max source (k+1)/(N+1), exact ties and no-hit null queries retained',
            'population_keys':['species','locus','epitope','mhc_a','mhc_b','mhc_class'] +
                (['normalization_stratum'] if 'normalization_stratum' in frame.columns else []),
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
