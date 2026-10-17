"""Fresh exhaustive single-chain reference calls and batched peptide-withholding nulls."""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('reference-calls', help='benchmark full-reference winning labels or peptide-withheld nulls')
    p.add_argument('--vdjdb', required=True)
    p.add_argument('--species', default='HomoSapiens')
    p.add_argument('--locus', choices=['TRA', 'TRB', 'paired'], required=True)
    p.add_argument('--scorer', choices=['original', 'stem4'], default='original')
    p.add_argument('--reference-population', choices=['all-valid', 'curated-benchmark'], default='curated-benchmark')
    p.add_argument('--control', help='matching raw repertoire; human TRB has a bundled input')
    p.add_argument('--scope', default='1,0,0,1')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--calibration-score', choices=['E', 'p_enrichment', 'score'], default='score',
                   help='optional recommendation statistic; never filters continuous outputs')
    p.add_argument('--min-target-records', type=int, default=30,
                   help='benchmark query targets only; the searchable reference remains complete')
    p.add_argument('--leave-epitope-out', action='store_true')
    p.add_argument('--exact-policy', choices=['exclude', 'retain'], default='exclude',
                   help='exclude query junctions for benchmarks; retain exact evidence only for ordinary-call calibration')
    p.add_argument('--output-prefix', required=True)
    p.set_defaults(func=main)
    return p


def main(args):
    import importlib.metadata
    import json
    import time
    import polars as pl
    from .. import db
    from ..api import Annotator, _prepare, _append_calls
    from ..aggregate.candidates import PMHC
    from ..evalue.control import raw_background
    from ..match.scope import search_params
    from ..match.scoring import load_vdjam, DEFAULT_SCALE

    started = time.monotonic()
    if args.threads < 1:
        raise ValueError('threads must be positive')
    if args.min_target_records < 1:
        raise ValueError('min-target-records must be positive')
    if args.scorer == 'stem4' and args.calibration_score != 'score':
        raise ValueError('stem4 calibration uses only its final kernel_score; E and p_enrichment are incompatible')
    if args.scorer == 'stem4':
        if args.control is not None:
            raise ValueError('stem4 uses the other-peptide reference background; --control is incompatible')
        if args.scope != '1,0,0,1':
            raise ValueError('stem4 uses fixed recorded geometry; nondefault --scope is incompatible')
        from ..evalue.control import _organism
        args.species = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[_organism(args.species)]
    if args.locus=='paired' and (args.scorer!='stem4' or not args.leave_epitope_out):
        raise ValueError('paired reference-calls requires --scorer stem4 --leave-epitope-out')
    reference = db.load(args.vdjdb, species=args.species, gene=None if args.locus=='paired' else args.locus)
    if args.reference_population == 'all-valid':
        reference = reference.with_columns(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).alias('reference_valid'))
    valid = reference.filter(pl.col('reference_valid'))
    if args.locus=='paired':
        return _stem4_paired(args,reference,valid,started)
    support = valid.group_by('epitope').len().rename({'len': 'reference_records'})
    eligible = support.filter(pl.col('reference_records') >= args.min_target_records)
    benchmark = valid.join(eligible.select('epitope'), on='epitope', how='semi')
    keys = ['cdr3', 'v', 'j']
    if args.leave_epitope_out:
        keys.append('epitope')
    raw = benchmark.select(keys).unique().sort(keys).with_row_index('query_id').with_columns(
        pl.lit(args.locus).alias('locus'))
    if args.scorer == 'stem4':
        return _stem4(args, reference, valid, support, eligible, benchmark, raw, keys, started)
    ann = Annotator.from_frame(valid)
    _, q = _prepare(raw)
    control, control_provenance = raw_background(args.locus, args.species, args.control)
    exclusions = raw.select('query_id', 'epitope') if args.leave_epitope_out else None
    params = search_params(args.scope, matrix=load_vdjam(), gap_open=DEFAULT_SCALE, gap_extend=DEFAULT_SCALE)
    _, candidates = ann._evidence(q, params, control=control,
        threads=args.threads, exclude_exact=getattr(args, 'exact_policy', 'exclude') == 'exclude', exclude_epitopes=exclusions)
    calls = _append_calls(raw, q, candidates, ann.loci, 'vdjmatch_')
    score = pl.when(pl.col('vdjmatch_status') == 'matched').then(
        pl.col('vdjmatch_' + args.calibration_score)).otherwise(None).alias('score')
    if args.leave_epitope_out:
        cases = calls.select('query_id', score,
            pl.lit('null').alias('outcome'), pl.lit('withheld').alias('null_source'))
    else:
        truth = raw.join(valid, on=['cdr3', 'v', 'j'], nulls_equal=True, validate='1:m').select(
            'query_id', *PMHC).unique()
        predicted = calls.select('query_id', *[pl.col('vdjmatch_'+x).alias(x) for x in PMHC])
        correct = predicted.join(truth, on=['query_id', *PMHC], how='semi', nulls_equal=True)['query_id']
        cases = calls.select('query_id', score,
            pl.when(pl.col('query_id').is_in(correct.implode()) & (pl.col('vdjmatch_status') == 'matched'))
            .then(pl.lit('correct')).otherwise(pl.lit('incorrect')).alias('outcome'),
            pl.lit(None, dtype=pl.String).alias('null_source'))
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    calls.write_csv(str(prefix)+'.calls.tsv', separator='\t')
    candidates.write_csv(str(prefix)+'.candidates.tsv', separator='\t')
    cases.write_csv(str(prefix)+'.calibration-cases.tsv', separator='\t')
    support.with_columns((pl.col('reference_records') >= args.min_target_records).alias('benchmark_eligible')).write_csv(
        str(prefix)+'.target-support.tsv', separator='\t')
    raw.join(valid, on=keys, nulls_equal=True, validate='1:m').select('query_id', *PMHC).unique().write_csv(
        str(prefix)+'.truth.tsv', separator='\t')
    provenance = db.provenance(args.vdjdb)
    manifest = {'reference': provenance, 'reference_sha256': provenance.get('input_sha256'),
        'control': control_provenance,
        'software': {x: importlib.metadata.version(x) for x in ['vdjmatch', 'seqtree', 'vdjtools']},
        'options': vars(args) | {'func': 'reference-calls'}, 'source_records': reference.height,
        'valid_reference_records': valid.height, 'invalid_reference_records': reference.height-valid.height,
        'eligible_benchmark_records': benchmark.height, 'ineligible_benchmark_records': valid.height-benchmark.height,
        'eligible_benchmark_peptides': eligible.height, 'all_reference_peptides': support.height,
        'target_support_rule': 'valid source records per species/locus/peptide, before self-exclusion; pool MHC restrictions',
        'search_reference': 'entire valid species/locus reference, including sparse peptides',
        'queries': raw.height, 'candidate_rows': candidates.height, 'wall_seconds': time.monotonic()-started,
        'self_exclusion': ('benchmark-only full-junction hit filtering; fixed model populations'
            if getattr(args, 'exact_policy', 'exclude') == 'exclude' else 'exact evidence retained for ordinary-call calibration'),
        'observation_unit': 'distinct full-junction/V/J'+('/truth peptide' if args.leave_epitope_out else ''),
        'peptide_withholding': 'all restrictions removed before reduction' if args.leave_epitope_out else None,
        'scoring': 'production default vdjam matrix, score_scale400, soft V weight0.25',
        'limitations': 'Fixed-ball enrichment is evidence, not a posterior for the winning label. No paired calibration.'}
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0


def _stem4(args, reference, valid, support, eligible, benchmark, raw, keys, started):
    """Recorded final scalar and selection, using one run-local grouped traversal."""
    import hashlib
    import importlib.metadata
    import json
    import time
    import polars as pl
    from .. import db
    from ..api import _prepare
    from ..aggregate.candidates import PMHC
    from ..evalue.kernel import grouped_kernel_selection, finish_kernel_selection, ladder_cutoffs
    from ..match.geometry import precompute_v_loop_distances, resolve_v_alleles

    species = {'HomoSapiens': 'human', 'MusMusculus': 'mouse'}.get(args.species)
    if species is None:
        raise ValueError('stem4 requires HomoSapiens or MusMusculus germline models')
    context = precompute_v_loop_distances(species=species, matrix='blosum62')
    matrix_scale = context['distance_units'][0]
    radius = 22.5 * matrix_scale
    model = context['model']
    allowed = [v for v, m in model.items() if m['locus'] == args.locus]
    def prepare(frame):
        return frame.with_columns(resolve_v_alleles(frame['v'], model, locus=args.locus).alias('_allele')).filter(
            pl.col('_allele').is_in(allowed) & (pl.col('cdr3').str.len_chars() >= 8))
    _, q = _prepare(raw)
    q = prepare(q)
    query_keys = ['cdr3', '_allele']
    withheld = None
    if args.leave_epitope_out:
        withheld = '_withheld_epitope'
        q = q.join(raw.select('query_id', pl.col('epitope').alias(withheld)), on='query_id', validate='1:1')
        query_keys.append(withheld)
    r = prepare(valid).unique(['cdr3', '_allele', *PMHC], maintain_order=True)
    groups = [(dict(zip(PMHC, task)), q, part) for task, part in r.group_by(PMHC, maintain_order=True)]
    geometry = dict(species=species, threads=args.threads, ctrim=3, v_loop_weight=1,
        alphabet='full', stem_depth=4, stem_weight=1.0, matrix='blosum62',
        exclude_exact=getattr(args, 'exact_policy', 'exclude') == 'exclude',
        v_regions='cdr12plus', position_weighting='uniform')
    state, hit_ids = grouped_kernel_selection(q, r, groups, query_keys=query_keys,
        search_keys=['cdr3', '_allele'], paired=False, cutoffs=ladder_cutoffs(matrix_scale),
        radius=radius, kernel_scale=radius, scorer='v-bounded', withheld_column=withheld,
        v_loop_distances=context, **geometry)
    selection = finish_kernel_selection(q.select('query_id'), state)
    calls = raw.join(selection, on='query_id', how='left', validate='1:1', maintain_order='left').with_columns(
        pl.col('query_id').is_in(q['query_id'].implode()).alias('available'))
    empty_fold_ids = []
    if withheld and r.height:
        exclusive = r.group_by('cdr3', '_allele').agg(pl.col('epitope').n_unique().alias('_labels'),
            pl.col('epitope').first()).filter(pl.col('_labels') == 1).group_by('epitope').len()
        empty_peptides = exclusive.filter(pl.col('len') == r.unique(['cdr3', '_allele']).height)['epitope']
        empty_fold_ids = raw.filter(pl.col('epitope').is_in(empty_peptides.implode()))['query_id'].to_list()
    calls = calls.with_columns(
        pl.when(~pl.col('available')).then(pl.lit('invalid_query_or_gene'))
        .when(pl.col('query_id').is_in(empty_fold_ids)).then(pl.lit('no_reference'))
        .when(pl.col('kernel_score').is_not_null()).then(pl.lit('searched'))
        .when(pl.col('query_id').is_in(list(hit_ids))).then(pl.lit('searched'))
        .otherwise(pl.lit('no_neighbours' if r.height else 'no_reference')).alias('status'),
        pl.col('selection_status').fill_null('no_support'))
    truth = raw.join(valid, on=keys, nulls_equal=True, validate='1:m').select('query_id', *PMHC).unique()
    if args.leave_epitope_out:
        outcome = pl.lit('null'); null_source = pl.lit('withheld')
    else:
        correct = calls.select('query_id', pl.col('selected_epitope').alias('epitope')).join(
            truth.select('query_id', 'epitope').unique(), on=['query_id', 'epitope'], how='semi')['query_id']
        outcome = pl.when(pl.col('query_id').is_in(correct.implode())).then(pl.lit('correct')).otherwise(pl.lit('incorrect'))
        null_source = pl.lit(None, dtype=pl.String)
    cases = calls.select('query_id', pl.col('kernel_score').alias('score'),
        outcome.alias('outcome'), null_source.alias('null_source'), 'available', 'status', 'selection_status')
    prefix = Path(args.output_prefix); prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix, table in [('calls', calls), ('calibration-cases', cases), ('truth', truth),
            ('target-support', support.with_columns((pl.col('reference_records') >= args.min_target_records).alias('benchmark_eligible')))]:
        table.write_csv(str(prefix)+'.'+suffix+'.tsv', separator='\t')
    provenance = db.provenance(args.vdjdb)
    from ..evalue import kernel
    from ..match import geometry as geometry_module
    import seqtree.gapblock as native_wrapper
    import seqtree._core as native_core
    source_hashes = {name: hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
        for name, module in [('kernel', kernel), ('geometry', geometry_module), ('seqtree_python', native_wrapper), ('seqtree_native', native_core)]}
    source_hashes['reference_calls'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    manifest = dict(reference=provenance, reference_sha256=provenance.get('input_sha256'),
        software={x: importlib.metadata.version(x) for x in ['vdjmatch', 'seqtree', 'vdjtools']},
        source_sha256=source_hashes, model_sha256=hashlib.sha256(json.dumps(model, sort_keys=True, default=str).encode()).hexdigest(),
        options=vars(args) | {'func': 'reference-calls'}, source_records=reference.height,
        valid_reference_records=valid.height, excluded_reference_records=reference.height-valid.height,
        searchable_reference_records=r.height, searchable_reference_keys=r.unique(['cdr3','_allele']).height,
        eligible_benchmark_records=benchmark.height, ineligible_benchmark_records=valid.height-benchmark.height,
        eligible_benchmark_peptides=eligible.height, all_reference_peptides=support.height, queries=raw.height,
        available_queries=q.height, geometry=geometry | {'radius':radius, 'kernel_scale':radius, 'cutoffs':ladder_cutoffs(matrix_scale),
            'cross_v':.25, 'floor':.01, 'scorer':'v-bounded', 'matrix_scale':matrix_scale},
        selection='maximum restriction per peptide, then maximum peptide; exact score ties remain ambiguous',
        reference_population=args.reference_population,
        target_support_rule='valid source records per species/locus/peptide before exclusion; MHC pooled only for eligibility',
        search_reference='all selected-policy valid species/locus observations, including sparse peptides; resolved V and junction length >=8',
        self_exclusion=('full-junction hit filtering across labels; fixed reference/control model populations'
            if geometry['exclude_exact'] else 'exact evidence retained for ordinary-call calibration'),
        peptide_withholding='per-query full peptide memberships removed before counts, sizes and exclusivity' if withheld else None,
        observation_unit='distinct full-junction/V/J'+('/truth peptide' if withheld else ''),
        calibration_statistic='final selected kernel_score; no component P/E transfer',
        limitations='Single-chain empirical calibration only; no paired calibration or binding posterior.',
        wall_seconds=time.monotonic()-started)
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0



def _stem4_paired(args,reference,valid,started):
    """Genuine peptide-withheld pairs; three independent native population batches."""
    import hashlib
    import importlib.metadata
    import json
    import time
    import polars as pl
    from .. import db
    from ..aggregate.candidates import PMHC
    from ..evalue.paired import build_paired_ref
    from ..evalue.kernel import grouped_kernel_selection,finish_kernel_selection,ladder_cutoffs
    from ..match.geometry import precompute_v_loop_distances,resolve_v_alleles
    from .calibration import PAIRED_STATISTIC,paired_call_event,paired_pipeline_sources

    valid=valid.filter(pl.col('gene').is_in(['TRA','TRB']))
    support=valid.group_by('gene','epitope').len().rename({'len':'reference_records'})
    eligible=support.filter(pl.col('reference_records')>=args.min_target_records).select('epitope').unique()
    pairs=build_paired_ref(valid)
    benchmark=pairs.join(eligible,on='epitope',how='semi')
    keys=['alpha','beta','v_alpha','j_alpha','v_beta','j_beta','epitope']
    raw=benchmark.select(keys).unique().sort(keys).with_row_index('query_id')
    truth=raw.join(benchmark,on=keys,nulls_equal=True,validate='1:m').select('query_id','complex_id',*PMHC).unique()
    species={'HomoSapiens':'human','MusMusculus':'mouse'}[args.species]
    context=precompute_v_loop_distances(species=species,matrix='blosum62');model=context['model']
    matrix_scale=context['distance_units'][0]
    radius=22.5*matrix_scale
    geometry=dict(species=species,threads=args.threads,ctrim=3,v_loop_weight=1,alphabet='full',
        stem_depth=4,stem_weight=1.,matrix='blosum62',exclude_exact=args.exact_policy=='exclude',
        v_regions='cdr12plus',position_weighting='uniform')
    def resolve(frame,seq,v,locus,allele):
        return frame.with_columns(resolve_v_alleles(frame[v],model,locus=locus).alias(allele)).filter(
            pl.col(allele).is_not_null() & pl.col(seq).str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]{8,}$').fill_null(False))
    def score(q,r,query_keys,search_keys,paired):
        groups=[(dict(zip(PMHC,task)),q,part) for task,part in r.group_by(PMHC,maintain_order=True)]
        state,hits=grouped_kernel_selection(q,r,groups,query_keys=query_keys+['_withheld_epitope'],
            search_keys=search_keys,paired=paired,cutoffs=ladder_cutoffs(matrix_scale),radius=radius,kernel_scale=radius,
            scorer='v-bounded',withheld_column='_withheld_epitope',v_loop_distances=context,**geometry)
        selected=finish_kernel_selection(q.select('query_id'),state)
        result=raw.select('query_id').join(selected,on='query_id',how='left',validate='1:1',maintain_order='left')
        empty=[]
        if r.height:
            exclusive=r.group_by(search_keys).agg(pl.col('epitope').n_unique().alias('_labels'),
                pl.col('epitope').first()).filter(pl.col('_labels')==1).group_by('epitope').len()
            empty_peptides=exclusive.filter(pl.col('len')==r.unique(search_keys).height)['epitope']
            empty=raw.filter(pl.col('epitope').is_in(empty_peptides.implode()))['query_id'].to_list()
        return result.with_columns(pl.col('query_id').is_in(q['query_id'].implode()).alias('available'),
            pl.when(~pl.col('query_id').is_in(q['query_id'].implode())).then(pl.lit('invalid_query_or_gene'))
            .when(pl.col('query_id').is_in(empty)).then(pl.lit('no_reference'))
            .when(pl.col('query_id').is_in(list(hits))).then(pl.lit('searched'))
            .otherwise(pl.lit('no_neighbours' if r.height else 'no_reference')).alias('status'),
            pl.col('selection_status').fill_null('no_support'))
    q=raw.rename({'alpha':'cdr3a','beta':'cdr3b'}).with_columns(pl.col('epitope').alias('_withheld_epitope'))
    linked_q=resolve(resolve(q,'cdr3a','v_alpha','TRA','_allelea'),'cdr3b','v_beta','TRB','_alleleb')
    linked_r=resolve(resolve(pairs,'alpha','v_alpha','TRA','_allele_alpha'),'beta','v_beta','TRB','_allele_beta')
    linked_keys=['alpha','beta','_allele_alpha','_allele_beta']
    linked_r=linked_r.unique([*linked_keys,*PMHC],maintain_order=True)
    calls=raw.join(score(linked_q,linked_r,['cdr3a','cdr3b','_allelea','_alleleb'],linked_keys,True),
        on='query_id',validate='1:1',maintain_order='left').rename({'selected_epitope':'selected_linked_epitope',
            'selection_status':'linked_selection_status','available':'linked_available','status':'linked_status'})
    axis_counts=dict(linked=dict(reference_keys=linked_r.unique(linked_keys).height,available_queries=linked_q.height))
    for suffix,locus,name in [('a','TRA','alpha'),('b','TRB','beta')]:
        chain_q=q.select('query_id',pl.col('cdr3'+suffix).alias('cdr3'),pl.col('v_'+name).alias('v'),'_withheld_epitope')
        chain_q=resolve(chain_q,'cdr3','v',locus,'_allele')
        chain_r=resolve(valid.filter(pl.col('gene')==locus),'cdr3','v',locus,'_allele').unique(
            ['cdr3','_allele',*PMHC],maintain_order=True)
        lane=score(chain_q,chain_r,['cdr3','_allele'],['cdr3','_allele'],False)
        lane=lane.rename({c:('kernel_'+name+'_score' if c=='kernel_score' else
            'kernel_'+name+'_competing_score' if c=='kernel_competing_score' else
            'selected_'+name+'_epitope' if c=='selected_epitope' else name+'_'+c)
            for c in lane.columns if c!='query_id'})
        calls=calls.join(lane,on='query_id',validate='1:1',maintain_order='left')
        axis_counts[name]=dict(reference_keys=chain_r.unique(['cdr3','_allele']).height,available_queries=chain_q.height)
    calls=calls.with_columns(pl.max_horizontal('kernel_score','kernel_alpha_score','kernel_beta_score').alias('kernel_any_axis_score'),
        pl.any_horizontal('linked_available','alpha_available','beta_available').alias('available'))
    cases=calls.select('query_id',pl.col('kernel_any_axis_score').alias('score'),'kernel_any_axis_score',
        pl.lit('null').alias('outcome'),pl.lit('withheld').alias('null_source'),'available',
        'linked_status','alpha_status','beta_status')
    prefix=Path(args.output_prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    for suffix,table in [('calls',calls),('calibration-cases',cases),('truth',truth),
            ('target-support',support.with_columns((pl.col('reference_records')>=args.min_target_records).alias('benchmark_eligible')))]:
        table.write_csv(str(prefix)+'.'+suffix+'.tsv',separator='\t')
    from ..evalue import kernel
    from ..match import geometry as geometry_module
    import seqtree.gapblock as wrapper
    import seqtree._core as native
    sources={name:hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
        for name,module in [('kernel',kernel),('geometry',geometry_module),('seqtree_python',wrapper),('seqtree_native',native)]}
    provenance=db.provenance(args.vdjdb)
    manifest=dict(reference=provenance,reference_sha256=provenance.get('input_sha256'),
        software={name:importlib.metadata.version(name) for name in ['vdjmatch','seqtree','vdjtools']},
        source_sha256=sources,paired_source_sha256=paired_pipeline_sources(),paired_call_event=paired_call_event(),
        model_sha256=hashlib.sha256(json.dumps(model,sort_keys=True,default=str).encode()).hexdigest(),
        options=vars(args)|{'func':'reference-calls'},source_records=reference.height,valid_reference_records=valid.height,
        excluded_reference_records=reference.height-valid.height,query_pairs=raw.height,axis_populations=axis_counts,
        eligible_benchmark_peptides=eligible.height,genuine_reference_pairs=pairs.height,
        native_population_batches=sum(x['reference_keys']>0 and x['available_queries']>0 for x in axis_counts.values()),
        geometry=geometry|dict(radius=radius,kernel_scale=radius,cutoffs=ladder_cutoffs(matrix_scale),cross_v=.25,floor=.01,
            scorer='v-bounded',matrix_scale=matrix_scale),reference_population=args.reference_population,
        selection='maximum restriction per peptide, then maximum peptide; exact score ties remain ambiguous',
        calibration_statistic=PAIRED_STATISTIC,
        target_support_rule='source records per species/chain/peptide before exclusion; either chain >= minimum; MHC pooled only for eligibility',
        observation_unit='distinct genuine alpha/beta full-junction/V/J pair and truth peptide; source complex IDs retained in truth',
        peptide_withholding='truth peptide removed from all linked and marginal memberships, counts, sizes and exclusivity',
        self_exclusion='both full junctions for linked; own full junction for each marginal' if geometry['exclude_exact'] else 'exact evidence retained',
        limitations='Empirical paired any-axis event; no joint label fusion or binding posterior.',wall_seconds=time.monotonic()-started)
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return 0
