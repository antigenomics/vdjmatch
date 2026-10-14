"""CLI for reproducible historical density baselines; no fitted combination."""
from collections import deque
from pathlib import Path


def _take_prepared(prepared, queries, reference, counts):
    if counts is None or not queries.height or not reference.height:
        return None
    if not prepared:
        raise ValueError('missing prepared density neighbour population')
    return prepared.popleft()


def register(subparsers):
    p=subparsers.add_parser('historical-density',help='reproduce the historical PSSM density component')
    p.add_argument('sample')
    p.add_argument('--vdjdb',required=True)
    p.add_argument('--reference-population',choices=['curated-benchmark','all-valid'],default='curated-benchmark')
    p.add_argument('--legacy-components',action='store_true',help='owning held-out NED, separate V/J/length odds and global first-hit gate; no empirical fusion')
    p.add_argument('--locus',choices=['TRA','TRB','paired'],required=True)
    p.add_argument('--epitope',action='append')
    p.add_argument('--mhc-a')
    p.add_argument('--mhc-b')
    p.add_argument('--mhc-match',choices=['exact','compatible'],default='exact')
    p.add_argument('--pool-reference',action='store_true',help='pool compatible restrictions per epitope before counting unique junctions')
    p.add_argument('--targets-from-sample',action='store_true',help='score each row against its declared assayed pMHC reference union')
    p.add_argument('--all-targets',action='store_true',help='emit original log evidence for every full pMHC; query target labels are ignored')
    p.add_argument('--species',default='HomoSapiens')
    p.add_argument('--control')
    p.add_argument('--cohort-column',help='paired mode: isolate rank fusion by this explicit raw metadata column and emit original single-chain scores')
    p.add_argument('--alpha-control',help='raw alpha junction control for original paired density')
    p.add_argument('--germline-background',help='raw prior table; include original sparse-reference unpaired score')
    p.add_argument('--likelihood-score',action='store_true',help='also report original components on a common log scale; paired adds chain log-density and original combined germline evidence')
    p.add_argument('--puncture-populations',action='store_true',
                   help='candidate: remove exact full junctions from density and prior population sizes/counts; preserves historical defaults')
    p.add_argument('--gapped-extension',action='store_true',help='experimental gapped density mode; matched-pssm rescores all accepted edges; no combined P-value')
    p.add_argument('--gap-geometry',choices=['clipped','historical-pssm','matched-pssm'],default='clipped',
                   help='geometry/historical-pssm add extension edges; matched-pssm uses one positional kernel and control CDF for all edges')
    p.add_argument('--gap-radius',type=int,default=90,help='clipped BLOSUM62 geometry extension radius only')
    p.add_argument('--gap-temperature',type=float,default=12,help='clipped BLOSUM62 geometry extension temperature only')
    p.add_argument('--pssm-kernel-scale',type=float,default=400,help='positive positional-kernel decay scale in integer penalty units; positional modes only')
    p.add_argument('--withhold-epitope',action='append',help='hide this peptide from candidate hits at every MHC restriction; keep other target models and controls fixed')
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--output-prefix',required=True)
    p.set_defaults(func=main)


def assayed_targets(queries, raw, reference, species):
    """Yield declared target, aligned query rows and compatible reference union.

    Query IDs are validated before joining; search-key deduplication never defines
    benchmark membership. Species aliases follow the existing control contract.
    """
    import polars as pl
    from ..aggregate.candidates import PMHC
    from ..db.schema import mhc_compatible
    from ..evalue.control import _organism
    required = ['query_id', 'species', *PMHC]
    missing = sorted(set(required) - set(raw.columns))
    if missing:
        raise ValueError('sample-declared targets require columns: '+', '.join(missing))
    if raw['query_id'].null_count() or raw['query_id'].n_unique() != raw.height:
        raise ValueError('query_id must be unique and non-null')
    organism = _organism(species)
    aliases = ['human', 'homosapiens', 'homo sapiens'] if organism == 'human' else ['mouse', 'musmusculus', 'mus musculus']
    metadata = raw.filter(pl.col('species').str.strip_chars().str.to_lowercase().is_in(aliases))
    joined = queries.join(metadata.select(required), on='query_id', how='inner',
                          validate='1:1', maintain_order='left')
    for key, part in joined.group_by(PMHC, maintain_order=True):
        task = dict(zip(PMHC, key))
        if any(not isinstance(value, str) or not value.strip() for value in key):
            selected = reference.head(0)
        else:
            selected = reference.filter(
                (pl.col('epitope') == task['epitope']) &
                (pl.col('mhc_class') == task['mhc_class']) &
                mhc_compatible('mhc_a', task['mhc_a']) &
                mhc_compatible('mhc_b', task['mhc_b']))
        yield task, part, selected


def _paired_columns(cohort_column):
    """Validate an optional cohort name against the fixed paired input/output schema."""
    import polars as pl
    from .. import io
    from ..aggregate.candidates import PMHC
    fields=['pair_id','cdr3a','va','ja','cdr3b','vb','jb','pair_status']
    score_columns={'status':pl.String,'score':pl.Float64,'germline_lr':pl.Float64,
                   'cohort_size':pl.UInt64,'rank_density_alpha':pl.Float64,
                   'rank_density_beta':pl.Float64,'rank_germline_prior':pl.Float64}
    for name in ['alpha','beta']:
        score_columns.update({c+'_'+name:dtype for c,dtype in
            [('density_score',pl.Float64),('germline_lr',pl.Float64),('p_enrichment',pl.Float64),
             ('n_reference',pl.UInt64),('control_size',pl.UInt64)]})
    if cohort_column:
        for name in ['alpha','beta']:
            score_columns.update({'unpaired_score_'+name:pl.Float64,'status_'+name:pl.String})
    if cohort_column:
        reserved=set(fields)|set(score_columns)|{'query_id','sequence_id','species',*PMHC,'binder','estimator','_locus',
            'query_id_alpha','query_id_beta','estimator_alpha','estimator_beta',
            *[c+suffix for c in ['cdr3','v','j','count','query_id','sequence_id'] for suffix in ['a','b']]}
        reserved.update(alias for aliases in io.columns.ALIASES.values() for alias in aliases)
        if cohort_column.lower() in reserved:
            raise ValueError('--cohort-column collides with a canonical input, parsed or score column')
    return fields,score_columns


def _paired_scores(a, reference, beta_control, prior_background):
    """Original independent-reference chain components and task cohort rank fusion."""
    import polars as pl
    from .. import io
    from ..aggregate.candidates import PMHC
    from ..io.airr import _read_table
    from ..evalue.control import raw_background, _organism
    from ..match.historical import density, germline_prior, paired_score, unpaired_score, shared_control_histograms, PreparedBudget
    from ..match.vgene import precompute_v_similarities
    fields,score_columns=_paired_columns(a.cohort_column)
    puncture=getattr(a,'puncture_populations',False)
    all_targets=getattr(a,'all_targets',False)
    raw=_read_table(a.sample)
    if 'query_id' not in raw.columns:
        raw=raw.with_row_index('query_id')
    link=io.columns._resolve(raw).get('pair_id')
    if link is None:
        raise ValueError('paired original scoring requires a pair_id or clone_id linkage')
    if all_targets and 'species' not in raw.columns:
        raw=raw.with_columns(pl.lit(a.species).alias('species'))
    required=['species']+([] if all_targets else PMHC)+([a.cohort_column] if a.cohort_column else [])
    missing=sorted(set(required)-set(raw.columns))
    if missing:
        raise ValueError('sample-declared targets require columns: '+', '.join(missing))
    if a.cohort_column:
        raw=raw.with_columns(pl.col(a.cohort_column).cast(pl.String))
        if raw[a.cohort_column].null_count() or raw.filter(pl.col(a.cohort_column).str.strip_chars()=='').height:
            raise ValueError('every paired row requires nonempty cohort metadata')
    if raw['query_id'].null_count() or raw['query_id'].n_unique()!=raw.height:
        raise ValueError('query_id must be unique and non-null')
    raw=raw.with_columns(pl.col(link).cast(pl.String).alias('pair_id'))
    if raw['pair_id'].null_count() or raw.filter(pl.col('pair_id')=='').height:
        raise ValueError('every paired row requires a nonempty link')
    species_aliases={s:_organism(s) for s in raw['species'].unique().to_list()}
    metadata_columns=required+(['binder'] if 'binder' in raw.columns and not all_targets else [])
    normalized=raw.with_columns(pl.col('species').replace_strict(species_aliases))
    conflicts=normalized.group_by('pair_id').agg(pl.col(metadata_columns).n_unique())
    if conflicts.filter(pl.any_horizontal(pl.col(metadata_columns)>1)).height:
        raise ValueError('linked paired rows have conflicting species, target or label metadata')
    # Taking a representative is safe only after all linked metadata agree.
    metadata=normalized.group_by('pair_id',maintain_order=True).agg(pl.col(metadata_columns).first())
    metadata=metadata.with_columns(pl.col('pair_id').alias('query_id'))
    # read_cell owns chain parsing and rejects ambiguous pairs; invalid amino acids
    # remain present for explicit chain status rather than disappearing silently.
    cells,ingestion=io.read_cell(a.sample,link=link,valid_aa=False,source='airr',return_report=True)
    paired=metadata.select('query_id','pair_id',*([a.cohort_column] if a.cohort_column else [])).join(cells.select(fields),on='pair_id',how='left',
                                                  validate='1:1',maintain_order='left')
    paired=paired.with_columns(pl.col('pair_status').fill_null('missing_both'))
    # Retain source chain IDs even when their missing junction was filtered by read_cell.
    locus=io.columns._resolve(raw).get('locus')
    if locus is None:
        v_column=io.columns._resolve(raw).get('v')
        if v_column is None:
            raise ValueError('paired rows require an explicit locus or V call')
        raw=raw.with_columns(pl.col(v_column).str.slice(0,3).str.to_uppercase().alias('_locus'))
        locus='_locus'
    multiplicity=raw.group_by('pair_id',locus).len()
    if multiplicity.filter(pl.col('len')>1).height:
        raise ValueError('ambiguous paired input: multiple raw rows for one cell/locus')
    for chain,name in [('TRA','alpha'),('TRB','beta')]:
        ids=raw.filter(pl.col(locus).str.to_uppercase()==chain).select(
            'pair_id',pl.col('query_id').alias('query_id_'+name))
        paired=paired.join(ids,on='pair_id',how='left',validate='1:1',maintain_order='left')
    alpha_control,alpha_provenance=raw_background('TRA',a.species,a.alpha_control)
    controls={'alpha':alpha_control,'beta':beta_control}
    v_similarities={name:precompute_v_similarities(paired['v'+suffix],
        reference.filter(pl.col('gene')==chain)['v'])
        for chain,name,suffix in [('TRA','alpha','a'),('TRB','beta','b')]}
    canonical=r'^[ACDEFGHIKLMNPQRSTVWY]+$'
    genes_present={suffix:pl.all_horizontal(pl.col('v'+suffix,'j'+suffix).str.strip_chars().ne('').fill_null(False))
                   for suffix in ['a','b']}
    def targets():
        if not all_targets:
            yield from assayed_targets(paired,metadata,reference,a.species)
        else:
            queries=paired.join(metadata.select('query_id','species'),on='query_id',validate='1:1').filter(pl.col('species')==_organism(a.species))
            for key,selected in reference.group_by(PMHC,maintain_order=True):
                task=dict(zip(PMHC,key))
                yield task,queries.with_columns(*[pl.lit(v,dtype=pl.String).alias(k) for k,v in task.items()]),selected
    shared_counts = {}
    prepared_batches = {name: deque() for name in ('alpha', 'beta')}
    prepared_budget = PreparedBudget()
    if (reference if all_targets else metadata).select(PMHC).unique().height > 1:
        for chain,name,suffix in [('TRA','alpha','a'),('TRB','beta','b')]:
            batches = ((part.filter(pl.col('cdr3'+suffix).str.contains(canonical).fill_null(False))
                .select(pl.col('cdr3'+suffix).alias('cdr3')), selected.filter(pl.col('gene')==chain))
                for _,part,selected in targets())
            shared_counts[name] = shared_control_histograms(controls[name], batches, threads=a.threads,
                prepared_batches=prepared_batches[name], prepared_budget=prepared_budget)
    parts=[]
    for task,part,selected in targets():
        out=part.drop('species',*PMHC)
        for chain,name,suffix in [('TRA','alpha','a'),('TRB','beta','b')]:
            ref=selected.filter(pl.col('gene')==chain)
            queries=part.filter(pl.col('cdr3'+suffix).str.contains(canonical).fill_null(False)).select(
                'query_id',pl.col('cdr3'+suffix).alias('cdr3'),
                pl.col('v'+suffix).alias('v'),pl.col('j'+suffix).alias('j'))
            if ref.height and queries.height:
                scores=density(queries,ref,controls[name],threads=a.threads,pool_reference=True,
                               control_counts=shared_counts.get(name),puncture_populations=puncture,
                               v_similarities=v_similarities[name],
                               prepared_neighbours=_take_prepared(prepared_batches[name],queries,ref,shared_counts.get(name)))
                prior=germline_prior(queries.filter(pl.all_horizontal(pl.col('v','j').str.strip_chars().ne('').fill_null(False))),
                                     ref,prior_background,include_length=chain=='TRB',exclude_exact=puncture)
                scores=scores.join(prior.select('query_id','germline_lr'),on='query_id',how='left',
                                   validate='1:1',maintain_order='left').drop(PMHC)
                if a.cohort_column:
                    valid=scores.filter(pl.col('germline_lr').is_not_null() & pl.col('score').is_not_null())
                    original=valid.select('query_id').with_columns(unpaired_score(
                        valid['score'],valid['germline_lr'],valid['n_reference'] if puncture else ref.unique('cdr3').height).alias('unpaired_score'))
                    scores=scores.join(original,on='query_id',how='left',validate='1:1',maintain_order='left')
                scores=scores.rename({c:('density_score_' if c=='score' else c+'_')+name
                                     for c in scores.columns if c!='query_id'})
                out=out.join(scores,on='query_id',how='left',validate='1:1',maintain_order='left')
            else:
                out=out.with_columns(pl.lit(None,dtype=pl.Float64).alias('density_score_'+name),
                    pl.lit(None,dtype=pl.Float64).alias('germline_lr_'+name),
                    pl.lit(None,dtype=pl.Float64).alias('p_enrichment_'+name))
            if not puncture or not (ref.height and queries.height):
                out=out.with_columns(pl.lit(ref.unique('cdr3').height,dtype=pl.UInt64).alias('n_reference_'+name),
                                     pl.lit(len(controls[name]),dtype=pl.UInt64).alias('control_size_'+name))
        out=out.with_columns(
            pl.when(pl.col('pair_status')!='paired').then(pl.col('pair_status'))
            .when(~pl.col('cdr3a').str.contains(canonical).fill_null(False)).then(pl.lit('invalid_alpha'))
            .when(~pl.col('cdr3b').str.contains(canonical).fill_null(False)).then(pl.lit('invalid_beta'))
            .when(pl.lit(any(not isinstance(v,str) or not v.strip() for v in task.values()))).then(pl.lit('invalid_target'))
            .when(~genes_present['a']).then(pl.lit('missing_v_or_j_alpha'))
            .when(~genes_present['b']).then(pl.lit('missing_v_or_j_beta'))
            .when(pl.col('n_reference_alpha')==0).then(pl.lit('no_reference_alpha'))
            .when(pl.col('n_reference_beta')==0).then(pl.lit('no_reference_beta'))
            .when(pl.col('control_size_alpha')==0).then(pl.lit('no_control_alpha'))
            .when(pl.col('control_size_beta')==0).then(pl.lit('no_control_beta'))
            .when(pl.col('germline_lr_alpha').is_null() | pl.col('germline_lr_beta').is_null()).then(pl.lit('no_prior_evidence'))
            .otherwise(pl.lit('scored')).alias('status'),
            (pl.col('germline_lr_alpha')+pl.col('germline_lr_beta')).alias('germline_lr'))
        if a.cohort_column:
            for name,suffix in [('alpha','a'),('beta','b')]:
                if 'unpaired_score_'+name not in out.columns:
                    out=out.with_columns(pl.lit(None,dtype=pl.Float64).alias('unpaired_score_'+name))
                out=out.with_columns(
                    pl.when(pl.col('query_id_'+name).is_null()).then(pl.lit('missing_'+name))
                    .when(~pl.col('cdr3'+suffix).str.contains(canonical).fill_null(False)).then(pl.lit('invalid_query'))
                    .when(pl.lit(any(not isinstance(v,str) or not v.strip() for v in task.values()))).then(pl.lit('invalid_target'))
                    .when(~genes_present[suffix]).then(pl.lit('missing_v_or_j'))
                    .when(pl.col('n_reference_'+name)==0).then(pl.lit('no_reference'))
                    .when(pl.col('control_size_'+name)==0).then(pl.lit('no_control'))
                    .when(pl.col('germline_lr_'+name).is_null()).then(pl.lit('no_prior_evidence'))
                    .otherwise(pl.lit('scored')).alias('status_'+name))
        valid=out.filter(pl.col('status')=='scored')
        out=out.with_columns((pl.col('status').eq('scored').sum().over(a.cohort_column)
            if a.cohort_column else pl.lit(valid.height)).cast(pl.UInt64).alias('cohort_size'))
        if valid.height:
            if a.cohort_column:
                fused=valid.with_columns(
                    (pl.col('density_score_beta').rank('average').over(a.cohort_column)-1).alias('rank_density_beta'),
                    (pl.col('density_score_alpha').rank('average').over(a.cohort_column)-1).alias('rank_density_alpha'),
                    (pl.col('germline_lr').rank('average').over(a.cohort_column)-1).alias('rank_germline_prior'))
                fused=fused.select('query_id','rank_density_beta','rank_density_alpha','rank_germline_prior',
                    (pl.col('rank_density_beta')+pl.col('rank_density_alpha')+pl.col('rank_germline_prior')).alias('score'))
            else:
                fused=valid.select('query_id').with_columns(
                    paired_score(valid['density_score_beta'],valid['density_score_alpha'],valid['germline_lr']),
                    (valid['density_score_beta'].rank('average')-1).alias('rank_density_beta'),
                    (valid['density_score_alpha'].rank('average')-1).alias('rank_density_alpha'),
                    (valid['germline_lr'].rank('average')-1).alias('rank_germline_prior'))
            out=out.join(fused,on='query_id',how='left',validate='1:1',maintain_order='left')
        else:
            out=out.with_columns(*[pl.lit(None,dtype=pl.Float64).alias(c) for c in
                                  ['score','rank_density_beta','rank_density_alpha','rank_germline_prior']])
        if all_targets:
            out=out.with_columns(*[pl.lit(v,dtype=pl.String).alias(k) for k,v in task.items()])
        parts.append(out)
    scores=pl.concat(parts,how='diagonal_relaxed') if parts else paired.head(0)
    # Preserve every original pair and both source chain IDs, even when its
    # species is not selected for this invocation's reference and controls.
    out=metadata.join(paired.drop('pair_id',*([a.cohort_column] if a.cohort_column else [])),on='query_id',how='left',validate='1:1',maintain_order='left')
    scoring=scores.drop([c for c in paired.columns if c!='query_id' and c in scores.columns])
    out=out.join(scoring,on='query_id',how='left',validate='1:m' if all_targets else '1:1',maintain_order='left')
    out=out.with_columns(*[pl.lit(None,dtype=dtype).alias(name) for name,dtype in score_columns.items()
                           if name not in out.columns])
    out=out.with_columns(
        pl.when(pl.col('species')!=_organism(a.species)).then(pl.lit('unselected_species'))
        .otherwise(pl.col('status').fill_null('no_reference' if all_targets and not reference.height else 'invalid_query')).alias('status'),
        pl.lit('original-paired-cohort-rank-sum-v1').alias('estimator'))
    if all_targets:
        out=out.with_columns(*[pl.lit(None,dtype=pl.String).alias(c) for c in PMHC if c not in out.columns])
    if a.cohort_column:
        out=out.with_columns(*[pl.when(pl.col('species')!=_organism(a.species)).then(pl.lit('unselected_species'))
            .otherwise(pl.col('status_'+name).fill_null('invalid_query')).alias('status_'+name) for name in ['alpha','beta']])
    return out,ingestion,alpha_provenance


def main(a):
    from importlib.metadata import version
    import json
    import math
    import time
    import polars as pl
    from .. import db,io
    from ..api import _prepare
    from ..db.cache import sha256
    from seqtree import _core
    from ..evalue.control import raw_background, _organism
    from ..match.historical import density,germline_prior,unpaired_score,shared_control_histograms
    from ..match.vgene import precompute_v_similarities
    start=time.perf_counter()
    legacy=getattr(a,'legacy_components',False)
    a.reference_population=getattr(a,'reference_population','curated-benchmark')
    if legacy and (a.locus=='paired' or not a.germline_background or a.gapped_extension
            or getattr(a,'puncture_populations',False) or a.likelihood_score or a.cohort_column
            or a.targets_from_sample or getattr(a,'all_targets',False) or getattr(a,'withhold_epitope',None)):
        raise ValueError('legacy components require TRA/TRB and germline-background without alternate reducers or withholding')
    a.puncture_populations=getattr(a,'puncture_populations',False)
    a.all_targets=getattr(a,'all_targets',False)
    a.withhold_epitope=getattr(a,'withhold_epitope',None)
    if a.withhold_epitope and not (a.all_targets or a.targets_from_sample):
        raise ValueError('--withhold-epitope requires --all-targets or --targets-from-sample')
    if a.all_targets and (a.targets_from_sample or a.epitope or a.mhc_a or a.mhc_b or a.pool_reference or a.cohort_column or a.gapped_extension):
        raise ValueError('--all-targets cannot select assayed targets, pooled restrictions, cohort ranks or an experimental gap extension')
    if a.all_targets and not a.likelihood_score:
        raise ValueError('--all-targets requires --likelihood-score')
    if a.puncture_populations and a.gapped_extension:
        raise ValueError('--puncture-populations does not support the experimental gapped extension')
    if a.threads<1:raise ValueError('threads must be positive')
    if not math.isfinite(a.pssm_kernel_scale) or a.pssm_kernel_scale<=0:
        raise ValueError('--pssm-kernel-scale must be finite and positive')
    if a.gap_geometry=='clipped' and a.pssm_kernel_scale!=400:
        raise ValueError('--pssm-kernel-scale requires a positional gap geometry')
    if a.gap_geometry in ('historical-pssm','matched-pssm'):
        if not a.gapped_extension:
            raise ValueError('--gap-geometry requires --gapped-extension')
        if a.gap_radius!=90 or a.gap_temperature!=12:
            raise ValueError('--gap-radius/--gap-temperature apply only to clipped BLOSUM62 geometry')
        from ..match.unified import require_historical_pssm_kernel
        require_historical_pssm_kernel()
    if a.cohort_column and a.locus!='paired':
        raise ValueError('--cohort-column requires --locus paired')
    if a.cohort_column:
        _paired_columns(a.cohort_column)
    if a.locus=='paired' and not ((a.targets_from_sample or a.all_targets) and a.germline_background and a.alpha_control):
        raise ValueError('--locus paired requires --targets-from-sample or --all-targets, --germline-background and --alpha-control')
    if a.locus=='paired' and a.gapped_extension:
        raise ValueError('experimental gapped extension currently supports single-chain scores only')
    if a.targets_from_sample:
        if a.epitope or a.mhc_a or a.mhc_b:
            raise ValueError('--targets-from-sample is mutually exclusive with --epitope/--mhc-a/--mhc-b')
    elif not legacy and not a.all_targets and (not a.epitope or not a.mhc_a):
        raise ValueError('explicit selection requires --epitope and --mhc-a')
    if a.locus!='paired':
        q,ingestion=io.read_rearrangement(a.sample,source='airr',return_report=True)
        _,q=_prepare(q)
        q=q.filter(pl.col('locus')==a.locus)
        if a.all_targets:
            from ..io.airr import _read_table
            raw=_read_table(a.sample)
            if 'query_id' not in raw.columns:raw=raw.with_row_index('query_id')
            if raw['query_id'].null_count() or raw['query_id'].n_unique()!=raw.height:raise ValueError('query_id must be unique and non-null')
            aliases=['human','homosapiens','homo sapiens'] if _organism(a.species)=='human' else ['mouse','musmusculus','mus musculus']
            selected=pl.col('species').str.strip_chars().str.to_lowercase().is_in(aliases).fill_null(False) if 'species' in raw.columns else pl.lit(True)
            raw=raw.with_columns(selected.alias('_species_selected'))
            q=q.join(raw.filter(pl.col('_species_selected')).select('query_id'),on='query_id',how='semi')
    organism=_organism(a.species)
    reference_species={'human':'HomoSapiens','mouse':'MusMusculus'}[organism]
    r=db.load(a.vdjdb,species=reference_species,gene=None if a.locus=='paired' else a.locus,epitope=a.epitope,
              mhc_a=a.mhc_a,mhc_b=a.mhc_b,mhc_match=a.mhc_match)
    if a.reference_population=='all-valid':
        r=r.with_columns(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).alias('reference_valid'))
    r=r.filter(pl.col('reference_valid'))
    if missing_withheld:=set(a.withhold_epitope or [])-set(r.filter(pl.col('reference_valid'))['epitope']):
        raise ValueError('withheld epitopes absent from valid species reference: '+', '.join(sorted(missing_withheld)))
    withheld_records=r.filter(pl.col('epitope').is_in(a.withhold_epitope or [])).height
    r=r.filter(~pl.col('epitope').is_in(a.withhold_epitope or []))
    missing=sorted(set(a.epitope or [])-set(r['epitope'].to_list()))
    if missing:
        raise ValueError('requested epitopes absent under the selected reference restriction: '+', '.join(missing))
    ctrl,provenance=raw_background('TRB' if a.locus=='paired' else a.locus,a.species,a.control)
    v_similarities=precompute_v_similarities(q['v'],r['v']) if a.locus!='paired' else None
    shared_counts = None
    prepared_batches = deque()
    def score_density(queries,reference,*,pool_reference=False):
        scores=density(queries,reference,ctrl,threads=a.threads,pool_reference=pool_reference,
                       control_counts=shared_counts,puncture_populations=a.puncture_populations,
                       v_similarities=v_similarities,
                       prepared_neighbours=_take_prepared(prepared_batches,queries,reference,shared_counts))
        if a.gapped_extension:
            if not pool_reference:
                raise ValueError('--gapped-extension requires pooled or sample-declared targets')
            from ..match.unified import gapped_extension,historical_pssm_extension
            if a.gap_geometry in ('historical-pssm','matched-pssm'):
                extension=historical_pssm_extension(queries,reference,ctrl,species=organism,threads=a.threads,
                    matched_background=a.gap_geometry=='matched-pssm',kernel_scale=a.pssm_kernel_scale,
                    v_similarities=v_similarities)
            else:
                extension=gapped_extension(queries,reference,ctrl,species=organism,threads=a.threads,
                                            radius=a.gap_radius,temperature=a.gap_temperature)
            scores=scores.rename({'score':'historical_density'}).join(extension,on='query_id',
                                    how='left',validate='1:1',maintain_order='left').with_columns(
                (pl.col('gapped_density') if a.gap_geometry=='matched-pssm' else
                 pl.col('historical_density')+pl.col('gapped_density')).alias('score'),
                pl.lit('matched-positional-kernel-apex6-v1' if a.gap_geometry=='matched-pssm' else
                       'pssm-plus-historical-pssm-apex6-extension-v2' if a.gap_geometry=='historical-pssm'
                       else 'pssm-plus-gapped-extension-v1').alias('estimator'))
        return scores
    prior_background=None
    if a.germline_background:
        prior_background=io.read_rearrangement(a.germline_background,source='legacy')
        if not (legacy or a.pool_reference or a.targets_from_sample or a.all_targets):
            raise ValueError('original unpaired scoring requires --pool-reference')
    if legacy:
        from ..match.historical import legacy_components
        from ..io.airr import _read_table
        out=legacy_components(q,r,ctrl,prior_background,threads=a.threads,v_similarities=v_similarities)
        raw=_read_table(a.sample)
        if 'query_id' not in raw.columns:raw=raw.with_row_index('query_id')
        names=io.columns._resolve(raw)
        if names.get('locus'):raw=raw.filter(pl.col(names['locus']).str.to_uppercase()==a.locus)
        elif names.get('v'):raw=raw.filter(pl.col(names['v']).str.slice(0,3).str.to_uppercase()==a.locus)
        identities=[c for c in ['query_id','sequence_id','clone_id'] if c in raw.columns]
        targets=r.unique('cdr3',keep='first',maintain_order=True)['epitope'].unique(maintain_order=True).to_list()
        axes=raw.select(identities).join(pl.DataFrame({'epitope':targets},schema={'epitope':pl.String}),how='cross')
        out=axes.join(out,on=['query_id','epitope'],how='left',validate='1:1',maintain_order='left').with_columns(
            pl.when(pl.col('ned_score').is_null()).then(pl.lit('invalid_query'))
                .when(~pl.col('feature_model_available')).then(pl.lit('no_feature_model'))
                .otherwise(pl.lit('scored')).alias('status'))
        prefix=Path(a.output_prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
        out.write_csv(str(prefix)+'.scores.tsv',separator='\t')
        Path(str(prefix)+'.manifest.json').write_text(json.dumps({
            'sample_sha256':sha256(Path(a.sample)),
            'software':{'vdjmatch':version('vdjmatch'),'seqtree':version('seqtree'),
                'seqtree_native_sha256':sha256(Path(_core.__file__)),
                'scorer_source_sha256':sha256(Path(legacy_components.__code__.co_filename)),
                'cli_source_sha256':sha256(Path(__file__))},
            'reference':db.provenance(a.vdjdb),'reference_population':a.reference_population,
            'control':provenance,'ingestion':ingestion,
            'germline_background':{'sha256':sha256(Path(a.germline_background)),
                'raw_rows':prior_background.height},
            'parameters':{k:v for k,v in vars(a).items() if k!='func'},
            'estimator':'legacy-heldout-components-v1',
            'representative':'first source-order peptide/V/J per union junction; input query projection and first J supplied by caller',
            'density':'sum(w*exp(-P/400)/(1+controlCDF)); w=1 same allele-stripped V, .25*directed V similarity otherwise; exact hits excluded; fixed populations',
            'features':'independent V/J/length Laplace1 log odds; denominator population+observed category count; unseen-in-both query category0; reference eligibility20 junctions',
            'gate':'global union first positive total-edit radius<=1; Poisson p_enrichment<.001 with existing rule-of-three empty-control convention',
            'scope':'positional5 substitutions, no indels; unit control5 total edits/2 insertions/2 deletions',
            'models':{name:sha256(Path(__file__).parents[1]/'resources'/folder/name) for folder,name in
                [('trimming','position_significance.tsv'),('vgene','human_v_cdr12.tsv')]},
            'not_included':['empirical negative-tail and four-channel Sidak fusion','ROC/AUC/confusion'],
            'wall_seconds':time.perf_counter()-start},indent=2)+'\n')
        return 0
    genes_present=pl.all_horizontal(pl.col('v','j').str.strip_chars().ne('').fill_null(False))
    alpha_provenance=None
    if a.locus=='paired':
        out,ingestion,alpha_provenance=_paired_scores(a,r,ctrl,prior_background)
    elif a.targets_from_sample:
        from ..aggregate.candidates import PMHC
        from ..io.airr import _read_table
        raw=_read_table(a.sample)
        if 'query_id' not in raw.columns:
            raw=raw.with_row_index('query_id')
        aliases=['human','homosapiens','homo sapiens'] if organism=='human' else ['mouse','musmusculus','mus musculus']
        species_selected=pl.col('species').str.strip_chars().str.to_lowercase().is_in(aliases).fill_null(False)
        q=q.join(raw.filter(species_selected).select('query_id'),on='query_id',how='semi',maintain_order='left')
        selected=raw
        if raw.filter(species_selected).select(PMHC).unique().height > 1:
            shared_counts = shared_control_histograms(ctrl,
                ((part, ref) for _,part,ref in assayed_targets(q,raw,r,a.species)), threads=a.threads,
                prepared_batches=prepared_batches)
        parts=[]
        for task,querypart,referencepart in assayed_targets(q,raw,r,a.species):
            if not referencepart.height:
                status='invalid_target' if any(not isinstance(v,str) or not v.strip() for v in task.values()) else 'no_reference'
                parts.append(querypart.select('query_id').with_columns(
                    *[pl.lit(v,dtype=pl.String).alias(k) for k,v in task.items()],
                    pl.lit(status).alias('status')))
                continue
            scores=score_density(querypart,referencepart,pool_reference=True)
            if prior_background is not None:
                prior=germline_prior(querypart.filter(genes_present),referencepart,prior_background,
                    include_length=a.locus=='TRB',exclude_exact=a.puncture_populations)
                valid=scores.join(prior.select('query_id','germline_lr'),on='query_id',validate='1:1',maintain_order='left')
                valid=valid.filter(pl.col('germline_lr').is_not_null() & pl.col('score').is_not_null())
                full=valid.select('query_id','germline_lr')
                if valid.height:
                    full=full.with_columns(unpaired_score(valid['score'],valid['germline_lr'],
                        valid['n_reference'] if a.puncture_populations else int(valid['n_reference'][0])))
                    if a.gapped_extension:
                        full=full.with_columns(unpaired_score(valid['historical_density'],valid['germline_lr'],
                                               int(valid['n_reference'][0])).alias('historical_score'))
                else:
                    full=full.with_columns(pl.lit(None,dtype=pl.Float64).alias('score'))
                    if a.gapped_extension:
                        full=full.with_columns(pl.lit(None,dtype=pl.Float64).alias('historical_score'))
                scores=scores.rename({'score':'density_score'}).join(full,on='query_id',how='left',
                                                                    validate='1:1',maintain_order='left')
            parts.append(scores.with_columns(
                *[pl.lit(v,dtype=pl.String).alias(k) for k,v in task.items()],
                (pl.when(pl.col('germline_lr').is_null()).then(pl.lit('missing_v_or_j'))
                 .otherwise(pl.lit('scored')) if prior_background is not None else pl.lit('scored')).alias('status')))
        # Raw metadata survives normalization's invalid/missing-sequence filtering.
        if 'locus' in selected.columns:
            declared=pl.col('locus').str.to_uppercase()
            selected=selected.filter((declared==a.locus)|~declared.is_in(['TRA','TRB']).fill_null(False))
        else:
            selected=selected.filter((pl.col('v_call').str.slice(0,3).str.to_uppercase()==a.locus)
                                     | pl.col('v_call').is_null() | (pl.col('v_call')==''))
        identifiers=['query_id']+[c for c in ('sequence_id','clone_id','binder') if c in selected.columns]
        out=selected.select(*identifiers,'species',*PMHC).with_row_index('_order')
        if parts:
            scores=pl.concat(parts,how='diagonal_relaxed').drop(PMHC)
            out=out.join(scores,on='query_id',how='left',validate='1:1',maintain_order='left')
        else:
            out=out.with_columns(pl.lit(None,dtype=pl.String).alias('status'))
        out=out.with_columns(pl.when(~species_selected).then(pl.lit('unselected_species'))
                            .otherwise(pl.col('status').fill_null('invalid_query')).alias('status')).sort('_order').drop('_order')
        if prior_background is not None:
            missing_gene_ids=q.filter(~genes_present)['query_id']
            out=out.with_columns(pl.when(pl.col('query_id').is_in(missing_gene_ids.implode()) &
                                        pl.col('status').is_in(['scored','no_reference']))
                                 .then(pl.lit('missing_v_or_j')).otherwise(pl.col('status')).alias('status'))
        score_columns={'score':pl.Float64,'p_enrichment':pl.Float64,'n_reference':pl.UInt64,
                       'control_size':pl.UInt64,'estimator':pl.String}
        if prior_background is not None:
            score_columns.update(density_score=pl.Float64,germline_lr=pl.Float64)
        if a.gapped_extension:
            score_columns.update(historical_density=pl.Float64,gapped_density=pl.Float64,
                                  availability=pl.Boolean,gapped_status=pl.String,M_gap=pl.UInt64,
                                  n_reference_unavailable=pl.UInt64,
                                  gapped_density_same_v_unweighted=pl.Float64,
                                  gapped_density_cross_v_unweighted=pl.Float64,
                                  gapped_edges_same_v=pl.UInt64,gapped_edges_cross_v=pl.UInt64,
                                  gapped_floor_density=pl.Float64,gapped_best_total_distance=pl.Int32,
                                  gapped_best_cdr3_distance=pl.Int32,gapped_best_vloop_distance=pl.Int32)
            if a.gap_geometry in ('historical-pssm','matched-pssm'):
                score_columns.update(gapped_best_kernel_penalty=pl.Int32,gapped_best_reference_junction=pl.String,
                    gapped_best_query_v=pl.String,
                    gapped_best_reference_v=pl.String,gapped_best_contribution=pl.Float64,
                    gapped_best_control_count=pl.UInt64,gapped_best_gap_length=pl.UInt32,
                    gapped_best_block_position=pl.Int32,gapped_best_alignment_status=pl.String)
            if prior_background is not None:
                score_columns['historical_score']=pl.Float64
        out=out.with_columns(*[pl.lit(None,dtype=dtype).alias(name) for name,dtype in score_columns.items()
                               if name not in out.columns])
    elif a.pool_reference or a.all_targets:
        from ..aggregate.candidates import PMHC
        parts=[]
        target_keys=[]
        groups=r.group_by(PMHC if a.all_targets else 'epitope',maintain_order=True)
        if a.all_targets and r.select(PMHC).unique().height > 1:
            shared_counts=shared_control_histograms(ctrl,((q,part) for _,part in groups),threads=a.threads,
                prepared_batches=prepared_batches)
        for _,part in groups:
            target_keys.append({c:part[c][0] for c in PMHC})
            scores=score_density(q,part,pool_reference=True)
            if prior_background is not None:
                prior=germline_prior(q.filter(genes_present),part,prior_background,
                    include_length=a.locus=='TRB',exclude_exact=a.puncture_populations)
                valid=scores.join(prior.select('query_id','germline_lr'),on='query_id',validate='1:1',maintain_order='left')
                valid=valid.filter(pl.col('germline_lr').is_not_null() & pl.col('score').is_not_null())
                full=valid.select('query_id','germline_lr')
                if valid.height:
                    full=full.with_columns(unpaired_score(valid['score'],valid['germline_lr'],
                        valid['n_reference'] if a.puncture_populations else int(valid['n_reference'][0])))
                    if a.gapped_extension:
                        full=full.with_columns(unpaired_score(valid['historical_density'],valid['germline_lr'],
                                               int(valid['n_reference'][0])).alias('historical_score'))
                else:
                    full=full.with_columns(pl.lit(None,dtype=pl.Float64).alias('score'))
                    if a.gapped_extension:
                        full=full.with_columns(pl.lit(None,dtype=pl.Float64).alias('historical_score'))
                scores=scores.rename({'score':'density_score'}).join(full,on='query_id',how='left',
                                                                    validate='1:1',maintain_order='left')
                scores=scores.with_columns(pl.when(pl.col('germline_lr').is_null()).then(pl.lit('missing_v_or_j'))
                                            .otherwise(pl.lit('scored')).alias('status'))
            parts.append(scores)
        out=pl.concat(parts,how='diagonal_relaxed') if parts else q.head(0)
        if a.all_targets:
            names=io.columns._resolve(raw)
            if names.get('locus'):raw=raw.filter(pl.col(names['locus']).str.to_uppercase()==a.locus)
            elif names.get('v'):raw=raw.filter(pl.col(names['v']).str.slice(0,3).str.to_uppercase()==a.locus)
            ids=raw.select('query_id','_species_selected')
            if parts:
                # Restore invalid/missing rows for every target without reading labels.
                out=pl.concat([ids.join(part,on='query_id',how='left',validate='1:1').with_columns(
                    *[pl.lit(target[c],dtype=pl.String).alias(c) for c in PMHC],
                    pl.when(~pl.col('_species_selected')).then(pl.lit('unselected_species')).otherwise(
                        pl.col('status').fill_null('invalid_query') if 'status' in part.columns else pl.when(pl.col('score').is_not_null()).then(pl.lit('scored')).otherwise(pl.lit('invalid_query'))).alias('status')).drop('_species_selected') for part,target in zip(parts,target_keys)],how='diagonal_relaxed')
            else:
                out=ids.with_columns(pl.when(pl.col('_species_selected')).then(pl.lit('no_reference')).otherwise(pl.lit('unselected_species')).alias('status'),pl.lit(None,dtype=pl.Float64).alias('score'),*[pl.lit(None,dtype=pl.String).alias(c) for c in PMHC]).drop('_species_selected')
            out=out.with_columns(*[pl.lit(None,dtype=t).alias(c) for c,t in [('density_score',pl.Float64),('germline_lr',pl.Float64),('n_reference',pl.UInt64),('control_size',pl.UInt64)] if c not in out.columns])
    else:
        out=score_density(q,r)
    if a.puncture_populations and a.locus!='paired':
        if 'status' not in out.columns:
            out=out.with_columns(pl.lit('scored').alias('status'))
        status=(pl.when(pl.col('n_reference')==0).then(pl.lit('no_reference'))
                .when(pl.col('control_size')==0).then(pl.lit('no_control')))
        if prior_background is not None:
            status=(status.when(pl.col('query_id').is_in(q.filter(~genes_present)['query_id'].implode()))
                    .then(pl.lit('missing_v_or_j'))
                    .when(pl.col('germline_lr').is_null()).then(pl.lit('no_prior_evidence')))
        status=status.otherwise(pl.col('status'))
        if a.all_targets:
            status=pl.when(pl.col('status').is_in(['invalid_query','unselected_species','no_reference'])).then(pl.col('status')).otherwise(status)
        out=out.with_columns(status.alias('status'))
    if a.likelihood_score:
        from ..match.historical import paired_score
        out=out.with_row_index('_likelihood_row')
        valid=out.filter(pl.col('status')=='scored' if 'status' in out.columns else pl.col('score').is_not_null())
        if a.locus=='paired':
            values=paired_score(valid['density_score_beta'],valid['density_score_alpha'],
                                valid['germline_lr'],log_scale=True)
        elif prior_background is not None:
            values=unpaired_score(valid['density_score'],valid['germline_lr'],
                                  valid['n_reference'],log_scale=True)
        else:
            values=(valid['score']+1e-6).log()
        # Diagnose the inherited prior activation without changing model counts,
        # matching, the historical score or the selected likelihood formula.
        log_density=((pl.col('density_score_beta')+1e-6).log()+(pl.col('density_score_alpha')+1e-6).log()
                     if a.locus=='paired' else (pl.col('density_score' if prior_background is not None else 'score')+1e-6).log())
        prior=pl.col('germline_lr') if prior_background is not None else pl.lit(None,dtype=pl.Float64)
        diagnostics=valid.select('_likelihood_row',log_density.alias('likelihood_density_only'),
                                 (log_density+prior).alias('likelihood_full_prior')).with_columns(values.alias('likelihood_score'))
        out=out.join(diagnostics,
                     on='_likelihood_row',how='left',validate='1:1',maintain_order='left').drop('_likelihood_row')
        if a.locus=='paired' and a.cohort_column:
            out=out.with_columns(*[pl.when(pl.col('status_'+side)=='scored').then(
                (pl.col('density_score_'+side)+1e-6).log()+(pl.col('germline_lr_'+side) if full else 0.))
                .otherwise(None).alias('likelihood_'+('full_prior_' if full else 'density_only_')+side)
                for side in ('alpha','beta') for full in (False,True)])
    prefix=Path(a.output_prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    out.write_csv(str(prefix)+'.scores.tsv',separator='\t')
    extension_manifest=None
    if a.gapped_extension:
        from ..match.unified import gapped_extension,historical_pssm_extension
        common={'geometry':a.gap_geometry,'native_sha256':sha256(Path(_core.__file__)),
                'significance':'p_enrichment remains the original component test; no combined P-value'}
        if a.gap_geometry in ('historical-pssm','matched-pssm'):
            from seqtree import SubstitutionMatrix
            from ..match.regions import significance_weights
            from ..match.vgene import vsim
            resource_root=Path(significance_weights.__code__.co_filename).parent.parent/'resources'
            scale=SubstitutionMatrix.blosum62().scale()
            extension_manifest={**common,
                'source_sha256':sha256(Path(historical_pssm_extension.__code__.co_filename)),
                'position_source_sha256':sha256(Path(significance_weights.__code__.co_filename)),
                'position_profile_sha256':sha256(resource_root/'trimming/position_significance.tsv'),
                'v_weight_source_sha256':sha256(Path(vsim.__wrapped__.__code__.co_filename)),
                'v_model_sha256':sha256(resource_root/'vgene/human_v_cdr12.tsv'),
                'position_frame':'longer full junction','matrix':'BLOSUM62 Gram','weight_scale':100,
                'weight_rounding':'max(1, round(100 * original significance_weights(length)))',
                'gap_open':2*scale*100,'gap_extend':scale*100,'gap_positions':[6],
                'gap_placement':'single Cys-relative block start min(6, shorter full junction length); equal-length prior ignored',
                'gap_charge':'2800 + (d - 1) * 1400 for d = abs(query length - reference length) > 0; zero at d = 0',
                'cutoff':5*scale*100,'kernel_scale':a.pssm_kernel_scale,
                'control_geometry':'same full-junction weighted single-gap kernel; unique full controls; full identity excluded',
                'score_composition':'matched_kernel_only' if a.gap_geometry=='matched-pssm' else 'original_plus_extension',
                'formula':('sum of all accepted positional kernel edges with the same weighted control CDF; ' if a.gap_geometry=='matched-pssm' else 'original density plus new positional kernel edges only; ') +
                          'V weight1 same allele-stripped gene/.25*vsim otherwise; original .01 expected-count denominator floor',
                'v_contract':'original raw V gene_family/vsim handling, including unresolved historical aliases; no alias-aware refinement',
                'retrieval':'all nonexact nonzero-V-weight reference representatives within cutoff; ' +
                    ('original edges rescored with matched CDF; no hit cap' if a.gap_geometry=='matched-pssm' else 'excluding original same-length<=5-substitution edges; no hit cap'),
                'diagnostics':'highest weighted-contribution reference junction/V, penalty/count/contribution/gap length; block position and alignment unavailable from scalar native output'}
        else:
            from ..match.geometry import distance_matrix, load_v_loops
            extension_manifest={**common,
                'source_sha256':sha256(Path(gapped_extension.__code__.co_filename)),
                'distance_source_sha256':sha256(Path(distance_matrix.__code__.co_filename)),
                'model_sha256':sha256(Path(load_v_loops.__code__.co_filename).parent.parent/'resources/vloops'/f'{organism}_v_loops.tsv'),
                'control_geometry':'trimmed restricted-gap CDR3; unique full controls; exact full identity punctured',
                'formula':'original density plus only new gapped/V-loop kernel edges; V weight1 same allele-stripped gene/.25 otherwise; no vsim multiplier',
                'diagnostics':'pre-prior unweighted same/cross V sums, edge counts, weighted floor contribution, and source-order best accepted edge total/CDR3/V-loop distances'}
    Path(str(prefix)+'.manifest.json').write_text(json.dumps({'sample_sha256':sha256(Path(a.sample)),
        'software':{'vdjmatch':version('vdjmatch'),'seqtree':version('seqtree'),
                    'seqtree_native_sha256':sha256(Path(_core.__file__)),
                    'scorer_source_sha256':sha256(Path(density.__code__.co_filename)),
                    'cli_source_sha256':sha256(Path(__file__))},
        'reference':db.provenance(a.vdjdb),'reference_population':a.reference_population,'control':provenance,'ingestion':ingestion,
        'alpha_control':alpha_provenance,
        'gapped_extension':extension_manifest,
        **({'likelihood':{'score_column':'likelihood_score',
            'formula':'log(D_alpha+1e-6)+log(D_beta+1e-6)+original combined germline term' if a.locus=='paired' else
                      'log(D+1e-6)+original germline term when sparse; density-only if no prior supplied',
            'components':'original density, V weight, control correction and germline counts unchanged',
            'cohort_dependent':False,'historical_score_preserved':True,'calibration':'requires identical complete selected-score null pipeline'}} if a.likelihood_score else {}),
        'parameters':{k:v for k,v in vars(a).items() if k!='func'},
        'withheld_reference_chain_records':withheld_records,
        'withholding':'candidate hits only; other target priors and raw controls retained',
        'target_selection':'every full pMHC; query labels ignored' if a.all_targets else 'assayed or explicit targets',
        'scope':'PSSM5 substitutions; unit controls5,2,2; significance nearest<=1; exact excluded',
        'representative':('first source-order V per declared task union junction' if a.targets_from_sample
                          else 'first source-order V per selected union junction' if a.pool_reference
                          else 'first source-order V per pMHC junction'),
        'germline_background':{'sha256':sha256(Path(a.germline_background)),
                               'raw_rows':prior_background.height} if prior_background is not None else None,
        'paired_contract':({'fusion':'task cohort average-rank sum of alpha/beta density and combined germline prior',
                            **({'cohort_column':a.cohort_column,'rank_partition':'species + full declared PMHC + explicit cohort',
                                'unpaired_components':'original sparse/dense alpha and beta scores; same beta raw prior background'} if a.cohort_column else {}),
                            'alpha_prior_background':'original beta raw background, V/J only',
                            'reference':'independent per-chain reference unions; exact junction excluded per chain',
                            'significance':'separate alpha/beta fixed P-values; no combined P-value'}
                           if a.locus=='paired' else None),
        'not_included':(['sparse-reference germline prior'] if prior_background is None else [])+
                       ([] if a.locus=='paired' else ['paired cohort rank fusion']),
        'wall_seconds':time.perf_counter()-start},indent=2)+'\n')
    return 0
