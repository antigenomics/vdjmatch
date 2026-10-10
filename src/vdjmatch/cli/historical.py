"""CLI for reproducible historical density baselines; no fitted combination."""
from pathlib import Path


def register(subparsers):
    p=subparsers.add_parser('historical-density',help='reproduce the historical PSSM density component')
    p.add_argument('sample')
    p.add_argument('--vdjdb',required=True)
    p.add_argument('--locus',choices=['TRA','TRB','paired'],required=True)
    p.add_argument('--epitope',action='append')
    p.add_argument('--mhc-a')
    p.add_argument('--mhc-b')
    p.add_argument('--mhc-match',choices=['exact','compatible'],default='exact')
    p.add_argument('--pool-reference',action='store_true',help='pool compatible restrictions per epitope before counting unique junctions')
    p.add_argument('--targets-from-sample',action='store_true',help='score each row against its declared assayed pMHC reference union')
    p.add_argument('--species',default='HomoSapiens')
    p.add_argument('--control')
    p.add_argument('--cohort-column',help='paired mode: isolate rank fusion by this explicit raw metadata column and emit original single-chain scores')
    p.add_argument('--alpha-control',help='raw alpha junction control for original paired density')
    p.add_argument('--germline-background',help='raw prior table; include original sparse-reference unpaired score')
    p.add_argument('--gapped-extension',action='store_true',help='experimental gapped density mode; matched-pssm rescores all accepted edges; no combined P-value')
    p.add_argument('--gap-geometry',choices=['tcrdist','historical-pssm','matched-pssm'],default='tcrdist',
                   help='tcrdist/historical-pssm add extension edges; matched-pssm uses one positional kernel and control CDF for all edges')
    p.add_argument('--gap-radius',type=int,default=90,help='TCRdist extension radius only')
    p.add_argument('--gap-temperature',type=float,default=12,help='TCRdist extension temperature only')
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
    from ..match.historical import density, germline_prior, paired_score, unpaired_score, control_histograms
    fields,score_columns=_paired_columns(a.cohort_column)
    raw=_read_table(a.sample)
    if 'query_id' not in raw.columns:
        raw=raw.with_row_index('query_id')
    link=io.columns._resolve(raw).get('pair_id')
    if link is None:
        raise ValueError('paired original scoring requires a pair_id or clone_id linkage')
    required=['species',*PMHC]+([a.cohort_column] if a.cohort_column else [])
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
    metadata_columns=required+(['binder'] if 'binder' in raw.columns else [])
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
    canonical=r'^[ACDEFGHIKLMNPQRSTVWY]+$'
    genes_present={suffix:pl.all_horizontal(pl.col('v'+suffix,'j'+suffix).str.strip_chars().ne('').fill_null(False))
                   for suffix in ['a','b']}
    parts=[]
    for task,part,selected in assayed_targets(paired,metadata,reference,a.species):
        out=part.drop('species',*PMHC)
        for chain,name,suffix in [('TRA','alpha','a'),('TRB','beta','b')]:
            ref=selected.filter(pl.col('gene')==chain)
            queries=part.filter(pl.col('cdr3'+suffix).str.contains(canonical).fill_null(False)).select(
                'query_id',pl.col('cdr3'+suffix).alias('cdr3'),
                pl.col('v'+suffix).alias('v'),pl.col('j'+suffix).alias('j'))
            if ref.height and queries.height:
                scores=density(queries,ref,controls[name],threads=a.threads,pool_reference=True)
                prior=germline_prior(queries.filter(pl.all_horizontal(pl.col('v','j').str.strip_chars().ne('').fill_null(False))),
                                     ref,prior_background,include_length=chain=='TRB')
                scores=scores.join(prior.select('query_id','germline_lr'),on='query_id',how='left',
                                   validate='1:1',maintain_order='left').drop(PMHC)
                if a.cohort_column:
                    valid=scores.filter(pl.col('germline_lr').is_not_null())
                    original=valid.select('query_id').with_columns(unpaired_score(
                        valid['score'],valid['germline_lr'],ref.unique('cdr3').height).alias('unpaired_score'))
                    scores=scores.join(original,on='query_id',how='left',validate='1:1',maintain_order='left')
                scores=scores.rename({c:('density_score_' if c=='score' else c+'_')+name
                                     for c in scores.columns if c!='query_id'})
                out=out.join(scores,on='query_id',how='left',validate='1:1',maintain_order='left')
            else:
                out=out.with_columns(pl.lit(None,dtype=pl.Float64).alias('density_score_'+name),
                    pl.lit(None,dtype=pl.Float64).alias('germline_lr_'+name),
                    pl.lit(None,dtype=pl.Float64).alias('p_enrichment_'+name))
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
        parts.append(out)
    scores=pl.concat(parts,how='diagonal_relaxed') if parts else paired.head(0)
    # Preserve every original pair and both source chain IDs, even when its
    # species is not selected for this invocation's reference and controls.
    out=metadata.join(paired.drop('pair_id',*([a.cohort_column] if a.cohort_column else [])),on='query_id',how='left',validate='1:1',maintain_order='left')
    scoring=scores.drop([c for c in paired.columns if c!='query_id' and c in scores.columns])
    out=out.join(scoring,on='query_id',how='left',validate='1:1',maintain_order='left')
    out=out.with_columns(*[pl.lit(None,dtype=dtype).alias(name) for name,dtype in score_columns.items()
                           if name not in out.columns])
    out=out.with_columns(
        pl.when(pl.col('species')!=_organism(a.species)).then(pl.lit('unselected_species'))
        .otherwise(pl.col('status').fill_null('invalid_query')).alias('status'),
        pl.lit('original-paired-cohort-rank-sum-v1').alias('estimator'))
    if a.cohort_column:
        out=out.with_columns(*[pl.when(pl.col('species')!=_organism(a.species)).then(pl.lit('unselected_species'))
            .otherwise(pl.col('status_'+name).fill_null('invalid_query')).alias('status_'+name) for name in ['alpha','beta']])
    return out,ingestion,alpha_provenance


def main(a):
    from importlib.metadata import version
    import json
    import time
    import polars as pl
    from .. import db,io
    from ..api import _prepare
    from ..db.cache import sha256
    from seqtree import _core
    from ..evalue.control import raw_background, _organism
    from ..match.historical import density,germline_prior,unpaired_score
    start=time.perf_counter()
    if a.threads<1:raise ValueError('threads must be positive')
    if a.gap_geometry in ('historical-pssm','matched-pssm'):
        if not a.gapped_extension:
            raise ValueError('--gap-geometry requires --gapped-extension')
        if a.gap_radius!=90 or a.gap_temperature!=12:
            raise ValueError('--gap-radius/--gap-temperature apply only to TCRdist geometry')
        from ..match.unified import require_historical_pssm_kernel
        require_historical_pssm_kernel()
    if a.cohort_column and a.locus!='paired':
        raise ValueError('--cohort-column requires --locus paired')
    if a.cohort_column:
        _paired_columns(a.cohort_column)
    if a.locus=='paired' and not (a.targets_from_sample and a.germline_background and a.alpha_control):
        raise ValueError('--locus paired requires --targets-from-sample, --germline-background and --alpha-control')
    if a.locus=='paired' and a.gapped_extension:
        raise ValueError('experimental gapped extension currently supports single-chain scores only')
    if a.targets_from_sample:
        if a.epitope or a.mhc_a or a.mhc_b:
            raise ValueError('--targets-from-sample is mutually exclusive with --epitope/--mhc-a/--mhc-b')
    elif not a.epitope or not a.mhc_a:
        raise ValueError('explicit selection requires --epitope and --mhc-a')
    if a.locus!='paired':
        q,ingestion=io.read_rearrangement(a.sample,source='airr',return_report=True)
        _,q=_prepare(q)
        q=q.filter(pl.col('locus')==a.locus)
    organism=_organism(a.species)
    reference_species={'human':'HomoSapiens','mouse':'MusMusculus'}[organism]
    r=db.load(a.vdjdb,species=reference_species,gene=None if a.locus=='paired' else a.locus,epitope=a.epitope,
              mhc_a=a.mhc_a,mhc_b=a.mhc_b,mhc_match=a.mhc_match)
    r=r.filter(pl.col('reference_valid'))
    missing=sorted(set(a.epitope or [])-set(r['epitope'].to_list()))
    if missing:
        raise ValueError('requested epitopes absent under the selected reference restriction: '+', '.join(missing))
    ctrl,provenance=raw_background('TRB' if a.locus=='paired' else a.locus,a.species,a.control)
    def score_density(queries,reference,*,pool_reference=False):
        scores=density(queries,reference,ctrl,threads=a.threads,pool_reference=pool_reference)
        if a.gapped_extension:
            if not pool_reference:
                raise ValueError('--gapped-extension requires pooled or sample-declared targets')
            from ..match.unified import gapped_extension,historical_pssm_extension
            if a.gap_geometry in ('historical-pssm','matched-pssm'):
                extension=historical_pssm_extension(queries,reference,ctrl,species=organism,threads=a.threads,
                    matched_background=a.gap_geometry=='matched-pssm')
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
        if not (a.pool_reference or a.targets_from_sample):
            raise ValueError('original unpaired scoring requires --pool-reference')
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
                prior=germline_prior(querypart.filter(genes_present),referencepart,prior_background,include_length=a.locus=='TRB')
                valid=scores.join(prior.select('query_id','germline_lr'),on='query_id',validate='1:1',maintain_order='left')
                full=valid.select('query_id','germline_lr')
                if valid.height:
                    full=full.with_columns(unpaired_score(valid['score'],valid['germline_lr'],int(valid['n_reference'][0])))
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
    elif a.pool_reference:
        parts=[]
        for _,part in r.group_by('epitope',maintain_order=True):
            scores=score_density(q,part,pool_reference=True)
            if prior_background is not None:
                prior=germline_prior(q.filter(genes_present),part,prior_background,include_length=a.locus=='TRB')
                valid=scores.join(prior.select('query_id','germline_lr'),on='query_id',validate='1:1',maintain_order='left')
                full=valid.select('query_id','germline_lr')
                if valid.height:
                    full=full.with_columns(unpaired_score(valid['score'],valid['germline_lr'],int(valid['n_reference'][0])))
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
        out=pl.concat(parts)
    else:
        out=score_density(q,r)
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
                'cutoff':5*scale*100,'kernel_scale':400,
                'control_geometry':'same full-junction weighted single-gap kernel; unique full controls; full identity excluded',
                'score_composition':'matched_kernel_only' if a.gap_geometry=='matched-pssm' else 'original_plus_extension',
                'formula':('sum of all accepted positional kernel edges with the same weighted control CDF; ' if a.gap_geometry=='matched-pssm' else 'original density plus new positional kernel edges only; ') +
                          'V weight1 same allele-stripped gene/.25*vsim otherwise; original .01 expected-count denominator floor',
                'v_contract':'original raw V gene_family/vsim handling, including unresolved historical aliases; no alias-aware refinement',
                'retrieval':'all nonexact nonzero-V-weight reference representatives within cutoff; ' +
                    ('original edges rescored with matched CDF; no hit cap' if a.gap_geometry=='matched-pssm' else 'excluding original same-length<=5-substitution edges; no hit cap'),
                'diagnostics':'highest weighted-contribution reference junction/V, penalty/count/contribution/gap length; block position and alignment unavailable from scalar native output'}
        else:
            from ..match.tcrdist import distance_matrix, load_v_loops
            extension_manifest={**common,
                'source_sha256':sha256(Path(gapped_extension.__code__.co_filename)),
                'distance_source_sha256':sha256(Path(distance_matrix.__code__.co_filename)),
                'model_sha256':sha256(Path(load_v_loops.__code__.co_filename).parent.parent/'resources/tcrdist'/f'{organism}_v_loops.tsv'),
                'control_geometry':'trimmed restricted-gap CDR3; unique full controls; exact full identity punctured',
                'formula':'original density plus only new gapped/V-loop kernel edges; V weight1 same allele-stripped gene/.25 otherwise; no vsim multiplier',
                'diagnostics':'pre-prior unweighted same/cross V sums, edge counts, weighted floor contribution, and source-order best accepted edge total/CDR3/V-loop distances'}
    Path(str(prefix)+'.manifest.json').write_text(json.dumps({'sample_sha256':sha256(Path(a.sample)),
        'software':{'vdjmatch':version('vdjmatch'),'seqtree':version('seqtree'),
                    'seqtree_native_sha256':sha256(Path(_core.__file__)),
                    'scorer_source_sha256':sha256(Path(density.__code__.co_filename)),
                    'cli_source_sha256':sha256(Path(__file__))},
        'reference':db.provenance(a.vdjdb),'control':provenance,'ingestion':ingestion,
        'alpha_control':alpha_provenance,
        'gapped_extension':extension_manifest,
        'parameters':{k:v for k,v in vars(a).items() if k!='func'},
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
