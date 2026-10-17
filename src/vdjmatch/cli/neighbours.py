"""Native junction/V-loop single-chain/linked-pair distance neighbourhoods."""
from pathlib import Path

GEOMETRY_EXTRAS = ('flank', 'framev', 'framemid', 'framej', 'freegap', 'jw',
                   'stem4j', 'jlen', 'vjlen', 'norm', 'jnorm', 'sw', 'nw')
J_GEOMETRY_EXTRAS = frozenset(('jw', 'stem4j', 'jlen', 'vjlen', 'jnorm'))


def _grouped_detailed_eligible(a,kernel,selection_only,groups,references):
    """Keep advanced radius/rank and null-label diagnostics on their existing route."""
    return (not selection_only and kernel and a.background=='vdjdb-other' and
        a.position_weighting=='uniform' and not (a.reference_radius or a.nearest_rank or
        a.conditional_j or a.epitope_rules or getattr(a,'chemistry',False)) and
        all(task['epitope'] is not None for task,_,_ in groups) and
        (len(groups)>1 or not references['epitope'].null_count()))


def _apply_geometry_extra(a):
    """Apply a retained fixed recipe before preparing its shared V matrices."""
    extra = getattr(a, 'geometry_extra', None)
    if extra is None:
        return
    if extra not in GEOMETRY_EXTRAS:
        raise ValueError('unknown geometry extra')
    if (a.neighbour_weighting != 'kernel' or not a.selection_only or
            a.background != 'vdjdb-other' or a.alphabet != 'full' or
            a.position_weighting != 'uniform' or a.v_regions != 'cdr12plus' or
            a.v_loop_weight != 1):
        raise ValueError('geometry extras require compact uniform full-alphabet kernel scoring with all V regions')
    from ..match.geometry import require_geometry_capabilities
    require_geometry_capabilities(extra,paired=a.locus in ('paired','auto'))
    a.distance_matrix = 'clipped' if extra in ('framev', 'framemid', 'framej', 'freegap') else 'blosum62'
    a.junction_ends = 'trim3'
    a.stem_depth = 4 if extra == 'stem4j' else 0
    a.stem_weight = 1.0


def register(subparsers, name='neighbours'):
    p=subparsers.add_parser(name,help='targeted stem4 kernel evidence' if name=='annotate' else 'native junction/V-loop distance balls')
    p.add_argument('sample')
    p.add_argument('--vdjdb',required=True)
    p.add_argument('--locus',choices=['TRA','TRB','paired'],required=True)
    p.add_argument('--species',default='HomoSapiens')
    p.add_argument('--reference-population',choices=['all-valid','curated-benchmark'],
                   default='all-valid' if name=='annotate' else 'curated-benchmark',
                   help='ordinary valid observations or the source-curated benchmark population')
    p.add_argument('--epitope',action='append')
    p.add_argument('--mhc-a')
    p.add_argument('--mhc-b')
    p.add_argument('--mhc-match',choices=['exact','compatible'],default='exact')
    p.add_argument('--targets-from-sample',action='store_true',help='score only each row\'s declared assayed pMHC with compatible MHC resolution')
    exports=p.add_mutually_exclusive_group()
    exports.add_argument('--selection-only',action='store_true',help='kernel selection with bounded table retention; omit detailed candidate and all-label evidence tables')
    if name!='annotate':
        exports.add_argument('--evidence-only',action='store_true',help='export all kernel evidence and query dispositions without candidate distance matrices')
    p.add_argument('--all-targets',action='store_true',help='emit separate evidence for every full pMHC; query target labels are ignored')
    p.add_argument('--withhold-epitope',action='append',help='hide this peptide from candidate hits at every MHC restriction; keep the full background model fixed')
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--verbose',action='store_true',help='report grouped query/reference workload and count progress')
    p.add_argument('--distance-matrix',choices=['clipped','blosum62','pam250','mj','vdjam'],default='clipped',
                   help='clipped BLOSUM62 costs or retained unclipped BLOSUM62/PAM250/MJ/VDJAM Gram penalties')
    p.add_argument('--alphabet',choices=['full','imgt'],default='full',
                   help='experimental ordered IMGT-class geometry; same projected matrix for CDR3/V loops, raw identities retained')
    p.add_argument('--junction-ends',choices=['trim32','trim3'],default='trim32',
                   help='clipped BLOSUM62 geometry trims 3 N/2 C residues; trim3 tests 3 residues at both ends')
    p.add_argument('--v-loop-weight',type=int,choices=[0,1],default=1,
                   help='fixed experimental CDR3-only(0) versus CDR3+V(1); identical V eligibility')
    p.add_argument('--v-regions',choices=['cdr12plus','cdr12'],default='cdr12plus',
                   help='aligned V regions: CDR1/CDR2/CDR2.5 (default), or retained vcdr12 extra omitting CDR2.5')
    p.add_argument('--stem-depth',type=int,default=0,
                   help='joint V-J hairpin stem levels added to the junction distance. '
                        'The CDR3 loop pairs position i with position L-1-i; a level counts '
                        'as matched only when both residues match. 0 reproduces '
                        'clipped BLOSUM62 geometry exactly; 4 uses four stem levels. '
                        'Defaults: 4 for annotate/evidence, 0 for neighbours.')
    p.add_argument('--stem-weight',type=float,default=1.0,
                   help='weight of the joint stem term, in shared substitution-matrix cost '
                        'units; 1.0 is the measured arm and spans the same range as one '
                        'aligned residue')
    p.add_argument('--background',choices=['real','generative','vdjdb-other'],
                   help='experimental assayed-target ball enrichment; no P-value')
    p.add_argument('--control',help='raw junction/V control table for real or generative background')
    p.add_argument('--position-weighting',choices=['uniform','significance','significance-trimmed','ndn'],default='uniform',
                   help='enrichment geometry; ndn uses fixed0.1V/J flanks and1NDN with verified ARDA boundaries')
    p.add_argument('--neighbour-weighting',choices=['ball','linear','local-rank','kernel'],default='ball',
                   help='ball/linear enrichment or experimental background-local rank evidence')
    p.add_argument('--kernel-scorer', choices=['v-bounded','all-bounded','vonly-bounded',
                   'v-rare','all-rare','vonly-rare'], default='v-bounded',
                   help='advanced ladder alternative; v-bounded is the targeted default')
    p.add_argument('--geometry-extra', choices=GEOMETRY_EXTRAS,
                   help='advanced retained fixed geometry; sets its matrix, trim3 ends and stem recipe; requires compact kernel scoring')
    p.add_argument('--rank-signal-prior',type=float,default=.5,
                   help='prior for the local-rank model only; not biological specificity prevalence')
    radius=p.add_mutually_exclusive_group()
    radius.add_argument('--radius',type=int,help='absolute maximum distance for historical reproduction; paired uses the linked alpha+beta sum')
    radius.add_argument('--radius-factor',type=float,help='radius / maximum effective matrix substitution cost; default22.5 (radius90 with clipped BLOSUM62)')
    p.add_argument('--reference-radius',action='store_true',help='reference-only Kneedle radius per compatible target; uniform nearest-rank mode')
    p.add_argument('--knee-floor-factor',type=float,default=.5,help='positive-knee median floor within this species/locus/geometry')
    p.add_argument('--exclude-exact',action='store_true')
    p.add_argument('--unbounded-nearest',action='store_true',help='exhaustive nearest reference for each sample-declared target; neighbour counts retain the radius')
    p.add_argument('--conditional-j',action='store_true',help='experimental J usage conditional on V and junction length; ranking only')
    p.add_argument('--epitope-rules',action='store_true',help='experimental joint V usage with V-conditioned CDR3/motif kernels; ranking only')
    p.add_argument('--rules-only',action='store_true',help='declared-target unbounded rule evidence and query statuses; omit ordinary counts and candidate distances')
    p.add_argument('--chemistry',action='store_true',help='experimental V-conditioned local KF4 mixture; requires bench extra; ranking only')
    p.add_argument('--chemistry-bandwidth-factor',type=float,default=1.0,help='Gaussian bandwidth / amino-acid-table KF4 standard deviation')
    p.add_argument('--nearest-rank',action='store_true',help='experimental background rank at the exhaustive nearest target distance; uniform geometry only')
    p.add_argument('--shortlist-audit',action='store_true',help='audit actual CDR3-only radius retrieval and unchanged V costs; uniform single-chain only')
    p.add_argument('--output-prefix',required=True)
    p.set_defaults(func=main)
    if name=='annotate':
        p.set_defaults(stem_depth=4,distance_matrix='blosum62',junction_ends='trim3',
                       neighbour_weighting='kernel',background='vdjdb-other')
    return p


def _paired_input(path, raw, model, organism):
    """Preserve raw pair/chain identities, including rows dropped during ingestion."""
    import polars as pl
    from .. import io
    from ..match.geometry import resolve_v_alleles
    from ..aggregate.candidates import PMHC
    if 'clone_id' not in raw.columns:
        raise ValueError('paired AIRR input requires clone_id linkage')
    locus=pl.col('locus').str.to_uppercase() if 'locus' in raw.columns else pl.col('v_call').str.slice(0,3).str.to_uppercase()
    raw=raw.with_columns(pl.col('clone_id').cast(pl.String).alias('pair_id'),locus.alias('_locus'))
    if raw.filter(pl.col('pair_id').is_null() | (pl.col('pair_id')=='')).height:
        raise ValueError('every paired row requires nonempty clone_id linkage')
    if raw.filter(~pl.col('_locus').is_in(['TRA','TRB']).fill_null(False)).height:
        raise ValueError('paired rows require explicit TRA/TRB roles when V calls do not establish locus')
    if raw.group_by('pair_id','_locus').len().filter(pl.col('len')>1).height:
        raise ValueError('ambiguous paired input: multiple rows for one clone/locus')
    metadata=[c for c in ['species',*PMHC] if c in raw.columns]
    if 'species' in metadata:
        raw=raw.with_columns(pl.col('species').str.strip_chars().str.to_lowercase().replace(
            {'human':'human','homosapiens':'human','homo sapiens':'human',
             'mouse':'mouse','musmusculus':'mouse','mus musculus':'mouse'}).alias('_species'))
    shared=['_species' if c=='species' else c for c in metadata]
    if shared and raw.group_by('pair_id').agg(*(pl.col(c).n_unique().alias(c) for c in shared)).filter(
            pl.any_horizontal(*(pl.col(c)>1 for c in shared))).height:
        raise ValueError('paired query has conflicting species or peptide/MHC annotations')
    v_column=io.columns._resolve(raw).get('v')
    pairs=raw.group_by('pair_id',maintain_order=True).agg(
        pl.col('clone_id').first(),*(pl.col(c).first() for c in metadata),
        *[pl.col(c).filter(pl.col('_locus')==locus).first().alias(c+suffix)
          for locus,suffix in [('TRA','a'),('TRB','b')]
          for c in ['query_id',*(['sequence_id'] if 'sequence_id' in raw.columns else [])]],
        *[(pl.col(v_column).filter(pl.col('_locus')==locus).first() if v_column else
           pl.lit(None,dtype=pl.String)).alias('v_call_'+name)
          for locus,name in [('TRA','alpha'),('TRB','beta')]])
    pairs=pairs.with_row_index('query_id')
    cells,ingestion=io.read_cell(path,link='clone_id',source='airr',valid_aa=False,return_report=True)
    # Chain IDs and V calls come from raw rows, including missing/invalid sequences.
    cells=cells.drop([c for c in cells.columns if c.startswith('query_id') or c.startswith('sequence_id')])
    pairs=pairs.join(cells,on='pair_id',how='left',validate='1:1',maintain_order='left')
    for suffix,locus,name in [('a','TRA','alpha'),('b','TRB','beta')]:
        allele='_allele'+suffix
        pairs=pairs.with_columns(resolve_v_alleles(pairs['v_call_'+name],model,locus=locus).alias(allele))
        allowed=[v for v,m in model.items() if m['locus']==locus]
        pairs=pairs.with_columns(
            pl.col('cdr3'+suffix).str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]{8,}$').fill_null(False).alias('_sequence_valid'+suffix),
            pl.col(allele).is_in(allowed).fill_null(False).alias('_gene_valid'+suffix))
        pairs=pairs.with_columns((pl.col('_sequence_valid'+suffix)&pl.col('_gene_valid'+suffix)).alias('_valid'+suffix))
    species_ok=(pl.col('species').str.strip_chars().str.to_lowercase().is_in(
        ['human','homosapiens','homo sapiens'] if organism=='human' else ['mouse','musmusculus','mus musculus'])
        if 'species' in pairs.columns else pl.lit(True))
    pairs=pairs.with_columns(
        (pl.col('_valida') & pl.col('_validb') & species_ok.fill_null(False)).alias('available'),
        pl.when(~species_ok.fill_null(False)).then(pl.lit('unselected_species'))
        .when(pl.col('query_ida').is_null()).then(pl.lit('missing_alpha'))
        .when(pl.col('query_idb').is_null()).then(pl.lit('missing_beta'))
        .when(~pl.col('_sequence_valida')).then(pl.lit('invalid_alpha_query'))
        .when(~pl.col('_sequence_validb')).then(pl.lit('invalid_beta_query'))
        .when(~pl.col('_gene_valida')).then(pl.lit('unknown_alpha_gene'))
        .when(~pl.col('_gene_validb')).then(pl.lit('unknown_beta_gene'))
        .otherwise(pl.lit('searched')).alias('status'))
    return pairs,ingestion


def main(a, *, marginal_tables=None, v_loop_distances=None):
    from importlib.metadata import version
    from importlib.resources import files
    import json
    import time
    import numpy as np
    import polars as pl
    from .. import db,io
    from ..api import _prepare
    from ..aggregate.candidates import PMHC
    from ..io.airr import _read_table
    from ..db.cache import sha256
    from seqtree import _core,gapblock
    from ..match.geometry import distance_matrix,paired_distance_matrix,load_v_loops,resolve_v_alleles,distance_units,distance_cutoff,_substitution_costs
    from ..evalue.control import _organism
    start=time.perf_counter()
    _apply_geometry_extra(a)
    geometry_extra=getattr(a,'geometry_extra',None)
    j_geometry=geometry_extra in J_GEOMETRY_EXTRAS
    kernel=a.neighbour_weighting=='kernel'
    selection_only=getattr(a,'selection_only',False)
    rules_only=getattr(a,'rules_only',False)
    evidence_only=getattr(a,'evidence_only',False)
    if evidence_only and (selection_only or not kernel or a.background!='vdjdb-other' or
            a.position_weighting!='uniform' or a.shortlist_audit or a.unbounded_nearest or
            a.reference_radius or a.nearest_rank or a.conditional_j or a.epitope_rules or
            a.chemistry or geometry_extra):
        raise ValueError('--evidence-only requires uniform vdjdb-other kernel scoring without other diagnostics')
    if rules_only and (not (a.epitope_rules and a.targets_from_sample and a.unbounded_nearest and
            a.background=='vdjdb-other' and a.neighbour_weighting=='linear' and
            a.position_weighting=='uniform' and a.distance_matrix=='clipped' and a.alphabet=='full') or
            a.all_targets or selection_only or a.nearest_rank or a.reference_radius or a.shortlist_audit or
            a.conditional_j or a.chemistry or a.stem_depth or geometry_extra or a.withhold_epitope):
        raise ValueError('--rules-only requires declared-target unbounded clipped epitope rules with linear vdjdb-other background and no other diagnostics')
    if marginal_tables is not None and (a.locus!='paired' or not kernel):
        raise ValueError('run-local marginal tables require paired kernel scoring')
    if selection_only and (not kernel or a.shortlist_audit or a.unbounded_nearest):
        raise ValueError('selection-only requires kernel scoring without candidate diagnostics')
    if kernel and not a.background:
        raise ValueError('kernel scoring requires a reference/control background')
    if a.position_weighting=='significance-trimmed' and (not kernel or not selection_only or a.junction_ends!='trim3'):
        raise ValueError('trimmed significance requires compact kernel scoring and trim3 junction ends')
    if kernel and (a.position_weighting not in ('uniform','significance-trimmed') or a.nearest_rank or a.reference_radius
                   or a.epitope_rules or a.conditional_j or a.chemistry):
        raise ValueError('kernel scoring requires uniform geometry and excludes other diagnostic reducers')
    if kernel and not (a.targets_from_sample or a.epitope or a.mhc_a or a.mhc_b):
        a.all_targets=True
    if a.stem_depth < 0:
        raise SystemExit('--stem-depth must be a nonnegative number of residues')
    if not np.isfinite(a.stem_weight) or a.stem_weight < 0:
        raise SystemExit('--stem-weight must be finite and nonnegative')
    if a.stem_depth and (a.epitope_rules or a.conditional_j or a.chemistry):
        raise ValueError('stem geometry is separate from J/epitope-rule/chemistry diagnostics')
    a.alphabet=getattr(a,'alphabet','full')
    a.v_regions=getattr(a,'v_regions','cdr12plus')
    matrix_scale,bin_width=distance_units(_substitution_costs(a.alphabet,a.distance_matrix))
    if a.alphabet!='full' and (a.epitope_rules or a.conditional_j or getattr(a,'chemistry',False)):
        raise ValueError('class alignment is separate from J/epitope-rule/chemistry diagnostics')
    if a.radius_factor is not None and (not np.isfinite(a.radius_factor) or a.radius_factor<=0):
        raise ValueError('radius factor must be finite and positive')
    absolute_radius=a.radius is not None
    if not absolute_radius:
        a.radius_factor=22.5 if a.radius_factor is None else a.radius_factor
        a.radius=a.radius_factor*matrix_scale
        if a.radius.is_integer():a.radius=int(a.radius)
    if not np.isfinite(a.radius) or a.radius*200>np.iinfo(np.int32).max:
        raise ValueError('radius must fit finite native distance units')
    requested_radius=a.radius
    native_radius=distance_cutoff(a.radius)
    a.radius=native_radius/200
    if a.radius.is_integer():a.radius=int(a.radius)
    if a.threads<1 or a.radius<0:raise ValueError('positive threads and nonnegative radius required')
    if a.reference_radius and (not a.nearest_rank or a.background!='vdjdb-other' or
            a.neighbour_weighting!='local-rank' or a.position_weighting!='uniform' or
            a.v_loop_weight!=1 or a.shortlist_audit or a.conditional_j or a.epitope_rules or a.radius<=0):
        raise ValueError('--reference-radius requires uniform nearest local-rank with V loops and vdjdb-other')
    if not np.isfinite(a.knee_floor_factor) or a.knee_floor_factor<=0:
        raise ValueError('knee floor factor must be finite and positive')
    a.all_targets=getattr(a,'all_targets',False)
    a.withhold_epitope=getattr(a,'withhold_epitope',None)
    if a.withhold_epitope and not (a.all_targets or a.targets_from_sample):
        raise ValueError('--withhold-epitope requires --all-targets or --targets-from-sample')
    if a.all_targets and (a.targets_from_sample or a.epitope or a.mhc_a or a.mhc_b):
        raise ValueError('--all-targets cannot be combined with assayed or explicit target selectors')
    if a.unbounded_nearest and not (a.targets_from_sample or a.all_targets):
        raise ValueError('--unbounded-nearest requires --targets-from-sample or --all-targets')
    if a.nearest_rank and (not a.background or not a.unbounded_nearest or a.position_weighting!='uniform'):
        raise ValueError('--nearest-rank requires a background, --unbounded-nearest and uniform geometry')
    if a.conditional_j and (not a.nearest_rank or a.neighbour_weighting!='local-rank' or a.background!='vdjdb-other'):
        raise ValueError('--conditional-j requires nearest local-rank evidence and vdjdb-other background')
    if a.epitope_rules and (a.background!='vdjdb-other' or a.position_weighting!='uniform' or a.radius<=0):
        raise ValueError('--epitope-rules requires vdjdb-other, uniform geometry and positive radius')
    if a.epitope_rules and a.conditional_j:
        raise ValueError('epitope rules and conditional J are separate diagnostic models')
    if getattr(a,'chemistry',False) and (a.background!='vdjdb-other' or a.epitope_rules or a.conditional_j):
        raise ValueError('--chemistry requires vdjdb-other and is separate from J/epitope-rule diagnostics')
    if a.shortlist_audit and (a.locus=='paired' or not a.targets_from_sample or
                             not a.unbounded_nearest or a.position_weighting!='uniform' or a.v_loop_weight!=1):
        raise ValueError('--shortlist-audit requires uniform single-chain target/nearest mode with V weight1')
    if a.background and not (a.targets_from_sample or a.all_targets or kernel):
        raise ValueError('--background requires --targets-from-sample or --all-targets')
    if a.locus=='paired' and a.background and (a.background!='vdjdb-other' or a.position_weighting not in ('uniform','ndn','significance-trimmed')):
        raise ValueError('paired background requires vdjdb-other and uniform or ndn positional weights; generative/real joint controls are not defined')
    if a.neighbour_weighting=='linear' and (not a.background or a.radius<=0):
        raise ValueError('linear neighbour weighting requires a background and positive radius')
    if a.neighbour_weighting=='local-rank' and (not a.background or not 0<a.rank_signal_prior<1):
        raise ValueError('local-rank requires a background and prior strictly between zero and one')
    if bool(a.control) != (a.background in {'real','generative'}):
        raise ValueError('--control is required only for --background real or generative')
    if a.position_weighting in ('significance','ndn') and (not a.background or a.junction_ends!='trim32'):
        raise ValueError('full positional weighting requires --background and excludes an additional trim3')
    raw=_read_table(a.sample)
    if 'query_id' not in raw.columns:raw=raw.with_row_index('query_id')
    if raw['query_id'].null_count() or raw['query_id'].n_unique()!=raw.height:
        raise ValueError('query_id must be non-null and unique')
    if a.all_targets:
        raw=raw.drop(PMHC,strict=False)
    organism=_organism(a.species)
    ctrim=2 if a.junction_ends=='trim32' else 3
    paired=a.locus=='paired'
    from ..match.geometry import precompute_v_loop_distances
    if v_loop_distances is None:
        v_loop_distances=precompute_v_loop_distances(species=organism,alphabet=a.alphabet,matrix=a.distance_matrix,v_regions=a.v_regions)
    elif any(v_loop_distances[key]!=value for key,value in
             [('species',organism),('alphabet',a.alphabet),('matrix',a.distance_matrix),('v_regions',a.v_regions)]):
        raise ValueError('run-local V tables must match the requested species and geometry')
    model=v_loop_distances['model']
    if paired:
        pair_rows,ingestion=_paired_input(a.sample,raw,model,organism)
        q=pair_rows.filter(pl.col('available'))
        metadata=pair_rows.select('query_id',*[c for c in ['species',*PMHC] if c in pair_rows.columns])
    else:
        data,ingestion=io.read_rearrangement(a.sample,source='airr',valid_aa=False,return_report=True)
        _,q=_prepare(data)
        metadata=raw
    reference_species={'human':'HomoSapiens','mouse':'MusMusculus'}[organism]
    from ..db.index import load_or_build_index
    load_or_build_index(a.vdjdb, threads=a.threads)
    r=db.load(a.vdjdb,species=reference_species,gene=None if paired else a.locus,epitope=None if kernel else a.epitope,
              mhc_a=None if kernel else a.mhc_a,mhc_b=None if kernel else a.mhc_b,mhc_match=a.mhc_match)
    original=r.height
    if a.reference_population=='all-valid':
        r=r.with_columns(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).alias('reference_valid'))
    if missing_withheld:=set(a.withhold_epitope or [])-set(r.filter(pl.col('reference_valid'))['epitope']):
        raise ValueError('withheld epitopes absent from valid species reference: '+', '.join(sorted(missing_withheld)))
    withheld_records=r.filter(pl.col('epitope').is_in(a.withhold_epitope or [])).height
    allowed=[v for v,m in model.items() if m['locus']==a.locus]
    def prepare(frame):
        return frame.with_columns(resolve_v_alleles(frame['v'],model,locus=a.locus).alias('_allele')).filter(
                pl.col('_allele').is_in(allowed)&(pl.col('cdr3').str.len_chars()>=8))
    if paired:
        from ..evalue.paired import build_paired_ref
        r=build_paired_ref(r)
        paired_observations=r.height
        for suffix,locus in [('alpha','TRA'),('beta','TRB')]:
            allowed_chain=[v for v,m in model.items() if m['locus']==locus]
            r=r.with_columns(resolve_v_alleles(r['v_'+suffix],model,locus=locus).alias('_allele_'+suffix)).filter(
                    pl.col('_allele_'+suffix).is_in(allowed_chain) & (pl.col(suffix).str.len_chars()>=8))
        search_keys=['alpha','beta','_allele_alpha','_allele_beta']
        if j_geometry:
            for column in ('j_alpha','j_beta'):
                r=r.with_columns(pl.col(column).fill_null('').str.split('*').list.first())
            for column in ('ja','jb'):
                q=q.with_columns(pl.col(column).fill_null('').str.split('*').list.first())
            search_keys += ['j_alpha','j_beta']
        r=r.unique(subset=[*search_keys,*PMHC,*(['j_alpha','j_beta'] if not j_geometry and (a.position_weighting=='ndn' or a.conditional_j) else [])],maintain_order=True)
    else:
        q=prepare(q.filter(pl.col('locus')==a.locus))
        r=prepare(r.filter(pl.col('reference_valid')))
        search_keys=['cdr3','_allele']
        if j_geometry:
            q=q.with_columns(pl.col('j').fill_null('').str.split('*').list.first())
            r=r.with_columns(pl.col('j').fill_null('').str.split('*').list.first())
            search_keys += ['j']
        r=r.unique(subset=[*search_keys,*PMHC,*(['j'] if not j_geometry and (a.position_weighting=='ndn' or a.conditional_j) else [])],maintain_order=True)
    if not paired and 'species' in metadata.columns:
        aliases=['human','homosapiens','homo sapiens'] if organism=='human' else ['mouse','musmusculus','mus musculus']
        q=q.join(metadata.filter(pl.col('species').str.strip_chars().str.to_lowercase().is_in(aliases)).select('query_id'),on='query_id',how='semi')
    if a.targets_from_sample:
        if a.epitope or a.mhc_a or a.mhc_b:
            raise ValueError('--targets-from-sample cannot be combined with explicit target selectors')
        from .historical import assayed_targets
        groups=assayed_targets(q.drop([c for c in ['species',*PMHC] if c in q.columns]),metadata,r,a.species)
    elif a.all_targets:
        if 'species' in metadata.columns:
            aliases=['human','homosapiens','homo sapiens'] if organism=='human' else ['mouse','musmusculus','mus musculus']
            q=q.join(metadata.filter(pl.col('species').str.strip_chars().str.to_lowercase().is_in(aliases)).select('query_id'),on='query_id',how='semi')
        groups=[(dict(zip(PMHC,key)),q,part) for key,part in r.group_by(PMHC,maintain_order=True)]
    elif kernel:
        selected=db.load(a.vdjdb,species=reference_species,gene=None if paired else a.locus,
                         epitope=a.epitope,mhc_a=a.mhc_a,mhc_b=a.mhc_b,mhc_match=a.mhc_match)
        target_keys=selected.select(PMHC).unique()
        targets=r.join(target_keys,on=PMHC,how='semi',nulls_equal=True)
        groups=[(dict(zip(PMHC,key)),q,part) for key,part in targets.group_by(PMHC,maintain_order=True)]
    else:
        groups=[(None,q,r)]
    radius=a.radius
    radius_rows=[]
    target_radii={}
    if a.reference_radius:
        from ..match.geometry import reference_nearest_distances,kneedle_radius,regularize_reference_radii
        radius_start=time.perf_counter()
        groups=list(groups)
        nearest_curves=[]
        for task,_,reference_part in groups:
            reference_keys=reference_part.unique(search_keys,maintain_order=True)
            fields=[]
            for i in range(2 if paired else 1):
                suffix=('alpha','beta')[i]
                fields.extend([pl.col(suffix if paired else 'cdr3').alias(f'seq{i}'),
                    pl.col('_allele_'+suffix if paired else '_allele').alias(f'v{i}')])
            distances=reference_nearest_distances(reference_keys.select(fields),paired=paired,
                species=organism,threads=a.threads,ctrim=ctrim,alphabet=a.alphabet,stem_depth=a.stem_depth,stem_weight=a.stem_weight,matrix=a.distance_matrix,v_loop_distances=v_loop_distances)
            finite=distances[np.isfinite(distances)]
            nearest_curves.append(finite)
            knee=kneedle_radius(finite)
            radius_rows.append({**task,'reference_keys':reference_keys.height,
                'raw_radius':knee,'nearest_keys':len(finite)})
        adjusted,floor,floor_source=regularize_reference_radii(
            [row['raw_radius'] if row['raw_radius'] is not None else np.nan for row in radius_rows],
            np.concatenate(nearest_curves) if nearest_curves else [],factor=a.knee_floor_factor)
        for row,value in zip(radius_rows,adjusted):
            # Round in matrix units, preserving the boundary under common rescaling.
            selected=int(np.ceil((value if np.isfinite(value) else a.radius)/bin_width))*bin_width
            row.update(radius=selected,radius_floor=floor,floor_source=floor_source,
                radius_status='no_knee_fixed_fallback' if row['raw_radius'] is None else
                    'degenerate_fixed_fallback' if not np.isfinite(value) else
                    'knee_floored' if value>row['raw_radius'] else 'knee')
            target_radii[tuple(row[c] for c in PMHC)]=selected
        radius_seconds=time.perf_counter()-radius_start
    if a.withhold_epitope:
        # Candidate masking must not refit complement populations or knee floors.
        groups=[(task,part,reference_part) for task,part,reference_part in groups
                if not a.all_targets or task['epitope'] not in a.withhold_epitope]
        if not a.all_targets:
            groups=[(task,part,reference_part.head(0) if task['epitope'] in a.withhold_epitope else reference_part)
                    for task,part,reference_part in groups]
    # Counts share one distance on references and controls; no dense control matrix.
    enrichment_rows=[]
    selection_state=selection_pending=None
    selection_peptide=None
    selection_hit_ids=set()
    rules_hit_ids=set()
    selection_dispositions={}
    background_provenance=None
    background_costs={'control_load':0.0,'native_count_batches':0.0}
    boundary_provenance=None
    boundary_map={}
    if a.background:
        background_start=time.perf_counter()
        groups=list(groups)  # assayed_targets is streamed; the new arm uses it twice
        def features(frame,query=False):
            fields=[pl.col('query_id')] if query else []
            for i in range(2 if paired else 1):
                suffix=('a','b')[i] if query else ('alpha','beta')[i]
                fields.extend([pl.col(('cdr3'+suffix if query else suffix) if paired else 'cdr3').alias(f'seq{i}'),
                    pl.col(('_allele'+suffix if query else '_allele_'+suffix) if paired else '_allele').alias(f'v{i}'),
                    pl.col(('j'+suffix if query else 'j_'+suffix) if paired else 'j').alias(f'j{i}')])
            return frame.select(fields)
        rule_batches=None
        if a.epitope_rules and (len(groups)>1 or rules_only):
            from ..evalue.epitope_rules import reference_epitope_rules,grouped_reference_epitope_rules
            rule_start=time.perf_counter()
            residuals=grouped_reference_epitope_rules(features(r).with_columns(r['epitope']),
                [(task,features(part,True),features(selected)) for task,part,selected in groups],
                paired=paired,exclude_exact=a.exclude_exact,species=organism,threads=a.threads,radius=a.radius,
                v_loop_distances=(v_loop_distances if a.distance_matrix=='clipped' and a.alphabet=='full' else None))
            rule_batches={tuple(task[c] for c in PMHC):residual for (task,_,_),residual in zip(groups,residuals)}
            background_costs['epitope_rules']=time.perf_counter()-rule_start
            if rules_only:
                for (task,part,_),residual in zip(groups,residuals):
                    if not part.height:continue
                    enrichment_rows.append(residual.with_columns(*[pl.lit(task[c],dtype=pl.String).alias(c) for c in PMHC]))
                    rules_hit_ids.update(residual.filter(pl.col('rule_reference_population')>
                        pl.col('rule_reference_excluded_hits'))['query_id'])
        from ..match.geometry import total_distance_count_batch,paired_total_distance_mass_batch
        from ..evalue.control import raw_gene_background
        linear=a.neighbour_weighting=='linear'
        local=a.neighbour_weighting=='local-rank' or kernel
        if kernel:
            from ..evalue.kernel import ladder_cutoffs,kernel_ladder
            kernel_cutoffs=ladder_cutoffs(matrix_scale,radius_factor=a.radius/matrix_scale)
        query_keys=['cdr3a','cdr3b','_allelea','_alleleb'] if paired else ['cdr3','_allele']
        if j_geometry:
            query_keys += ['ja','jb'] if paired else ['j']
        options=dict(species=organism,threads=a.threads,ctrim=ctrim,
                     position_weighting=a.position_weighting,exclude_exact=a.exclude_exact,linear_mass=linear,v_loop_weight=a.v_loop_weight,alphabet=a.alphabet,stem_depth=a.stem_depth,stem_weight=a.stem_weight,matrix=a.distance_matrix,v_loop_distances=v_loop_distances)
        def boundaries(sequences,genes):
            return [boundary_map[(s,v)] for s,v in zip(sequences,genes)]
        def counts(part,population,thresholds=None,*,same_v_gene=False):
            shared_thresholds=thresholds is None and kernel
            if thresholds is None:
                bound=int(round(200*radius))
                cutoffs=[200*t for t in kernel_cutoffs] if kernel else list(range(0,bound+1,200*bin_width))
                if not kernel and cutoffs[-1]!=bound:cutoffs.append(bound)
                thresholds=cutoffs if shared_thresholds else ([cutoffs]*part.height if local else [[bound]]*part.height)
            width=len(thresholds) if shared_thresholds else (len(thresholds[0]) if thresholds else int(np.floor(radius/bin_width))+1)
            if not part.height:
                return [],[],np.zeros((0,width),dtype=np.uint64) if local else None
            if not population.height:
                return [0]*part.height,[0.0]*part.height,np.zeros((part.height,width),dtype=np.uint64) if local else None
            count_start=time.perf_counter()
            if paired:
                native=paired_total_distance_mass_batch(
                    part['cdr3a'].to_list(),part['cdr3b'].to_list(),
                    population['alpha'].to_list(),population['beta'].to_list(),
                    part['_allelea'].to_list(),part['_alleleb'].to_list(),
                    population['_allele_alpha'].to_list(),population['_allele_beta'].to_list(),
                    thresholds,species=organism,threads=a.threads,same_v_gene=same_v_gene,shared_thresholds=shared_thresholds,v_loop_distances=v_loop_distances,
                    ctrim=ctrim,exclude_exact=a.exclude_exact,position_weighting=a.position_weighting,v_loop_weight=a.v_loop_weight,alphabet=a.alphabet,stem_depth=a.stem_depth,stem_weight=a.stem_weight,matrix=a.distance_matrix,
                    **({
                        'query_boundaries_alpha':boundaries(part['cdr3a'],part['_allelea']),
                        'query_boundaries_beta':boundaries(part['cdr3b'],part['_alleleb']),
                        'reference_boundaries_alpha':boundaries(population['alpha'],population['_allele_alpha']),
                        'reference_boundaries_beta':boundaries(population['beta'],population['_allele_beta'])
                    } if a.position_weighting=='ndn' else {}))
                curve=np.asarray(native[0],dtype=np.uint64) if local else None
                result=[row[-1] for row in native[0]]
                mass=[row[0]/native_radius for row in native[1]] if linear else result
            else:
                native=total_distance_count_batch(part['cdr3'].to_list(),population['cdr3'].to_list(),
                    part['_allele'].to_list(),population['_allele'].to_list(),
                    thresholds,**options,same_v_gene=same_v_gene,shared_thresholds=shared_thresholds,**({
                        'query_boundaries':boundaries(part['cdr3'],part['_allele']),
                        'reference_boundaries':boundaries(population['cdr3'],population['_allele'])
                    } if a.position_weighting=='ndn' else {}))
                curve=np.asarray(native,dtype=np.uint64) if local else None
                result=[row[-1] for row in (native[0] if linear else native)]
                mass=[row[0]/native_radius for row in native[1]] if linear else result
            background_costs['native_count_batches']+=time.perf_counter()-count_start
            return result,mass,curve
        def populations(part,population):
            # Exact exclusion changes hit counts only, never the reference model.
            return [population.height] * part.height
        shared=None
        controls=None
        if a.background!='vdjdb-other':
            load_start=time.perf_counter()
            controls,background_provenance=raw_gene_background(a.control,a.locus,organism,deduplicate=a.background!='generative',include_j=a.position_weighting=='ndn')
            background_costs['control_load']=time.perf_counter()-load_start
        if a.position_weighting=='ndn':
            from ..match.regions import annotate_ndn_boundaries
            from arda.cdr3fix import markup_batch
            from arda.paths import vdj_dir
            anchor_path=vdj_dir(organism)/'cdr3_anchors.tsv'
            if not anchor_path.is_file():raise ValueError('NDN geometry requires the ARDA anchor table')
            annotation_start=time.perf_counter()
            if paired:
                frames=[q.select(pl.col('cdr3'+s).alias('cdr3'),pl.col('_allele'+s).alias('_allele'),pl.col('j'+s).alias('j')) for s in ('a','b')]
                frames += [r.select(pl.col(s).alias('cdr3'),pl.col('_allele_'+s).alias('_allele'),pl.col('j_'+s).alias('j')) for s in ('alpha','beta')]
            else:
                frames=[q.select('cdr3','_allele','j'),r.select('cdr3','_allele','j')]
                if controls is not None:frames.append(controls.select('cdr3','_allele','j'))
            annotated=annotate_ndn_boundaries(pl.concat(frames,how='vertical_relaxed'),species=organism).unique(['cdr3','_allele'],maintain_order=True)
            boundary_map={(s,v):(ve,js) if status=='annotated' else None for s,v,ve,js,status in
                annotated.select('cdr3','_allele','v_end','j_start','boundary_status').iter_rows()}
            boundary_path=Path(str(a.output_prefix)+'.boundaries.tsv');boundary_path.parent.mkdir(parents=True,exist_ok=True)
            annotated.select('cdr3','_allele','v_end','j_start','boundary_status').write_csv(boundary_path,separator='\t')
            boundary_provenance={'producer':'arda.cdr3fix.markup_batch','version':version('arda-mapper'),
                'wrapper_source_sha256':sha256(Path(annotate_ndn_boundaries.__code__.co_filename)),
                'producer_source_sha256':sha256(Path(markup_batch.__code__.co_filename)),
                'anchor_table_sha256':sha256(anchor_path),
                'coordinate':'0-based half-open unchanged full junction; v_end count, j_start first J residue',
                'pair_substitution_weight':'max(aligned query weight, aligned reference weight)',
                'weights':{'V_J':0.1,'NDN':1.0},'fallback':'uniform1 on missing/conflicting/ambiguous/changed/unusable annotation; no row removed',
                'status_counts':dict(annotated.group_by('boundary_status').len().iter_rows()),
                'final_table_sha256':sha256(boundary_path),'annotation_seconds':time.perf_counter()-annotation_start}
            if controls is not None and a.background!='generative':
                controls=controls.unique(['cdr3','_allele'],maintain_order=True)
                background_provenance['population_size']=controls.height
        if a.background!='vdjdb-other':
            # One native control batch, shared transiently across declared targets.
            parts=[part.select('cdr3','_allele') for _,part,_ in groups if part.height]
            usable=pl.concat(parts).unique(maintain_order=True) if parts else q.select('cdr3','_allele').head(0)
            m,w,mc=counts(usable,controls)
            shared=usable.with_columns(pl.Series('_m',m,dtype=pl.UInt64),pl.Series('control_weight',w,dtype=pl.Float64),
                                       pl.Series('_M',populations(usable,controls),dtype=pl.UInt64))
            if local:
                shared=shared.with_columns(pl.Series('_rank_control',mc))
            del mc
        else:
            background_provenance={'kind':'vdjdb_other_epitopes','reference':db.provenance(a.vdjdb),
                'target_removal':'all records of target peptide, before key deduplication',
                'shared_key_policy':'other epitope annotations retained; no arbitrary label selection',
                'observation_unit':('distinct linked full alpha/beta junctions + resolved V alleles' if paired else
                                    'distinct full junction + resolved V allele')+
                                    (' + literal J families' if j_geometry else '')+'; no abundance weights',
                'search_keys':search_keys,
                'restriction':'same species/locus; other epitopes may have other MHC restrictions'}
        grouped_kernel=(selection_only and kernel and a.background=='vdjdb-other' and
            a.position_weighting in ('uniform','significance-trimmed') and not (a.conditional_j or a.epitope_rules or getattr(a,'chemistry',False)))
        grouped_detailed=_grouped_detailed_eligible(a,kernel,selection_only,groups,r)
        if evidence_only and not grouped_detailed:
            raise ValueError('--evidence-only requires grouped detailed kernel eligibility and non-null target epitopes')
        global_curves=None
        if not (grouped_kernel or grouped_detailed) and kernel and a.background=='vdjdb-other' and len(groups)>1:
            # All-label targets share the same query population. Concatenating
            # that population once per target creates a query × target temporary.
            if not a.targets_from_sample:
                global_queries=q.select(query_keys).unique(maintain_order=True)
            else:
                parts=[part.select(query_keys) for _,part,_ in groups if part.height]
                global_queries=(pl.concat(parts).unique(maintain_order=True) if parts else
                                q.select(query_keys).head(0))
            global_queries=global_queries.with_row_index('_global_row')
            whole=r.unique(search_keys,maintain_order=True)
            _,_,global_curves=counts(global_queries,whole)
            global_sizes=np.asarray(populations(global_queries,whole),dtype=np.uint64)
            spans=r.group_by(search_keys).agg(pl.col('epitope').n_unique().alias('_epitopes'))
            background_provenance['reduction']='whole search-key counts minus keys exclusive to target peptide; shared-label keys retained'
        if selection_only:
            # Preserve the detailed route's final disposition order before
            # sorting restrictions so one peptide can be completed at a time.
            for task,part,reference_part in (groups[-1:] if a.all_targets else groups):
                disposition=('invalid_target' if any(not isinstance(v,str) or not v.strip() for v in task.values())
                             else 'searched' if reference_part.height else 'no_reference')
                selection_dispositions.update(dict.fromkeys(part['query_id'],disposition))
            if a.all_targets and not groups:
                selection_dispositions=dict.fromkeys(q['query_id'],'no_reference')
            groups.sort(key=lambda group: group[0]['epitope'] or '')
        if grouped_kernel:
            from ..evalue.kernel import grouped_kernel_selection
            count_start=time.perf_counter()
            selection_state,selection_hit_ids=grouped_kernel_selection(q,r,groups,
                query_keys=query_keys,search_keys=search_keys,paired=paired,
                cutoffs=kernel_cutoffs,radius=a.radius,kernel_scale=22.5*matrix_scale,
                scorer=a.kernel_scorer,assayed=a.targets_from_sample,
                geometry_preset=geometry_extra,
                progress=getattr(a,'verbose',False),
                **{k:v for k,v in options.items() if k!='linear_mass'})
            background_costs['native_count_batches']+=time.perf_counter()-count_start
            background_provenance['reduction']='one grouped reference traversal; whole search-key counts minus keys exclusive to target peptide; shared-label keys retained'
        if grouped_detailed:
            from ..evalue.kernel import grouped_kernel_enrichment,grouped_kernel_selection
            count_start=time.perf_counter()
            detailed=(grouped_kernel_selection if evidence_only else grouped_kernel_enrichment)
            result=detailed(q,r,groups,query_keys=query_keys,search_keys=search_keys,
                paired=paired,cutoffs=kernel_cutoffs,radius=a.radius,kernel_scale=22.5*matrix_scale,
                scorer=a.kernel_scorer,assayed=a.targets_from_sample,
                progress=getattr(a,'verbose',False),
                **({'_detailed':True,'_return_hits':True} if evidence_only else {}),
                **{k:v for k,v in options.items() if k!='linear_mass'})
            if evidence_only:
                evidence,selection_hit_ids=result
            else:
                evidence=result
            if evidence is not None:
                enrichment_rows.append(evidence)
            background_costs['native_count_batches']+=time.perf_counter()-count_start
            background_provenance['reduction']='one grouped reference traversal; whole search-key counts minus keys exclusive to target peptide; shared-label keys retained'
        population_queries=(q.select(query_keys).unique(maintain_order=True)
                            if kernel and not a.targets_from_sample else None)
        for task,part,reference_part in ([] if rules_only or grouped_kernel or grouped_detailed else groups):
            radius=target_radii[tuple(task[c] for c in PMHC)] if a.reference_radius else a.radius
            unique=(population_queries if population_queries is not None else
                    part.select(query_keys).unique(maintain_order=True))
            if not unique.height:
                continue
            reference_keys=reference_part.unique(search_keys,maintain_order=True)
            n,reference_weight,nc=counts(unique,reference_keys)
            N=populations(unique,reference_keys)
            if shared is not None:
                evidence=unique.join(shared,on=query_keys,how='left',validate='1:1',maintain_order='left')
            else:
                if global_curves is not None:
                    only=(r.filter(pl.lit(False) if task['epitope'] is None else pl.col('epitope')==task['epitope']).join(spans,on=search_keys,how='left')
                          .filter(pl.col('_epitopes')==1).unique(search_keys,maintain_order=True))
                    _,_,removed=counts(unique,only)
                    indices=unique.join(global_queries,on=query_keys,how='left',validate='1:1',
                                        maintain_order='left')['_global_row'].to_numpy()
                    mc=global_curves[indices]-removed
                    M=global_sizes[indices]-np.asarray(populations(unique,only),dtype=np.uint64)
                    m=mc[:,-1];w=m
                else:
                    other=r.filter(pl.col('epitope')!=task['epitope']).unique(search_keys,maintain_order=True)
                    m,w,mc=counts(unique,other)
                    M=populations(unique,other)
                evidence=unique.with_columns(pl.Series('_m',m,dtype=pl.UInt64),pl.Series('control_weight',w,dtype=pl.Float64),
                                            pl.Series('_M',M,dtype=pl.UInt64))
                if local:
                    evidence=evidence.with_columns(pl.Series('_rank_control',mc))
                del mc
            evidence=evidence.with_columns(pl.Series('n_reference',n,dtype=pl.UInt64),
                                          pl.Series('reference_weight',reference_weight,dtype=pl.Float64),
                                          pl.Series('reference_population',N,dtype=pl.UInt64))
            if a.reference_radius:
                evidence=evidence.with_columns(pl.lit(radius,dtype=pl.UInt32).alias('selected_radius'))
            evidence=evidence.rename({'_m':'n_control','_M':'control_population'})
            available=(pl.col('reference_population')>0)&(pl.col('control_population')>0)
            evidence=evidence.with_columns(
                pl.when(available).then(pl.col('reference_population')*(pl.col('control_weight')+.5)/(pl.col('control_population')+1)).alias('expected_count'),
                pl.when(available).then(pl.lit('available')).otherwise(pl.lit('empty_population')).alias('enrichment_status'))
            evidence=evidence.with_columns((pl.col('reference_weight')/pl.col('expected_count')).alias('enrichment'),
                (pl.col('n_control')==0).alias('zero_control_hits'),(pl.col('control_weight')==0).alias('zero_control_weight'))
            if kernel:
                _,_,same=counts(unique,reference_keys,same_v_gene=True)
                scores=kernel_ladder(nc,same,evidence['_rank_control'].to_numpy(),N,
                    evidence['control_population'].to_numpy(),kernel_cutoffs,kernel_scale=22.5*matrix_scale,
                    scorer=a.kernel_scorer)
                seq_cols=['cdr3a','cdr3b'] if paired else ['cdr3']
                ref_seq_cols=['alpha','beta'] if paired else ['cdr3']
                exact=reference_keys.select(*[pl.col(r).alias(q) for q,r in zip(seq_cols,ref_seq_cols)]).group_by(seq_cols).len()
                exact_counts=unique.join(exact,on=seq_cols,how='left',maintain_order='left')['len'].fill_null(0)
                if selection_only:
                    radius_counts=(nc[:,kernel_cutoffs.index(radius)] if radius in kernel_cutoffs else
                        counts(unique,reference_keys,[[distance_cutoff(radius)]]*unique.height)[2][:,0])
                    evidence=evidence.with_columns(pl.Series('_radius_hit',radius_counts>0))
                evidence=evidence.drop('_rank_control').with_columns(
                    (pl.lit(0,dtype=pl.UInt64) if a.exclude_exact else exact_counts.cast(pl.UInt64)).alias('exact_reference_count'),
                    *(pl.Series(k,v,nan_to_null=True) for k,v in scores.items()),
                    pl.lit('relative_evidence_uncalibrated').alias('kernel_status'))
            elif local:
                from ..evalue.local_rank import local_rank_evidence
                rank=local_rank_evidence(nc,evidence['_rank_control'].to_numpy(),N,
                                         evidence['control_population'].to_numpy(),prior=a.rank_signal_prior)
                evidence=evidence.drop('_rank_control').with_columns(
                    *(pl.Series(k,v) for k,v in rank.items()),
                    pl.lit(a.rank_signal_prior).alias('rank_signal_prior'),
                    pl.when(available).then(pl.lit('conditional_rank_model_only')).otherwise(pl.lit('empty_population')).alias('rank_model_status'))
                evidence=evidence.with_columns(
                    pl.when(available).then(pl.col('rank_bayes_factor')).alias('rank_bayes_factor'),
                    *(pl.when(available).then(pl.col(c)).alias(c) for c in ['p_rank_bound','p_nearest_rank']))
            del nc
            expanded=part.select('query_id',*query_keys).join(evidence,
                on=query_keys,how='left',validate='m:1',maintain_order='left').drop(query_keys)
            if a.conditional_j or a.epitope_rules or getattr(a,'chemistry',False):
                j_start=time.perf_counter()
                if a.conditional_j:
                    from ..evalue.germline import conditional_j_usage
                    residual=conditional_j_usage(features(part,True),features(reference_part),
                        features(r.filter(pl.col('epitope')!=task['epitope'])),paired=paired,exclude_exact=a.exclude_exact)
                elif a.epitope_rules:
                    from ..evalue.epitope_rules import reference_epitope_rules
                    residual=(rule_batches[tuple(task[c] for c in PMHC)] if rule_batches is not None else
                        reference_epitope_rules(features(part,True),features(reference_part),
                        features(r.filter(pl.col('epitope')!=task['epitope'])),paired=paired,exclude_exact=a.exclude_exact,
                        species=organism,threads=a.threads,radius=a.radius,
                        v_loop_distances=(v_loop_distances if a.distance_matrix=='clipped' and a.alphabet=='full' else None)))
                else:
                    from ..match.chemistry import local_chemistry_evidence
                    residual=local_chemistry_evidence(features(part,True),features(reference_part),
                        features(r.filter(pl.col('epitope')!=task['epitope'])),paired=paired,
                        exclude_exact=a.exclude_exact,bandwidth_factor=a.chemistry_bandwidth_factor)
                expanded=expanded.join(residual,on='query_id',how='left',validate='1:1',maintain_order='left')
                field='conditional_j_usage' if a.conditional_j else 'epitope_rules' if a.epitope_rules else 'chemistry'
                background_costs[field]=background_costs.get(field,0)+time.perf_counter()-j_start
            if selection_only:
                selection_hit_ids.update(expanded.filter(pl.col('_radius_hit'))['query_id'])
                from ..evalue.kernel import kernel_peptides,update_kernel_selection
                if task['epitope'] not in (a.withhold_epitope or []):
                    peptide=task['epitope']
                    if selection_pending is not None and peptide!=selection_peptide:
                        selection_state=update_kernel_selection(selection_state,selection_pending)
                        selection_pending=None
                    current=expanded.select('query_id','kernel_score','exact_reference_count').with_columns(
                        pl.lit(peptide,dtype=pl.String).alias('epitope')).select(
                            'query_id','epitope','kernel_score','exact_reference_count')
                    selection_pending=kernel_peptides(current if selection_pending is None else
                        pl.concat([selection_pending,current],how='vertical_relaxed'))
                    selection_peptide=peptide
            else:
                enrichment_rows.append(expanded.with_columns(*(pl.lit(task[c],dtype=pl.String).alias(c) for c in PMHC)))
        if selection_only and selection_pending is not None:
            selection_state=update_kernel_selection(selection_state,selection_pending)
            selection_pending=None
        background_costs['total']=time.perf_counter()-background_start
    rows=[]
    shortlist_rows=[];shortlist_costs=[]
    dispositions=selection_dispositions if selection_only else ({query_id:'no_reference' for query_id in q['query_id']} if a.all_targets and not groups else {})
    for task,query_part,reference_part in ([] if selection_only else groups):
        radius=target_radii[tuple(task[c] for c in PMHC)] if a.reference_radius else a.radius
        if task is not None:
            reference_part=reference_part.unique(search_keys,maintain_order=True)
            for query_id in query_part['query_id']:
                dispositions[query_id]=('invalid_target' if any(not isinstance(v,str) or not v.strip() for v in task.values())
                                        else 'searched' if reference_part.height else 'no_reference')
        if rules_only or evidence_only:continue
        if not reference_part.height or not query_part.height:
            continue
        if a.shortlist_audit:
            from ..match.geometry import shortlist_audit_batch
            computation=query_part.unique(['cdr3','_allele'],maintain_order=True)
            audit=shortlist_audit_batch(computation['cdr3'],reference_part['cdr3'],
                computation['_allele'],reference_part['_allele'],radius=a.radius,
                species=organism,threads=a.threads,ctrim=ctrim,exclude_exact=a.exclude_exact,alphabet=a.alphabet,stem_depth=a.stem_depth,stem_weight=a.stem_weight,matrix=a.distance_matrix,v_loop_distances=v_loop_distances)
            shared=computation.select('cdr3','_allele').with_columns(
                pl.Series('cdr3_shortlist_keys',audit['cdr3_shortlist_keys'],dtype=pl.UInt64),
                pl.Series('post_v_radius_keys',audit['post_v_radius_keys'],dtype=pl.UInt64),
                pl.Series('nearest_within_shortlist_distance',audit['nearest_within_shortlist_distance'],dtype=pl.Float64))
            shortlist_rows.append(query_part.select('query_id','cdr3','_allele').join(
                shared,on=['cdr3','_allele'],validate='m:1',maintain_order='left').drop('cdr3','_allele'))
            shortlist_costs.append(dict(target=task,query_appearances=query_part.height,
                unique_computations=computation.height,**{k:audit[k] for k in ['stage_seconds','projected_sparse_payload_bytes','route','native_output_bytes']}))
        # Paired calls hold two matrices plus native/gather temporaries; bound them together.
        matrix_budget=32*1024**2 if paired else 64*1024**2
        bytes_per_pair=16 if paired else 4
        if bytes_per_pair*reference_part.height>matrix_budget:
            raise ValueError('reference exceeds the matrix budget for one query; use a narrower target')
        batch=max(1,matrix_budget//(bytes_per_pair*reference_part.height))
        if not paired:
            rs,rv=reference_part['cdr3'].to_list(),reference_part['_allele'].to_list()
        for part in query_part.iter_slices(batch):
            if paired:
                qa,qb=part['cdr3a'].to_list(),part['cdr3b'].to_list()
                ra,rb=reference_part['alpha'].to_list(),reference_part['beta'].to_list()
                distances=paired_distance_matrix(qa,qb,ra,rb,part['_allelea'].to_list(),part['_alleleb'].to_list(),
                    reference_part['_allele_alpha'].to_list(),reference_part['_allele_beta'].to_list(),species=organism,threads=a.threads,ctrim=ctrim,v_loop_weight=a.v_loop_weight,alphabet=a.alphabet,stem_depth=a.stem_depth,stem_weight=a.stem_weight,matrix=a.distance_matrix,v_loop_distances=v_loop_distances)
            else:
                seqs=part['cdr3'].to_list()
                distances=distance_matrix(seqs,rs,part['_allele'].to_list(),rv,species=organism,threads=a.threads,ctrim=ctrim,v_loop_weight=a.v_loop_weight,alphabet=a.alphabet,stem_depth=a.stem_depth,stem_weight=a.stem_weight,matrix=a.distance_matrix,v_loop_distances=v_loop_distances)
            keep=distances<=radius
            nonexact=True
            if a.exclude_exact:
                if paired:
                    nonexact=~((np.asarray(qa)[:,None]==np.asarray(ra)[None,:]) &
                               (np.asarray(qb)[:,None]==np.asarray(rb)[None,:]))
                else:
                    nonexact=np.asarray(seqs)[:,None]!=np.asarray(rs)[None,:]
                keep&=nonexact
            if a.unbounded_nearest:
                sentinel=np.iinfo(np.int32).max
                nearest=np.min(distances,axis=1,where=nonexact,initial=sentinel)
                available=nearest!=sentinel
                h=part.filter(pl.Series(available)).select('query_id').with_columns(
                    pl.Series('distance',nearest[available]),
                    pl.Series('n_neighbours',keep.sum(axis=1,dtype=np.uint32)[available]))
                rows.append(h.with_columns(*(pl.lit(task[c],dtype=pl.String).alias(c) for c in PMHC)))
                continue
            qi,ri=np.nonzero(keep)
            h=pl.DataFrame({'_qi':qi,'_ri':ri,'distance':distances[qi,ri]})
            h=h.join(part.with_row_index('_qi').select('_qi','query_id'),on='_qi')
            if task is None:
                h=h.join(reference_part.with_row_index('_ri').select('_ri',*PMHC),on='_ri')
            else:
                h=h.with_columns(*(pl.lit(task[c],dtype=pl.String).alias(c) for c in PMHC))
            rows.append(h.group_by('query_id',*PMHC).agg(pl.col('distance').min(),pl.len().alias('n_neighbours')))
    schema={'query_id':q.schema['query_id'],**{c:pl.String for c in PMHC},'distance':pl.Int32,'n_neighbours':pl.UInt32}
    out=pl.concat(rows) if rows else pl.DataFrame(schema=schema)
    out=out.sort('query_id','distance',*PMHC,nulls_last=True)
    if a.nearest_rank:
        from scipy.stats import hypergeom
        rank_start=time.perf_counter();nearest_rows=[]
        for task,part,_ in groups:
            unique=part.unique(query_keys,maintain_order=True)
            target_distances=out.filter(pl.all_horizontal([pl.col(c).eq_missing(pl.lit(task[c],dtype=pl.String)) for c in PMHC]))
            nearest=unique.join(target_distances.select('query_id','distance'),on='query_id',how='inner',validate='1:1')
            if not nearest.height:continue
            population=(r.filter(pl.col('epitope')!=task['epitope']).unique(search_keys,maintain_order=True)
                        if a.background=='vdjdb-other' else controls)
            thresholds=[[200*int(d)] for d in nearest['distance']]
            m,_,_=counts(nearest,population,thresholds=thresholds)
            values=nearest.select(*query_keys,'query_id').with_columns(pl.Series('nearest_rank_control_count',m,dtype=pl.UInt64))
            nearest_rows.append(part.select('query_id',*query_keys).join(values.drop('query_id'),on=query_keys,how='inner',validate='m:1').drop(query_keys).with_columns(*[pl.lit(task[c],dtype=pl.String).alias(c) for c in PMHC]))
        nearest_fields=(pl.concat(nearest_rows) if nearest_rows else pl.DataFrame(schema={
            'query_id':q.schema['query_id'],'nearest_rank_control_count':pl.UInt64,**{c:pl.String for c in PMHC}}))
        for i,evidence in enumerate(enrichment_rows):
            evidence=evidence.join(nearest_fields,on=['query_id',*PMHC],how='left',validate='1:1',nulls_equal=True)
            usable=(evidence['reference_population']>0)&(evidence['control_population']>0)&evidence['nearest_rank_control_count'].is_not_null()
            tail=np.full(evidence.height,np.nan)
            if usable.any():
                # SciPy support arithmetic must not underflow on unsigned counts.
                n=evidence.filter(usable);N=n['reference_population'].cast(pl.Float64).to_numpy();M=n['control_population'].cast(pl.Float64).to_numpy()
                tail[usable.to_numpy()]=hypergeom.sf(0,N+M,N,n['nearest_rank_control_count'].cast(pl.Float64).to_numpy()+1)
            if a.neighbour_weighting=='local-rank':
                from ..evalue.local_rank import nearest_conditioned_rank
                conditional=np.full(evidence.height,np.nan)
                if usable.any():
                    conditional[usable.to_numpy()]=nearest_conditioned_rank(
                        n['rank_bayes_factor'].to_numpy(),n['n_reference'].to_numpy(),N,M,
                        n['nearest_rank_control_count'].to_numpy())
                evidence=evidence.with_columns(pl.Series('conditional_rank_bayes_factor',conditional,nan_to_null=True))
                if a.conditional_j:
                    evidence=evidence.with_columns((pl.col('conditional_rank_bayes_factor')*pl.col('conditional_j_usage_ratio')).alias('j_conditioned_rank_score'))
            enrichment_rows[i]=evidence.with_columns(pl.Series('p_nearest_unbounded',tail,nan_to_null=True)).with_columns(
                (-pl.col('p_nearest_unbounded')).alias('nearest_rank_score'),
                pl.when(pl.col('p_nearest_unbounded').is_not_null()).then(pl.lit('conditional_rank_model_only'))
                .when(pl.col('enrichment_status')=='empty_population').then(pl.lit('empty_population'))
                .otherwise(pl.lit('no_nonexcluded_reference')).alias('nearest_rank_status'))
        background_costs['unbounded_nearest_rank']=time.perf_counter()-rank_start
        background_costs['total']+=background_costs['unbounded_nearest_rank']
    if paired:
        status=pair_rows.select(*[c for c in ['query_id','pair_id','clone_id','query_ida','query_idb',
            'sequence_ida','sequence_idb','species',*PMHC] if c in pair_rows.columns],'available','status',
            'v_call_alpha','v_call_beta',
            pl.col('_allelea').alias('resolved_v_allele_alpha'),pl.col('_alleleb').alias('resolved_v_allele_beta'))
        hit_ids=selection_hit_ids if selection_only or evidence_only else rules_hit_ids if rules_only else set(out['query_id'])
        status=status.with_columns(pl.Series('status',[
            old if old!='searched' else dispositions.get(query_id) if dispositions.get(query_id) in {'invalid_target','no_reference'} else
            'matched' if query_id in hit_ids else 'no_neighbours'
            for query_id,old in zip(status['query_id'],status['status'])]))
    else:
        status=raw.select([c for c in ('query_id','sequence_id') if c in raw.columns]).join(q.select('query_id').with_columns(pl.lit(True).alias('available')),on='query_id',how='left')
    if not paired:
        v_column=io.columns._resolve(raw).get('v')
        raw_calls=raw[v_column] if v_column else pl.Series([None]*raw.height,dtype=pl.String)
        status=status.with_columns(raw_calls.alias('v_call'),
            resolve_v_alleles(raw_calls,model,locus=a.locus).alias('resolved_v_allele'))
    status=status.with_columns(pl.col('available').fill_null(False))
    if (a.targets_from_sample or a.all_targets) and not paired:
        valid_target_ids=set(dispositions)
        hit_ids=selection_hit_ids if selection_only or evidence_only else rules_hit_ids if rules_only else set(out['query_id'])
        species_ids=set(raw.filter(pl.col('species').str.strip_chars().str.to_lowercase().is_in(
            ['human','homosapiens','homo sapiens'] if organism=='human' else ['mouse','musmusculus','mus musculus']))['query_id']) if 'species' in raw.columns else set(raw['query_id'])
        status=status.with_columns(pl.Series('status',[
            'unselected_species' if query_id not in species_ids else
            'invalid_query_or_gene' if query_id not in valid_target_ids else
            dispositions[query_id] if dispositions[query_id] in {'no_reference','invalid_target'} else
            'matched' if query_id in hit_ids else 'no_neighbours'
            for query_id in status['query_id']]))
        status=status.with_columns((pl.col('available') & (pl.col('status')!='unselected_species')).alias('available'))
    prefix=Path(a.output_prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    if a.reference_radius:
        radius_schema={**{c:pl.String for c in PMHC},'reference_keys':pl.UInt64,
            'raw_radius':pl.Float64,'nearest_keys':pl.UInt64,'radius':pl.UInt32,
            'radius_floor':pl.Float64,'floor_source':pl.String,'radius_status':pl.String}
        pl.DataFrame(radius_rows,schema=radius_schema).write_csv(str(prefix)+'.radii.tsv',separator='\t')
    if not (selection_only or rules_only or evidence_only):
        out.write_csv(str(prefix)+'.candidates.tsv',separator='\t')
    status.write_csv(str(prefix)+'.queries.tsv',separator='\t')
    if a.shortlist_audit:
        diagnostic=pl.concat(shortlist_rows) if shortlist_rows else pl.DataFrame(schema={
            'query_id':q.schema['query_id'],'cdr3_shortlist_keys':pl.UInt64,
            'post_v_radius_keys':pl.UInt64,'nearest_within_shortlist_distance':pl.Float64})
        status.join(raw.select('query_id',*PMHC),on='query_id',validate='1:1').join(
            diagnostic,on='query_id',how='left',validate='1:1').write_csv(str(prefix)+'.shortlist.tsv',separator='\t')
    if selection_only:
        from ..evalue.kernel import finish_kernel_selection
        selected=finish_kernel_selection(status,selection_state)
        selected.write_csv(str(prefix)+'.selection.tsv',separator='\t')
    if a.background and not selection_only:
        schema={'query_id':q.schema['query_id'],**{c:pl.String for c in PMHC},
                **{c:pl.UInt64 for c in ['n_reference','reference_population','n_control','control_population']},
                'reference_weight':pl.Float64,'control_weight':pl.Float64,
                'expected_count':pl.Float64,'enrichment_status':pl.String,'enrichment':pl.Float64,
                'zero_control_hits':pl.Boolean,'zero_control_weight':pl.Boolean}
        enrichment=pl.concat(enrichment_rows) if enrichment_rows else pl.DataFrame(schema=schema)
        if kernel and not enrichment_rows:
            enrichment=enrichment.with_columns(*[pl.lit(None,dtype=pl.Float64).alias('kernel_'+w+'_'+d)
                for w in ('all','v','vonly') for d in ('bounded','rare')],
                pl.lit(None,dtype=pl.Float64).alias('kernel_score'),pl.lit(None,dtype=pl.String).alias('kernel_status'),
                pl.lit(None,dtype=pl.UInt64).alias('exact_reference_count'))
        if a.neighbour_weighting=='local-rank' and not enrichment_rows:
            enrichment=enrichment.with_columns(
                *(pl.lit(None,dtype=pl.Float64).alias(k) for k in ['rank_bayes_factor','p_rank_bound','p_nearest_rank','posterior_rank_signal','rank_signal_prior']),
                pl.lit(None,dtype=pl.String).alias('rank_model_status'))
        if a.nearest_rank and not enrichment_rows:
            enrichment=enrichment.with_columns(
                pl.lit(None,dtype=pl.UInt64).alias('nearest_rank_control_count'),
                *(pl.lit(None,dtype=pl.Float64).alias(k) for k in ['p_nearest_unbounded','nearest_rank_score']),
                pl.lit(None,dtype=pl.String).alias('nearest_rank_status'))
            if a.neighbour_weighting=='local-rank':
                enrichment=enrichment.with_columns(pl.lit(None,dtype=pl.Float64).alias('conditional_rank_bayes_factor'))
                if a.conditional_j:
                    enrichment=enrichment.with_columns(pl.lit(None,dtype=pl.Float64).alias('j_conditioned_rank_score'),
                        pl.lit(None,dtype=pl.Float64).alias('conditional_j_usage_ratio'),pl.lit(None,dtype=pl.String).alias('conditional_j_status'),
                        *[pl.lit(None,dtype=pl.Float64).alias('j_'+side+'_probability') for side in ['reference','control']],
                        *[pl.lit(None,dtype=pl.Int64).alias('j_'+side+'_'+field) for side in ['reference','control'] for field in ['population','stratum','joint','global']])
        if a.epitope_rules and not enrichment.height:
            from ..evalue.epitope_rules import reference_epitope_rules
            schema={'seq'+str(i):pl.String for i in range(2 if paired else 1)}
            schema.update({'v'+str(i):pl.String for i in range(2 if paired else 1)})
            pop=pl.DataFrame(schema=schema)
            empty=reference_epitope_rules(pop.with_row_index('query_id'),pop,pop,paired=paired,radius=a.radius)
            enrichment=enrichment.with_columns(*[pl.Series(c,[],dtype=t) for c,t in empty.schema.items() if c!='query_id'])
        if getattr(a,'chemistry',False) and not enrichment.height:
            from ..match.chemistry import local_chemistry_evidence
            pop=pl.DataFrame(schema={c+str(i):pl.String for c in ('seq','v') for i in range(2 if paired else 1)})
            empty=local_chemistry_evidence(pop.with_row_index('query_id'),pop,pop,paired=paired,
                bandwidth_factor=a.chemistry_bandwidth_factor)
            enrichment=enrichment.with_columns(*[pl.Series(c,[],dtype=t) for c,t in empty.schema.items() if c!='query_id'])
        if rules_only:
            enrichment=enrichment.select('query_id',*PMHC,*[c for c in enrichment.columns if c.startswith('rule_')])
        enrichment.sort('query_id',*PMHC).write_csv(str(prefix)+'.enrichment.tsv',separator='\t')
        if kernel:
            from ..evalue.kernel import select_kernel
            selected=select_kernel(status,enrichment.filter(~pl.col('epitope').is_in(a.withhold_epitope or [])))
            selected.write_csv(str(prefix)+'.selection.tsv',separator='\t')
    if paired and kernel:
        from copy import copy
        # Independent full-reference chain axes retain unlinked support; no joint P-value.
        combined=None if selection_only else enrichment.select('query_id',*PMHC,pl.col('kernel_score').alias('kernel_linked_score'))
        selected=selected.rename({'selected_epitope':'selected_linked_epitope',
                                  'selection_status':'linked_selection_status'})
        for suffix,locus,name in [('a','TRA','alpha'),('b','TRB','beta')]:
            marginal=copy(a)
            marginal.locus=locus
            marginal.output_prefix=str(prefix)+'.'+locus
            if marginal_tables is None:
                main(marginal,v_loop_distances=v_loop_distances)
            else:
                lane=marginal_tables[locus]
                required={'queries','selection'} if selection_only else {'queries','selection','enrichment'}
                if (selection_only and set(lane)!=required) or not required.issubset(lane):
                    raise ValueError('union marginals require query/selection and detailed evidence tables')
                for kind,frame in lane.items():
                    frame.write_csv(marginal.output_prefix+'.'+kind+'.tsv',separator='\t')
            mapping=pair_rows.select('query_id',pl.col('query_id'+suffix).alias('_chain_id'))
            if not selection_only:
                lane=pl.read_csv(marginal.output_prefix+'.enrichment.tsv',separator='\t',
                                 schema_overrides={'query_id':raw.schema['query_id']})
                lane=mapping.join(lane.rename({'query_id':'_chain_id'}),on='_chain_id',how='inner').select(
                    'query_id',*PMHC,pl.col('kernel_score').alias('kernel_'+name+'_score'))
                combined=combined.join(lane,on=['query_id',*PMHC],how='full',coalesce=True,nulls_equal=True)
            calls=pl.read_csv(marginal.output_prefix+'.selection.tsv',separator='\t',
                              schema_overrides={'query_id':raw.schema['query_id']})
            calls=mapping.join(calls.rename({'query_id':'_chain_id'}),on='_chain_id',how='left').select(
                'query_id',pl.col('selected_epitope').alias('selected_'+name+'_epitope'),
                pl.col('selection_status').alias(name+'_selection_status'),
                pl.col('kernel_score').cast(pl.Float64).alias('kernel_'+name+'_score'),
                pl.col('kernel_competing_score').cast(pl.Float64).alias('kernel_'+name+'_competing_score'),
                pl.col('tied_peptides').alias(name+'_tied_peptides'),
                pl.col('exact_peptides').alias(name+'_exact_peptides'))
            selected=selected.join(calls,on='query_id',how='left',validate='1:1',maintain_order='left')
        if not selection_only:
            combined.sort('query_id',*PMHC).write_csv(str(prefix)+'.paired-evidence.tsv',separator='\t')
        selected.with_columns(
            pl.max_horizontal('kernel_score','kernel_alpha_score','kernel_beta_score').alias('kernel_any_axis_score'),
            pl.any_horizontal(*[
            pl.col(c).is_in(['ranked_evidence','ambiguous']).fill_null(False) for c in
            ['linked_selection_status','alpha_selection_status','beta_selection_status']]).alias('has_support')).write_csv(
                str(prefix)+'.selection.tsv',separator='\t')
    if paired and kernel:
        from .calibration import paired_call_event,paired_pipeline_sources
    manifest={'sample_sha256':sha256(Path(a.sample)),
        'selection_only':selection_only,
        **({'evidence_only':True} if evidence_only else {}),
        **({'rules_only':True,'status_derivation':'preserve ingestion/target precedence; matched iff fixed reference population exceeds excluded exact hits; availability unchanged'} if rules_only else {}),
        **({'paired_call_event':paired_call_event(),'paired_source_sha256':paired_pipeline_sources()} if paired and kernel else {}),
        **({'geometry_extra':geometry_extra,
            'geometry_extra_status':'advanced retained diagnostic',
            'normalization':('fixed V-weighted reference population'
                if geometry_extra in ('norm','jnorm') else 'fixed reference population'),
            **({'paired_normalization':'extension using the joint resolved V-family tuple; no historical joint validation'}
               if paired and geometry_extra in ('norm','jnorm') else {}),
            'calibration_transfer':'requires its own geometry and all-label selection calibration'} if geometry_extra else {}),
        **({'marginal_execution':'run-local union chain population services'} if marginal_tables is not None else {}),
        'output_tables':(['queries','enrichment'] if rules_only else ['queries','selection'] if selection_only else ['queries']+([] if evidence_only else ['candidates'])+
            (['enrichment','selection'] if kernel else ['enrichment'] if a.background else [])+
            (['paired-evidence'] if paired and kernel else [])),
        'software':{'vdjmatch':version('vdjmatch'),'seqtree':version('seqtree'),
                    'seqtree_native_sha256':sha256(Path(_core.__file__)),
                    'seqtree_python_source_sha256':sha256(Path(gapblock.__file__)),
                    'scorer_source_sha256':sha256(Path(distance_matrix.__code__.co_filename)),
                    **({'kernel_source_sha256':sha256(Path(kernel_ladder.__code__.co_filename))} if kernel else {}),
                    'cli_source_sha256':sha256(Path(__file__))},
        'reference':db.provenance(a.vdjdb),'reference_population':a.reference_population,'ingestion':ingestion,'reference_observations':original,'withheld_reference_chain_records':withheld_records,
        'v_loop_model':{'sha256':sha256(Path(str(files('vdjmatch.resources')/'vloops'/f'{organism}_v_loops.tsv'))),
                        'precomputation':'one exact aligned-loop table per locus and invocation; reused across population batches'},
        'reference_usable_keys':r.height,'query_rows':raw.height,'query_usable_rows':2*q.height if paired else q.height,
        'withholding':'candidate hits only; full background populations and geometry parameters retained',
        **({'query_pairs':pair_rows.height,'query_usable_pairs':q.height,'reference_complete_pairs':paired_observations,
            'query_status_counts':dict(status.group_by('status').len().iter_rows()),
            'reference_key_unit':'distinct linked alpha/beta junction+V-allele+pMHC keys',
            'distance_unit':'alpha+beta sum to the same linked reference pair',
            'exact_exclusion':'hit-list filter when both full junctions equal; model unchanged',
            'matrix_work_budget_bytes':32*1024**2,
            'independent_marginal_baseline':('full compatible TRA/TRB reference axes in separate outputs' if kernel else 'not emitted')} if paired else {}),
        **({'kernel':{'score_column':'kernel_score','scorer':a.kernel_scorer,
            'source_score_column':'kernel_'+a.kernel_scorer.replace('-','_'),
            'matrix':a.distance_matrix,'cutoffs':kernel_cutoffs,
            'count_cutoff':kernel_cutoffs[-1],'count_columns':'n_reference and n_control count through the ladder reach',
            'kernel_scale':22.5*matrix_scale,
            'cross_v':.25 if a.kernel_scorer.startswith('v-') else 1. if a.kernel_scorer.startswith('all-') else 0.,
            'denominator':'1+control_count' if a.kernel_scorer.endswith('bounded') else 'max(control_count,0.01)',
            'normalization':'fixed target population; exact hits filtered only','interpretation':'relative evidence; not a calibrated P/E or binding probability',
            'paired':'linked alpha+beta sum plus separate full-reference chain scores; no joint probability or fitted fusion'}} if kernel else {}),
        'parameters':{k:v for k,v in vars(a).items() if k!='func'},
        **({'chemistry':{'source_sha256':sha256(Path(local_chemistry_evidence.__code__.co_filename)),
            'property_table_sha256':sha256(Path(str(files('vdjtools.resources').joinpath('aa_property_table.txt')))),
            'sklearn_version':version('scikit-learn'),
            'definition':'equal-weight receptor mixture of ordered four-residue KF4 core windows; linked pairs use joint eight-dimensional vectors; query windows averaged',
            'core':'full junction[3:-3]','bandwidth_factor':a.chemistry_bandwidth_factor,
            'scale':'population standard deviation of published vdjtools KF4 amino-acid table; same geometry in reference/control',
            'exact_exclusion':'all full-junction identities across V categories, before both conditional and V-usage denominators',
            'interpretation':'experimental local density/V ratios; no P/E-value or binding probability'}} if getattr(a,'chemistry',False) else {}),
        'distance_scale':{'matrix':a.distance_matrix,'alphabet':a.alphabet,'matrix_unit':matrix_scale,'radius':a.radius,'radius_factor':a.radius_factor,
            'requested_radius':requested_radius,
            'radius_mode':'absolute' if absolute_radius else 'matrix_factor','bin_width':bin_width,
            'gap_cost':3*matrix_scale,'cdr3_weight':3,'v_loop_weight':a.v_loop_weight,'v_regions':a.v_regions,
            'stem_depth':a.stem_depth,'stem_weight':a.stem_weight,
            **({'stem':'joint V-J hairpin; level i pairs junction positions i and L-1-i on '
                       'untrimmed junctions; a level costs max(end,start) so it breaks when '
                       'either residue differs; mean over levels times stem_weight, in shared '
                       +a.distance_matrix+' cost units; one term per chain on paired distances'}
               if a.stem_depth else {}),
            'native_units_per_distance_unit':200,'native_radius':native_radius,
            'definition':'maximum effective '+a.distance_matrix+' substitution cost; shared CDR3/V-loop matrix' if a.alphabet=='full' else
                'zero within IMGT11 class; minimum clipped BLOSUM62 cost across classes; shared CDR3/V-loop matrix; raw unknown symbols remain singletons',
            'boundary_rounding':'floor to native integer units'},
        **({'reference_radius':{'selection':'reference-only nonexact leave-one-junction-out nearest-distance strongest confirmed Kneedle S=1 bend; no smoothing',
            'floor_factor':a.knee_floor_factor,'floor':floor,'floor_source':floor_source,
            'scope':'one species/locus/uniform geometry; compatible full-pMHC target unions',
            'fallback_radius':a.radius,'rounding':'ceil to matrix-quantum distance bins',
            'stage_seconds':radius_seconds,'table_sha256':sha256(Path(str(prefix)+'.radii.tsv')),
            'probability_interpretation':'none; radius selection does not transfer fixed-radius model calibration'}} if a.reference_radius else {}),
        'selection':('exhaustive nearest nonexcluded reference; n_neighbours counts only within radius'
                     if a.unbounded_nearest else 'radius ball'),
        'distance':(('clipped BLOSUM62 geometry default' if ctrim==2 else 'experimental symmetric-trim')+
                    (' 3*CDR3+CDR1+CDR2'+('+CDR2.5' if a.v_regions=='cdr12plus' else '')) if a.v_loop_weight else
                    'experimental CDR3-only:3*CDR3')+
                   (f'+joint stem depth {a.stem_depth}' if a.stem_depth else '')+'; '+
                   organism+' combo_xcr_2024-03-05',
        'junction_end_trim':{'n_terminal':3,'c_terminal':ctrim},
        **({'epitope_rules':{'germline':'joint V-gene predictive probabilities; one shared pooled pseudo-observation; no J/length marginal or added V-loop distance',
            'source_sha256':sha256(Path(reference_epitope_rules.__code__.co_filename)),
            'execution':('one conditional-V grouped traversal plus unsupported-category global fallback; scalar prior/motif reductions unchanged'
                if rule_batches is not None else 'single-task scalar traversal'),
            'sequence':'3/3-trimmed CDR3-only exhaustive multi-reference reciprocal rank, exact V-gene conditioning when both populations support the category; otherwise global CDR3-only rank and neutral V; fixed local kernels retained as diagnostics',
            'motifs':'exact3, IMGT11-class3, internal-wildcard4, IMGT11-class internal-wildcard4; unique per-receptor tokens; paired overlap fractions multiplied within the same reference pair',
            'combination':'V ratio counted once times CDR3 rank times (1+mean positive conditional motif frequency difference); bounded correlated sequence bonus, not independent likelihood evidence; unsupported V has no motif bonus',
            'rank_bins':'exact distance multiples of three, through 12*sum(maximum full junction length per axis); exhaustive population coverage checked',
            'prior':'p_g=(n_g+m_g+.5)/(N+M+.5K), K population joint V categories plus one; conditional kernel pseudo-observation c=(full_model_mass_ref+full_model_mass_bg+.5)/(n_g+m_g+1)',
            'exact_policy':'hit-only exclusion; full model N/M/n/m, pooled priors and motif enrichment frequencies stay fixed; excluded ranks contribute zero, not a censored tail',
            'source':'https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html',
            'interpretation':'empirical predictive ranking only; no P/E-value, biological posterior or benchmark-label fitting'}} if a.epitope_rules else {}),
        **({'shortlist_audit':{'stages':shortlist_costs,'scope':'CDR3 radius shortlist then unchanged nonnegative V costs; same eligibility/full-junction puncture. Original exhaustive outputs unchanged; no global K cap, biological confidence or paired claim.'}} if a.shortlist_audit else {}),
        **({'background_enrichment':{'background':background_provenance,'stage_wall_seconds':background_costs,
            'v_loop_weight':a.v_loop_weight,'v_regions':a.v_regions,
            'geometry':'fixed0.1VJ/1NDN' if a.position_weighting=='ndn' else 'full-profile decay' if a.position_weighting=='significance' else 'trim3 mean-one significance; junction floor before V/stem addition' if a.position_weighting=='significance-trimmed' else f'trim3/{ctrim}',
            **({'boundaries':boundary_provenance} if boundary_provenance else {}),
            'formula':'reference_weight/[reference_population*(control_weight+0.5)/(control_population+1)]',
            'kernel':'max(0,1-distance/radius)' if a.neighbour_weighting=='linear' else 'distance<=radius',
            'denominators':'reference distinct receptor/V-allele keys; controls use declared population units; both fixed under optional exact-hit filtering',
            'calibration':('none; posterior mean kernel mass under unit Dirichlet prior split equally at distance0/outside-radius; ranking only'
                           if a.neighbour_weighting=='linear' else 'none; Jeffreys background predictive smoothing, not a Bayes factor or P-value'),
            'nearest_output':'distance and native neighbour counts retain their unweighted comparator geometry'}} if a.background else {}),
        **({'local_rank_model':{
            'score_column':'rank_bayes_factor','distance_bins':'ceil(native distance / 200), fixed integer cutoffs 0..radius',
            'rank_weight':'(M+1)/(H_(M+1)*rank); ties averaged; outside-radius references have conditional null mean',
            'null':'exchangeable reference/control distance labels conditional on pooled distances',
            'alternative':'null label law tilted by the mean normalized reciprocal reference rank',
            'p_rank_bound':'min(1,1/rank_bayes_factor); conservative Markov bound, not exact tail',
            'p_nearest_rank':'hypergeometric first-target rank tail; upper-bin ties conservative; separate statistic',
            'posterior_rank_signal':'prior*BF/(1-prior+prior*BF); model probability, not biological Prob(TP)',
            'prior':a.rank_signal_prior,'sampling_status':'exchangeability is not established by source declaration; generated draws versus distinct references may violate it',
            'selection':'fixed radius, geometry and target; no selected-route or multiple-target calibration'}}
           if a.neighbour_weighting=='local-rank' else {}),
        **({'nearest_rank_model':{'score_column':'nearest_rank_score','formula':'-hypergeom.sf(0,N+M,N,m_nearest+1)',
            'distance':'exhaustive nearest full CDR3+V distance to the same linked reference; controls counted at that exact distance',
            'ties':'all control ties included; conservative first-reference tail',
            'sampling_status':'conditional exchangeable-label model only; no biological confidence, E-value or selected-label calibration',
            'selection':'no radius scan or fitted blend; fixed-radius fields unchanged'}} if a.nearest_rank else {}),
        **({'conditional_j_model':{'score_column':'j_conditioned_rank_score',
            'formula':'conditional_rank_bayes_factor * P_ref(J|V,length)/P_control(J|V,length)',
            'predictive':'(n_stratum_J + (n_J+.5)/(N+.5*K))/(n_stratum+1); K observed population J categories plus one unseen category',
            'pairing':'joint linked alpha/beta V,length and J tuples, no chain independence assumption',
            'missing':'missing query J neutral; missing/conflicting population J explicit categories',
            'exact_policy':'model-only frequencies and pooled priors stay fixed; exact-hit filtering belongs to the separate sequence statistic',
            'interpretation':'empirical ranking residual; no V or length marginal factor; no calibrated BF/P/E-value/posterior or biological confidence; J and retained junction may remain dependent'}} if a.conditional_j else {}),
        **({'conditional_rank_model':{'score_column':'conditional_rank_bayes_factor',
            'formula':'existing BF when n_radius>0; otherwise (M+1)/(N*H_(M+1))*[1/(k+1)+(N-1)*mean(1/r,r=k+1..M+1)]',
            'symbols':'N target population; M control population; k controls at exhaustive nearest target distance',
            'interpretation':'conditional all-reference rank mean with refined all-censored imputation; control ties conservative; no biological confidence or transferred P/E-value',
            'selection':'identical rule for every task; no fitted weights or task routing'}}
           if a.nearest_rank and a.neighbour_weighting=='local-rank' else {}),
        'calibration':'none','wall_seconds':time.perf_counter()-start}
    if rules_only:
        manifest.update(selection='fixed-model unbounded hit eligibility; nearest distances omitted',
            distance='rule-only clipped 3*CDR3, trim3/3; V usage counted once, no V-loop distance',
            junction_end_trim={'n_terminal':3,'c_terminal':3})
        manifest['background_enrichment'].update(v_loop_weight=0,v_regions='unused',
            geometry='rule-only clipped CDR3, trim3/3',
            formula='rule_* predictive ratios and conditional multi-reference ranks',
            kernel='spatial and motif kernels declared in epitope_rules; ordinary comparator omitted',
            calibration='none; empirical reference ranking only',nearest_output='omitted')
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return 0
