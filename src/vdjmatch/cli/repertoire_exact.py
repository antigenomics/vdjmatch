"""Exact reference-panel enrichment preserving source-row and read counts."""
from pathlib import Path
from ..config import GRADED_FAMILIES


def register(sub):
    p=sub.add_parser('repertoire-exact',help=__doc__)
    p.add_argument('metadata',help='TSV sample_id,path and optional FDR-group columns')
    p.add_argument('--vdjdb',required=True)
    p.add_argument('--control',required=True,help='Distinct junction_aa,v,j probe keys; optional joint_pgen column')
    p.add_argument('--abundance-metrics',action='store_true',help='Compare full-depth curve-average richness, Shannon and Simpson Hill diversity against reference mocks (includes Chao1)')
    p.add_argument('--chao1',action='store_true',help='Compare matched source-clonotype Chao1 from vdjtools against the same mock panels')
    p.add_argument('--pgen-features',help='Verified marginal Pgen features; use raw probabilities and bins without corrections')
    p.add_argument('--model-vj-counts',help='Generated-model v,j,records counts from the selected source')
    p.add_argument('--model-vj-probabilities',help='Exact model v,j,p probabilities; alternative to generated V/J counts')
    p.add_argument('--include-nonfunctional',action='store_true',help='Retain nonfunctional/unsupported genes and zero-mass probes in the lowest Pgen bin; modern default filters them')
    p.add_argument('--model-informed-prior',action='store_true',help='Use one total model-distributed prior observation per V group; default preserves historical add-one')
    p.add_argument('--epitope',action='append')
    p.add_argument('--pool-mhc',action='store_true',help='Explicit original epitope-union reproduction')
    p.add_argument('--locus',choices=['TRA','TRB'],default='TRB')
    p.add_argument('--species',choices=['human','mouse'],default='human')
    p.add_argument('--selection-q',type=float,default=None,help='Q override; default is calibrated for source/species/locus')
    p.add_argument('--selection-profile',help='Opt-in length | vj-length | vdjtools calibration JSON')
    p.add_argument('--source',default=None,help='Model source; default isalgo (legacy human: olga-legacy)')
    p.add_argument('--legacy',action='store_true',help='Compiled reproduction with historical source and Q=1')
    p.add_argument('--pgen-mismatches',type=int,choices=[0,1],default=0,help='0 exact joint Pgen; 1 inclusive Hamming-1 joint mass')
    p.add_argument('--allow-missing-model-genes',action='store_true',help='Retain absent model genes with zero joint mass and the existing add-one V/J-count smoothing; original calls are unchanged')
    p.add_argument('--matcher',choices=['exact','hamming','stem4','blosum62','positional-blosum62',*GRADED_FAMILIES],default='exact')
    p.add_argument('--min-score',type=float,help='Explicit dimensionless stem4 score threshold; fixed geometry')
    p.add_argument('--radius-factor',type=float,help='Normalized distance bound for BLOSUM or graded family matching')
    p.add_argument('--mismatches',type=int,default=1,help='Hamming substitutions or graded total edits; bulk retains exact J')
    p.add_argument('--sampler',choices=['numpy','r'],default='numpy')
    p.add_argument('--mock-replicates',type=int,default=200)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--fdr-group',action='append',default=[],help='Metadata columns defining FDR families within each target')
    p.add_argument('--threads',type=int,default=0,help='Native CPU budget; 0 selects the allocation')
    p.add_argument('--output-prefix',required=True)
    p.set_defaults(func=main)



def _selection_options(args, *, external_pgen=False):
    """Resolve model and Q together; external probability scales need explicit Q."""
    from vdjtools.model.selection import selection_factor, selection_mass
    source = args.source or ('olga-legacy' if args.legacy and args.species == 'human' else 'isalgo')
    if args.selection_q is not None:
        return source, selection_mass(1., args.selection_q), 'explicit'
    if args.legacy:
        return source, 1., 'legacy'
    if external_pgen:
        raise ValueError('external joint_pgen requires explicit --selection-q or --legacy')
    return source, selection_factor(args.locus, organism=args.species, source=source), 'calibrated'


def _source(path, *, locus, functional_genes=None, repair_v=False, species="human"):
    import polars as pl
    from .. import io
    frame=io.read_rearrangement(path,valid_aa=True).filter(pl.col("locus")==locus)
    frame=frame.select(pl.col('cdr3').alias('junction_aa'),
        pl.col('v','j').str.split('*').list.first(),'count')
    if repair_v:
        from ..io.columns import repair_unknown_v
        frame=repair_unknown_v(frame,locus=locus,species=species)
    if functional_genes is not None:
        frame=frame.filter(pl.all_horizontal(pl.col(k).is_in(v) for k,v in functional_genes.items()))
    return frame


def _functional_genes(model):
    """Gene-level functional, positive-choice support for the fixed model."""
    import polars as pl
    allowed={}
    for kind in ['v','j']:
        allele=kind+'_allele'
        positive=model.tables[kind+'_choice'].filter(pl.col('p')>0).select(allele).unique()
        allowed[kind]=model.genomic['genes_'+kind].filter(pl.col('functional')).join(
            positive,on=allele,how='semi')['gene'].unique().to_list()
    return allowed


def _joint_pgen(probes,*,species,locus,threads,source="isalgo",profile=None,mismatches=0,allow_missing_model_genes=False):
    import polars as pl
    from vdjtools.model import load_bundled,native
    model=load_bundled(locus,source,organism=species,collapse=False)
    if allow_missing_model_genes:
        import numpy as np
        if profile is not None:raise ValueError('missing model genes require model counts, not a selection profile')
        _,v_index,j_index=native.pack(model)
        known=set(v_index)|set(j_index)
        known.update(a.split('*')[0] for a in tuple(known))
        supported=probes.select(pl.all_horizontal(pl.col('v','j').is_in(known))).to_series()
        valid=probes.filter(supported)
        values=np.zeros(probes.height)
        if valid.height:
            values[supported.to_numpy()]=native.pgen_aa_gene_batch(model,valid['junction_aa'].to_list(),
                valid['v'].to_list(),valid['j'].to_list(),mismatches=mismatches,threads=threads)
        return probes.with_columns(pl.Series('pgen',values),supported.alias('pgen_model_supported'))
    values=native.pgen_aa_gene_batch(model,probes['junction_aa'].to_list(),probes['v'].to_list(),
                                   probes['j'].to_list(),mismatches=mismatches,threads=threads)
    if profile is not None:
        from vdjtools.model.selection import selection_profile_factors
        values=selection_profile_factors(profile,probes['junction_aa'].str.len_chars().to_numpy(),
            probes['v'].to_list(),probes['j'].to_list()) * values
    return probes.with_columns(pl.Series('pgen',values))


def _summaries(counts):
    import numpy as np
    import polars as pl
    rows=[]
    for (sample,target),part in counts.partition_by('sample_id','group_id',as_dict=True).items():
        observed=part.filter(pl.col('replicate')==-1).row(0,named=True)
        mock=part.filter(pl.col('replicate')>=0)
        present=mock.filter(pl.col('matched_records')>0)
        row=dict(sample_id=sample,group_id=target,matched_records=observed['matched_records'],matched_reads=observed['matched_reads'],mock_replicates=mock.height,original_present_mocks=present.height)
        for endpoint in ['records','reads']+[e for e in ['chao1','curve_mean','shannon','simpson','hill0','hill1','hill2'] if 'matched_'+e in counts.columns]:
            direction=-1 if endpoint in {'shannon','simpson','hill0','hill1','hill2'} else 1
            actual=observed['matched_'+endpoint];values=mock['matched_'+endpoint].to_numpy()
            original=present['matched_'+endpoint].to_numpy()
            transformed=np.log2(original+1.) if endpoint=='reads' else original
            effect_actual=np.log2(actual+1.) if endpoint=='reads' else actual
            mean=float(transformed.mean()) if len(original) else None
            sd=float(transformed.std(ddof=1)) if len(original)>1 else None
            with np.errstate(divide='ignore',invalid='ignore'):
                effect=float(np.divide(direction*(effect_actual-mean),sd)) if sd is not None else None
            if effect is not None and np.isnan(effect):effect=None
            row.update({'matched_'+endpoint:actual,'mock_mean_'+endpoint:float(values.mean()),'p_'+endpoint:float((1+(direction*values>=direction*actual).sum())/(len(values)+1)),
                'original_p_'+endpoint:float((1+(direction*original>=direction*actual).sum())/(len(original)+1)) if observed['matched_records']>0 else None,
                'original_mean_'+endpoint:mean,'original_sd_'+endpoint:sd,
                'original_effect_'+endpoint:effect,
                'excess_'+endpoint:float(direction*(actual-values.mean()))})
        if 'sample_clonotypes' in part.columns:
            for denominator in ['sample_clonotypes','epitope_clonotypes','sample_reads']:
                if part[denominator].n_unique()!=1:raise ValueError('real/mock normalization denominators must agree')
                row[denominator]=observed[denominator]
            for measure in ['normalized_clonotype_matches','matched_read_fraction']:
                row[measure]=observed[measure]
                row['mock_mean_'+measure]=float(mock[measure].mean())
        rows.append(row)
    return pl.DataFrame(rows,infer_schema_length=None).with_columns(pl.col('^original_(p|mean|sd|effect)_.+$').cast(pl.Float64))


def _adjust_summary(summary,endpoints,fdr_groups):
    import polars as pl
    from vdjtools.biomarker.stats import fdr_bh
    parts=[]
    for part in summary.partition_by(['group_id',*fdr_groups]):
        for endpoint in endpoints:
            for stem in ['', 'original_']:
                col=stem+'p_'+endpoint
                values=part[col].to_list()
                adjusted=iter(fdr_bh(part[col].drop_nulls().to_numpy()))
                part=part.with_columns(pl.Series(stem+'q_'+endpoint,[next(adjusted) if value is not None else None for value in values],dtype=pl.Float64))
            part=part.with_columns(((pl.col('matched_records')>0)&(pl.col('p_'+endpoint)<.05)&(pl.col('excess_'+endpoint)>0)).alias('significant_nominal_'+endpoint),
                ((pl.col('matched_records')>0)&(pl.col('q_'+endpoint)<.1)&(pl.col('excess_'+endpoint)>0)).alias('significant_bh_'+endpoint),
                ((pl.col('original_q_'+endpoint)<.1)&(pl.col('original_effect_'+endpoint)>0)).alias('original_significant_bh_'+endpoint))
        parts.append(part)
    return pl.concat(parts)


def main(a):
    import json,time
    from importlib.metadata import version
    import polars as pl
    from .. import db
    from ..db.cache import sha256
    from ..io.airr import _read_table
    from ..aggregate.candidates import PMHC
    from ..aggregate.repertoire_enrichment import exact_reference_counts,kernel_reference_counts
    from ..evalue.mock_reference import compiled_reference_strata,sample_stratified_panels,_pgen_model_identity
    from ..runner.annotation import _validate_output_paths
    from vdjtools.cores import kernel_threads
    a.threads = kernel_threads(a.threads)
    abundance_metrics=getattr(a,'abundance_metrics',False)
    chao1=getattr(a,'chao1',False) or abundance_metrics
    if chao1 and a.matcher=='stem4':raise ValueError('--chao1 requires a membership-count matcher')
    if a.mock_replicates<2 or a.seed<0 or a.threads<1:raise ValueError('at least two mocks, nonnegative seed and positive threads required')
    if not 0<=a.mismatches<=3:raise ValueError('zero to three mismatches required')
    if a.matcher!='stem4' and a.min_score is not None:raise ValueError('--min-score requires --matcher stem4')
    if a.matcher=='stem4' and (a.min_score is None or not 0<a.min_score<=1):
        raise ValueError('stem4 requires an explicit dimensionless threshold in (0,1]')
    factor=getattr(a,'radius_factor',None)
    if a.matcher in {'blosum62','positional-blosum62'}:
        from ..match.scope import blosum_radius_params
        blosum_radius_params(a.matcher,factor,1)
    elif a.matcher in GRADED_FAMILIES:
        from ..match.scope import graded_family_params
        graded_family_params(a.matcher,a.mismatches,factor)
    elif factor is not None:raise ValueError('--radius-factor requires a BLOSUM or graded family matcher')
    start=time.monotonic();meta_path=Path(a.metadata).expanduser().resolve();meta=_read_table(meta_path)
    stages={};sample_timings=[];mark=start
    needed=['sample_id','path',*a.fdr_group]
    if not set(needed)<=set(meta.columns) or not meta.height or meta['sample_id'].n_unique()!=meta.height:
        raise ValueError('unique sample_id,path and declared FDR columns required')
    if any(meta[c].null_count() or meta[c].str.strip_chars().eq('').any() for c in needed):raise ValueError('sample identities, paths and FDR groups cannot be missing')
    paths=[(meta_path.parent/Path(p).expanduser()).resolve() for p in meta['path']]
    if len(set(paths))!=len(paths):raise ValueError('sample paths must be unique')
    marginal=getattr(a,'pgen_features',None)
    if marginal and any([a.legacy,a.model_vj_counts,a.model_vj_probabilities,a.selection_profile,a.model_informed_prior,a.selection_q is not None,a.allow_missing_model_genes]):raise ValueError('--pgen-features uses raw marginal Pgen without model/selection corrections')
    model_probabilities=getattr(a,'model_vj_probabilities',None)
    if model_probabilities and a.model_vj_counts:
        raise ValueError('supply model V/J counts or probabilities, not both')
    if a.selection_profile and (a.legacy or a.model_vj_counts or model_probabilities):
        raise ValueError('selection profiles require modern mode and replace model V/J counts')
    if a.allow_missing_model_genes and a.selection_profile:
        raise ValueError('--allow-missing-model-genes requires model counts, not a selection profile')
    if not marginal and not a.selection_profile and not a.model_vj_counts and not model_probabilities:
        raise ValueError('--model-vj-counts is required without a selection profile')
    inputs=[meta_path,Path(a.vdjdb).expanduser().resolve(),Path(a.control).expanduser().resolve(),*paths]
    if a.model_vj_counts:
        inputs.append(Path(a.model_vj_counts).expanduser().resolve())
    if model_probabilities:
        inputs.append(Path(model_probabilities).expanduser().resolve())
    if a.selection_profile and a.selection_profile not in {'length','vj-length'}:
        inputs.append(Path(a.selection_profile).expanduser().resolve())
    if marginal:
        feature_path=Path(marginal).expanduser().resolve();inputs.extend([feature_path,feature_path.with_name(feature_path.name.removesuffix('.features.parquet')+'.manifest.json')])
    prefix=_validate_output_paths(a.output_prefix,[str(p) for p in inputs])
    suffixes=['v-repairs.tsv','source-audit.tsv','summary.tsv','counts.tsv','vj-audit.tsv','reference-bins.tsv','manifest.json']
    if a.allow_missing_model_genes or a.include_nonfunctional:suffixes.append('missing-model-genes.tsv')
    if chao1:suffixes.append('abundance-spectra.tsv')
    if abundance_metrics:suffixes.extend(['discovery-curves.tsv','hill-rarefaction.tsv'])
    if any(Path(str(prefix)+'.'+s).exists() for s in suffixes):raise FileExistsError('fresh output prefix required')
    before={str(p):sha256(p) for p in inputs}
    import seqtree._core as search_native
    import vdjtools._core as native_core
    import vdjtools.model.io as model_io
    import vdjtools.model.native as model_native
    import vdjtools.model.selection as model_selection
    source_files=[*Path(__file__).resolve().parents[1].rglob('*.py'),Path(native_core.__file__),Path(search_native.__file__),Path(model_io.__file__),Path(model_native.__file__),Path(model_selection.__file__),Path(model_selection.__file__).with_name('selection_factors.json'),Path(model_selection.__file__).with_name('selection_profiles.json')]
    if chao1:
        from vdjtools.stats.inext import asymptotic_diversity
        source_files.append(Path(asymptotic_diversity.__code__.co_filename))
    if abundance_metrics:
        from vdjtools.stats import subset_rarefaction_batch,shannon_wiener,inverse_simpson,inext_batch
        source_files.extend(Path(f.__code__.co_filename) for f in [subset_rarefaction_batch,shannon_wiener,inverse_simpson,inext_batch])
    source_before={str(p):sha256(p) for p in source_files}
    now=time.monotonic();stages['input_validation_and_freeze']=now-mark;mark=now
    reference=db.load(inputs[1],species={'human':'HomoSapiens','mouse':'MusMusculus'}[a.species]).filter(pl.col('gene')==a.locus)
    full_reference=reference.clone()
    if a.epitope:
        if set(a.epitope)-set(reference['epitope']):raise ValueError('requested epitope absent from selected reference population')
        reference=reference.filter(pl.col('epitope').is_in(a.epitope))
    if not reference.height:raise ValueError('empty selected reference population')
    labels=['epitope'] if a.pool_mhc else PMHC
    inventory=reference.select(labels).unique().sort(labels).with_columns(pl.struct(labels).struct.json_encode().alias('group_id'))
    keys=['junction_aa','v','j']
    reference=reference.join(inventory,on=labels,nulls_equal=True).select('group_id',pl.col('cdr3').alias('junction_aa'),pl.col('v','j').str.split('*').list.first()).unique()
    if marginal:
        from .repertoire_reference import _read_control
        probes=_read_control(inputs[2],species=a.species,locus=a.locus)
    else:probes=_read_table(inputs[2])
    now=time.monotonic();stages['reference_and_probe_loading']=now-mark;mark=now
    model_identity=None
    shape=None
    if marginal:
        a.source=a.source or 'isalgo';a.selection_q=1.;q_origin='unused'
    elif a.selection_profile:
        if 'joint_pgen' in probes.columns:
            raise ValueError('selection profiles require fresh native probabilities, not external joint_pgen')
        from vdjtools.model import load_bundled
        from vdjtools.model.selection import load_selection_profile, selection_mass
        a.source=a.source or 'isalgo'
        shape,profile_q=load_selection_profile(a.selection_profile,
            load_bundled(a.locus,a.source,organism=a.species,collapse=False))
        a.selection_q=profile_q if a.selection_q is None else selection_mass(1.,a.selection_q)
        q_origin='profile'
    else:
        a.source,a.selection_q,q_origin=_selection_options(a,external_pgen='joint_pgen' in probes.columns)
    if not set(keys)<=set(probes.columns):raise ValueError('control requires junction_aa,v,j')
    from vdjtools.model import load_bundled
    from ..io.columns import repair_unknown_v
    functional_genes=None
    filtering=not (a.legacy or a.include_nonfunctional)
    if filtering:
        working_model=load_bundled(a.locus,a.source,organism=a.species,collapse=False)
        functional_genes=_functional_genes(working_model)
    if not a.legacy:
        reference=repair_unknown_v(reference,locus=a.locus,species=a.species).select('group_id',*keys).unique()
        probes=repair_unknown_v(probes,locus=a.locus,species=a.species)
    filter_audit={'enabled':filtering,'reference_keys_before':reference.height,'control_keys_before':probes.height}
    if filtering:
        predicate=pl.all_horizontal(pl.col(k).is_in(v) for k,v in functional_genes.items())
        reference=reference.filter(predicate)
        probes=probes.filter(predicate)
    filter_audit.update(reference_keys_after_gene_filter=reference.height,control_keys_after_gene_filter=probes.height)
    if not reference.height or not probes.height:raise ValueError('no functional reference/control keys remain')

    if 'event_weight' in probes.columns:
        from ..io.columns import _integer_count
        control_weights=probes.select(*keys,_integer_count('event_weight')).group_by(keys).agg(pl.col('event_weight').sum())
    else:control_weights=probes.group_by(keys).len(name='event_weight')
    if a.allow_missing_model_genes and 'joint_pgen' in probes.columns:
        raise ValueError('--allow-missing-model-genes requires fresh native probabilities')
    if 'joint_pgen' in probes.columns and a.pgen_mismatches:
        raise ValueError('neighbourhood Pgen requires fresh native probabilities')
    if marginal:
        from ..evalue.mock_reference import load_pgen_features
        model_identity=_pgen_model_identity(a.species,a.locus,source=a.source)
        probabilities,feature_identity=load_pgen_features(marginal,reference,probes,species=a.species,locus=a.locus,mismatches=a.pgen_mismatches,source=a.source)
        probes=pl.concat([probes.select(keys),reference.select(keys)]).unique(maintain_order=True).join(probabilities,on='junction_aa',how='left',validate='m:1')
        pgen_origin='verified raw marginal Pgen features'
    elif 'joint_pgen' in probes.columns:
        probes=probes.select(*keys,pl.col('joint_pgen').cast(pl.Float64).alias('pgen'))
        if reference.select(keys).join(probes.select(keys),on=keys,how='anti').height:raise ValueError('external joint_pgen probes must include every reference key')
        pgen_origin='external gene-masked probability asset'
    else:
        model_identity=_pgen_model_identity(a.species,a.locus,source=a.source)
        probes=pl.concat([probes.select(keys),reference.select(keys)]).unique(maintain_order=True)
        probes=_joint_pgen(probes,species=a.species,locus=a.locus,threads=a.threads,source=a.source,profile=shape,mismatches=a.pgen_mismatches,allow_missing_model_genes=a.allow_missing_model_genes or (a.include_nonfunctional and shape is None))
        pgen_origin='native gene-allele-summed joint probability'
    now=time.monotonic();stages['pgen_and_selection']=now-mark;mark=now
    if filtering:
        probes=probes.filter(pl.col('pgen')>0)
        reference=reference.join(probes.select(keys),on=keys,how='semi')
    filter_audit.update(reference_keys_after_mass_filter=reference.height,control_keys_after_mass_filter=probes.height)
    if not reference.height or not probes.height:raise ValueError('no positive-mass reference/control keys remain')
    vj=[];source_audit=[];v_repairs=[]
    for sample,path in zip(meta['sample_id'],paths):
        frame=_source(path,locus=a.locus,repair_v=not a.legacy,species=a.species)
        records_before=frame.height;reads_before=int(frame['count'].sum())
        repaired=int(frame['v_repaired'].sum()) if 'v_repaired' in frame.columns else 0
        if repaired:
            v_repairs.append(frame.filter(pl.col('v_repaired')).select(*keys,'original_v','count').with_columns(pl.lit(sample).alias('sample_id')))
        if filtering:frame=frame.filter(predicate)
        source_audit.append(dict(sample_id=sample,source_records=records_before,retained_records=frame.height,
            excluded_records=records_before-frame.height,source_reads=reads_before,retained_reads=int(frame['count'].sum()),
            excluded_reads=reads_before-int(frame['count'].sum()),repaired_v_records=repaired))
        vj.append(frame.group_by('v','j').len(name='records'))
    now=time.monotonic();stages['source_vj_loading']=now-mark;mark=now
    if marginal:
        from ..evalue.mock_reference import pgen_strata
        probes=probes.join(pgen_strata(probes.select('junction_aa','pgen').unique()),on='junction_aa',validate='m:1').with_columns((pl.col('pgen')<=0).alias('zero_model_mass'))
        audit=pl.concat(vj).group_by('v','j').agg(pl.col('records').sum())
    elif model_probabilities:
        probes,audit=compiled_reference_strata(probes,pl.concat(vj),selection_q=a.selection_q,
            model_vj_probabilities=_read_table(Path(model_probabilities).expanduser().resolve()),model_informed_prior=a.model_informed_prior)
    elif shape is None:
        model=_read_table(Path(a.model_vj_counts).expanduser().resolve())
        probes,audit=compiled_reference_strata(probes,pl.concat(vj),model,selection_q=a.selection_q,model_informed_prior=a.model_informed_prior)
    else:
        selected=pl.DataFrame(shape['selected_vj']).rename({'v_call':'v','j_call':'j'})
        probes,audit=compiled_reference_strata(probes,pl.concat(vj),selection_q=a.selection_q,
            model_vj_probabilities=selected,model_informed_prior=a.model_informed_prior)
    if not a.legacy:
        supported=probes.join(control_weights,on=keys,how='inner').filter((pl.col('event_weight')>0)&~pl.col('zero_model_mass'))
        lowest=int(supported['pgen_bin'].min()) if supported.height else int(probes['pgen_bin'].min())
        probes=probes.with_columns(pl.when(pl.col('zero_model_mass')).then(lowest).otherwise(pl.col('pgen_bin')).alias('pgen_bin'))
    reference=reference.join(probes.select(*keys,'pgen_bin'),on=keys,validate='m:1')
    now=time.monotonic();stages['strata_and_reference_bins']=now-mark;mark=now
    control=probes.select(*keys,'pgen_bin').join(control_weights,on=keys,validate='1:1').filter(pl.col('event_weight')>0)
    mocks=sample_stratified_panels(reference,control,replicates=a.mock_replicates,seed=a.seed,sampler=a.sampler)
    panels=pl.concat([reference.select('group_id',*keys).with_columns(pl.lit(-1,dtype=pl.Int64).alias('replicate')),mocks.select('group_id',*keys,pl.col('replicate').cast(pl.Int64))])
    panel_sizes=panels.group_by('group_id','replicate').len(name='epitope_clonotypes')
    now=time.monotonic();stages['mock_panel_generation']=now-mark;mark=now
    counts=[];cohort=[];totals=[]
    for i,(sample,path) in enumerate(zip(meta['sample_id'],paths),1):
        sample_start=time.monotonic()
        observations=_source(path,locus=a.locus,functional_genes=functional_genes,repair_v=not a.legacy,species=a.species)
        loaded=time.monotonic()
        totals.append(dict(sample_id=sample,sample_clonotypes=observations.height,sample_reads=int(observations['count'].sum())))
        if a.matcher=='stem4':cohort.append(observations.with_columns(pl.lit(sample).alias('sample_id')))
        else:counts.append(exact_reference_counts(observations,panels,matcher=a.matcher,mismatches=a.mismatches,threads=a.threads,radius_factor=factor,locus=a.locus,species=a.species,with_abundances=chao1).with_columns(pl.lit(sample).alias('sample_id')))
        finished=time.monotonic()
        sample_timings.append(dict(sample_id=sample,records=observations.height,load_seconds=loaded-sample_start,match_or_stage_seconds=finished-loaded))
        print(f'loaded sample {i}/{len(paths)}; load={loaded-sample_start:.3f}s match_or_stage={finished-loaded:.3f}s elapsed={finished-start:.3f}s',flush=True)
    now=time.monotonic();stages['sample_load_and_matching']=now-mark;mark=now
    if a.matcher=='stem4':
        score_inventory=inventory
        if a.pool_mhc:
            classes=full_reference.filter(pl.col('epitope').is_in(inventory['epitope'])).select('epitope','mhc_class').unique()
            if classes['epitope'].n_unique()!=classes.height:raise ValueError('pooled kernel targets require one MHC class per peptide')
            score_inventory=inventory.join(classes,on='epitope',validate='1:1').with_columns(pl.lit(None,dtype=pl.String).alias('mhc_a'),pl.lit(None,dtype=pl.String).alias('mhc_b'))
        counts=kernel_reference_counts(pl.concat(cohort),panels,score_inventory,full_reference.rename({'gene':'locus'}),locus=a.locus,species=a.species,min_score=a.min_score,threads=a.threads,progress=lambda replicate:print(f'scored reference panel {replicate}',flush=True))
    else:counts=pl.concat(counts)
    if chao1:
        vectors=counts['matched_abundances'].to_list()
        counts=counts.with_columns(pl.Series('matched_chao1',[
            float(asymptotic_diversity(v,q=0,se=False)['estimator'][0]) if any(v) else 0. for v in vectors]))
        spectra=counts.select('sample_id','group_id','replicate','matched_abundances').explode('matched_abundances',empty_as_null=True).drop_nulls('matched_abundances').rename({'matched_abundances':'reads_per_source_row'}).group_by('sample_id','group_id','replicate','reads_per_source_row').len(name='source_rows').sort('sample_id','group_id','replicate','reads_per_source_row')
        counts=counts.with_columns(pl.col('matched_abundances').list.count_matches(1).alias('matched_singletons'),pl.col('matched_abundances').list.count_matches(2).alias('matched_doubletons')).drop('matched_abundances')
    counts=counts.join(panel_sizes,on=['group_id','replicate'],validate='m:1',maintain_order='left').join(pl.DataFrame(totals),on='sample_id',validate='m:1',maintain_order='left').with_columns(
        (pl.col('matched_records')/pl.col('sample_clonotypes')/pl.col('epitope_clonotypes')).alias('normalized_clonotype_matches'),
        (pl.col('matched_reads')/pl.col('sample_reads')).alias('matched_read_fraction'))
    if abundance_metrics:
        import numpy as np
        from scipy.integrate import trapezoid
        depth=int(counts['sample_reads'].min())
        grid=np.unique(np.linspace(0,depth,513).astype(np.int64))
        curves=subset_rarefaction_batch(vectors,counts['sample_reads'].to_numpy(),sizes=grid,threads=a.threads)
        matrix=curves['qD'].to_numpy().reshape(len(vectors),len(grid))
        means=trapezoid(matrix,x=grid,axis=1)/depth
        coarse_indices=np.unique(np.append(np.arange(0,len(grid),2),len(grid)-1))
        coarse=trapezoid(matrix[:,coarse_indices],x=grid[coarse_indices],axis=1)/depth
        counts=counts.with_columns(pl.Series('matched_curve_mean',means),
            pl.Series('curve_mean_grid_delta',means-coarse),
            pl.Series('matched_shannon',[shannon_wiener(np.asarray(v)) if any(v) else 0. for v in vectors]),
            pl.Series('matched_simpson',[inverse_simpson(np.asarray(v)) if any(v) else 0. for v in vectors]))
        curves=curves.join(counts.with_row_index('sample').select('sample','sample_id','group_id','replicate'),on='sample',validate='m:1')
        curve_summary=curves.group_by('sample_id','group_id','m').agg(
            pl.col('qD').filter(pl.col('replicate')==-1).first().alias('observed_richness'),
            pl.col('qD').filter(pl.col('replicate')>=0).mean().alias('mock_mean_richness'),
            pl.col('qD').filter(pl.col('replicate')>=0).std().alias('mock_sd_richness')).sort('sample_id','group_id','m')
        del curves,matrix
        limits=counts.group_by('sample_id','group_id').agg(pl.col('matched_reads').min().alias('hill_max_depth'))
        max_hill_depth=int(limits['hill_max_depth'].max())
        hill_keys=counts.with_row_index('sample').select('sample','sample_id','group_id','replicate','matched_records','matched_reads').join(limits,on=['sample_id','group_id'],validate='m:1',maintain_order='left').filter(pl.col('hill_max_depth')>0)
        if hill_keys.height:
            hill=inext_batch([vectors[i] for i in hill_keys['sample']],q=(0,1,2),sizes=range(1,max_hill_depth+1),se=False,threads=a.threads)
            hill=hill.join(hill_keys.drop('sample').with_row_index('sample'),on='sample',validate='m:1').filter(pl.col('m')<=pl.col('hill_max_depth'))
        else:hill=pl.DataFrame(schema={'m':pl.Int64})
        hill_parts=[]
        for (m,),part in hill.partition_by('m',as_dict=True).items():
            wide=part.pivot(on='order_q',index=['sample_id','group_id','replicate','matched_records','matched_reads'],values='qD').rename({'0':'matched_hill0','1':'matched_hill1','2':'matched_hill2'})
            h=_summaries(wide).join(meta.select('sample_id',*a.fdr_group),on='sample_id',validate='m:1')
            h=_adjust_summary(h,['hill0','hill1','hill2'],a.fdr_group).with_columns(pl.lit(m).alias('matched_read_depth'),(1-1/pl.col('matched_hill2')).alias('gini_simpson'))
            hill_parts.append(h)
        hill_summary=(pl.concat(hill_parts).join(limits,on=['sample_id','group_id'],validate='m:1') if hill_parts else limits.clear().with_columns(pl.lit(None,dtype=pl.Int64).alias('matched_read_depth')))
        unavailable=limits.filter(pl.col('hill_max_depth')==0).with_columns(pl.lit(None,dtype=pl.Int64).alias('matched_read_depth'))
        hill_summary=pl.concat([hill_summary,unavailable],how='diagonal_relaxed').sort('sample_id','group_id','matched_read_depth')
        counts=counts.with_columns(pl.when(pl.col('matched_simpson')>0).then(1-1/pl.col('matched_simpson')).otherwise(0.).alias('gini_simpson'))
    now=time.monotonic();stages['cohort_kernel_and_count_assembly']=now-mark;mark=now
    summary=_summaries(counts).join(meta.select('sample_id',*a.fdr_group),on='sample_id',validate='m:1')
    summary=_adjust_summary(summary,['records','reads']+(['chao1'] if chao1 else [])+(['curve_mean','shannon','simpson'] if abundance_metrics else []),a.fdr_group).join(inventory,on='group_id',validate='m:1')
    now=time.monotonic();stages['statistics_and_fdr']=now-mark;mark=now
    outputs={'source-audit.tsv':pl.DataFrame(source_audit),'summary.tsv':summary,'counts.tsv':counts,'vj-audit.tsv':audit,'reference-bins.tsv':reference.group_by('group_id','pgen_bin').len(name='reference_keys')}
    if chao1:outputs['abundance-spectra.tsv']=spectra
    if abundance_metrics:outputs.update({'discovery-curves.tsv':curve_summary,'hill-rarefaction.tsv':hill_summary})
    if v_repairs:outputs['v-repairs.tsv']=pl.concat(v_repairs)
    if 'pgen_model_supported' in probes.columns:
        outputs['missing-model-genes.tsv']=probes.filter(~pl.col('pgen_model_supported'))
    prefix.parent.mkdir(parents=True,exist_ok=True)
    for suffix,frame in outputs.items():frame.write_csv(str(prefix)+'.'+suffix,separator='\t')
    now=time.monotonic();stages['table_export']=now-mark;mark=now
    if before!={str(p):sha256(p) for p in inputs} or source_before!={str(p):sha256(p) for p in source_files}:raise RuntimeError('inputs or source changed during calculation')
    if model_identity is not None and model_identity!=_pgen_model_identity(a.species,a.locus,source=a.source):raise RuntimeError('recombination model changed during calculation')
    stages['final_identity_checks']=time.monotonic()-mark
    manifest=dict(stage_seconds=stages,sample_timings=sample_timings,timing_scope='CLI body after argument validation; external process-tree profile includes startup',inputs_frozen=True,source_frozen=True,inputs=before,source=source_before,options=vars(a)|{'func':None},pgen_origin=pgen_origin,probability='selection_q times AA-max joint Pgen times within-V add-one J usage ratio; floor2^-80',counting='source rows and reads; '+('fixed stem4 score; production V geometry imputation; no hard V/J predicate' if a.matcher=='stem4' else 'exact V/J; '+a.matcher+' junction membership'),fdr_family=['group_id',*a.fdr_group],seconds=time.monotonic()-start,software={n:version(n) for n in ['vdjmatch','vdjtools','seqtree','polars','numpy']},outputs={Path(str(prefix)+'.'+s).name:sha256(Path(str(prefix)+'.'+s)) for s in outputs})
    if a.sampler=='r':
        import subprocess,shutil
        manifest['r_version']=subprocess.run([shutil.which('Rscript'),'--version'],capture_output=True,text=True,check=True).stdout.strip()
    if a.model_informed_prior:
        manifest['probability']='selection_q times AA-max joint Pgen times within-V model-informed J usage ratio; one total prior observation per supported V; missing/zero mass goes to lowest control bin'
    if marginal:manifest.update(probability='raw marginal junction Pgen; no corrections',prepared_pgen=feature_identity)
    if chao1:manifest['chao1']=dict(api='vdjtools.stats.inext.asymptotic_diversity(q=0,se=False)',unit='matched source-row clonotypes; read abundances',empty_matches=0,internal_bootstrap=False,mock_replicates=a.mock_replicates)
    if abundance_metrics:manifest['abundance_metrics']=dict(curve_mean=dict(api='vdjtools.stats.subset_rarefaction_batch',maximum_depth=depth,grid=grid.tolist(),integration='scipy.integrate.trapezoid / maximum_depth',coarse_grid_max_absolute_delta=float(np.max(np.abs(means-coarse)))),shannon='vdjtools.stats.shannon_wiener: exp(H), lower-tail concentration',simpson='vdjtools.stats.inverse_simpson: lower-tail concentration',hill_rarefaction='vdjtools.stats.inext_batch(q=(0,1,2),se=False); integer matched-read depths 1 through each library minimum real/mock matched mass; lower tails; no extrapolated values used',empty_matches=0,unit='original source-row clonotypes with producer read abundances',tails=dict(chao1='upper',curve_mean='upper',shannon='lower',simpson='lower',hill0='lower',hill1='lower',hill2='lower'))
    manifest['functional_filter']=filter_audit
    manifest['model_identity']=model_identity
    manifest['selection_q_origin']=q_origin
    manifest['selection_q']=a.selection_q
    manifest['model_source']=a.source
    manifest['legacy']=a.legacy
    if a.allow_missing_model_genes:
        manifest['missing_model_genes']=dict(policy='Zero joint mass for genes absent from the fixed model; AA-max join and 2^-80 floor unchanged; add-one V/J counts over all probe keys; source calls retained',
            keys=int((~probes['pgen_model_supported']).sum()))
    if a.matcher in {'blosum62','positional-blosum62'}:
        manifest['matching_model']=dict(preset=a.matcher,radius_factor=factor,scope='3,0,0,3',normalization='maximum matrix substitution penalty; positional weights use fixed scale100')
    if a.matcher in GRADED_FAMILIES:
        _,scale=graded_family_params(a.matcher,a.mismatches,factor)
        manifest['matching_model']=dict(preset=a.matcher,radius_factor=factor,max_edits=a.mismatches,normalization=scale,match_v=a.matcher.startswith('v-'),match_j=True)
        manifest['counting']='source rows and reads; individual alignment distance threshold; exact J; hard V only for v-prefixed families'
    if shape is not None:
        manifest['selection_profile_kind']=shape['kind']
        manifest['probability']='global Q times AA-max selected joint mass times residual within-V J usage ratio; floor2^-80'
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return 0
