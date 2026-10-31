"""Exact reference-panel enrichment preserving source-row and read counts."""
from pathlib import Path


def register(sub):
    p=sub.add_parser('repertoire-exact',help=__doc__)
    p.add_argument('metadata',help='TSV sample_id,path and optional FDR-group columns')
    p.add_argument('--vdjdb',required=True)
    p.add_argument('--control',required=True,help='Distinct junction_aa,v,j probe keys; optional joint_pgen column')
    p.add_argument('--model-vj-counts',help='Generated-model v,j,records counts from the selected source')
    p.add_argument('--epitope',action='append')
    p.add_argument('--pool-mhc',action='store_true',help='Explicit original epitope-union reproduction')
    p.add_argument('--locus',choices=['TRA','TRB'],default='TRB')
    p.add_argument('--species',choices=['human','mouse'],default='human')
    p.add_argument('--selection-q',type=float,default=None,help='Q override; default is calibrated for source/species/locus')
    p.add_argument('--selection-profile',help='Opt-in length | vj-length | vdjtools calibration JSON')
    p.add_argument('--source',default=None,help='Model source; default isalgo (legacy human: olga-legacy)')
    p.add_argument('--legacy',action='store_true',help='Compiled reproduction with historical source and Q=1')
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


def _source(path, *, locus):
    import polars as pl
    from .. import io
    frame=io.read_rearrangement(path,valid_aa=False).filter(pl.col("locus")==locus)
    return frame.select(pl.col('cdr3').alias('junction_aa'),
        pl.col('v','j').str.split('*').list.first(),'count')


def _joint_pgen(probes,*,species,locus,threads,source="isalgo",profile=None):
    import polars as pl
    from vdjtools.model import load_bundled,native
    model=load_bundled(locus,source,organism=species,collapse=False)
    values=native.pgen_aa_gene_batch(model,probes['junction_aa'].to_list(),probes['v'].to_list(),
                                   probes['j'].to_list(),mismatches=0,threads=threads)
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
        for endpoint in ['records','reads']:
            actual=observed['matched_'+endpoint];values=mock['matched_'+endpoint].to_numpy()
            original=present['matched_'+endpoint].to_numpy()
            transformed=np.log2(original+1.) if endpoint=='reads' else original
            effect_actual=np.log2(actual+1.) if endpoint=='reads' else actual
            mean=float(transformed.mean()) if len(original) else None
            sd=float(transformed.std(ddof=1)) if len(original)>1 else None
            with np.errstate(divide='ignore',invalid='ignore'):
                effect=float(np.divide(effect_actual-mean,sd)) if sd is not None else None
            if effect is not None and np.isnan(effect):effect=None
            row.update({'mock_mean_'+endpoint:float(values.mean()),'p_'+endpoint:float((1+(values>=actual).sum())/(len(values)+1)),
                'original_p_'+endpoint:float((1+(original>=actual).sum())/(len(original)+1)) if observed['matched_records']>0 else None,
                'original_mean_'+endpoint:mean,'original_sd_'+endpoint:sd,
                'original_effect_'+endpoint:effect,
                'excess_'+endpoint:float(actual-values.mean())})
        rows.append(row)
    return pl.DataFrame(rows,infer_schema_length=None)


def main(a):
    import json,time
    from importlib.metadata import version
    import polars as pl
    from vdjtools.biomarker.stats import fdr_bh
    from .. import db
    from ..db.cache import sha256
    from ..io.airr import _read_table
    from ..aggregate.candidates import PMHC
    from ..aggregate.repertoire_enrichment import exact_reference_counts
    from ..evalue.mock_reference import compiled_reference_strata,sample_stratified_panels,_pgen_model_identity
    from ..runner.annotation import _validate_output_paths
    from vdjtools.cores import kernel_threads
    a.threads = kernel_threads(a.threads)
    if a.mock_replicates<2 or a.seed<0 or a.threads<1:raise ValueError('at least two mocks, nonnegative seed and positive threads required')
    start=time.monotonic();meta_path=Path(a.metadata).expanduser().resolve();meta=_read_table(meta_path)
    needed=['sample_id','path',*a.fdr_group]
    if not set(needed)<=set(meta.columns) or not meta.height or meta['sample_id'].n_unique()!=meta.height:
        raise ValueError('unique sample_id,path and declared FDR columns required')
    if any(meta[c].null_count() or meta[c].str.strip_chars().eq('').any() for c in needed):raise ValueError('sample identities, paths and FDR groups cannot be missing')
    paths=[(meta_path.parent/Path(p).expanduser()).resolve() for p in meta['path']]
    if len(set(paths))!=len(paths):raise ValueError('sample paths must be unique')
    if a.selection_profile and (a.legacy or a.model_vj_counts):
        raise ValueError('selection profiles require modern mode and replace model V/J counts')
    if not a.selection_profile and not a.model_vj_counts:
        raise ValueError('--model-vj-counts is required without a selection profile')
    inputs=[meta_path,Path(a.vdjdb).expanduser().resolve(),Path(a.control).expanduser().resolve(),*paths]
    if a.model_vj_counts:
        inputs.append(Path(a.model_vj_counts).expanduser().resolve())
    if a.selection_profile and a.selection_profile not in {'length','vj-length'}:
        inputs.append(Path(a.selection_profile).expanduser().resolve())
    prefix=_validate_output_paths(a.output_prefix,[str(p) for p in inputs])
    suffixes=['summary.tsv','counts.tsv','vj-audit.tsv','reference-bins.tsv','manifest.json']
    if any(Path(str(prefix)+'.'+s).exists() for s in suffixes):raise FileExistsError('fresh output prefix required')
    before={str(p):sha256(p) for p in inputs}
    import vdjtools._core as native_core
    import vdjtools.model.io as model_io
    import vdjtools.model.native as model_native
    import vdjtools.model.selection as model_selection
    source_files=[*Path(__file__).resolve().parents[1].rglob('*.py'),Path(native_core.__file__),Path(model_io.__file__),Path(model_native.__file__),Path(model_selection.__file__),Path(model_selection.__file__).with_name('selection_factors.json'),Path(model_selection.__file__).with_name('selection_profiles.json')]
    source_before={str(p):sha256(p) for p in source_files}
    reference=db.load(inputs[1],species={'human':'HomoSapiens','mouse':'MusMusculus'}[a.species]).filter(pl.col('gene')==a.locus)
    if a.epitope:
        if set(a.epitope)-set(reference['epitope']):raise ValueError('requested epitope absent from selected reference population')
        reference=reference.filter(pl.col('epitope').is_in(a.epitope))
    if not reference.height:raise ValueError('empty selected reference population')
    labels=['epitope'] if a.pool_mhc else PMHC
    inventory=reference.select(labels).unique().sort(labels).with_columns(pl.struct(labels).struct.json_encode().alias('group_id'))
    keys=['junction_aa','v','j']
    reference=reference.join(inventory,on=labels,nulls_equal=True).select('group_id',pl.col('cdr3').alias('junction_aa'),pl.col('v','j').str.split('*').list.first()).unique()
    probes=_read_table(inputs[2])
    model_identity=None
    shape=None
    if a.selection_profile:
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
    if 'joint_pgen' in probes.columns:
        probes=probes.select(*keys,pl.col('joint_pgen').cast(pl.Float64).alias('pgen'))
        if reference.select(keys).join(probes.select(keys),on=keys,how='anti').height:raise ValueError('external joint_pgen probes must include every reference key')
        pgen_origin='external gene-masked probability asset'
    else:
        model_identity=_pgen_model_identity(a.species,a.locus,source=a.source)
        probes=pl.concat([probes.select(keys),reference.select(keys)]).unique(maintain_order=True)
        probes=_joint_pgen(probes,species=a.species,locus=a.locus,threads=a.threads,source=a.source,profile=shape)
        pgen_origin='native gene-allele-summed joint probability'
    vj=[]
    for path in paths:
        frame=_source(path,locus=a.locus);vj.append(frame.group_by('v','j').len(name='records'))
    if shape is None:
        model=_read_table(Path(a.model_vj_counts).expanduser().resolve())
        probes,audit=compiled_reference_strata(probes,pl.concat(vj),model,selection_q=a.selection_q)
    else:
        selected=pl.DataFrame(shape['selected_vj']).rename({'v_call':'v','j_call':'j'})
        probes,audit=compiled_reference_strata(probes,pl.concat(vj),selection_q=a.selection_q,
            model_vj_probabilities=selected)
    reference=reference.join(probes.select(*keys,'pgen_bin'),on=keys,validate='m:1')
    mocks=sample_stratified_panels(reference,probes.select(*keys,'pgen_bin'),replicates=a.mock_replicates,seed=a.seed,sampler=a.sampler)
    panels=pl.concat([reference.select('group_id',*keys).with_columns(pl.lit(-1,dtype=pl.Int64).alias('replicate')),mocks.select('group_id',*keys,pl.col('replicate').cast(pl.Int64))])
    counts=[]
    for i,(sample,path) in enumerate(zip(meta['sample_id'],paths),1):
        counts.append(exact_reference_counts(_source(path,locus=a.locus),panels).with_columns(pl.lit(sample).alias('sample_id')))
        print(f'counted sample {i}/{len(paths)}',flush=True)
    counts=pl.concat(counts);summary=_summaries(counts).join(meta.select('sample_id',*a.fdr_group),on='sample_id',validate='m:1')
    parts=[]
    for part in summary.partition_by(['group_id',*a.fdr_group]):
        for endpoint in ['records','reads']:
            for stem in ['', 'original_']:
                col=stem+'p_'+endpoint
                values=part[col].to_list()
                adjusted=iter(fdr_bh(part[col].drop_nulls().to_numpy()))
                part=part.with_columns(pl.Series(stem+'q_'+endpoint,[next(adjusted) if value is not None else None for value in values],dtype=pl.Float64))
            part=part.with_columns(((pl.col('p_'+endpoint)<.05)&(pl.col('excess_'+endpoint)>0)).alias('significant_nominal_'+endpoint),
                ((pl.col('q_'+endpoint)<.1)&(pl.col('excess_'+endpoint)>0)).alias('significant_bh_'+endpoint),
                ((pl.col('original_q_'+endpoint)<.1)&(pl.col('original_effect_'+endpoint)>0)).alias('original_significant_bh_'+endpoint))
        parts.append(part)
    summary=pl.concat(parts).join(inventory,on='group_id',validate='m:1')
    outputs={'summary.tsv':summary,'counts.tsv':counts,'vj-audit.tsv':audit,'reference-bins.tsv':reference.group_by('group_id','pgen_bin').len(name='reference_keys')}
    prefix.parent.mkdir(parents=True,exist_ok=True)
    for suffix,frame in outputs.items():frame.write_csv(str(prefix)+'.'+suffix,separator='\t')
    if before!={str(p):sha256(p) for p in inputs} or source_before!={str(p):sha256(p) for p in source_files}:raise RuntimeError('inputs or source changed during calculation')
    if model_identity is not None and model_identity!=_pgen_model_identity(a.species,a.locus,source=a.source):raise RuntimeError('recombination model changed during calculation')
    manifest=dict(inputs_frozen=True,source_frozen=True,inputs=before,source=source_before,options=vars(a)|{'func':None},pgen_origin=pgen_origin,probability='selection_q times AA-max joint Pgen times within-V add-one J usage ratio; floor2^-80',counting='source rows and reads; exact V/J/junction',fdr_family=['group_id',*a.fdr_group],seconds=time.monotonic()-start,software={n:version(n) for n in ['vdjmatch','vdjtools','polars','numpy']},outputs={Path(str(prefix)+'.'+s).name:sha256(Path(str(prefix)+'.'+s)) for s in outputs})
    if a.sampler=='r':
        import subprocess,shutil
        manifest['r_version']=subprocess.run([shutil.which('Rscript'),'--version'],capture_output=True,text=True,check=True).stdout.strip()
    manifest['model_identity']=model_identity
    manifest['selection_q_origin']=q_origin
    manifest['selection_q']=a.selection_q
    manifest['model_source']=a.source
    manifest['legacy']=a.legacy
    if shape is not None:
        manifest['selection_profile_kind']=shape['kind']
        manifest['probability']='global Q times AA-max selected joint mass times residual within-V J usage ratio; floor2^-80'
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return 0
