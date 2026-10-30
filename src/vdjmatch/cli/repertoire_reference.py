"""Reference-null repertoire enrichment adapter; historical resampling stays separate."""
from pathlib import Path


def _genes(frame,locus,*,species='human'):
    import polars as pl
    from .repertoire_enrichment import _single_v
    from ..evalue.mock_reference import _modal_gene
    frame=_single_v(frame,'v',species=species,locus=locus)
    frame=frame.with_columns(pl.col('j').str.replace(r'\*.*$','').alias('j'))
    good=pl.col('j').str.contains('^'+locus+'J[A-Z0-9/.-]+$').fill_null(False)
    if frame.filter(~good).height:
        modal=_modal_gene(frame.filter(good),'j',species=species,locus=locus,kind='j').split('*')[0]
        frame=frame.with_columns(pl.when(good).then(pl.col('j')).otherwise(pl.lit(modal)).alias('j'))
    return frame

def _read_control(path,*,species,locus,batch_size=65536):
    import polars as pl
    import polars.selectors as cs
    from ..evalue.mock_reference import control_events
    from ..io.airr import _read_table
    control_columns=cs.matches(r'(?i)^(junction_aa|cdr3aa|cdr3|cdr3_aa|v|j|v_call|j_call|v_gene|j_gene|v\.segm|j\.segm|locus|chain|species|organism|junction|cdr3nt|junction_nt|cdr3_nt|sequence_id|sample_id|sample|repertoire_id|occurrences|convergence|incidence|original_v|original_j|v_imputed|j_imputed)$')
    p=Path(path)
    suffix=p.with_suffix('').suffix if p.suffix.lower()=='.gz' else p.suffix
    scan=(pl.scan_parquet(p) if suffix.lower() in {'.parquet','.pq'} else
          pl.scan_csv(p,separator=',' if suffix.lower()=='.csv' else '\t',infer_schema_length=0))
    # Explicit occurrence counts are authoritative in control_events. NT/event
    # identity columns are then unused; do not materialize them for huge controls.
    if 'occurrences' in {name.lower() for name in scan.collect_schema().names()}:
        control_columns=control_columns-cs.matches(r'(?i)^(junction|cdr3nt|junction_nt|cdr3_nt|sample_id|sample|repertoire_id|convergence|incidence)$')
        from vdjtools.io import iter_table_batches
        parts=[];matched_population=False
        for chunk in iter_table_batches(path,batch_size=batch_size,columns=control_columns):
            try:
                part=control_events(chunk,species=species,locus=locus)
            except ValueError as error:
                if str(error)=='control has no rows matching requested species/locus':
                    continue
                if str(error)=='control has no positive-weight junctions':
                    matched_population=True
                    continue
                raise
            matched_population=True
            # Chunk-local modes are not whole-control modes. Unresolved genes
            # require the existing global path; never silently change imputation.
            if any(part[c].any() for c in ('v_imputed','j_imputed') if c in part.columns):
                break
            canonical_species,canonical_locus=part['species'][0],part['locus'][0]
            parts.append(part.select('junction_aa','v','j','event_weight'))
        else:
            if not parts:
                raise ValueError('control has no positive-weight junctions' if matched_population else
                                 'control has no rows matching requested species/locus')
            return pl.concat(parts).group_by('junction_aa','v','j').agg(pl.col('event_weight').sum()).with_columns(
                pl.lit(canonical_species).alias('species'),pl.lit(canonical_locus).alias('locus')).sort('junction_aa','v','j')
        del parts

    return control_events(_read_table(path,columns=control_columns),species=species,locus=locus)

def _observations(frame,*,locus,species):
    import polars as pl
    loci=['TRA','TRB'] if locus=='paired' else [locus]
    if locus!='paired':frame=frame.with_columns(pl.col('locus').fill_null(locus))
    valid=pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)&pl.col('locus').is_in(loci).fill_null(False)
    selected=frame.filter(valid);parts=[]
    for (chain,),part in selected.partition_by('locus',as_dict=True).items():parts.append(_genes(part,chain,species=species))
    if not parts:raise ValueError('sample has no valid selected-locus receptors')
    selected=pl.concat(parts)
    if locus=='paired':
        if selected['pair_id'].null_count() or selected['pair_id'].str.strip_chars().eq('').any():raise ValueError('paired mode requires non-null pair_id for every chain, including orphans')
        obs=selected.select(pl.col('pair_id').alias('receptor_id'),'locus',pl.col('cdr3').alias('junction_aa'),'v','j','count')
    else:
        obs=selected.group_by('locus','cdr3','v','j').agg(pl.col('count').sum()).sort('locus','cdr3','v','j').with_row_index('_id').select(pl.col('_id').cast(pl.String).alias('receptor_id'),'locus',pl.col('cdr3').alias('junction_aa'),'v','j','count')
    return selected,obs

def main(a):
    import json,time
    from importlib.metadata import version
    import polars as pl
    from vdjtools.biomarker.stats import fdr_bh
    from .. import db,io
    from ..aggregate.candidates import PMHC
    from ..aggregate.repertoire_enrichment import reference_enrichment
    from ..evalue.mock_reference import compute_pgen,load_pgen_features,_pgen_model_identity
    from ..io.airr import _read_table
    from ..runner.annotation import _validate_output_paths
    from .repertoire_enrichment import _repertoire_strata

    started=time.monotonic()
    if a.mock_replicates<200 or a.seed<0 or a.threads<1:
        raise ValueError('at least 200 mocks, nonnegative seed and positive threads required')
    if not 0<a.nominal_alpha<1 or not 0<a.fdr_alpha<1:
        raise ValueError('significance thresholds must lie between zero and one')
    for mode,alpha in [('nominal',a.nominal_alpha),('bh',a.fdr_alpha)]:
        if a.significance in (mode,'both') and 1/(a.mock_replicates+1)>=alpha:
            raise ValueError(f'mock resolution cannot reach {mode} threshold; increase --mock-replicates')
    metadata_path=Path(a.metadata).expanduser().resolve()
    metadata=_read_table(metadata_path)
    if not {'sample_id','path'}<=set(metadata.columns) or not metadata.height or metadata['sample_id'].n_unique()!=metadata.height:
        raise ValueError('nonempty unique sample_id/path metadata required')
    if any(metadata[c].null_count() or metadata[c].cast(pl.String).str.strip_chars().eq('').any() for c in ['sample_id','path']):
        raise ValueError('sample_id/path cannot be missing')
    paths=[(metadata_path.parent/Path(p).expanduser()).resolve() for p in metadata['path']]
    if len(set(paths))!=len(paths):raise ValueError('each raw repertoire path must appear once')
    loci=['TRA','TRB'] if a.locus=='paired' else [a.locus]
    control_paths={locus:Path(a.control_tra if locus=='TRA' and a.locus=='paired' else a.control).expanduser().resolve() for locus in loci} if a.locus!='paired' or a.control_tra else {}
    if not control_paths:raise ValueError('paired mode requires --control-tra and --control (TRB)')
    features={locus:(a.pgen_features_tra if locus=='TRA' and a.locus=='paired' else a.pgen_features) for locus in loci}
    inputs=[metadata_path,Path(a.vdjdb).expanduser().resolve(),*paths,*control_paths.values()]
    for p in features.values():
        if p:
            p=Path(p).expanduser().resolve();inputs.extend([p,p.with_name(p.name.removesuffix('.features.parquet')+'.manifest.json')])
    prefix=_validate_output_paths(a.output_prefix,[str(p) for p in inputs])
    suffixes=['summary.tsv','mocks.parquet','hits.tsv','ingestion.tsv','reference-audit.tsv','manifest.json']
    if any(Path(str(prefix)+'.'+s).exists() for s in suffixes):raise FileExistsError('fresh output files must not exist')
    from ..db.cache import sha256
    def identities():return {str(p):db.provenance(p) if p.is_dir() else sha256(p) for p in inputs}
    import seqtree._core as seq_native
    import vdjtools._core as pgen_native
    source_files=[*Path(__file__).resolve().parents[1].rglob('*.py'),Path(seq_native.__file__),Path(pgen_native.__file__)]
    def freeze():return {str(p):sha256(p) for p in source_files}
    before,source_before=identities(),freeze()
    models={locus:_pgen_model_identity(a.species,locus) for locus in loci}
    reference=db.load(inputs[1],species={'human':'HomoSapiens','mouse':'MusMusculus'}[a.species])
    valid=pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)&pl.col('gene').is_in(loci)
    population=a.reference_population
    if population=='curated-benchmark':valid &= pl.col('reference_valid')
    reference=reference.filter(valid)
    if a.epitope:
        missing=set(a.epitope)-set(reference['epitope'])
        if missing:raise ValueError(f'requested epitopes absent from reference: {sorted(missing)}')
        reference=reference.filter(pl.col('epitope').is_in(a.epitope))
    if not reference.height:raise ValueError('no eligible reference records')
    inventory=reference.select(PMHC).unique().sort(PMHC,nulls_last=True).with_columns(pl.struct(PMHC).struct.json_encode().alias('group_id'))
    reference=reference.join(inventory,on=PMHC,validate='m:1',nulls_equal=True)
    pools={};ref_parts=[];audits=[];prepared={};stages={}
    for locus in loci:
        t=time.monotonic()
        raw=_read_control(control_paths[locus],species=a.species,locus=locus)
        stages[locus+'_control_ingestion_seconds']=time.monotonic()-t
        print(f'{locus}: control ingestion complete ({raw.height} keys)',flush=True)
        control=_genes(raw,locus,species=a.species).select('junction_aa','v','j').unique().sort('junction_aa','v','j')
        part=reference.filter(pl.col('gene')==locus).rename({'cdr3':'junction_aa','gene':'locus'})
        if not part.height:continue
        part=_genes(part,locus,species=a.species)
        if features[locus]:
            pgen,prepared[locus]=load_pgen_features(features[locus],part,control,species=a.species,locus=locus,mismatches=a.pgen_mismatches)
        else:
            pgen=compute_pgen(part,control,species=a.species,locus=locus,threads=a.threads,mismatches=a.pgen_mismatches)
        bins=_repertoire_strata(pgen)
        part=part.join(bins,on='junction_aa',validate='m:1')
        pools[locus]=control.join(bins,on='junction_aa',validate='m:1')
        audits.append(part.group_by('group_id','locus','pgen_bin').agg(pl.len().alias('source_records'),pl.struct('junction_aa','v','j').n_unique().alias('reference_keys')))
        ref_parts.append(part.select('group_id','locus','junction_aa','v','j','pgen_bin'))
        stages[locus+'_prepare_seconds']=time.monotonic()-t
        print(f'{locus}: {part.height} reference records, {control.height} distinct control keys',flush=True)
    ref=pl.concat(ref_parts)
    summaries=[];mocks=[];hits=[];reports=[]
    for number,(sample,path) in enumerate(zip(metadata['sample_id'],paths)):
        t=time.monotonic();frame,report=io.read_rearrangement(path,valid_aa=False,return_report=True)
        selected,obs=_observations(frame,locus=a.locus,species=a.species)
        summary,mock,hit=reference_enrichment(obs,ref,pools,replicates=a.mock_replicates,seed=a.seed,policy=a.mock_policy,threads=a.threads,mismatches=a.mismatches,match_j=a.match_j)
        summaries.append(summary.with_columns(pl.lit(sample).alias('sample_id')))
        mocks.append(mock.with_columns(pl.lit(sample).alias('sample_id')))
        hits.append(hit.with_columns(pl.lit(sample).alias('sample_id')))
        reports.append(dict(sample_id=sample,**report,selected_chain_rows=selected.height,excluded_rows=frame.height-selected.height,imputed_v_rows=int(selected['v_imputed'].sum()),unique_receptors=obs['receptor_id'].n_unique()))
        stages[f'sample_{number+1}_seconds']=time.monotonic()-t
        print(f'sample {number+1}/{len(paths)} complete in {stages[f"sample_{number+1}_seconds"]:.3f}s',flush=True)
    summary=pl.concat(summaries)
    if a.significance in ('bh','both'):
        if a.fdr_family=='target':
            summary=pl.concat([part.with_columns(pl.Series('q_receptors',fdr_bh(part['p_receptors'].fill_null(1.).to_numpy()))) for part in summary.partition_by('group_id')])
        else:summary=summary.with_columns(pl.Series('q_receptors',fdr_bh(summary['p_receptors'].fill_null(1.).to_numpy())))
        summary=summary.with_columns(pl.when(pl.col('status')=='ok').then(pl.col('q_receptors')).otherwise(None).alias('q_receptors'))
        summary=summary.with_columns(((pl.col('q_receptors')<a.fdr_alpha)&(pl.col('excess_receptors')>0)).alias('significant_bh'))
    if a.significance in ('nominal','both'):
        summary=summary.with_columns(((pl.col('p_receptors')<a.nominal_alpha)&(pl.col('excess_receptors')>0)).alias('significant_nominal'))
    after,source_after=identities(),freeze()
    if before!=after or source_before!=source_after or models!={locus:_pgen_model_identity(a.species,locus) for locus in loci}:
        raise ValueError('raw inputs, model or software changed during the run')
    prefix.parent.mkdir(parents=True,exist_ok=True);outputs={}
    for suffix,table in [('summary.tsv',summary),('mocks.parquet',pl.concat(mocks)),('hits.tsv',pl.concat(hits)),('reference-audit.tsv',pl.concat(audits))]:
        table=table.join(inventory,on='group_id',validate='m:1').drop('group_id')
        path=Path(str(prefix)+'.'+suffix)
        if suffix.endswith('parquet'):table.write_parquet(path)
        else:table.write_csv(path,separator='\t')
        outputs[path.name]=sha256(path)
    path=Path(str(prefix)+'.ingestion.tsv');pl.DataFrame(reports,infer_schema_length=None).write_csv(path,separator='\t');outputs[path.name]=sha256(path)
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(dict(null='reference',mock_policy=a.mock_policy,inputs=before,inputs_frozen=True,source=source_before,source_frozen=True,
        models=models,prepared_pgen=prepared,software={n:version(n) for n in ['vdjmatch','vdjtools','seqtree','polars','numpy','scipy']},
        options=vars(a)|{'func':a.func.__name__},reference_population=population,counting_unit='unique receptor identity',
        paired_rule='union of available chain supports, counted once per receptor/target; independently calibrated by the same mock procedure',
        pgen_definition=f'marginal full-junction Hamming<={a.pgen_mismatches}; rounded log2(Pgen+2**-60)',
        control_law='uniform unique junction/V/J keys, without replacement per target/bin; source event/read weights are not sampling weights',
        archive_differences='marginal model Pgen without cohort-fitted VJ factor; complete zero-hit rows; explicit matching and conditional mock options; independent locus streams',
        significance_family=a.fdr_family,significance=a.significance,unavailable_null='null-valued inference; retained as nonrejectable hypotheses in BH family',stages=stages,seconds=time.monotonic()-started,outputs=outputs),indent=2)+'\n')
    return 0
