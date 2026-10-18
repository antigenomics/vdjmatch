"""Fresh strict Pgen/count-matched repertoire enrichment one repertoire and one mock at a time."""
from pathlib import Path


def _single_v(frame,column):
    """Pool same-gene alleles; refuse unresolved multi-gene calls without dropping rows."""
    import polars as pl
    genes=pl.col(column).cast(pl.String).str.to_uppercase().str.split(',').list.eval(pl.element().str.split('*').list.first().str.strip_chars()).list.unique()
    ambiguous=frame.select((genes.list.len()>1).sum()).item()
    if ambiguous:raise ValueError(f'{ambiguous}unresolved multi-gene V calls; exact V matching cannot choose a gene or discard these records')
    return frame.with_columns(genes.list.first().alias(column))


def register(subparsers):
    p=subparsers.add_parser('repertoire-enrichment',help=__doc__)
    p.add_argument('metadata',help='TSV with unique sample_id,path; paths relative to this file')
    p.add_argument('--vdjdb',required=True);p.add_argument('--control',required=True)
    p.add_argument('--pgen-mismatches',type=int,choices=[0,1],default=1,help='0 exact sequence or1 inclusive Hamming ball; prepared and fresh probabilities must agree')
    p.add_argument('--pgen-features',help='validated durable control Pgen bank; missing query masses are freshly calculated')
    p.add_argument('--species',choices=['human','mouse'],default='human')
    p.add_argument('--locus',choices=['TRA','TRB'],default='TRB')
    p.add_argument('--reference-population',choices=['all-valid','curated-benchmark'],default='all-valid')
    p.add_argument('--mock-replicates',type=int,default=200);p.add_argument('--seed',type=int,default=0)
    p.add_argument('--threads',type=int,default=4);p.add_argument('--output-prefix',required=True)
    p.set_defaults(func=main)


def main(a):
    import json,time
    from importlib.metadata import version
    import numpy as np
    import polars as pl
    import seqtree._core as search_native
    import vdjtools._core as pgen_native
    import vdjtools.io.convert as converter
    import vdjtools.model.io as model_io
    from .. import db,io
    from ..aggregate.candidates import PMHC
    from ..aggregate.repertoire_enrichment import matched_enrichment
    from ..db.cache import sha256
    from ..evalue.mock_reference import control_events,compute_pgen,pgen_strata,load_pgen_features,_pgen_model_identity,EVENT_SAMPLER
    from ..io.airr import _read_table
    from ..io.columns import _resolve
    from ..match.engine import VdjdbIndex
    from ..match.scope import search_params
    from ..runner.annotation import _validate_output_paths
    start=time.monotonic();stages={}
    if type(a.threads) is not int or a.threads<1 or a.mock_replicates<200 or a.seed<0:
        raise ValueError('positive threads, at least200mocks and nonnegative seed required')
    metadata_path=Path(a.metadata).expanduser().resolve();metadata=_read_table(metadata_path)
    if not {'sample_id','path'}<=set(metadata.columns) or metadata['sample_id'].n_unique()!=metadata.height or not metadata.height:
        raise ValueError('metadata requires nonempty unique sample_id and path')
    if metadata.select(pl.any_horizontal(pl.col(c).is_null()|(pl.col(c).str.strip_chars()=='') for c in ['sample_id','path']).any()).item():
        raise ValueError('metadata sample identities and paths cannot be missing')
    paths=[(metadata_path.parent/Path(p).expanduser()).resolve() for p in metadata['path']]
    if len(set(paths))!=len(paths):raise ValueError('each raw repertoire path must appear once')
    inputs=[metadata_path,Path(a.control).expanduser().resolve(),Path(a.vdjdb).expanduser().resolve(),*paths]
    if a.pgen_features:
        p=Path(a.pgen_features).expanduser().resolve();inputs.extend([p,p.with_name(p.name.removesuffix('.features.parquet')+'.manifest.json')])
    prefix=_validate_output_paths(a.output_prefix,[str(p) for p in inputs])
    suffixes=['summary.tsv','mocks.parquet','matching-audit.tsv','ingestion.tsv','manifest.json']
    if any(Path(str(prefix)+'.'+s).exists() for s in suffixes):raise FileExistsError('fresh output files must not exist')
    def identities():return {str(p):db.provenance(p) if p.is_dir() else sha256(p) for p in inputs}
    source_files=[*Path(__file__).resolve().parents[1].rglob('*.py'),Path(search_native.__file__),Path(pgen_native.__file__),Path(converter.__file__),*Path(model_io.__file__).parent.rglob('*.py')]
    def freeze():return {str(p):sha256(p) for p in source_files}
    before,source_before=identities(),freeze();model_before=_pgen_model_identity(a.species,a.locus)
    raw_control=_read_table(inputs[1]);v_column=_resolve(raw_control).get('v')
    if v_column:raw_control=_single_v(raw_control,v_column)
    controls=control_events(raw_control,species=a.species,locus=a.locus);del raw_control
    stages['control_ingestion_seconds']=time.monotonic()-start
    print(f'control ingestion: {controls.height} event keys',flush=True)
    t=time.monotonic();prepared=None
    if a.pgen_features:
        pgen,prepared=load_pgen_features(a.pgen_features,controls,controls,species=a.species,locus=a.locus,mismatches=a.pgen_mismatches)
    else:
        pgen=compute_pgen(controls,controls.head(0),species=a.species,locus=a.locus,threads=a.threads,mismatches=a.pgen_mismatches)
    controls=controls.join(pgen_strata(pgen),on='junction_aa',validate='m:1').select('junction_aa','v','event_weight','pgen_bin').with_row_index('query_id')
    stages['control_pgen_seconds']=time.monotonic()-t
    print('control Pgen bins ready',flush=True)
    t=time.monotonic()
    reference=db.load(inputs[2],species={'human':'HomoSapiens','mouse':'MusMusculus'}[a.species],gene=a.locus)
    valid=pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
    reference=reference.filter(valid & (pl.col('reference_valid') if a.reference_population=='curated-benchmark' else pl.lit(True))).with_columns(pl.lit(True).alias('reference_valid'))
    inventory=reference.select(PMHC).unique().sort(PMHC,nulls_last=True)
    index=VdjdbIndex.build(reference);del reference
    stages['reference_seconds']=time.monotonic()-t
    target_index=inventory.with_row_index('_target')
    reports=[];summaries=[];mock_tables=[];audits=[];search_batches=0;search_keys=0;source_records=0
    for sample_number,(sample,path) in enumerate(zip(metadata['sample_id'],paths)):
        t=time.monotonic()
        frame,report=io.read_rearrangement(path,valid_aa=False,return_report=True)
        valid=pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
        locus=pl.col('locus').eq(a.locus).fill_null(False)
        called=(pl.col('v').is_not_null()&(pl.col('v').str.strip_chars()!='')).fill_null(False)
        selected=_single_v(frame.filter(valid&locus&called),'v')
        if not selected.height:raise ValueError('repertoire has no valid selected-locus V-labelled records')
        reports.append(dict(sample_id=sample,**report,selected_records=selected.height,selected_reads=int(selected['count'].sum()),
            invalid_aa_records=frame.filter(~valid).height,other_or_missing_locus_records=frame.filter(valid&~locus).height,
            missing_v_records=frame.filter(valid&locus&~called).height,
            retained_noncanonical_anchor_records=selected.filter(~pl.col('cdr3').str.contains(r'^C[ACDEFGHIKLMNPQRSTVWY]*[FW]$')).height))
        del frame
        missing=selected.select(pl.col('cdr3').alias('junction_aa')).unique().join(pgen.select('junction_aa'),on='junction_aa',how='anti')
        sample_pgen=pgen.join(selected.select(pl.col('cdr3').alias('junction_aa')).unique(),on='junction_aa',how='semi')
        if missing.height:
            sample_pgen=pl.concat([sample_pgen,compute_pgen(missing,controls.head(0),species=a.species,locus=a.locus,threads=a.threads,mismatches=a.pgen_mismatches)])
        observations=selected.select('cdr3','count').with_row_index('query_id',offset=controls.height).join(pgen_strata(sample_pgen),left_on='cdr3',right_on='junction_aa',validate='m:1').with_columns(pl.lit(sample).alias('sample_id')).select('sample_id','query_id','count','pgen_bin')
        real_keys=selected.select('cdr3','v');source_records+=selected.height
        del selected,sample_pgen,missing
        def measure(rows,counts):
            nonlocal real_keys,search_batches,search_keys
            if real_keys is not None:
                batch=real_keys[rows-controls.height];real_keys=None
            else:
                batch=controls[rows].select(pl.col('junction_aa').alias('cdr3'),'v')
            batch=batch.with_row_index('query_id').with_columns(pl.Series('count',counts))
            hits=index.annotate(batch,search_params('1,0,0,1'),gene=a.locus,match_v=True,exclude_exact=False,threads=a.threads)
            member=hits.select('query_id',*PMHC).unique();del hits
            totals=member.join(batch.select('query_id','count'),on='query_id',validate='m:1').group_by(PMHC).agg(pl.len().alias('records'),pl.col('count').sum().alias('reads')).join(target_index,on=PMHC,validate='1:1',nulls_equal=True)
            records=np.zeros(inventory.height,dtype=np.int64);reads=np.zeros_like(records)
            records[totals['_target'].to_numpy()]=totals['records'].to_numpy()
            reads[totals['_target'].to_numpy()]=totals['reads'].to_numpy()
            search_batches+=1;search_keys+=batch.height
            if search_batches%(a.mock_replicates+1) in (1,0) or search_batches%20==0:
                print(f'sample {sample_number+1}/{len(paths)}, matched batch {search_batches}, rows {batch.height}',flush=True)
            return records,reads
        summary,mocks,audit=matched_enrichment(observations,controls,None,inventory,replicates=a.mock_replicates,seed=a.seed+sample_number,measure=measure)
        summaries.append(summary);mock_tables.append(mocks);audits.append(audit)
        stages[f'sample_{sample_number+1}_seconds']=time.monotonic()-t
        del observations,measure
    summary=pl.concat(summaries);mocks=pl.concat(mock_tables);audit=pl.concat(audits)
    after,source_after=identities(),freeze()
    if before!=after or source_before!=source_after or model_before!=_pgen_model_identity(a.species,a.locus):
        raise ValueError('raw, model or source/native changed during fresh repertoire enrichment')
    prefix.parent.mkdir(parents=True,exist_ok=True);tables={}
    for suffix,table in [('summary.tsv',summary),('mocks.parquet',mocks),('matching-audit.tsv',audit),('ingestion.tsv',pl.DataFrame(reports,infer_schema_length=None))]:
        p=Path(str(prefix)+'.'+suffix)
        if suffix.endswith('parquet'):table.write_parquet(p)
        else:table.write_csv(p,separator='\t')
        tables[p.name]=sha256(p)
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(dict(inputs=before,input_after=after,inputs_frozen=True,
        source_before=source_before,source_after=source_after,source_frozen=True,pgen_model=model_before,prepared_pgen=prepared,
        software={n:version(n) for n in ['vdjmatch','vdjtools','seqtree','polars','numpy','scipy']},species=a.species,locus=a.locus,
        source_records=source_records,control_keys=controls.height,matched_query_rows=search_keys,native_search_batches=search_batches,
        processing="one real repertoire or one mock matched at a time; fixed reference index reused; only aggregate mock endpoints retained",sample_seeds=[a.seed+i for i in range(len(paths))],
        scope='full-junction Hamming<=1 and exact V gene, allele-pooled; exact reference matches retained',reference_population=a.reference_population,
        mock_replicates=a.mock_replicates,seed=a.seed,threads=a.threads,event_sampler=EVENT_SAMPLER,
        pgen_mismatches=a.pgen_mismatches,pgen_definition=f'inclusive Hamming<={a.pgen_mismatches} marginal full-junction probability, summed over model V/J',
        matching='exact rounded-log2Pgen bins; no missing-bin remapping or dropped source records; counts fixed within bins; event slots without replacement',
        inference='independent per-sample/per-full-pMHC one-sided empirical (1+mock>=real)/(B+1); central95percent mock intervals, no donor/timepoint contrasts or correction',
        automatic_detection_threshold=None,calculation='fresh raw queries/mocks/matches; explicitly supplied immutable Pgen assets only; new query masses remain in-memory',
        stages=stages,seconds=time.monotonic()-start,tables=tables),indent=2)+'\n')
    return 0
