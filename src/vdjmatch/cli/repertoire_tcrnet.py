"""TCRnet-filtered repertoire matching with counts and auditable filter evidence."""
from pathlib import Path


def register(subparsers):
    p=subparsers.add_parser('repertoire-tcrnet',help=__doc__)
    p.add_argument('metadata',help='sample_id/path manifest')
    p.add_argument('--vdjdb',required=True);p.add_argument('--control',required=True)
    p.add_argument('--control-tra',help='TRA control in paired mode; --control supplies TRB')
    p.add_argument('--species',choices=['human','mouse'],default='human')
    p.add_argument('--locus',choices=['TRA','TRB','paired'],default='TRB')
    p.add_argument('--epitope',action='append')
    p.add_argument('--tcrnet-fdr',type=float,default=.1)
    p.add_argument('--scope',default='3,1,1,3',help='flexible matching subs,ins,dels,total; TCRnet itself uses Hamming1')
    p.add_argument('--match-v',action='store_true')
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--output-prefix',required=True)
    p.set_defaults(func=main)


def main(a):
    import json,time
    from importlib.metadata import version
    import polars as pl
    import vdjtools,seqtree._core as seq_native
    from .. import db,io
    from ..db.cache import sha256
    from ..aggregate.repertoire_tcrnet import TcrnetMatcher
    from ..evalue.mock_reference import _pgen_model_identity
    from ..runner.annotation import _input_rows,_validate_output_paths
    from .repertoire_reference import _genes,_read_control,_observations
    start=time.monotonic();stages={}
    rows=_input_rows(None,a.metadata)
    loci=['TRA','TRB'] if a.locus=='paired' else [a.locus]
    if a.locus=='paired' and not a.control_tra:raise ValueError('paired mode requires --control-tra')
    control_paths={locus:Path(a.control_tra if locus=='TRA' and a.locus=='paired' else a.control).expanduser().resolve() for locus in loci}
    reference_path=Path(a.vdjdb).expanduser().resolve()
    inputs=[Path(a.metadata).expanduser().resolve(),reference_path,*[p for _,p,_ in rows],*control_paths.values()]
    prefix=_validate_output_paths(a.output_prefix,inputs)
    suffixes=['summary.tsv','hits.tsv','sample-filter.tsv','reference-filter.tsv','ingestion.tsv','manifest.json']
    if any(Path(str(prefix)+'.'+suffix).exists() for suffix in suffixes):raise FileExistsError('fresh output files must not exist')
    def identities():return {str(p):db.provenance(p) if p.is_dir() else sha256(p) for p in inputs}
    sources=[*Path(__file__).resolve().parents[1].rglob('*.py'),*Path(vdjtools.__file__).parent.rglob('*.py'),Path(seq_native.__file__)]
    def freeze():return {str(p):sha256(p) for p in sources}
    before,source_before=identities(),freeze()
    models={locus:_pgen_model_identity(a.species,locus) for locus in loci}
    reference=db.load(reference_path,species={'human':'HomoSapiens','mouse':'MusMusculus'}[a.species]).filter(pl.col('reference_valid')&pl.col('gene').is_in(loci))
    if a.epitope:
        missing=set(a.epitope)-set(reference['epitope'])
        if missing:raise ValueError(f'requested epitopes absent from reference: {sorted(missing)}')
        reference=reference.filter(pl.col('epitope').is_in(a.epitope))
    if not reference.height:raise ValueError('no eligible reference records')
    reference=pl.concat([_genes(part,locus,species=a.species) for (locus,),part in reference.partition_by('gene',as_dict=True).items()])
    controls={locus:_read_control(path,species=a.species,locus=locus).select('junction_aa').unique() for locus,path in control_paths.items()}
    matcher=TcrnetMatcher(reference,controls,species=a.species,filter_fdr=a.tcrnet_fdr,scope=a.scope,match_v=a.match_v,threads=a.threads)
    stages['reference_and_controls_seconds']=time.monotonic()-start
    print('reference/control filtering complete',flush=True)
    summaries=[];hits=[];filters=[];reports=[]
    for number,(sample,path,_) in enumerate(rows):
        t=time.monotonic();frame,report=io.read_rearrangement(path,valid_aa=False,return_report=True)
        selected,observations=_observations(frame,locus=a.locus,species=a.species)
        summary,hit,filtering=matcher.match(observations)
        for table,container in [(summary,summaries),(hit,hits),(filtering,filters)]:
            container.append(table.with_columns(pl.lit(sample).alias('sample_id')))
        reports.append(dict(sample_id=sample,**report,selected_chain_rows=selected.height,excluded_rows=frame.height-selected.height,
                            imputed_v_rows=int(selected['v_imputed'].sum()),unique_receptors=observations['receptor_id'].n_unique()))
        stages[f'sample_{number+1}_seconds']=time.monotonic()-t
        print(f'sample {number+1}/{len(rows)} complete',flush=True)
    if before!=identities() or source_before!=freeze() or models!={locus:_pgen_model_identity(a.species,locus) for locus in loci}:
        raise ValueError('inputs, gene-imputation models or software changed during run')
    prefix.parent.mkdir(parents=True,exist_ok=True);outputs={}
    for suffix,table in [('summary.tsv',pl.concat(summaries)),('hits.tsv',pl.concat(hits,how='diagonal_relaxed')),
            ('sample-filter.tsv',pl.concat(filters)),('reference-filter.tsv',matcher.reference_filtering),
            ('ingestion.tsv',pl.DataFrame(reports,infer_schema_length=None))]:
        path=Path(str(prefix)+'.'+suffix);table.write_csv(path,separator='\t');outputs[path.name]=sha256(path)
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(dict(method='tcrnet-filtered-flexible-match',
        inputs=before,inputs_frozen=True,source=source_before,source_frozen=True,gene_imputation_models=models,
        software={n:version(n) for n in ['vdjmatch','vdjtools','seqtree','polars','numpy','scipy']},options=vars(a)|{'func':a.func.__name__},
        repertoire_p_values=False,pgen_probabilities_used=False,reference_population='curated-benchmark',
        filter_families='sample/locus and full-pMHC/locus separately; existing TCRnet BH',
        filter_scope='Hamming<=1; target self-hit removed, exact background hits retained',
        counting_unit='unique original receptor; any retained matching chain; raw denominator preserved',
        control_unique_junctions={locus:len(index) for locus,index in matcher.controls.items()},
        stages=stages,seconds=time.monotonic()-start,outputs=outputs),indent=2)+'\n')
    return 0
