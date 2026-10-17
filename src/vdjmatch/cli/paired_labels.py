"""Assign linked-first labels without changing kernel evidence or any-call detection."""
from pathlib import Path


def register(subparsers):
    p=subparsers.add_parser('paired-labels',help=__doc__)
    p.add_argument('evidence',help='native wide paired-evidence TSV, optionally gzip compressed')
    p.add_argument('--selection',required=True,help='matching native paired selection TSV')
    p.add_argument('--input-manifest',required=True,help='native paired all-target kernel manifest')
    p.add_argument('--output-prefix',required=True)
    p.add_argument('--threads',type=int,default=4)
    p.set_defaults(func=main)


def main(a):
    import json
    import time
    import polars as pl
    from ..aggregate.paired_labels import AXES,paired_labels
    from ..db.cache import sha256
    from ..runner.annotation import _validate_output_paths

    started=time.monotonic()
    paths=[Path(a.evidence),Path(a.selection),Path(a.input_manifest)]
    prefix=_validate_output_paths(a.output_prefix,paths)
    identities={str(p):sha256(p) for p in paths}
    source_paths=[Path(__file__),Path(paired_labels.__code__.co_filename),
        *[Path(__file__).with_name(n+'.py') for n in ('neighbours','reference_calls','calibration')],
        Path(__file__).parent.parent/'evalue/kernel.py',Path(__file__).parent.parent/'match/geometry.py']
    sources={str(p):sha256(p) for p in source_paths}
    native=json.loads(paths[2].read_text());options=native.get('parameters',{})
    required={'locus':'paired','all_targets':True,'background':'vdjdb-other','neighbour_weighting':'kernel',
        'stem_depth':4,'kernel_scorer':'v-bounded'}
    if any(options.get(k)!=v for k,v in required.items()):
        raise ValueError('paired-labels requires native all-target paired stem4 v-bounded VDJdb-other kernel evidence')
    species={'human':'HomoSapiens','mouse':'MusMusculus','HomoSapiens':'HomoSapiens','MusMusculus':'MusMusculus'}.get(options.get('species'))
    if species is None: raise ValueError('native manifest requires explicit human or mouse species')
    # Preserve lexical identities and metadata exactly, including leading zeros.
    selection=pl.read_csv(paths[1],separator='\t',infer_schema=False)
    if selection['query_id'].null_count() or selection['query_id'].n_unique()!=selection.height:
        raise ValueError('native selection requires unique non-null query IDs')
    maxima=selection.select('query_id',pl.lit(True).alias('_known'),
        *[pl.col(maximum).cast(pl.Float64).alias('_max_'+axis) for axis,_,maximum,_,_ in AXES])
    wide=pl.scan_csv(paths[0],separator='\t',schema_overrides={'query_id':pl.String,'epitope':pl.String},
        infer_schema_length=10000).select('query_id','epitope',*[column for _,column,_,_,_ in AXES])
    joined=wide.join(maxima.lazy(),on='query_id',how='left')
    # Valid exports retain only exact native maxima. Include values above a
    # declared maximum to reject score drift rather than silently discard it.
    retain=pl.col('_known').is_null()
    for axis,column,_,_,_ in AXES:
        value=pl.col(column).cast(pl.Float64);maximum=pl.col('_max_'+axis)
        retain=retain | (value.is_finite() & (value>0) &
            ((value>=maximum) | maximum.is_null() | ~maximum.is_finite() | (maximum<=0)))
    compact=joined.filter(retain).select('query_id','epitope',*[column for _,column,_,_,_ in AXES]).collect(engine='streaming')
    calls,scores=paired_labels(compact,selection,species=species)
    if identities!={str(p):sha256(p) for p in paths} or sources!={str(p):sha256(p) for p in source_paths}:
        raise ValueError('paired label inputs or source changed during reduction')
    prefix.parent.mkdir(parents=True,exist_ok=True)
    outputs={}
    for kind,frame in [('calls',calls),('scores',scores)]:
        path=Path(str(prefix)+'.'+kind+'.tsv');frame.write_csv(path,separator='\t');outputs[path.name]=sha256(path)
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(dict(inputs=identities,source_sha256=sources,tables=outputs,
        species=species,definition='positive finite linked axis: all linked maximum peptides; otherwise union of all alpha/beta maximum peptides; no supported axis: explicit no-support',
        scorer='stem4_linked_first',score_kind='selected decision indicator only: 1 for selected labels, null for no support; not continuous AP/ROC evidence or probability',
        detection_changed=False,calibration_transfer='none; native any-call scalar and bindings unchanged; label policy has no calibrated confidence',
        queries=selection.height,label_rows=calls.height,seconds=time.monotonic()-started),indent=2)+'\n')
    return 0
