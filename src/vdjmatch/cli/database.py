"""Whole-database junction presence and an independent nearest-match null tail."""
from pathlib import Path


def register(subparsers):
    p=subparsers.add_parser('database-evidence',help=__doc__)
    p.add_argument('sample')
    p.add_argument('--vdjdb',required=True)
    p.add_argument('--null-sample',action='append',required=True,help='matching raw null junction draws; repetitions and no-hits remain in the denominator')
    p.add_argument('--species',default='human')
    p.add_argument('--locus',choices=['TRA','TRB'],required=True)
    p.add_argument('--max-edits',type=int,default=5,help='maximum full-junction Levenshtein distance; applied identically to queries and nulls')
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--input-format',choices=['auto','airr','legacy','custom'],default='auto')
    p.add_argument('--sequence-convention',choices=['junction'])
    p.add_argument('--output-prefix',required=True)
    p.set_defaults(func=main)


def prepare(path,a):
    """Preserve every source row; normalization never fabricates anchors."""
    import polars as pl
    from ..io.airr import normalize_annotation_input
    from ..evalue.control import _organism
    raw,report=normalize_annotation_input(path,source=a.input_format,sequence_convention=a.sequence_convention)
    if any(c.startswith('database_') for c in raw.columns):
        raise ValueError('input contains reserved database_ diagnostic columns')
    if any(c.startswith('_database_') for c in raw.columns):
        raise ValueError('input contains reserved _database_ diagnostic columns')
    raw=raw.with_row_index('_database_row').with_columns(
        pl.col('junction_aa').str.strip_chars().str.to_uppercase().alias('database_junction'))
    species_ok=pl.lit(True)
    if 'species' in raw.columns:
        aliases={}
        for value in raw['species'].drop_nulls().unique():
            try: aliases[value]=_organism(value)==_organism(a.species)
            except ValueError: aliases[value]=False
        species_ok=pl.col('species').replace_strict(aliases,default=False,return_dtype=pl.Boolean)
    locus_ok=pl.col('locus').is_null() | pl.col('locus').eq(a.locus)
    seq=pl.col('database_junction')
    raw=raw.with_columns(pl.when(seq.is_null() | seq.eq('')).then(pl.lit('missing_sequence'))
        .when(~seq.str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$')).then(pl.lit('invalid_sequence'))
        .when(~locus_ok).then(pl.lit('locus_mismatch'))
        .when(~species_ok).then(pl.lit('species_mismatch')).otherwise(pl.lit('ready')).alias('database_status'))
    return raw,report


def exact_clonotypes(query,reference,locus):
    """Literal junction/V-gene/J-gene lookup; ambiguous calls stay unavailable."""
    import polars as pl
    gene=lambda c: pl.col(c).cast(pl.String).str.strip_chars().str.to_uppercase().str.split('*').list.first()
    single=lambda c,kind: pl.col(c).str.strip_chars().str.to_uppercase().str.contains(
        '^'+locus+kind+r'[A-Z0-9/.-]+(?:\*[A-Z0-9]+)?$').fill_null(False) & pl.col(c).str.to_uppercase().str.count_matches(r'TR[A-Z][VJ]').eq(1)
    keys=reference.filter(single('v','V') & single('j','J')).select(
        pl.col('cdr3').alias('database_junction'),gene('v').alias('_database_v'),gene('j').alias('_database_j')).unique()
    query=query.with_columns(gene('v_call').alias('_database_v'),gene('j_call').alias('_database_j'))
    available=pl.col('database_status').is_in(['ready','matched','no_hit'])
    for col,kind in [('v_call','V'),('j_call','J')]:
        available=available & single(col,kind)
    query=query.join(keys.with_columns(pl.lit(True).alias('_database_present')),
        on=['database_junction','_database_v','_database_j'],how='left',validate='m:1',maintain_order='left')
    return query.with_columns(pl.when(available).then(pl.col('_database_present').fill_null(False)).otherwise(None)
        .alias('database_exact_gene_clonotype')).drop('_database_v','_database_j','_database_present')


def main(a):
    import json
    import time
    from importlib.metadata import version
    import polars as pl
    from seqtree import Index
    import seqtree._core as native
    from .. import db
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..match.scope import search_params
    from ..runner.annotation import _validate_output_paths
    from .calibration import null_tail_probabilities

    start=time.monotonic()
    if a.threads<1: raise ValueError('positive threads required')
    if a.max_edits<0: raise ValueError('max-edits must be nonnegative')
    if len({Path(p).resolve() for p in a.null_sample})!=len(a.null_sample):
        raise ValueError('null sources must be distinct files')
    if Path(a.sample).resolve() in {Path(p).resolve() for p in a.null_sample}:
        raise ValueError('query sample cannot also define its null')
    prefix=_validate_output_paths(a.output_prefix,[a.sample,a.vdjdb,*a.null_sample])
    species={'human':'HomoSapiens','mouse':'MusMusculus'}[_organism(a.species)]
    reference=db.load(a.vdjdb,species=species,gene=a.locus)
    source_records=reference.height
    # Ordinary presence searches all valid sequences, without benchmark curation,
    # study-support cutoffs, epitope selection or exact-query exclusion.
    reference=reference.filter(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False))
    targets=reference['cdr3'].unique().sort().to_list()
    if not targets: raise ValueError('reference has no valid junctions')
    query,ingestion=prepare(a.sample,a)
    nulls=[];null_reports=[]
    for i,path in enumerate(a.null_sample):
        frame,report=prepare(path,a)
        if not frame.height or frame['database_status'].ne('ready').any():
            raise ValueError('null samples must contain only valid matching-species/locus junction rows; no rows may be silently removed')
        nulls.append(frame.select('database_junction').with_row_index('draw').with_columns(
            pl.lit(f'null_{i+1}').alias('null_source')))
        null_reports.append(dict(source=f'null_{i+1}',sha256=sha256(Path(path)),draws=frame.height,ingestion=report))
    null=pl.concat(nulls).with_row_index('query_id')
    seqs=sorted(set(query.filter(pl.col('database_status')=='ready')['database_junction']) | set(null['database_junction']))
    index=Index.build(targets,'aa')
    params=search_params(','.join([str(a.max_edits)]*4),engine='seqtrie',mode='top')
    hits=index.search_batch(seqs,params,a.threads)
    distance={seq:(h[0].score if h else None) for seq,h in zip(seqs,hits)}
    lookup=pl.DataFrame({'database_junction':seqs,'database_nearest_edits':pl.Series([distance[s] for s in seqs],dtype=pl.Int64)})
    query=query.join(lookup,on='database_junction',how='left',validate='m:1',maintain_order='left').with_columns(
        pl.when(pl.col('database_status')=='ready').then(pl.col('database_nearest_edits')).otherwise(None).alias('database_nearest_edits'))
    null=null.join(lookup,on='database_junction',how='left',validate='m:1',maintain_order='left').with_columns(
        (-pl.col('database_nearest_edits')).alias('score'),pl.lit('null').alias('outcome'))
    scored,tails=null_tail_probabilities(query.filter(pl.col('database_status')=='ready').select(
        pl.col('_database_row').alias('query_id'),(-pl.col('database_nearest_edits')).alias('score')),null)
    query=query.join(scored.select(pl.col('query_id').alias('_database_row'),pl.col('p_null').alias('database_p')),
        on='_database_row',how='left',validate='1:1',maintain_order='left')
    query=query.with_columns(pl.when(pl.col('database_status')=='ready').then(pl.col('database_junction').is_in(targets)).otherwise(None).alias('database_exact_junction'),
        pl.when(pl.col('database_status')=='ready').then(pl.when(pl.col('database_nearest_edits').is_not_null())
        .then(pl.lit('matched')).otherwise(pl.lit('no_hit'))).otherwise(pl.col('database_status')).alias('database_status'))
    query=exact_clonotypes(query,reference,a.locus)
    tails=tails.join(query.select(pl.col('_database_row').alias('query_id'),
        pl.col('database_nearest_edits').is_not_null().alias('_distance_available')),
        on='query_id',how='left',validate='m:1',maintain_order='left').with_columns(
        pl.when(pl.col('_distance_available')).then(pl.col('null_tail_count')).otherwise(None).alias('null_tail_count')
        ).drop('_distance_available','E_null')
    prefix.parent.mkdir(parents=True,exist_ok=True)
    paths=[Path(str(prefix)+suffix) for suffix in ('.database.tsv','.null.tsv','.null-tails.tsv')]
    query.rename({'_database_row':'database_input_row'}).write_csv(paths[0],separator='\t')
    null.drop('score','outcome').write_csv(paths[1],separator='\t')
    tails.rename({'query_id':'database_input_row','p_null':'database_p'}).write_csv(paths[2],separator='\t')
    manifest=dict(calculation='database-nearest-junction-null-tail-v1',reference=db.provenance(a.vdjdb),
        sample_sha256=sha256(Path(a.sample)),null_sources=null_reports,species=species,locus=a.locus,max_edits=a.max_edits,threads=a.threads,
        search='one native batch; nearest full-junction unit edit distance to any valid reference; exact matches retained',
        reference_source_records=source_records,valid_reference_records=reference.height,unique_reference_junctions=len(targets),query_rows=query.height,
        cpu_budget=dict(native_workers=a.threads,table_threads=pl.thread_pool_size(),blas_threads_environment=1,outer_pools=0),
        query_status_counts=dict(query.group_by('database_status').len().iter_rows()),unique_searched_junctions=len(seqs),native_batches=1,
        statistic='negative nearest full-junction edit distance; no-hit scores unavailable',
        null_tail='max source (1 + null draws at least as close)/(1 + all null draws); ties and no-hits retained; no-hit query P=1',
        interpretation='junction-level database resemblance under the supplied null; exact presence is a lookup, not a probability; no epitope or paired-binding claim',
        exact_gene_clonotype_policy='junction plus single supplied V/J gene calls, alleles pooled; ambiguous reference calls never choose a first gene',
        assumptions='fixed reference and identical query/null procedure; query and null draws exchangeable under the declared null',
        ingestion=ingestion,software={p:version(p) for p in ('vdjmatch','seqtree','polars')},
        source_sha256=sha256(Path(__file__)),native_sha256=sha256(Path(native.__file__)),
        tables={p.name:sha256(p) for p in paths},seconds=time.monotonic()-start)
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return 0
