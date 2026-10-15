"""Naive unit-edit fixed-scope comparison, without gene or background weighting."""
from pathlib import Path


def register(subparsers, *, trained=False):
    name = 'trained-score' if trained else 'fixed-scope'
    description='frozen manuscript same-antigen model and hit informativeness' if trained else __doc__
    p=subparsers.add_parser(name,help=description,description=description)
    p.add_argument('sample',help='raw junction table with declared assayed full-pMHC columns')
    p.add_argument('--vdjdb',required=True)
    p.add_argument('--species',default='human')
    p.add_argument('--paired',action='store_true',help=('retain declared pairing and trained marginal scores; no trained joint model' if trained else 'use declared pair_id; retain incomplete pairs and both marginals'))
    p.add_argument('--scope',default='4,2,2,6' if trained else '1,0,0,1',
        help='substitution,insertion,deletion,total operation caps')
    if trained:
        p.add_argument('--max-indels',type=int,default=2,help='aggregate insertion plus deletion budget (default: 2)')
    p.add_argument('--exclude-exact',action='store_true')
    p.add_argument('--reference-population',choices=['all-valid','curated-benchmark'],default='all-valid')
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--output-prefix',required=True)
    p.set_defaults(func=main,score_model='trained' if trained else 'naive')


def main(a):
    import json,time
    import polars as pl
    from importlib.metadata import version
    from .. import db,io
    from ..aggregate.candidates import PMHC
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..io.airr import _read_table
    from ..io.columns import _resolve
    from ..runner.annotation import _validate_output_paths
    from ..runner.fixed_scope import fixed_scope_scores
    from ..match.scope import parse_scope
    start=time.monotonic();limits=parse_scope(a.scope)
    trained=getattr(a,'score_model','naive') == 'trained'
    if a.threads<1:raise ValueError('positive threads required')
    organism=_organism(a.species)
    prefix=_validate_output_paths(a.output_prefix,[a.sample,a.vdjdb])
    raw=_read_table(a.sample)
    if trained and any(
            name.startswith(('classical_','_classical_')) for name in raw.columns):
        raise ValueError('query fields collide with reserved classical scorer outputs')
    if not trained and 'species' in raw.columns and any(_organism(value)!=organism for value in raw['species'].unique()):
        raise ValueError('query species must match the declared species; use separate population invocations')
    if 'query_id' not in raw.columns:raw=raw.with_row_index('query_id')
    queries,ingestion=io.read_rearrangement(a.sample,valid_aa=False,return_report=True)
    queries=queries.with_columns(pl.col('query_id').cast(pl.String))
    raw=raw.with_columns(pl.col('query_id').cast(pl.String))
    if missing:=set(PMHC)-set(raw.columns):raise ValueError('declared full pMHC required: '+', '.join(sorted(missing)))
    if queries.select('query_id').join(raw.select('query_id'),on='query_id',how='anti').height:
        raise ValueError('ingester query identities do not match source rows')
    fields=_resolve(raw)
    base=raw.select('query_id',*PMHC,*(['species'] if trained and 'species' in raw.columns else []),*[pl.col(fields[c]).cast(pl.String).alias('_raw_'+c) if c in fields else pl.lit(None,dtype=pl.String).alias('_raw_'+c) for c in ('v','j','locus','pair_id')])
    queries=base.join(queries,on='query_id',how='left',validate='1:1',maintain_order='left').with_columns(
        *[pl.coalesce(c,'_raw_'+c).alias(c) for c in ('v','j','pair_id')],
        pl.coalesce('locus',pl.col('_raw_locus').str.to_uppercase(),pl.col('_raw_v').str.slice(0,3).str.to_uppercase()).alias('locus')).drop(*['_raw_'+c for c in ('v','j','locus','pair_id')])
    reference=db.load(a.vdjdb,species={'human':'HomoSapiens','mouse':'MusMusculus'}[organism])
    valid=pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
    reference=reference.filter(pl.col('reference_valid') if a.reference_population=='curated-benchmark' else valid)
    # Pair validation uses the chosen population, including ordinary valid exacts.
    reference=reference.with_columns(valid.alias('reference_valid'))
    if trained:
        from ..runner.classical import classical_scores
        import hashlib
        import seqtree._core as native
        scores,statistics=classical_scores(queries,reference,species=organism,paired=a.paired,scope=a.scope,
            max_indels=a.max_indels,exclude_exact=a.exclude_exact,threads=a.threads)
        extras=[c for c in raw.columns if c not in scores.columns]
        if extras:
            scores=scores.join(raw.select('query_id',*extras),on='query_id',how='left',validate='1:1',maintain_order='left')
        prefix.parent.mkdir(parents=True,exist_ok=True)
        scores.write_csv(str(prefix)+'.scores.tsv',separator='\t')
        source=db.provenance(a.vdjdb)
        reference_hash=source.get('input_sha256')
        if reference_hash is None:
            digest=hashlib.sha256()
            for name,value in sorted(source['files'].items()):
                digest.update(name.encode()+b'\0'+bytes.fromhex(value))
            reference_hash=digest.hexdigest()
        common=dict(sample_sha256=sha256(Path(a.sample)),reference_sha256=reference_hash,
            reference_hash_kind='file_sha256' if Path(a.vdjdb).is_file() else 'relative_name_and_file_sha256_digest',
            source_query_rows=raw.height,output_chain_rows=scores.height,reference_population=a.reference_population,
            scope=a.scope,max_indels=a.max_indels,exclude_exact=a.exclude_exact,native_threads=a.threads,
            seqtree_native_sha256=sha256(Path(native.__file__)),
            scorer_source_sha256=sha256(Path(classical_scores.__code__.co_filename)),
            vdjmatch_version=version('vdjmatch'),seqtree_version=version('seqtree'),
            interpretation='frozen manuscript ranking; no calibrated match P-value or binding posterior',
            paired_rule='separate trained marginals; no trained joint model',whole_seconds=time.monotonic()-start)
        populations=statistics['populations'] or [dict(status='no_scorable_population')]
        pl.DataFrame([common|population for population in populations],infer_schema_length=None).write_csv(
            str(prefix)+'.provenance.tsv',separator='\t')
        return 0
    scores,statistics=fixed_scope_scores(queries,reference,species=organism,paired=a.paired,scope=a.scope,
        exclude_exact=a.exclude_exact,threads=a.threads)
    prefix.parent.mkdir(parents=True,exist_ok=True);scores.write_csv(str(prefix)+'.scores.tsv',separator='\t')
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(dict(reference=db.provenance(a.vdjdb),
        sample_sha256=sha256(Path(a.sample)),reference_population=a.reference_population,source_query_rows=raw.height,
        output_cases=scores.height,scope=a.scope,operation_budgets=dict(zip(['substitutions','insertions','deletions','total'],limits)),
        matrix='none',VJ_weighting='none',background_weighting='none',prior='none',exclude_exact=a.exclude_exact,
        sequence_contract='unchanged full junction including conserved anchors',
        count_units='marginal unique junctions and joint junction/V/J keys; linked unique junction pairs and joint paired V/J keys; source observations counted separately',
        paired_reduction='linked and both marginal rankings separate; any_hit is their support OR, with no availability veto',
        nearest='marginal minimum unit edits; linked minimum max(alpha_edits,beta_edits); lower distance is stronger',
        nohits='all cases retained, nearest distance null and neighbor counts zero',
        interpretation='naive fixed-scope rankings; no binding posterior or calibrated P-value',
        native_allocation_guard='native list-of-lists is materialized before hit-count validation; outer process-tree RSS guard required for wide scopes',
        statistics=statistics,ingestion=ingestion,software={p:version(p) for p in ('vdjmatch','seqtree','polars')},
        source_sha256=sha256(Path(fixed_scope_scores.__code__.co_filename)),seconds=time.monotonic()-start),indent=2)+'\n')
    return 0
