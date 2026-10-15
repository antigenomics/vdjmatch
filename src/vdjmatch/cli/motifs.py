"""Raw-input, source-curated V+K-mer evidence for declared assay targets."""
from pathlib import Path


def register(subparsers):
    p=subparsers.add_parser('motif-evidence',help='experimental masked-motif V+K-mer evidence for sample-declared pMHC')
    p.add_argument('sample')
    p.add_argument('--vdjdb',required=True)
    p.add_argument('--locus',choices=['TRA','TRB'],required=True)
    p.add_argument('--species',default='HomoSapiens')
    p.add_argument('--exclude-exact',action='store_true',help='benchmark-only exact-hit filtering; motif models remain fixed')
    p.add_argument('--output-prefix',required=True)
    p.add_argument('--threads',type=int,default=1,help='positive table thread budget (default: 1); BLAS stays serial')
    p.set_defaults(func=main)


def main(a):
    import json
    import polars as pl
    from .. import db,io
    from ..api import _LOCUS
    from ..aggregate.candidates import PMHC
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..io.airr import _read_table
    from ..match.motifs import masked_motif_scores
    from .historical import assayed_targets
    from ..runner.annotation import _validate_output_paths
    prefix=_validate_output_paths(a.output_prefix, [a.sample,a.vdjdb])
    raw=_read_table(a.sample)
    if 'query_id' not in raw.columns:raw=raw.with_row_index('query_id')
    queries,ingestion=io.read_rearrangement(a.sample,source='airr',valid_aa=False,return_report=True)
    aliases=io.columns._resolve(raw)
    locus=pl.col(aliases['locus']) if 'locus' in aliases else pl.col(aliases['v']).str.slice(0,3)
    base=raw.select('query_id',locus.cast(pl.String).str.to_uppercase().replace_strict(_LOCUS,default=None).alias('locus'),
        (pl.col(aliases['pair_id']).cast(pl.String) if 'pair_id' in aliases else pl.lit(None,dtype=pl.String)).alias('pair_id'))
    queries=base.join(queries.drop('locus','pair_id'),on='query_id',how='left',validate='1:1')
    queries=queries.filter(pl.col('locus')==a.locus)
    organism=_organism(a.species)
    ref=db.load(a.vdjdb,species={'human':'HomoSapiens','mouse':'MusMusculus'}[organism],gene=a.locus)
    ref=ref.filter(pl.col('reference_valid'))
    parts=[];targets=[]
    for task,part,selected in assayed_targets(queries,raw,ref,a.species):
        part=part.with_columns(pl.lit(selected.height>0).alias('_has_reference'))
        parts.append(part)
        if selected.height:
            targets.append(selected.select('cdr3','v').with_columns(
                *[pl.lit(value,dtype=pl.String).alias(key) for key,value in task.items()]))
    if not parts:raise ValueError('no queries for selected species/locus')
    q=pl.concat(parts,how='diagonal_relaxed')
    valid=pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)&pl.col('v').str.strip_chars().ne('').fill_null(False)
    q=q.with_columns(valid.alias('_valid'))
    if targets:
        keys=['cdr3','v',*PMHC]
        representatives=q.filter(pl.col('_valid')).unique(keys,maintain_order=True)
        evidence=masked_motif_scores(representatives,ref,exclude_exact=a.exclude_exact,
                                   target_reference=pl.concat(targets))
        q=q.join(representatives.select(*keys,pl.col('query_id').alias('_motif_query')),
            on=keys,how='left',validate='m:1',nulls_equal=True).join(
            evidence.rename({'query_id':'_motif_query'}),on=['_motif_query',*PMHC],
            how='left',validate='m:1',nulls_equal=True).drop('_motif_query')
    else:q=q.with_columns(pl.lit(None,dtype=pl.Float64).alias('motif_score'),pl.lit(None,dtype=pl.UInt32).alias('supported_motifs'))
    q=q.with_columns(pl.when(~pl.col('_valid')).then(pl.lit('invalid_query'))
        .when(~pl.col('_has_reference')).then(pl.lit('no_reference'))
        .when(pl.col('supported_motifs').is_null()).then(pl.lit('no_supported_motif'))
        .otherwise(pl.lit('scored')).alias('status'),pl.col('supported_motifs').fill_null(0))
    q=q.with_columns(pl.when(pl.col('status').is_in(['scored','no_supported_motif']))
        .then(pl.col('motif_score').fill_null(0)).otherwise(None).alias('motif_score')).drop('_valid','_has_reference')
    prefix.parent.mkdir(parents=True,exist_ok=True)
    q.write_csv(str(prefix)+'.scores.tsv',separator='\t')
    provenance=db.provenance(a.vdjdb)
    report={'sample_sha256':sha256(Path(a.sample)),'reference':provenance,
        'reference_sha256':provenance.get('input_sha256'),
        'species':organism,'locus':a.locus,'exclude_exact':a.exclude_exact,'queries':q.height,
        'ingestion':ingestion,'counts':q.group_by('status').len().to_dicts(),
        'parameters':{'k':4,'masks':[0,1,2,3],'minimum_length_v_target':2,'minimum_v_target':3,'minimum_target':5},
        'target_selection':'compatible full-pMHC reference union; global frequencies from all valid same-species/locus reference',
        'interpretation':'MF times IKF ranking evidence; no switched score, P-value or binding probability'}
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(report,indent=2,default=str)+'\n')
