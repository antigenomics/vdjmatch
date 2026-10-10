"""Native TCRdist3-compatible single-chain neighbourhoods, without confidence claims."""
from pathlib import Path


def register(subparsers):
    p=subparsers.add_parser('tcrdist-neighbours',help='native human TCRdist3-compatible distance balls')
    p.add_argument('sample')
    p.add_argument('--vdjdb',required=True)
    p.add_argument('--locus',choices=['TRA','TRB'],required=True)
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--radius',type=int,default=90)
    p.add_argument('--exclude-exact',action='store_true')
    p.add_argument('--output-prefix',required=True)
    p.set_defaults(func=main)


def main(a):
    import json
    import time
    import numpy as np
    import polars as pl
    from .. import db,io
    from ..api import _prepare
    from ..aggregate.candidates import PMHC
    from ..io.airr import _read_table
    from ..db.cache import sha256
    from ..match.tcrdist import distance_matrix,load_v_loops
    start=time.perf_counter()
    if a.threads<1 or a.radius<0:raise ValueError('positive threads and nonnegative radius required')
    raw=_read_table(a.sample)
    if 'query_id' not in raw.columns:raw=raw.with_row_index('query_id')
    if raw['query_id'].null_count() or raw['query_id'].n_unique()!=raw.height:
        raise ValueError('query_id must be non-null and unique')
    data,ingestion=io.read_rearrangement(a.sample,source='airr',valid_aa=False,return_report=True)
    _,q=_prepare(data)
    r=db.load(a.vdjdb,species='HomoSapiens',gene=a.locus).filter(pl.col('reference_valid'))
    original=r.height
    model=load_v_loops()
    allowed=[v for v,m in model.items() if m['locus']==a.locus]
    def prepare(frame):
        return frame.with_columns(pl.when(pl.col('v').str.contains(r'\*')).then(pl.col('v'))
            .otherwise(pl.col('v')+'*01').alias('_allele')).filter(
                pl.col('_allele').is_in(allowed)&(pl.col('cdr3').str.len_chars()>=8))
    q=prepare(q.filter(pl.col('locus')==a.locus))
    r=prepare(r).unique(subset=['cdr3','_allele',*PMHC],maintain_order=True)
    rows=[]
    if r.height and q.height:
        # The kernel's64MiB output budget determines the batch size; no all-cohort matrix.
        batch=max(1,64*1024**2//(4*r.height))
        rs,rv=r['cdr3'].to_list(),r['_allele'].to_list()
        for part in q.iter_slices(batch):
            seqs=part['cdr3'].to_list()
            distances=distance_matrix(seqs,rs,part['_allele'].to_list(),rv,threads=a.threads)
            keep=distances<=a.radius
            if a.exclude_exact:keep&=np.asarray(seqs)[:,None]!=np.asarray(rs)[None,:]
            qi,ri=np.nonzero(keep)
            h=pl.DataFrame({'_qi':qi,'_ri':ri,'distance':distances[qi,ri]})
            h=h.join(part.with_row_index('_qi').select('_qi','query_id'),on='_qi')
            h=h.join(r.with_row_index('_ri').select('_ri',*PMHC),on='_ri')
            rows.append(h.group_by('query_id',*PMHC).agg(pl.col('distance').min(),pl.len().alias('n_neighbours')))
    schema={'query_id':raw.schema['query_id'],**{c:pl.String for c in PMHC},'distance':pl.Int32,'n_neighbours':pl.UInt32}
    out=pl.concat(rows) if rows else pl.DataFrame(schema=schema)
    out=out.sort('query_id','distance',*PMHC,nulls_last=True)
    status=raw.select([c for c in ('query_id','sequence_id') if c in raw.columns]).join(q.select('query_id').with_columns(pl.lit(True).alias('available')),on='query_id',how='left')
    status=status.with_columns(pl.col('available').fill_null(False))
    prefix=Path(a.output_prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    out.write_csv(str(prefix)+'.candidates.tsv',separator='\t')
    status.write_csv(str(prefix)+'.queries.tsv',separator='\t')
    Path(str(prefix)+'.manifest.json').write_text(json.dumps({'sample_sha256':sha256(Path(a.sample)),
        'reference':db.provenance(a.vdjdb),'ingestion':ingestion,'reference_observations':original,
        'reference_usable_keys':r.height,'query_rows':raw.height,'query_usable_rows':q.height,
        'parameters':{k:v for k,v in vars(a).items() if k!='func'},
        'distance':'TCRdist3 default3*CDR3+CDR1+CDR2+CDR2.5; human combo_xcr_2024-03-05',
        'calibration':'none','wall_seconds':time.perf_counter()-start},indent=2)+'\n')
    return 0
