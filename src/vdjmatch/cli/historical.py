"""CLI for reproducible historical density baselines; no fitted combination."""
from pathlib import Path


def register(subparsers):
    p=subparsers.add_parser('historical-density',help='reproduce the historical PSSM density component')
    p.add_argument('sample')
    p.add_argument('--vdjdb',required=True)
    p.add_argument('--locus',choices=['TRA','TRB'],required=True)
    p.add_argument('--epitope',action='append',required=True)
    p.add_argument('--mhc-a',required=True)
    p.add_argument('--species',default='HomoSapiens')
    p.add_argument('--control')
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--output-prefix',required=True)
    p.set_defaults(func=main)


def main(a):
    import json
    import time
    import polars as pl
    from .. import db,io
    from ..api import _prepare
    from ..db.cache import sha256
    from ..evalue.control import raw_background
    from ..match.historical import density
    start=time.perf_counter()
    if a.threads<1:raise ValueError('threads must be positive')
    q,ingestion=io.read_rearrangement(a.sample,source='airr',return_report=True)
    _,q=_prepare(q)
    q=q.filter(pl.col('locus')==a.locus)
    r=db.load(a.vdjdb,species=a.species,gene=a.locus,epitope=a.epitope,mhc_a=[a.mhc_a])
    r=r.filter(pl.col('reference_valid'))
    missing=sorted(set(a.epitope)-set(r['epitope'].to_list()))
    if missing:
        raise ValueError('requested epitopes absent under the selected reference restriction: '+', '.join(missing))
    ctrl,provenance=raw_background(a.locus,a.species,a.control)
    out=density(q,r,ctrl,threads=a.threads)
    prefix=Path(a.output_prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    out.write_csv(str(prefix)+'.scores.tsv',separator='\t')
    Path(str(prefix)+'.manifest.json').write_text(json.dumps({'sample_sha256':sha256(Path(a.sample)),
        'reference':db.provenance(a.vdjdb),'control':provenance,'ingestion':ingestion,
        'parameters':{k:v for k,v in vars(a).items() if k!='func'},
        'scope':'PSSM5 substitutions; unit controls5,2,2; significance nearest<=1; exact excluded',
        'representative':'first source-order V per pMHC junction; baseline only',
        'not_included':['sparse-reference germline prior','paired cohort rank fusion'],
        'wall_seconds':time.perf_counter()-start},indent=2)+'\n')
    return 0
