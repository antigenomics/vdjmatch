"""Source-reference support for curated exact clonotype/pMHC identities."""
from pathlib import Path


def register(subparsers):
    parser = subparsers.add_parser('reference-support', help='count distinct sources per curated clonotype/pMHC')
    parser.add_argument('--vdjdb', required=True)
    parser.add_argument('--min-refs', type=int, default=2)
    parser.add_argument('--gene-level', action='store_true', help='strip V/J allele suffixes explicitly')
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=1)
    parser.set_defaults(func=main)
    return parser


def main(args):
    import json
    import time
    import polars as pl
    from .. import db
    from ..db.cache import sha256

    started = time.monotonic()
    path = Path(args.output)
    manifest = path.with_suffix(path.suffix+'.json')
    source = Path(args.vdjdb).expanduser().resolve()
    if any(p.resolve()==source or (source.is_dir() and source in p.resolve().parents)
           for p in (path,manifest)):
        raise ValueError('support outputs must not overwrite the reference or enter its input directory')
    provenance = db.provenance(args.vdjdb)
    reference = db.load(args.vdjdb)
    out = db.replicated(reference, min_refs=args.min_refs, gene_level=args.gene_level)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(path, separator='\t')
    manifest.write_text(json.dumps({
        'reference':provenance, 'min_refs':args.min_refs,
        'producer_sha256':sha256(Path(__file__)),
        'reducer_sha256':sha256(Path(db.replicated.__code__.co_filename)),
        'source_records':reference.height,
        'valid_source_records':reference.filter(pl.col('reference_valid')).height,
        'retained_clonotype_pmhc_keys':out.height, 'table_sha256':sha256(path),
        'identity':'full junction/V/J, locus, species and full pMHC', 'gene_level':args.gene_level,
        'interpretation':'distinct source IDs; not guaranteed independent experiments or binding',
        'seconds':time.monotonic()-started,
    }, indent=2)+'\n')
    return 0
