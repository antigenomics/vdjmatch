"""Reproduce accepted legacy matrix and epitope-withheld retrieval tables."""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('scoring-audit', help='reproduce legacy matrix/LOO tables from explicit raw VDJdb inputs')
    p.add_argument('--task', choices=['matrix', 'loo'], required=True)
    p.add_argument('--pmhc', required=True, help='raw VDJdb export; legacy cdr3 denotes anchored junction')
    p.add_argument('--species', choices=['HomoSapiens', 'MusMusculus'], default='HomoSapiens')
    p.add_argument('--locus', choices=['TRA', 'TRB'], required=True)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--seed', type=int, default=0, help='legacy composition-cap seed (LOO only)')
    p.add_argument('--min-epi', type=int, default=30, help='legacy unique-junction eligibility; not the full-VDJdb source-record gate')
    p.add_argument('--vlen', type=int, default=4, help='matrix task V prefix length')
    p.add_argument('--jlen', type=int, default=6, help='matrix task J suffix length')
    p.add_argument('--legacy', help='explicit legacy similarity TSV; default packaged vdjam.txt')
    p.add_argument('--retention', help='explicit V/J retention model TSV; default packaged human model')
    p.add_argument('--significance', help='explicit positional significance model TSV; default packaged model')
    p.add_argument('--structural-inc', help='explicit seqtree structural.inc; omit to skip structural row, never alias BLOSUM')
    p.add_argument('--top', type=int, default=8, help='LOO held-out epitopes by descending unique-junction count')
    p.add_argument('--max-queries', type=int, default=600, help='LOO lexicographically first unique queries per epitope')
    p.add_argument('--subs', type=int, default=2, help='LOO substitution-only candidate radius')
    p.add_argument('--output-prefix', required=True)
    p.set_defaults(func=main)
    return p


def main(a):
    import json
    import math
    import time
    from importlib import metadata, resources
    import polars as pl
    from .. import db
    from ..match import regions, scoring_audit as legacy
    from ..db.cache import sha256

    for name in ('threads', 'min_epi', 'top', 'max_queries'):
        if getattr(a, name) < 1:
            raise ValueError(f'{name.replace("_", "-")} must be positive')
    if min(a.vlen, a.jlen, a.subs) < 0:
        raise ValueError('region lengths and substitution radius must be nonnegative')
    if a.task == 'loo' and a.species == 'MusMusculus' and not a.retention:
        raise ValueError('mouse LOO requires an explicit species-matched --retention model')
    start = time.perf_counter()
    root = resources.files('vdjmatch.resources')
    paths = {'legacy': Path(a.legacy or root / 'vdjam.txt')} if a.task == 'matrix' else {}
    if a.task == 'loo':
        paths.update(retention=Path(a.retention or root / 'trimming' / 'human_vj_retention.tsv'),
                     significance=Path(a.significance or root / 'trimming' / 'position_significance.tsv'))
        if a.structural_inc:
            paths['structural'] = Path(a.structural_inc)
    models = {name: {'path': str(path.resolve()), 'sha256': sha256(path)} for name, path in paths.items()}
    raw = db.load(a.pmhc, species=a.species, gene=a.locus)
    # The legacy estimators use standard-residue sanitation, not modern curation.
    valid = legacy.valid_cdr3(raw)
    prefix = Path(a.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if a.task == 'matrix':
        matrices, correlations = legacy.matrix_audit(raw, min_epi=a.min_epi, vlen=a.vlen,
            jlen=a.jlen, legacy=legacy.legacy_vdjam_offdiag(paths['legacy']), threads=a.threads)
        tables = {'matrix': matrices, 'correlations': correlations}
    else:
        structural = None
        if a.structural_inc:
            structural = legacy.structural_dissim(a.structural_inc)
            if structural is None or len(structural) != len(legacy.AA)*(len(legacy.AA)-1):
                raise ValueError('structural-inc does not supply the complete standard-AA structural matrix')
        retention = regions.load_retention(paths['retention'])
        significance = regions.load_significance(paths['significance'])
        if not retention or not all(significance):
            raise ValueError('retention and significance models must be nonempty')
        if any(not math.isfinite(p) or not 0 <= p <= 1
               for profile in [*retention.values(), *significance] for p in profile):
            raise ValueError('retention/significance probabilities must be finite and in [0,1]')
        tasks, summary, pairs, queries, refs = legacy.loo_audit(raw, locus=a.locus, retention=retention,
            significance=significance, structural=structural, min_epi=a.min_epi, top=a.top,
            max_queries=a.max_queries, subs=a.subs, seed=a.seed, threads=a.threads)
        tables = {'tasks': tasks, 'summary': summary, 'pairs': pairs, 'queries': queries,
                  'references': refs.with_row_index('reference_key_id')}
    for name, table in tables.items():
        table.with_columns(pl.lit(a.species).alias('species'), pl.lit(a.locus).alias('locus')).write_csv(
            str(prefix)+'.'+name+'.tsv', separator='\t')
    manifest = {
        'reference': db.provenance(a.pmhc), 'models': models,
        'options': vars(a) | {'func': 'scoring-audit'},
        'versions': {name: metadata.version(name) for name in ('vdjmatch', 'seqtree', 'polars', 'numpy')},
        'source_chain_records': raw.height, 'valid_aa_chain_records': valid.height,
        'invalid_aa_chain_records': raw.height-valid.height,
        'contract': 'legacy accepted evidence; not current annotation or full-VDJdb benchmark eligibility',
        'sequence_convention': 'producer-normalized anchored junction; no fabricated anchors',
        'matrix_estimator': 'same-epitope unordered Hamming-1 pairs; position background; alpha=1 log2 odds',
        'loo_contract': 'spike-study records removed; unique(cdr3,v,j,epitope,mhc_class); cap=3000; seeded canonical input; first epitope per junction; exact sequence pairs excluded; NDN training positions [4,L-6)',
        'metric': 'legacy tied-score-group trapezoidal balanced PR-AUC at reference prevalence 0.5; not average precision',
        'structural_status': 'explicit input' if a.structural_inc else 'omitted; no BLOSUM fallback',
        'output_rows': {name: table.height for name, table in tables.items()},
        'seconds': time.perf_counter()-start,
    }
    Path(str(prefix)+'.provenance.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0
