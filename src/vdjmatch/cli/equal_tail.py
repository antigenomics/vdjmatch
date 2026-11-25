"""Experimental declared-target ranking using distance and kernel reference tails."""
from . import cpu_budget


def register(subparsers):
    p = subparsers.add_parser('equal-tail', help=__doc__, description=(
        'Rank receptor evidence for explicitly supplied candidate pMHCs. The input must '
        'provide epitope, mhc_a, mhc_b and mhc_class for every query. Use distinct query '
        'or pair IDs when evaluating the same receptor against multiple candidates.'))
    p.add_argument('sample', help='junction-bearing receptor table; linked TRA/TRB rows require a pair or clone ID')
    p.add_argument('--vdjdb', required=True, help='local raw VDJdb release, directory or table')
    p.add_argument('--locus', required=True, choices=['TRA', 'TRB', 'paired'])
    p.add_argument('--species', default='human', help='human or mouse')
    p.add_argument('--background-scope', choices=['all-hla', 'compatible-hla'], default='all-hla',
                   help='other-peptide calibration keys from all restrictions or compatible pMHCs only')
    p.add_argument('--exclude-exact', action='store_true', help='exclude exact full-junction hits from both components')
    p.add_argument('--threads', type=cpu_budget, default='0')
    p.add_argument('--output-prefix', required=True)
    p.set_defaults(func=main)


def _queries(path, locus, species):
    import polars as pl
    from ..aggregate.candidates import PMHC
    from ..db.schema import normalize_mhc_class
    from ..evalue.control import _organism
    from ..io.airr import normalize_annotation_input
    raw, report = normalize_annotation_input(path)
    missing = [c for c in PMHC if c not in raw]
    if missing:
        raise ValueError('declared candidate pMHC requires columns: ' + ', '.join(missing))
    if 'species' not in raw:
        raw = raw.with_columns(pl.lit(species).alias('species'))
    species_map = {v:_organism(v) for v in raw['species'].drop_nulls().unique()}
    raw = raw.with_columns(pl.col('species').replace_strict(species_map,default=None).alias('species'),
                           pl.col('locus').str.strip_chars().str.to_uppercase().replace({'ALPHA':'TRA','BETA':'TRB'}),
                           pl.col('junction_aa').str.strip_chars().str.to_uppercase())
    raw = raw.with_columns(pl.col(*PMHC).str.strip_chars()).with_columns(
        normalize_mhc_class(pl.col('mhc_class')))
    if locus == 'paired':
        if 'pair_id' not in raw:
            raise ValueError('linked queries require pair_id, cell_id or clone_id')
        if raw.filter(pl.col('pair_id').is_null() | (pl.col('pair_id').str.strip_chars()=='')).height:
            raise ValueError('missing paired query identity')
        if raw.group_by('pair_id').agg(pl.col('species', *PMHC).n_unique()).filter(
            pl.any_horizontal(pl.col('species', *PMHC)>1)).height:
            raise ValueError('paired query has conflicting species or candidate pMHC')
    raw = raw.filter(pl.col('species') == species)
    if locus in ('TRA','TRB'):
        raw = raw.filter(pl.col('locus') == locus)
    if raw.filter(pl.any_horizontal([pl.col(c).str.strip_chars().eq('').fill_null(True) for c in PMHC])).height:
        raise ValueError('each query requires nonempty epitope, mhc_a, mhc_b and mhc_class')
    if locus != 'paired':
        c = 'a' if locus == 'TRA' else 'b'
        return raw.select('query_id',*[k for k in ('sequence_id','pair_id','source_query_id') if k in raw],
                          'species',*PMHC,pl.col('junction_aa').alias('cdr3'+c),pl.col('v_call').alias('v'+c)), report
    if raw.filter(~pl.col('locus').is_in(['TRA','TRB']).fill_null(False)).height:
        raise ValueError('paired input requires unambiguous TRA/TRB chain roles')
    if raw.group_by('pair_id','locus').len().filter(pl.col('len')>1).height:
        raise ValueError('multiple chains per pair and locus')
    q = raw.group_by('pair_id', maintain_order=True).agg(
        pl.col('species',*PMHC).first(),
        *[pl.col(column).filter(pl.col('locus')==axis).first().alias(name+c)
          for axis,c in [('TRA','a'),('TRB','b')]
          for column,name in [('junction_aa','cdr3'),('v_call','v'),('query_id','chain_query_id')]],
    )
    return q.with_columns(pl.col('pair_id').alias('query_id')), report


def main(a):
    from pathlib import Path
    from importlib.metadata import version
    from importlib.resources import files
    import hashlib
    import json
    from .. import db
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..evalue.equal_tail import equal_tail_scores
    from ..runner.annotation import _validate_output_paths
    prefix = _validate_output_paths(a.output_prefix, [a.sample,a.vdjdb])
    species = _organism(a.species)
    q, ingestion = _queries(a.sample, a.locus, species)
    if q.is_empty():
        raise ValueError('no queries for the selected species and locus')
    ref = db.load(a.vdjdb,species={'human':'HomoSapiens','mouse':'MusMusculus'}[species]).filter('reference_valid')
    scores = equal_tail_scores(q, ref, locus=a.locus, species=species, threads=a.threads,
                              exclude_exact=a.exclude_exact, background_scope=a.background_scope)
    prefix.parent.mkdir(parents=True,exist_ok=True)
    scores.write_csv(str(prefix)+'.scores.tsv',separator='\t')
    template = files('vdjmatch.resources')/'vloops'/(species+'_v_loops.tsv')
    report = dict(command='equal-tail', reference=db.provenance(a.vdjdb), sample_sha256=sha256(Path(a.sample)),
                  v_loops_sha256=hashlib.sha256(template.read_bytes()).hexdigest(), ingestion=ingestion,
                  species=species, locus=a.locus, exclude_exact=a.exclude_exact, threads=a.threads,
                  background_scope=a.background_scope, reference_population='curated-benchmark',
                  versions={p:version(p) for p in ('vdjmatch','vdjtools','seqtree')},
                  target_selection='declared compatible full pMHC', queries=len(scores),
                  status_counts=scores.group_by('status').len().to_dicts(),
                  score='negative log product of two plus-one empirical marginal tails',
                  nearest='clipped BLOSUM62, resolved V, junction trim 3/2, minimum length 5',
                  kernel='v-bounded BLOSUM62, trim 3/3, stem depth 4, minimum length 8')
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(report,indent=2)+'\n')
