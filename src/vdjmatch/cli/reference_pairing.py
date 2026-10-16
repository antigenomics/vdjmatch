"""Count genuine observed TRA/TRB links across fixed epitope-specific communities.

Validate real community exports against their graph manifests and raw reference,
then join exact full-clonotype/pMHC identities. Preserve unlinked marginals and
missing memberships; report linked observations and distinct paired clonotypes
separately, without fitting graphs or inferring pairs from marginal frequencies.
"""
from pathlib import Path

PMHC = ['epitope', 'mhc_a', 'mhc_b', 'mhc_class']
CONTEXT = ['species', *PMHC]
KEY = [*CONTEXT, 'locus', 'cdr3', 'v', 'j']
PAIR_KEY = ['alpha', 'v_alpha', 'j_alpha', 'beta', 'v_beta', 'j_beta']


def register(subparsers):
    from .epitope import add_options

    p = subparsers.add_parser('reference-pairing', help=__doc__.splitlines()[0])
    p.add_argument('--vdjdb', required=True)
    p.add_argument('--species', default='human')
    p.add_argument('--epitope', action='append')
    add_options(p)
    p.add_argument('--reference-population', choices=['curated-benchmark', 'all-valid'], default='curated-benchmark')
    p.add_argument('--alpha-communities', action='append', required=True)
    p.add_argument('--beta-communities', action='append', required=True)
    p.add_argument('--geometry', choices=['levenshtein', 'sw'], required=True)
    p.add_argument('--output-prefix', required=True)
    p.set_defaults(func=main)


def _communities(path, *, locus, species, geometry, policy, reference,
        epitope_profiles=False, min_records=30):
    """Read only real memberships after checking their owning graph manifest."""
    import json
    import polars as pl
    from ..db.cache import sha256
    from ..evalue.control import _organism

    path = Path(path).resolve()
    if not path.name.endswith('.communities.tsv'):
        raise ValueError('community input must be PREFIX.communities.tsv')
    manifest_path = path.with_name(path.name.removesuffix('.communities.tsv')+'.manifest.json')
    manifest = json.loads(manifest_path.read_text())
    if (not isinstance(manifest, dict) or manifest.get('locus') != locus or
            _organism(manifest.get('species')) != _organism(species) or
            manifest.get('geometry') != geometry or manifest.get('reference_population') != policy or
            'targets_sha256' not in manifest or manifest['targets_sha256'] is not None):
        raise ValueError('community manifest must match species/locus/geometry/population and use exact pMHC groups')
    if (manifest.get('epitope_profiles', False) != epitope_profiles or
            (epitope_profiles and manifest.get('min_records') != min_records)):
        raise ValueError('community manifest grouping must match --epitope-profiles and --min-records')
    raw_identity = manifest.get('reference', {})
    identity_key = 'input_sha256' if 'input_sha256' in reference else 'files'
    if not isinstance(raw_identity, dict) or not reference.get(identity_key) or raw_identity.get(identity_key) != reference[identity_key]:
        raise ValueError('community manifest raw reference identity does not match --vdjdb')
    digest = sha256(path)
    if not isinstance(manifest.get('tables'), dict) or manifest['tables'].get(path.name) != digest:
        raise ValueError('community table SHA256 does not match its companion manifest')
    frame = pl.scan_csv(path, separator='\t', schema_overrides={
        **{c: pl.String for c in KEY}, 'population': pl.String, 'replicate': pl.Int64,
        'community': pl.Int64, 'geometry_available': pl.Boolean})
    if not {*KEY, 'population', 'replicate', 'community', 'geometry_available'} <= set(frame.collect_schema().names()):
        raise ValueError('community table lacks full clonotype/pMHC keys or membership dispositions')
    real = frame.filter((pl.col('population') == 'reference') & (pl.col('replicate') == -1)).select(
        *KEY, 'community', 'geometry_available').collect(engine='streaming')
    if (real.filter((pl.col('species') != species).fill_null(True) |
            (pl.col('locus') != locus).fill_null(True)).height or
            real.select(KEY).n_unique() != real.height or
            real['geometry_available'].null_count() or
            real.filter((pl.col('geometry_available') != pl.col('community').is_not_null()) |
                (pl.col('community') < 0).fill_null(False)).height):
        raise ValueError('real community rows require unique exact population keys and consistent modeled dispositions')
    return real.select(*KEY, 'community', 'geometry_available'), dict(
        path=str(path), sha256=digest, manifest_path=str(manifest_path), manifest_sha256=sha256(manifest_path))


def pairing_tables(reference, alpha, beta):
    """Join fixed memberships by full keys; preserve all unlinked chain marginals."""
    import polars as pl
    from ..evalue.paired import build_paired_ref

    # Caller has applied the requested policy. The paired helper otherwise applies
    # curated validity again, which would incorrectly shrink an all-valid request.
    pairs = build_paired_ref(reference.with_columns(pl.lit(True).alias('reference_valid')))
    memberships = pl.concat([alpha, beta])
    marginals = reference.rename({'gene': 'locus'}).group_by(KEY).len(name='source_records')
    marginals = marginals.join(memberships, on=KEY, how='left', validate='1:1', nulls_equal=True).with_columns(
        pl.when(pl.col('geometry_available').is_null()).then(pl.lit('missing_community'))
        .when(~pl.col('geometry_available')).then(pl.lit('geometry_unavailable'))
        .otherwise(pl.lit('assigned')).alias('community_status'))
    chain_links = []
    for locus, side, sequence in [('TRA', 'alpha', 'alpha'), ('TRB', 'beta', 'beta')]:
        chain = pairs.select(*CONTEXT, 'complex_id', pl.lit(locus).alias('locus'),
            pl.col(sequence).alias('cdr3'), pl.col('v_'+side).alias('v'), pl.col('j_'+side).alias('j'))
        chain_links.append(chain)
        lookup = marginals.filter(pl.col('locus') == locus).select(*CONTEXT,
            pl.col('cdr3').alias(sequence), pl.col('v').alias('v_'+side), pl.col('j').alias('j_'+side),
            pl.col('community').alias('community_'+side),
            pl.col('community_status').alias('community_status_'+side))
        pairs = pairs.join(lookup, on=[*CONTEXT, sequence, 'v_'+side, 'j_'+side],
            how='left', validate='m:1', nulls_equal=True)
    links = pl.concat(chain_links)
    linked_keys = links.group_by(KEY).agg(pl.col('complex_id').n_unique().alias('complete_complex_links'))
    marginals = marginals.join(linked_keys, on=KEY, how='left', validate='1:1', nulls_equal=True).with_columns(
        pl.col('complete_complex_links').fill_null(0))
    pair_groups = [*CONTEXT, 'community_alpha', 'community_beta', 'community_status_alpha', 'community_status_beta']
    community_pairs = pairs.group_by(pair_groups).agg(
        pl.col('complex_id').n_unique().alias('complete_complex_links'),
        pl.struct(PAIR_KEY).n_unique().alias('distinct_paired_clonotypes'))
    pair_totals = pairs.group_by(CONTEXT).agg(
        pl.col('complex_id').n_unique().alias('complete_complex_links'),
        pl.struct(PAIR_KEY).n_unique().alias('distinct_paired_clonotypes'))
    coverage = marginals.group_by(*CONTEXT, 'locus').agg(
        pl.col('source_records').sum().alias('chain_source_records'),
        pl.len().alias('chain_full_clonotypes'),
        (pl.col('complete_complex_links') > 0).sum().alias('paired_chain_full_clonotypes'),
        (pl.col('community_status') == 'missing_community').sum().alias('missing_community_clonotypes'),
        (pl.col('community_status') == 'geometry_unavailable').sum().alias('geometry_unavailable_clonotypes'))
    coverage = coverage.join(pair_totals, on=CONTEXT, how='left', validate='m:1', nulls_equal=True).with_columns(
        pl.col('complete_complex_links', 'distinct_paired_clonotypes').fill_null(0),
        (pl.col('paired_chain_full_clonotypes')/pl.col('chain_full_clonotypes')).alias('paired_chain_clonotype_fraction'))
    return {'.linked.tsv': pairs.sort([*CONTEXT, 'complex_id'], nulls_last=True),
        '.community-pairs.tsv': community_pairs.sort(pair_groups, nulls_last=True),
        '.pairing-coverage.tsv': coverage.sort([*CONTEXT, 'locus'], nulls_last=True),
        '.marginals.tsv': marginals.sort(KEY, nulls_last=True)}


def main(a):
    import json
    import time
    import polars as pl
    from .. import db
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..runner.annotation import _validate_output_paths

    started = time.monotonic()
    if a.min_records < 1:
        raise ValueError('--min-records must be positive')
    species = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[_organism(a.species)]
    paths = [[Path(p).resolve() for p in chain] for chain in (a.alpha_communities, a.beta_communities)]
    all_paths = [p for chain in paths for p in chain]
    protected = [a.vdjdb, *all_paths, *[p.with_name(p.name.removesuffix('.communities.tsv')+'.manifest.json') for p in all_paths]]
    prefix = _validate_output_paths(a.output_prefix, protected)
    provenance = db.provenance(a.vdjdb)
    frames, identities = [], []
    for chain_paths, locus in zip(paths, ['TRA', 'TRB']):
        batches = []
        for path in chain_paths:
            frame, identity = _communities(path, locus=locus, species=species, geometry=a.geometry,
                policy=a.reference_population, reference=provenance,
                epitope_profiles=a.epitope_profiles, min_records=a.min_records)
            batches.append(frame)
            identities.append(identity)
        frame = pl.concat(batches)
        if frame.select(KEY).n_unique() != frame.height:
            raise ValueError('community batches must have unique full clonotype keys across files')
        frames.append(frame)
    grouping = dict(epitope_profiles=False, min_records=a.min_records, modal_mhc_rule=None)
    if a.epitope_profiles:
        from .epitope import load_reference
        raw, grouping = load_reference(a, species, keep_tail_chains=True)
    else:
        raw = db.load(a.vdjdb, species=species, epitope=a.epitope)
    valid = pl.col('gene').is_in(['TRA', 'TRB']) & pl.col('cdr3').str.contains(
        r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
    if a.reference_population == 'curated-benchmark':
        valid &= pl.col('reference_valid').fill_null(False)
    selected = raw.filter(valid)
    tables = pairing_tables(selected, *frames)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for suffix, frame in tables.items():
        path = Path(str(prefix)+suffix)
        frame.with_columns(pl.lit(a.geometry).alias('geometry')).write_csv(path, separator='\t')
        outputs[path.name] = sha256(path)
    manifest = dict(reference=provenance, species=species, epitope=a.epitope,
        reference_population=a.reference_population, geometry=a.geometry, communities=identities,
        input_chain_records=raw.height, selected_chain_records=selected.height,
        excluded_chain_records=raw.height-selected.height,
        pairing='existing build_paired_ref; complete unambiguous linked TRA/TRB observations only; no inferred marginal pairing',
        joining='exact species/full-pMHC/full-junction/V/J; missing calls remain explicit; null keys compare equal',
        unit='complex links and distinct full paired clonotypes reported separately; chain source records are not independent assays',
        tables=outputs, source_sha256=sha256(Path(__file__)), seconds=time.monotonic()-started)
    manifest.update(grouping)
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0
