"""Epitope-level diagnostic populations; MHC is a descriptive modal label."""

MODAL_MHC_RULE = ('most frequent observed MHC-A/MHC-B pair among valid reference source records '
    'pooled across TRA/TRB within species/epitope; ties lexical, nulls last; '
    'MHC class is the mode within that pair; all MHC contexts contribute')
ELIGIBILITY_RULE = 'at least min_records valid source records per species/chain/epitope; tail excluded from classification'


def pool_epitopes(raw, *, reference_population='curated-benchmark', min_records=30,
                  keep_tail_chains=False):
    """Pool observations without deduplication; retain original MHC in source columns."""
    import polars as pl
    if isinstance(min_records, bool) or not isinstance(min_records, int) or min_records < 1:
        raise ValueError('min-records must be a positive integer')
    context = ['species', 'epitope']
    key = [*context, 'gene']
    mhc = ['mhc_a', 'mhc_b', 'mhc_class']
    raw = raw.filter(pl.col('gene').is_in(['TRA', 'TRB']))
    valid = pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
    if reference_population == 'curated-benchmark':
        valid &= pl.col('reference_valid').fill_null(False)
    elif reference_population != 'all-valid':
        raise ValueError('unknown reference population')
    selected = raw.filter(valid)
    counts = selected.group_by(key).len(name='valid_source_records')
    inventory = raw.group_by(key).len(name='source_records').join(counts, on=key,
        how='left', nulls_equal=True).with_columns(pl.col('valid_source_records').fill_null(0)).with_columns(
        (pl.col('valid_source_records') >= min_records).alias('classification_eligible')).sort(key)
    pairs = selected.group_by([*context, *mhc[:2]]).len(name='modal_mhc_records').sort(
        [*context, 'modal_mhc_records', *mhc[:2]], descending=[False, False, True, False, False],
        nulls_last=True).unique(context, keep='first', maintain_order=True)
    classes = selected.join(pairs.select(*context, *mhc[:2]), on=[*context, *mhc[:2]],
        nulls_equal=True).group_by([*context, 'mhc_class']).len().sort(
        [*context, 'len', 'mhc_class'], descending=[False, False, True, False], nulls_last=True).unique(
            context, keep='first', maintain_order=True).drop('len')
    labels = pairs.join(classes, on=context, nulls_equal=True, validate='1:1').sort(context)
    eligible = inventory.filter(pl.col('classification_eligible')).select(key)
    if keep_tail_chains:
        eligible = eligible.select(context).unique()
        join_key = context
    else:
        join_key = key
    pooled = raw.join(eligible, on=join_key, how='semi', nulls_equal=True).rename(
        {c: 'source_'+c for c in mhc}).join(labels.select(*context, *mhc), on=context,
            how='left', nulls_equal=True, validate='m:1')
    info = dict(epitope_profiles=True, min_records=min_records, modal_mhc_rule=MODAL_MHC_RULE,
        eligibility_rule=ELIGIBILITY_RULE, epitope_labels=labels.to_dicts(),
        epitope_inventory=inventory.to_dicts())
    return pooled, info


def load_reference(a, species, *, locus=None, keep_tail_chains=False):
    """Use the existing loader; pool before selecting a chain so labels agree."""
    from .. import db
    locus = locus if locus is not None else getattr(a, 'locus', None)
    if not getattr(a, 'epitope_profiles', False):
        return db.load(a.vdjdb, species=species, gene=locus, epitope=a.epitope), dict(epitope_profiles=False)
    if getattr(a, 'targets_from_sample', None):
        raise ValueError('epitope profiles cannot use compatible-MHC target unions')
    raw = db.load(a.vdjdb, species=species, epitope=a.epitope)
    pooled, info = pool_epitopes(raw, reference_population=a.reference_population,
        min_records=a.min_records, keep_tail_chains=keep_tail_chains)
    if locus is not None:
        import polars as pl
        pooled = pooled.filter(pl.col('gene') == locus)
    return pooled, info


def add_options(parser, *, default=False):
    parser.add_argument('--epitope-profiles', action='store_true', default=default,
        help='pool MHC contexts per epitope; use a shared modal MHC-A/B label')
    parser.add_argument('--min-records', type=int, default=30,
        help='minimum valid source records per species/chain/epitope for profiles (default: 30)')


def write_tables(prefix, info):
    """Flat cohort/label deliverables, including the tail reserved for later studies."""
    if not info.get('epitope_profiles'):
        return []
    import polars as pl
    from pathlib import Path
    paths = []
    for key, suffix in [('epitope_labels', 'epitopes'), ('epitope_inventory', 'eligibility')]:
        path = Path(str(prefix)+'.'+suffix+'.tsv')
        path.parent.mkdir(parents=True, exist_ok=True)
        pl.DataFrame(info[key]).write_csv(path, separator='\t')
        paths.append(path)
    return paths
