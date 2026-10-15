"""Fresh Pgen-stratified mock reference memberships from raw control events."""
from __future__ import annotations

import numpy as np
import polars as pl

from .control import _organism
from ..io.columns import _integer_count, _resolve


def normalize_control(frame: pl.DataFrame, *, species='human', locus='TRB') -> pl.DataFrame:
    """Collapse raw controls to unique full junctions and recombination-event weights.

    Read abundance is ignored. Prefer pooled ``occurrences``, otherwise distinct
    sample/NT/V/J events (pooled NT incidence if sample IDs are absent),
    supplied convergence, then one event per unique AA. AA-only convergence
    represents available convergence counts, not reconstructed pooled events.
    NT and supplied gene keys are normalized to uppercase; supplied sample IDs
    remain exact. Incomplete NT event metadata requires explicit event counts
    or convergence rather than silently collapsing the population to AA keys.
    Supplied population metadata selects matching rows; unknown/missing metadata
    never inherits the requested population.
    """
    organism = _organism(species)
    locus = locus.upper()
    if locus not in {'TRA', 'TRB'}:
        raise ValueError('control locus must be TRA or TRB')
    lower = {c.lower(): c for c in frame.columns}
    fields = _resolve(frame)
    sample = next((lower[c] for c in ('sample_id', 'sample', 'repertoire_id') if c in lower), None)
    nt = next((lower[c] for c in ('junction', 'cdr3nt', 'junction_nt', 'cdr3_nt') if c in lower), None)
    sequence = next((lower[c] for c in ('junction_aa', 'cdr3aa', 'cdr3', 'cdr3_aa') if c in lower), None)
    if sequence is None:
        raise ValueError('control requires a full junction_aa or legacy cdr3aa column')
    if sequence.lower() == 'cdr3_aa' or (sequence.lower() == 'cdr3' and
            any(c in lower for c in ('v_call', 'j_call', 'sequence_id'))):
        raise ValueError('AIRR control requires junction_aa; bare CDR3 excludes anchors')
    for name in ('species', 'organism'):
        if name in lower:
            aliases = {}
            for value in frame[lower[name]].drop_nulls().unique():
                try:
                    aliases[value] = _organism(value)
                except ValueError:
                    aliases[value] = None
            frame = frame.filter(pl.col(lower[name]).replace_strict(
                aliases, default=None, return_dtype=pl.String).eq(organism).fill_null(False))
    gene_columns = [fields[k] for k in ('v', 'j') if k in fields]
    frame = frame.with_columns(*[pl.col(c).cast(pl.String).str.strip_chars().str.to_uppercase()
        for c in gene_columns])
    if gene_columns:
        alpha = pl.any_horizontal(*[pl.col(c).str.contains(r'TRA[VJ]').fill_null(False) for c in gene_columns])
        beta = pl.any_horizontal(*[pl.col(c).str.contains(r'TRB[VJ]').fill_null(False) for c in gene_columns])
        if frame.filter(alpha & beta).height:
            raise ValueError('control has conflicting V/J locus annotations')
        inferred = pl.when(alpha).then(pl.lit('TRA')).when(beta).then(pl.lit('TRB')).otherwise(None)
    else:
        inferred = pl.lit(None, dtype=pl.String)
    if 'locus' in fields:
        declared = pl.col(fields['locus']).cast(pl.String).str.strip_chars().str.to_uppercase().replace(
            {'ALPHA': 'TRA', 'BETA': 'TRB'})
        if frame.filter(inferred.is_not_null() & declared.is_not_null() & (inferred != declared)).height:
            raise ValueError('control V/J locus conflicts with explicit locus metadata')
        frame = frame.filter(declared.eq(locus).fill_null(False))
    elif gene_columns:
        if nt and 'occurrences' not in lower and 'convergence' not in lower and frame.filter(inferred.is_null()).height:
            raise ValueError('control has incomplete NT event metadata or unresolved V/J locus')
        frame = frame.filter(inferred.eq(locus).fill_null(False))
    if not frame.height:
        raise ValueError('control has no rows matching requested species/locus')
    frame = frame.with_columns(pl.col(sequence).cast(pl.String).str.strip_chars().str.to_uppercase().alias('junction_aa'))
    if frame.select((pl.col('junction_aa').is_null() |
            ~pl.col('junction_aa').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$')).fill_null(True).any()).item():
        raise ValueError('control contains missing or invalid amino-acid junctions')
    if nt:
        frame = frame.with_columns(pl.col(nt).cast(pl.String).str.strip_chars().str.to_uppercase())
    event_columns = [c for c in (nt, *gene_columns, sample) if c]
    incomplete_events = bool(nt) and any(frame[c].null_count() or
        frame[c].cast(pl.String).str.strip_chars().eq('').any() for c in event_columns)
    if incomplete_events and 'occurrences' not in lower and 'convergence' not in lower:
        raise ValueError('control has incomplete NT event metadata; supply complete event keys or occurrences/convergence')
    if 'occurrences' in lower:
        weight = lower['occurrences']
        frame = frame.with_columns(_integer_count(weight).alias('event_weight'))
        if frame['event_weight'].null_count():
            raise ValueError('occurrences must contain non-negative integer event counts')
        result = frame.group_by('junction_aa').agg(pl.col('event_weight').sum())
    elif nt and not incomplete_events:
        keys = ['junction_aa', *event_columns]
        if sample is None and 'incidence' in lower:
            frame = frame.with_columns(_integer_count(lower['incidence']).alias('event_weight'))
            if frame['event_weight'].null_count():
                raise ValueError('incidence must contain non-negative integer event counts')
            result = frame.group_by(keys).agg(pl.col('event_weight').max()).group_by('junction_aa').agg(
                pl.col('event_weight').sum())
        else:
            result = frame.unique(subset=keys).group_by('junction_aa').len(name='event_weight')
    elif 'convergence' in lower:
        frame = frame.with_columns(_integer_count(lower['convergence']).alias('event_weight'))
        if frame['event_weight'].null_count():
            raise ValueError('convergence must contain non-negative integer event counts')
        result = frame.group_by('junction_aa', *([sample] if sample else [])).agg(
            pl.col('event_weight').max()).group_by('junction_aa').agg(pl.col('event_weight').sum())
    else:
        result = frame.select('junction_aa').unique().with_columns(pl.lit(1).alias('event_weight'))
    result = result.filter(pl.col('event_weight') > 0).with_columns(
        pl.col('event_weight').cast(pl.Int64), pl.lit(organism).alias('species'), pl.lit(locus).alias('locus'))
    if not result.height:
        raise ValueError('control has no positive-weight junctions')
    return result.sort('junction_aa')


def compute_pgen(reference: pl.DataFrame, control: pl.DataFrame, *, species='human', locus='TRB', threads=1) -> pl.DataFrame:
    """Compute marginal full-junction Pgen once over the reference/control union.

    Bundled human OLGA or mouse ARDA models are selected by species/locus.
    Supplied Pgen columns are ignored; no calculated values are persisted.
    """
    from vdjtools.model import load_bundled, native

    organism = _organism(species)
    if locus not in {'TRA', 'TRB'} or threads < 1:
        raise ValueError('Pgen requires TRA/TRB and positive threads')
    sequences = pl.concat([reference.select('junction_aa'), control.select('junction_aa')]).unique().sort('junction_aa')
    model = load_bundled(locus, 'olga' if organism == 'human' else 'arda', organism=organism, collapse=False)
    probabilities = np.asarray(native.pgen_aa_batch(model, sequences['junction_aa'].to_list(),
        mismatches=0, threads=threads), dtype=np.float64)
    if probabilities.shape != (sequences.height,) or not np.all(
            np.isfinite(probabilities) & (probabilities >= 0) & (probabilities <= 1)):
        raise RuntimeError('native Pgen returned invalid probabilities')
    return sequences.with_columns(pl.Series('pgen', probabilities))


def pgen_strata(pgen: pl.DataFrame) -> pl.DataFrame:
    """Return unique junction/Pgen-bin keys: rounded log2(Pgen), zero=-2147483648.

    Fixed one-bit bins use nearest-integer, ties-to-even rounding. Zero mass
    has its own sentinel rather than being mixed with small positive values.
    """
    if pgen['junction_aa'].n_unique() != pgen.height:
        raise ValueError('Pgen map must have unique junctions')
    values = pgen['pgen'].cast(pl.Float64).to_numpy()
    if not np.all(np.isfinite(values) & (values >= 0) & (values <= 1)):
        raise ValueError('Pgen map must contain finite probabilities from zero to one')
    bins = np.full(len(values), -2147483648, dtype=np.int64)
    positive = values > 0
    bins[positive] = np.rint(np.log2(values[positive])).astype(np.int64)
    return pgen.select('junction_aa').with_columns(pl.Series('pgen_bin', bins))


def iter_mock_references(reference: pl.DataFrame, control: pl.DataFrame, *, replicates: int,
        seed: int = 0, pgen: pl.DataFrame | None = None, species='human', locus='TRB', threads=1):
    """Yield one replicate of null event memberships matching every reference source row.

    ``reference`` has ``group_id,junction_aa``; repeated source rows remain in
    the target counts. ``control`` is normalized by :func:`normalize_control`.
    Match rounded log2(Pgen) in fixed one-bit bins, with a separate zero bin.
    Sample integer event slots without replacement within each group/replicate;
    independent events sharing an AA junction retain repeated membership rows.
    Search-key deduplication belongs to the downstream scorer. Incidental
    reference/control overlaps remain eligible random matches. Insufficient
    strata fail; sampled controls are not asserted experimental nonbinders.
    ``pgen`` is an optional in-memory computed map for composition/testing.
    NumPy's bounded-memory hypergeometric backend requires fewer than one
    billion events per stratum; larger pools fail explicitly.
    """
    if replicates < 1:
        raise ValueError('replicates must be positive')
    reference = reference.select('group_id', 'junction_aa')
    if not reference.height or any(reference[c].null_count() for c in reference.columns):
        raise ValueError('reference requires nonempty non-null group_id and junction_aa')
    if control['junction_aa'].n_unique() != control.height:
        raise ValueError('control must have unique normalized junctions')
    # Always discard any source-supplied Pgen; only a fresh explicit map is used.
    pgen = compute_pgen(reference, control, species=species, locus=locus, threads=threads) if pgen is None else pgen
    strata = pgen_strata(pgen).rename({'pgen_bin': '_bin'})
    targets = reference.join(strata, on='junction_aa', how='left', validate='m:1')
    pool = control.select('junction_aa', 'event_weight').join(strata, on='junction_aa', how='left', validate='1:1')
    if targets['_bin'].null_count() or pool['_bin'].null_count():
        raise ValueError('Pgen map must cover every reference and control junction')
    pool = pool.with_columns(_integer_count('event_weight').alias('event_weight'))
    if pool['event_weight'].null_count() or (pool['event_weight'] <= 0).any():
        raise ValueError('control event weights must be positive integers')
    pool = pool.sort('junction_aa')
    pools = {key[0]: (part['junction_aa'].to_numpy(), part['event_weight'].to_numpy())
             for key, part in pool.partition_by('_bin', as_dict=True).items()}
    counts = targets.group_by('group_id', '_bin').len().sort('group_id', '_bin')
    for group, bin_id, count in counts.iter_rows():
        available = int(pools[bin_id][1].sum()) if bin_id in pools else 0
        if available < count:
            label = 'zero_mass' if bin_id == -2147483648 else str(bin_id)
            raise ValueError(f'insufficient control events for group {group!r}, Pgen bin {label}: need {count}, have {available}')
        if available >= 10**9:
            raise ValueError('event stratum exceeds NumPy hypergeometric backend limit: fewer than 1000000000 events required')
    rng = np.random.default_rng(seed)
    for replicate in range(replicates):
        results = []
        for group, bin_id, count in counts.iter_rows():
            sequences, weight = pools[bin_id]
            multiplicity = rng.multivariate_hypergeometric(weight, count)
            selected = np.repeat(sequences, multiplicity)
            results.append(pl.DataFrame({'replicate': replicate, 'group_id': group, 'junction_aa': selected}))
        yield pl.concat(results).with_columns(pl.col('group_id').cast(reference.schema['group_id'])).sort('group_id', 'junction_aa')


def sample_mock_references(reference: pl.DataFrame, control: pl.DataFrame, *, replicates: int,
        seed: int = 0, pgen: pl.DataFrame | None = None, species='human', locus='TRB', threads=1) -> pl.DataFrame:
    """Collect mock memberships; use ``iter_mock_references`` for bounded memory."""
    return pl.concat(iter_mock_references(reference, control, replicates=replicates,
        seed=seed, pgen=pgen, species=species, locus=locus, threads=threads))


def iter_mock_paired_references(reference: pl.DataFrame, control: pl.DataFrame, *,
        species='human', replicates: int, seed: int = 0, threads: int = 1,
        pgen_alpha: pl.DataFrame | None = None, pgen_beta: pl.DataFrame | None = None):
    """Yield Pgen-matched pairs under an explicitly independent-chain null.

    Genuine wide reference pair rows have ``group_id,junction_aa_alpha,
    junction_aa_beta``. Every source pair is retained in each mock group.
    Control rows must declare their TRA/TRB locus; supplied species metadata
    selects only the requested species. Each chain's event slots are sampled
    independently with its own Pgen strata, then independently permuted within
    each group/replicate before pairing. Sorted AA lanes are never rank-zipped.

    One fresh native Pgen batch per chain covers its reference/control union,
    unless an explicit transient computed map is supplied. This sampler does
    not calibrate linked/all-marginal axes or selected full scores; their frozen
    statistic and selection event remain the scorer adapter's responsibility.
    """
    if replicates < 1 or threads < 1:
        raise ValueError('paired sampler requires positive replicates and threads')
    fields = _resolve(control)
    if 'locus' not in fields:
        raise ValueError('paired control requires supplied TRA/TRB locus metadata')
    required = ['group_id', 'junction_aa_alpha', 'junction_aa_beta']
    if not set(required).issubset(reference.columns) or not reference.height:
        raise ValueError('paired reference requires nonempty genuine paired junction rows')
    if reference['group_id'].null_count():
        raise ValueError('paired reference group_id cannot be missing')
    for lane in ('alpha', 'beta'):
        name = f'junction_aa_{lane}'
        if reference.filter(pl.col(name).is_null() |
                ~pl.col(name).cast(pl.String).str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$')).height:
            raise ValueError('paired reference requires both valid full junction chains on every row')
    for column in ('species', 'organism'):
        if column in reference.columns and (reference[column].null_count() or any(
                _organism(value) != _organism(species) for value in reference[column].unique())):
            raise ValueError('paired reference species does not match requested species')
    lanes = []
    streams = np.random.SeedSequence(seed).spawn(4)
    for lane, locus, probabilities, stream in zip(('alpha', 'beta'), ('TRA', 'TRB'),
            (pgen_alpha, pgen_beta), streams[:2]):
        raw = reference.select('group_id', pl.col(f'junction_aa_{lane}').alias('junction_aa'))
        pool = normalize_control(control, species=species, locus=locus)
        if probabilities is None:
            probabilities = compute_pgen(raw, pool, species=species, locus=locus, threads=threads)
        lanes.append(iter_mock_references(raw, pool, replicates=replicates, seed=stream,
            pgen=probabilities, species=species, locus=locus, threads=threads))
    alpha_rng, beta_rng = map(np.random.default_rng, streams[2:])
    for alpha, beta in zip(*lanes, strict=True):
        paired = []
        beta_groups = beta.partition_by('group_id', as_dict=True)
        for key, part in alpha.partition_by('group_id', as_dict=True).items():
            other = beta_groups[key]
            if part.height != other.height:
                raise RuntimeError('paired mock marginal source counts disagree')
            paired.append(part.select('group_id', 'replicate').with_columns(
                pl.Series('draw_id', np.arange(part.height, dtype=np.int64)),
                pl.Series('junction_aa_alpha', alpha_rng.permutation(part['junction_aa'].to_numpy())),
                pl.Series('junction_aa_beta', beta_rng.permutation(other['junction_aa'].to_numpy()))))
        yield pl.concat(paired).sort('group_id', 'draw_id')
