"""Fresh Pgen-stratified mock reference memberships from raw control events."""
from __future__ import annotations

import numpy as np
import polars as pl

from .control import _organism
from ..io.columns import _integer_count, _resolve


EVENT_SAMPLER = "uniform-event-slots-v1"
PGEN_BINNING = "rounded log2(Pgen), ties-to-even; absent reference bins, including zero mass, share the lowest occupied control bin with enough event slots for all mapped source rows plus its existing occupants"


def _pgen_model_identity(species, locus):
    """Hash the existing immutable bundled model inputs and native implementation."""
    import hashlib
    import json
    from ..precursor.model import bundled_model_path
    from pathlib import Path
    import vdjtools._core as native_core
    from ..db.cache import sha256

    organism = _organism(species)
    source = 'olga' if organism == 'human' else 'arda'
    root = bundled_model_path(locus, source, organism)
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(root.iterdir(), key=lambda p: p.name) if p.is_file()}
    return dict(source=source, organism=organism, locus=locus, files=hashes,
        sha256=hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest(),
        metadata=json.loads((root/'manifest.json').read_text()),
        native_sha256=sha256(Path(native_core.__file__)))


def load_pgen_features(path, reference, control, *, species='human', locus='TRB', mismatches=1):
    """Validate a durable feature bank with the requested probability and model preprocessing.

    These requested durable model features are distinct from resampled controls.
    Require the complete reference/control junction union; never fill missing mass
    or trust a Pgen column supplied in an ordinary rearrangement table.
    """
    import hashlib
    import json
    from pathlib import Path
    from ..db.cache import sha256

    path = Path(path).resolve()
    if not path.name.endswith('.features.parquet'):
        raise ValueError('Pgen features require PREFIX.features.parquet and PREFIX.manifest.json')
    manifest_path = path.with_name(path.name.removesuffix('.features.parquet')+'.manifest.json')
    manifest = json.loads(manifest_path.read_text())
    if not isinstance(manifest, dict) or any(not isinstance(manifest.get(key), dict)
            for key in ('features', 'options', 'model')):
        raise ValueError('Pgen feature manifest requires features, options and model identities')
    feature, options, model = (manifest[key] for key in ('features', 'options', 'model'))
    organism = _organism(species)
    if (type(options.get('mismatches')) is not int or options['mismatches'] != mismatches or
            _organism(options.get('species')) != organism or options.get('locus') != locus):
        raise ValueError(f'Pgen features must match species/locus and Hamming mismatches={mismatches}')
    recorded_path = feature.get('path')
    if not isinstance(recorded_path, str) or Path(recorded_path).name != path.name:
        raise ValueError('Pgen feature manifest basename does not match the requested table')
    relative = Path(recorded_path).name if Path(recorded_path).is_absolute() else recorded_path
    if (manifest_path.parent/relative).resolve() != path:
        raise ValueError('Pgen feature manifest path does not resolve beside its table')
    table_hash = sha256(path)
    if feature.get('sha256') != table_hash:
        raise ValueError('Pgen feature table SHA256 does not match its manifest')
    files, metadata = model.get('files'), model.get('metadata')
    expected = _pgen_model_identity(organism, locus)
    if (model.get('collapse_alleles') is not False or not isinstance(files, dict) or not files or
            any(expected['files'].get(name) != digest for name, digest in files.items()) or
            not isinstance(metadata, dict) or _organism(metadata.get('organism')) != organism or metadata.get('locus') != locus):
        raise ValueError('Pgen feature model files/metadata do not match the uncollapsed species/locus model')
    expected_metadata = expected['metadata']
    if any(metadata.get(key) != expected_metadata[key] for key in ('events', 'chain_type')):
        raise ValueError('Pgen feature model event/chain metadata does not match the bundled model')
    names = [*expected_metadata['events'], 'genes_v', 'genes_j']
    if expected_metadata['chain_type'] == 'VDJ':
        names.append('genes_d')
    consumed = {'manifest.json', *[next(name+suffix for suffix in ('.parquet', '.tsv', '.csv')
        if name+suffix in expected['files']) for name in names]}
    if 'training.json' in expected['files']:
        consumed.add('training.json')
    if set(files) != consumed:
        raise ValueError('Pgen feature model must identify every consumed event and gene file')
    digest = hashlib.sha256(''.join(f'{name}\t{files[name]}\n' for name in sorted(files)).encode()).hexdigest()
    if model.get('sha256') != digest:
        raise ValueError('Pgen feature model SHA256 does not match its consumed files')
    from vdjtools.model import io as model_io, native as native_bridge
    software = manifest.get('software', {})
    implementation = {'model_io_sha256': sha256(Path(model_io.__file__)),
        'native_bridge_sha256': sha256(Path(native_bridge.__file__))}
    if model.get('anchor_repair') != 'vdjtools load_model' or any(
            software.get(name) != digest for name, digest in implementation.items()):
        raise ValueError('Pgen feature model preprocessing/anchor repair differs from current implementation')
    frame = pl.read_parquet(path)
    required = {'junction_aa', 'pgen', 'species', 'locus'}
    if (not required <= set(frame.columns) or not frame.schema['pgen'].is_numeric() or
            any(frame.schema[c] != pl.String for c in ('junction_aa', 'species', 'locus'))):
        raise ValueError('Pgen feature table requires junction_aa, numeric pgen, species and locus')
    frame = frame.with_columns(pl.col('pgen').cast(pl.Float64))
    if (any(frame[c].null_count() for c in required) or
            frame['junction_aa'].n_unique() != frame.height or not frame['junction_aa'].is_sorted() or
            frame.filter(~pl.col('junction_aa').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$') |
                ~pl.col('pgen').is_finite() | ~pl.col('pgen').is_between(0, 1) |
                (pl.col('species') != organism) | (pl.col('locus') != locus)).height):
        raise ValueError('Pgen features require sorted unique canonical junctions, matching species/locus and finite probabilities in [0,1]')
    needed = pl.concat([reference.select('junction_aa'), control.select('junction_aa')]).unique()
    missing = needed.join(frame.select('junction_aa'), on='junction_aa', how='anti').height
    if missing:
        raise ValueError(f'Pgen features do not cover {missing} required reference/control junctions')
    identity = dict(path=str(path), sha256=table_hash, manifest_path=str(manifest_path),
        manifest_sha256=sha256(manifest_path), options=options, model=model, rows=frame.height)
    return frame.select('junction_aa', pl.col('pgen').cast(pl.Float64)), identity


def _sample_event_indices(cumulative, count, rng):
    """Draw uniform event slots, then map to keys without expanding the pool.

    Subsets have the multivariate hypergeometric law. Slot order is not uniform;
    callers assigning distinct source roles must shuffle the returned indices.
    ``cumulative`` contains positive integer event counts cumulatively summed.
    """
    slots = rng.choice(int(cumulative[-1]), count, replace=False, shuffle=False)
    return np.searchsorted(cumulative, slots, side="right")


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
    return _normalize_control(frame, species=species, locus=locus, gene_keys=False)


def _modal_gene(frame, column, *, species, locus, kind='v'):
    """Choose a label-independent observed mode, falling back to the bundled model."""
    observed=frame.group_by(column).len().sort('len',column,descending=[True,False])
    if observed.height:return observed[column][0]
    from vdjtools.model import load_bundled
    model=load_bundled(locus,'olga' if _organism(species)=='human' else 'arda',organism=_organism(species),collapse=False)
    prior=model.tables[kind+'_choice'];allele=kind+'_allele'
    return prior.filter(pl.col(allele).str.contains('^'+locus+kind.upper())).sort('p',allele,descending=[True,False])[allele][0]


def control_events(frame: pl.DataFrame, *, species='human', locus='TRB') -> pl.DataFrame:
    """Retain the joint junction/V/J event law for full-score randomization.

    Read counts do not define events. Missing V/J is imputed from the selected population's modal supplied call,
    or the bundled model mode when absent. Source event identities stay distinct
    and imputation flags accompany affected outputs. Convergence
    totals shared between conflicting V/J keys cannot reconstruct this law.
    """
    return _normalize_control(frame, species=species, locus=locus, gene_keys=True)


def _normalize_control(frame, *, species, locus, gene_keys):
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
        if not gene_keys and frame.filter(alpha & beta).height:
            raise ValueError('control has conflicting V/J locus annotations')
        inferred = pl.when(alpha & beta).then(None).when(alpha).then(pl.lit('TRA')).when(beta).then(pl.lit('TRB')).otherwise(None)
    else:
        inferred = pl.lit(None, dtype=pl.String)
    if 'locus' in fields:
        declared = pl.col(fields['locus']).cast(pl.String).str.strip_chars().str.to_uppercase().replace(
            {'ALPHA': 'TRA', 'BETA': 'TRB'})
        if gene_keys:
            declared=pl.coalesce(pl.when(declared.eq('')).then(None).otherwise(declared),inferred,pl.lit(locus))
        if frame.filter(inferred.is_not_null() & declared.is_not_null() & (inferred != declared)).height:
            raise ValueError('control V/J locus conflicts with explicit locus metadata')
        frame = frame.filter(declared.eq(locus).fill_null(False))
    elif gene_columns:
        if not gene_keys and nt and 'occurrences' not in lower and 'convergence' not in lower and frame.filter(inferred.is_null()).height:
            raise ValueError('control has incomplete NT event metadata or unresolved V/J locus')
        frame = frame.filter(inferred.fill_null(locus if gene_keys else '').eq(locus))
    if not frame.height:
        raise ValueError('control has no rows matching requested species/locus')
    frame = frame.with_columns(pl.col(sequence).cast(pl.String).str.strip_chars().str.to_uppercase().alias('junction_aa'))
    if frame.select((pl.col('junction_aa').is_null() |
            ~pl.col('junction_aa').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$')).fill_null(True).any()).item():
        raise ValueError('control contains missing or invalid amino-acid junctions')
    if nt:
        frame = frame.with_columns(pl.col(nt).cast(pl.String).str.strip_chars().str.to_uppercase())
    keys = ['junction_aa']
    if gene_keys:
        fields=fields.copy();original_gene_columns=[];imputed=False
        for kind in ('v','j'):
            column=fields.get(kind,'_missing_'+kind)
            if kind not in fields:
                frame=frame.with_columns(pl.lit(None,dtype=pl.String).alias(column));fields[kind]=column
            value=pl.col(column).cast(pl.String).str.strip_chars().str.to_uppercase()
            valid=value.str.contains('^'+locus+kind.upper()+r'[A-Z0-9/.-]+(?:\*[A-Z0-9]+)?$').fill_null(False)
            missing=frame.select((~valid).any()).item()
            prior=pl.col(kind+'_imputed') if kind+'_imputed' in frame.columns else pl.lit(False)
            original=pl.when(prior).then(pl.col('original_'+kind)).otherwise(value) if 'original_'+kind in frame.columns else value
            imputed=imputed or frame.select(prior.any()).item()
            frame=frame.with_columns(pl.when(original.is_null()|original.eq('')).then(pl.lit('__MISSING_'+kind.upper()+'__')).otherwise(original).alias('_original_control_'+kind),(prior|~valid).alias(kind+'_imputed'))
            if missing:
                imputed=True
                replacement=_modal_gene(frame.filter(valid),column,species=organism,locus=locus,kind=kind)
                frame=frame.with_columns(pl.when(valid).then(value).otherwise(pl.lit(replacement)).alias(column))
            original_gene_columns.append('_original_control_'+kind)
            frame=frame.with_columns(pl.col(column).alias('_control_'+kind))
        gene_columns=[fields[k] for k in ('v','j')]
        keys += ['_control_v', '_control_j']
        if not nt and 'occurrences' not in lower and 'convergence' in lower:
            ambiguous = frame.select(*keys, *([sample] if sample else [])).unique().group_by(
                'junction_aa', *([sample] if sample else [])).len().filter(pl.col('len') > 1)
            if ambiguous.height:
                raise ValueError('AA convergence totals cannot resolve conflicting V/J event keys')
    event_columns = [c for c in (nt, *gene_columns, sample) if c]
    if gene_keys and imputed:event_columns += original_gene_columns
    incomplete_events = bool(nt) and any(frame[c].null_count() or
        frame[c].cast(pl.String).str.strip_chars().eq('').any() for c in event_columns)
    if incomplete_events and 'occurrences' not in lower and 'convergence' not in lower:
        raise ValueError('control has incomplete NT event metadata; supply complete event keys or occurrences/convergence')
    if 'occurrences' in lower:
        weight = lower['occurrences']
        frame = frame.with_columns(_integer_count(weight).alias('event_weight'))
        if frame['event_weight'].null_count():
            raise ValueError('occurrences must contain non-negative integer event counts')
        result = frame.group_by(keys).agg(pl.col('event_weight').sum())
    elif nt and not incomplete_events:
        event_keys = list(dict.fromkeys([*keys, *event_columns]))
        if sample is None and 'incidence' in lower:
            frame = frame.with_columns(_integer_count(lower['incidence']).alias('event_weight'))
            if frame['event_weight'].null_count():
                raise ValueError('incidence must contain non-negative integer event counts')
            result = frame.group_by(event_keys).agg(pl.col('event_weight').max()).group_by(keys).agg(
                pl.col('event_weight').sum())
        else:
            result = frame.unique(subset=event_keys).group_by(keys).len(name='event_weight')
    elif 'convergence' in lower:
        frame = frame.with_columns(_integer_count(lower['convergence']).alias('event_weight'))
        if frame['event_weight'].null_count():
            raise ValueError('convergence must contain non-negative integer event counts')
        result = frame.group_by(*keys, *([sample] if sample else [])).agg(
            pl.col('event_weight').max()).group_by(keys).agg(pl.col('event_weight').sum())
    else:
        if gene_keys and imputed:
            result=frame.select(*keys,*original_gene_columns).unique().group_by(keys).len(name='event_weight')
        else:result = frame.select(keys).unique().with_columns(pl.lit(1).alias('event_weight'))
    result = result.filter(pl.col('event_weight') > 0).with_columns(
        pl.col('event_weight').cast(pl.Int64), pl.lit(organism).alias('species'), pl.lit(locus).alias('locus'))
    if not result.height:
        raise ValueError('control has no positive-weight junctions')
    if gene_keys:
        if imputed:
            flags=frame.group_by(keys).agg(pl.col('v_imputed','j_imputed').any())
            result=result.join(flags,on=keys,validate='1:1',nulls_equal=True)
        result = result.rename({'_control_v': 'v', '_control_j': 'j'})
    return result.sort(['junction_aa', 'v', 'j'] if gene_keys else 'junction_aa')


def compute_pgen(reference: pl.DataFrame, control: pl.DataFrame, *, species='human', locus='TRB', threads=1, mismatches=0) -> pl.DataFrame:
    """Compute marginal full-junction Pgen once over the reference/control union.

    Bundled human OLGA or mouse ARDA models are selected by species/locus.
    Supplied Pgen columns are ignored; no calculated values are persisted.
    ``mismatches`` selects exact probability (0) or inclusive Hamming<=1 mass (1).
    """
    from vdjtools.model import load_bundled, native

    organism = _organism(species)
    if locus not in {'TRA', 'TRB'} or threads < 1 or type(mismatches) is not int or mismatches not in (0,1):
        raise ValueError('Pgen requires TRA/TRB, positive threads and mismatches0or1')
    sequences = pl.concat([reference.select('junction_aa'), control.select('junction_aa')]).unique().sort('junction_aa')
    model = load_bundled(locus, 'olga' if organism == 'human' else 'arda', organism=organism, collapse=False)
    probabilities = np.asarray(native.pgen_aa_batch(model, sequences['junction_aa'].to_list(),
        mismatches=mismatches, threads=threads), dtype=np.float64)
    if probabilities.shape != (sequences.height,) or not np.all(
            np.isfinite(probabilities) & (probabilities >= 0) & (probabilities <= 1)):
        raise RuntimeError('native Pgen returned invalid probabilities')
    return sequences.with_columns(pl.Series('pgen', probabilities))


def vdjbet_strata(probes: pl.DataFrame, sample_vj: pl.DataFrame,
                  model_vj: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Reproduce the compiled vdjbet probability strata and V/J audit.

    Probes contain unique junction_aa/v/j keys and joint gene-masked ``pgen``.
    The original AA-only probability join takes the maximum over keys sharing
    a junction. Count tables contain v/j/records (source rows, not read weights).
    Add-one V/J frequencies are normalized within V, preserving the compiled
    R notebook's retained dplyr grouping. The compiled floor is 2**-80.
    """
    keys=['junction_aa','v','j']
    if not {*keys,'pgen'}<=set(probes.columns) or not probes.height:
        raise ValueError('nonempty probe keys and gene-masked pgen required')
    if probes.select(keys).n_unique()!=probes.height:
        raise ValueError('probe receptor keys must be unique')
    if any(probes[c].null_count() for c in [*keys,'pgen']):
        raise ValueError('non-null probe keys and probabilities required')
    probabilities=probes['pgen'].to_numpy()
    if not np.all(np.isfinite(probabilities)&(probabilities>=0)&(probabilities<=1)):
        raise ValueError('probe probabilities must be finite and between zero and one')
    counts=[]
    for name,frame in [('sample',sample_vj),('model',model_vj)]:
        if not {'v','j','records'}<=set(frame.columns) or not frame.height:
            raise ValueError('nonempty V/J record-count tables required')
        frame=frame.select('v','j',_integer_count('records').alias(name+'_records'))
        if frame[name+'_records'].null_count() or (frame[name+'_records']<0).any():
            raise ValueError('nonnegative integer V/J record counts required')
        counts.append(frame.group_by('v','j').agg(pl.col(name+'_records').sum()))
    audit=counts[0].join(counts[1],on=['v','j'],how='full',coalesce=True,nulls_equal=True)
    audit=audit.with_columns(pl.col('sample_records','model_records').fill_null(0)+1)
    audit=audit.with_columns((pl.col(name+'_records')/pl.col(name+'_records').sum().over('v')).alias(name+'_probability') for name in ['sample','model'])
    audit=audit.with_columns((pl.col('sample_probability')/pl.col('model_probability')).alias('vj_factor')).sort('v','j',nulls_last=True)
    result=probes.with_columns(pl.col('pgen').max().over('junction_aa').alias('pcdr')).join(
        audit.select('v','j','vj_factor'),on=['v','j'],how='left',validate='m:1',maintain_order='left')
    if result['vj_factor'].null_count():
        raise ValueError('probe V/J keys absent from both count populations')
    result=result.with_columns((pl.col('pcdr')*pl.col('vj_factor')+2.**-80).alias('adjusted_pgen'))
    bins=np.rint(np.log2(result['adjusted_pgen'].to_numpy())).astype(np.int64)
    return result.with_columns(pl.Series('pgen_bin',bins)),audit


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


def _prepare_mock_control(reference, control, *, pgen=None, species='human', locus='TRB', threads=1):
    """Prepare immutable event pools for this execution, never persisted or RNG state.

    The tuple binds the exact control/Pgen objects and species/locus. Callers
    must keep those frames unchanged and pass the same objects when reusing it.
    """
    event_keys = ['junction_aa', *[c for c in ('v', 'j') if c in control.columns]]
    if control.select(event_keys).n_unique() != control.height:
        raise ValueError('control must have unique normalized event keys')
    # Always discard any source-supplied Pgen; only a fresh explicit map is used.
    pgen = compute_pgen(reference, control, species=species, locus=locus, threads=threads) if pgen is None else pgen
    strata = pgen_strata(pgen).rename({'pgen_bin': '_bin'})
    pool = control.select(*event_keys, 'event_weight').join(strata, on='junction_aa', how='left', validate='m:1')
    if pool['_bin'].null_count():
        raise ValueError('Pgen map must cover every reference and control junction')
    pool = pool.with_columns(_integer_count('event_weight').alias('event_weight'))
    if pool['event_weight'].null_count() or (pool['event_weight'] <= 0).any():
        raise ValueError('control event weights must be positive integers')
    pool = pool.sort(event_keys)
    pools = {key[0]: (part.select(event_keys), np.cumsum(part['event_weight'].to_numpy(), dtype=np.int64))
             for key, part in pool.partition_by('_bin', as_dict=True).items()}
    if pools:
        needed = reference.select('junction_aa').join(strata.filter(pl.col('junction_aa').is_in(
            reference['junction_aa'].unique().implode())), on='junction_aa', how='left').group_by('_bin').len()
        counts = dict(needed.iter_rows())
        absent = [key for key in counts if key is not None and key not in pools]
        if absent:
            required = sum(counts[key] for key in absent)
            floor = min((key for key, (_, cumulative) in pools.items()
                if cumulative[-1] >= required + counts.get(key, 0)), default=None)
            if floor is None:
                raise ValueError(f'insufficient control events for absent Pgen bins: need {required} mapped rows plus existing floor occupants')
            strata = strata.with_columns(pl.when(pl.col('_bin').is_in(absent))
                .then(pl.lit(floor)).otherwise(pl.col('_bin')).alias('_bin'))
    for _, cumulative in pools.values():
        cumulative.flags.writeable = False
    return control, pgen, _organism(species), locus, strata, event_keys, pools


def iter_mock_references(reference: pl.DataFrame, control: pl.DataFrame, *, replicates: int,
        seed: int = 0, pgen: pl.DataFrame | None = None, species='human', locus='TRB', threads=1,
        prepared_control=None):
    """Yield one replicate of null event memberships matching every reference source row.

    ``reference`` has ``group_id,junction_aa``; repeated source rows remain in
    the target counts. ``control`` is normalized by :func:`normalize_control`.
    Match rounded log2(Pgen) in fixed one-bit bins. Absent reference bins,
    including zero mass, share the lowest occupied control bin with enough
    event slots for every mapped source row plus its existing occupants.
    Sample integer event slots without replacement within each group/replicate;
    independent events sharing an AA junction retain repeated membership rows.
    Search-key deduplication belongs to the downstream scorer. Incidental
    reference/control overlaps remain eligible random matches. Insufficient
    strata fail; sampled controls are not asserted experimental nonbinders.
    ``pgen`` is an optional in-memory computed map for composition/testing.
    Internal ``prepared_control`` reuses only this execution's validated pools;
    target counts, RNG streams and sampled memberships remain fresh per call.
    The supported pool limit remains fewer than one billion events per stratum.
    Uniform event-slot sampling preserves the previous hypergeometric law but
    changes seeded realizations; provenance records ``EVENT_SAMPLER``.
    """
    if replicates < 1:
        raise ValueError('replicates must be positive')
    reference = reference.select('group_id', 'junction_aa')
    if not reference.height or any(reference[c].null_count() for c in reference.columns):
        raise ValueError('reference requires nonempty non-null group_id and junction_aa')
    if prepared_control is None:
        prepared_control = _prepare_mock_control(reference, control, pgen=pgen,
            species=species, locus=locus, threads=threads)
    elif (prepared_control[0] is not control or prepared_control[1] is not pgen or
            prepared_control[2:4] != (_organism(species), locus)):
        raise ValueError('prepared control requires the same control/Pgen objects and species/locus')
    _, _, _, _, strata, event_keys, pools = prepared_control
    requested = strata.filter(pl.col('junction_aa').is_in(reference['junction_aa'].unique().implode()))
    targets = reference.join(requested, on='junction_aa', how='left', validate='m:1')
    if targets['_bin'].null_count():
        raise ValueError('Pgen map must cover every reference and control junction')
    counts = targets.group_by('group_id', '_bin').len().sort('group_id', '_bin')
    for group, bin_id, count in counts.iter_rows():
        available = int(pools[bin_id][1][-1]) if bin_id in pools else 0
        if available < count:
            label = 'zero_mass' if bin_id == -2147483648 else str(bin_id)
            raise ValueError(f'insufficient control events for group {group!r}, Pgen bin {label}: need {count}, have {available}')
        if available >= 10**9:
            raise ValueError('event stratum exceeds supported limit: fewer than 1000000000 events required')
    rng = np.random.default_rng(seed)
    for replicate in range(replicates):
        results = []
        for group, bin_id, count in counts.iter_rows():
            events, cumulative = pools[bin_id]
            selected = events[_sample_event_indices(cumulative, count, rng)]
            results.append(selected.with_columns(pl.lit(replicate).alias('replicate'), pl.lit(group).alias('group_id')))
        yield pl.concat(results).with_columns(pl.col('group_id').cast(reference.schema['group_id'])).sort('group_id', *event_keys)


def sample_mock_references(reference: pl.DataFrame, control: pl.DataFrame, *, replicates: int,
        seed: int = 0, pgen: pl.DataFrame | None = None, species='human', locus='TRB', threads=1) -> pl.DataFrame:
    """Collect mock memberships; use ``iter_mock_references`` for bounded memory."""
    return pl.concat(iter_mock_references(reference, control, replicates=replicates,
        seed=seed, pgen=pgen, species=species, locus=locus, threads=threads))


def iter_mock_models(reference: pl.DataFrame, controls: dict[str, pl.DataFrame], *,
                     pgen: dict[str, pl.DataFrame], replicates: int, seed: int = 0):
    """Randomize a complete long reference, preserving every source-row role.

    Canonical rows have ``gene,cdr3,v,j`` and retain their epitope, MHC and pair
    linkage. Within each chain/Pgen stratum, sample event slots once without
    replacement across the entire database and uniformly assign them to source
    draw roles. Identical linked ``complex_id/gene/cdr3/v/j`` observations
    share one sampled role and are restored with their original multiplicity;
    unlinked source rows each remain independent roles. Consequently both total
    source N and every epitope's source M_i/Pgen histogram are
    preserved without reusing a control event between distinct roles. The number
    of event draws is the number of roles, not duplicate evidence source rows. TRA/TRB
    streams are independent: retained pair linkage joins random chains. Only
    sequence/gene fields change; evidence filtering must precede randomization.
    """
    required = {'gene', 'cdr3', 'v', 'j'}
    if not required.issubset(reference.columns) or not reference.height or replicates < 1:
        raise ValueError('mock models require nonempty canonical reference and positive replicates')
    if reference.filter(~pl.col('gene').is_in(['TRA', 'TRB']).fill_null(False)).height:
        raise ValueError('mock models require explicit TRA/TRB source roles')
    indexed = reference.with_row_index('_mock_row')
    if 'complex_id' in indexed.columns:
        linked = ~pl.col('complex_id').cast(pl.String).fill_null('0').is_in(['', '0'])
        role_keys = ['complex_id', 'gene', 'cdr3', 'v', 'j']
        roles = indexed.filter(linked).group_by(role_keys, maintain_order=True).agg(
            pl.col('_mock_row').min().alias('_mock_role'))
        indexed = indexed.join(roles, on=role_keys, how='left', nulls_equal=True, validate='m:1',
            maintain_order='left').with_columns(pl.col('_mock_role').fill_null(pl.col('_mock_row')))
    else:
        indexed = indexed.with_columns(pl.col('_mock_row').alias('_mock_role'))
    # Identical linked chain observations share one biological draw role. Restore
    # all original source observations after sampling; unlinked rows stay independent.
    role_source = indexed.unique('_mock_role', keep='first', maintain_order=True)
    plans = []
    streams = np.random.SeedSequence(seed).spawn(2)
    for locus, stream in zip(('TRA', 'TRB'), streams):
        source = role_source.filter(pl.col('gene') == locus)
        if not source.height:
            continue
        if locus not in controls or locus not in pgen:
            raise ValueError(f'missing {locus} control/Pgen population')
        strata = pgen_strata(pgen[locus])
        source = source.join(strata, left_on='cdr3', right_on='junction_aa', how='left', validate='m:1')
        events = controls[locus].join(strata, on='junction_aa', how='left', validate='m:1')
        if source['pgen_bin'].null_count() or events['pgen_bin'].null_count():
            raise ValueError('Pgen map must cover every reference/control event')
        pools = events.partition_by('pgen_bin', as_dict=True)
        bins = []
        for key, part in source.partition_by('pgen_bin', as_dict=True).items():
            pool = pools.get(key)
            available = int(pool['event_weight'].sum()) if pool is not None else 0
            if available < part.height:
                raise ValueError(f'insufficient {locus} control events in Pgen bin {key[0]}: need {part.height}, have {available}')
            if available >= 10**9:
                raise ValueError('event stratum exceeds supported limit: fewer than 1000000000 events required')
            pool = pool.sort('junction_aa', 'v', 'j')
            bins.append((part['_mock_role'].to_numpy(), pool,
                np.cumsum(pool['event_weight'].to_numpy(), dtype=np.int64)))
        plans.append((np.random.default_rng(stream), bins))
    for _ in range(replicates):
        assigned = []
        for rng, bins in plans:
            for rows, pool, cumulative in bins:
                selection = _sample_event_indices(cumulative, len(rows), rng)
                rng.shuffle(selection)
                assigned.append(pool[selection].select(pl.col('junction_aa').alias('cdr3'), 'v', 'j')
                    .with_columns(pl.Series('_mock_role', rows)))
        replacement = pl.concat(assigned)
        yield indexed.drop('cdr3', 'v', 'j').join(replacement, on='_mock_role', validate='m:1',
            maintain_order='left').select(reference.columns)


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
        pool = (control_events if {'v', 'j'}.issubset(fields) else normalize_control)(
            control, species=species, locus=locus)
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
            a = part[alpha_rng.permutation(part.height)]
            b = other[beta_rng.permutation(other.height)]
            columns = [c for c in ('junction_aa', 'v', 'j') if c in a.columns]
            paired.append(pl.concat([a.select('group_id', 'replicate').with_columns(
                pl.Series('draw_id', np.arange(part.height, dtype=np.int64))),
                a.select(*[pl.col(c).alias(c + '_alpha') for c in columns]),
                b.select(*[pl.col(c).alias(c + '_beta') for c in columns])], how='horizontal'))
        yield pl.concat(paired).sort('group_id', 'draw_id')


def iter_epitope_mocks(reference: pl.DataFrame, control: pl.DataFrame, *,
                       replicates: int, seed: int = 0, policy: str = 'pgen', unavailable: list | None = None):
    """Yield Pgen-matched replacement reference corpora, retaining source roles.

    Inputs contain junction_aa/v/j/pgen_bin; reference also has group_id. Rows
    represent distinct receptor keys, not reads or weighted recombination events.
    Known target keys are excluded. ``pgen-vj`` additionally preserves V/J;
    ``pgen-different-v`` excludes every target junction and changes each source V.
    The latter samples source-V strata in fresh random order without replacement
    within a target/bin; it is a sequential conditional law, not a uniform draw
    over all feasible constrained assignments. Infeasible draws fail explicitly.
    No bin widening, motif-neighbour rejection, or persistent memberships.
    By default unavailable targets raise. With an ``unavailable`` audit list,
    record their group/reason and omit them from draws; callers must mark their
    entire null inference unavailable, including any earlier successful draws.
    Extra reference columns (including pairing roles) follow the replaced row.
    """
    import hashlib

    if policy not in {'pgen', 'pgen-vj', 'pgen-different-v'}:
        raise ValueError('unknown mock reference policy')
    if type(replicates) is not int or replicates < 1 or type(seed) is not int or seed < 0:
        raise ValueError('positive replicates and nonnegative seed required')
    keys = ['junction_aa', 'v', 'j']
    required = [*keys, 'pgen_bin']
    for frame, columns in [(reference, ['group_id', *required]), (control, required)]:
        if not set(columns) <= set(frame.columns) or not frame.height:
            raise ValueError('nonempty reference/control with receptor keys and Pgen bins required')
        if any(frame[c].null_count() for c in columns) or not frame.schema['pgen_bin'].is_integer():
            raise ValueError('non-null receptor keys and integer Pgen bins required')
        if frame.select(pl.any_horizontal(pl.col(c).cast(pl.String).str.strip_chars() == '' for c in keys).any()).item():
            raise ValueError('receptor keys cannot be empty')
    if control.select(keys).n_unique() != control.height:
        raise ValueError('control receptor keys must be unique')
    # Source roles may repeat a chain key across different linked receptors.
    source = reference.with_row_index('_source_row')
    pool = control.select(required).sort(keys).with_row_index('_pool_row')
    pool_index_dtype = pool['_pool_row'].to_numpy().dtype
    def partitions(frame,columns):
        grouped=frame.group_by(columns).agg(pl.col('_pool_row').sort())
        return {key:values.to_numpy() for key,values in zip(grouped.select(columns).iter_rows(),grouped['_pool_row'])}
    def excluded_positions(eligible,blocked):
        positions=np.searchsorted(eligible,blocked)
        valid=positions<len(eligible)
        positions,blocked=positions[valid],blocked[valid]
        return positions[eligible[positions]==blocked]
    # Partition raw key indices once, rather than scanning the full pool per target.
    bin_pools=partitions(pool,['pgen_bin'])
    strata=['pgen_bin','v','j'] if policy=='pgen-vj' else ['pgen_bin','v']
    conditional=partitions(pool,strata) if policy!='pgen' else {}
    empty_indices=np.empty(0,dtype=np.uint32)
    if policy=='pgen-different-v':
        conditional={role:np.setdiff1d(bin_pools.get((role[0],),empty_indices),
            conditional.get(role,empty_indices),assume_unique=True)
            for role in source.select('pgen_bin','v').unique().iter_rows()}
    exclusion_keys=['junction_aa'] if policy=='pgen-different-v' else keys
    excluded=partitions(source.select('group_id',*exclusion_keys).unique().join(
        pool.select(*exclusion_keys,'pgen_bin','_pool_row'),on=exclusion_keys,how='inner'),['group_id','pgen_bin'])
    plans=[];rejected=set()
    def reject(group,reason):
        if unavailable is None:raise ValueError(reason)
        rejected.add(group)
        unavailable.append(dict(group_id=group,reason=reason))
    for (group,),target in source.partition_by('group_id',as_dict=True).items():
        groups=[]
        for (bin_id,),part in target.partition_by('pgen_bin',as_dict=True).items():
            candidates=bin_pools.get((bin_id,),empty_indices)
            blocked=excluded.get((group,bin_id),empty_indices)
            if len(candidates)-len(excluded_positions(candidates,blocked))<part.height:
                reject(group,f'insufficient control keys for group {group!r}, bin {bin_id}')
                break
            roles=[]
            strata=['v','j'] if policy=='pgen-vj' else ['v'] if policy=='pgen-different-v' else ['pgen_bin']
            for role,members in part.partition_by(strata,as_dict=True).items():
                eligible=candidates if policy=='pgen' else conditional.get((bin_id,*role),empty_indices)
                blocked_positions=excluded_positions(eligible,blocked)
                available=len(eligible)-len(blocked_positions)
                if available<members.height:
                    reject(group,f'insufficient control keys for group {group!r}, bin {bin_id}, role {role}: need {members.height}, have {available}')
                    break
                roles.append((members['_source_row'].to_numpy(),eligible,blocked_positions))
            if group in rejected:break
            groups.append(roles)
        identity=int.from_bytes(hashlib.sha256(str(group).encode()).digest()[:8],'little')
        if group not in rejected:plans.append((group,identity,groups))
    for replicate in range(replicates):
        selected = np.zeros(source.height, dtype=np.int64)
        for group,identity,groups in plans:
            if group in rejected:continue
            rng = np.random.default_rng(np.random.SeedSequence([seed, identity, replicate]))
            for roles in groups:
                if group in rejected:break
                # Match pool indices: mixed uint32/int64 searchsorted copies
                # the complete control stratum even for an empty needle array.
                used = np.empty(0, dtype=pool_index_dtype)
                order = rng.permutation(len(roles)) if policy == 'pgen-different-v' else range(len(roles))
                for i in order:
                    positions,eligible,blocked_positions=roles[i]
                    blocked=np.union1d(blocked_positions,excluded_positions(eligible,used))
                    remaining=len(eligible)-len(blocked)
                    if remaining < len(positions):
                        reject(group,f'insufficient eligible keys in sequential different-V draw {replicate} for group {group!r}; enlarge the control pool')
                        break
                    ranks=rng.choice(remaining,len(positions),replace=False)
                    # Map uniform ranks in the complement without materialising it.
                    offsets=np.searchsorted(blocked-np.arange(len(blocked)),ranks,side='right')
                    draw=eligible[ranks+offsets]
                    selected[positions] = draw
                    used = np.concatenate([used, draw])
        replacement = pool[selected].select(keys)
        yield reference.with_columns(*(replacement[c] for c in keys), pl.lit(replicate).alias('replicate')).filter(~pl.col('group_id').is_in(list(rejected)))
