"""Marginal full-junction generation probability from an explicit raw model."""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('pgen', help='marginal full-junction generation probability')
    p.add_argument('sample', help='AIRR junction table, or a local reference with --vdjdb')
    p.add_argument('--vdjdb', action='store_true', help='load SAMPLE as a local VDJdb reference table, ZIP or directory')
    p.add_argument('--control', help='native vdjtools control AA table; score the distinct sample/control union and write .features.parquet')
    p.add_argument('--model-path', required=True, help='explicit native vdjtools raw model directory')
    p.add_argument('--species', required=True, help='human or mouse (standard species aliases accepted)')
    p.add_argument('--locus', choices=['TRA', 'TRB'], required=True)
    p.add_argument('--mismatches', type=int, choices=[0, 1], default=0,
                   help='inclusive AA Hamming radius, without indels (default: 0)')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--output-prefix', required=True)
    p.set_defaults(func=main)


def main(a):
    import hashlib
    from importlib.metadata import version
    import json
    import time
    import numpy as np
    import polars as pl
    from vdjtools import _core
    from vdjtools.model import load_model, native
    from vdjtools.model import io as model_io
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..io.airr import _read_table
    from ..io.columns import normalize_query, _resolve

    start = time.perf_counter()
    if a.threads < 1:
        raise ValueError('positive threads required')
    from ..runner.annotation import _validate_output_paths
    prefix = _validate_output_paths(a.output_prefix, [a.sample, a.model_path] + ([a.control] if a.control else []))
    organism = _organism(a.species)
    model_path = Path(a.model_path)
    if not model_path.is_dir():
        raise ValueError('--model-path must be an explicit raw model directory')
    # load_model retains all alleles; load_bundled would collapse them by default.
    model = load_model(model_path, validate=True)
    if _organism(model.manifest.organism) != organism:
        raise ValueError('model species does not match --species')
    if model.manifest.locus != a.locus:
        raise ValueError('model locus does not match --locus')
    # Pin the files actually consumed by load_model, including format priority.
    names = [*model.manifest.events, 'genes_v', 'genes_j']
    if model.manifest.chain_type == 'VDJ':
        names.append('genes_d')
    paths = [model_path / 'manifest.json']
    for name in names:
        paths.append(next(model_path / (name + suffix) for suffix in ('.parquet', '.tsv', '.csv')
                          if (model_path / (name + suffix)).is_file()))
    if (model_path / 'training.json').is_file():
        paths.append(model_path / 'training.json')
    hashes = {path.name: sha256(path) for path in sorted(paths)}
    model_hash = hashlib.sha256(''.join(f'{name}\t{digest}\n' for name, digest in hashes.items()).encode()).hexdigest()

    reference = None
    if a.vdjdb:
        from .. import db
        raw = db.load(a.sample, species={'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[organism], gene=a.locus)
        reference = {**db.provenance(a.sample),
                     'loader_sha256': sha256(Path(db.load.__code__.co_filename)),
                     'normalizer_sha256': sha256(Path(db.schema.normalize.__code__.co_filename))}
        raw = raw.with_columns(pl.col('cdr3').alias('junction_aa'), pl.col('v').alias('v_call'),
                               pl.col('j').alias('j_call'), pl.col('gene').alias('locus'))
    else:
        raw = _read_table(a.sample)
    if 'query_id' not in raw.columns:
        raw = raw.with_row_index('query_id')
    if raw['query_id'].null_count() or raw['query_id'].n_unique() != raw.height:
        raise ValueError('query_id must be non-null and unique')
    generated = {'junction_aa_normalized', 'pgen', 'log10_pgen', 'available', 'status', '_pgen_locus'}
    if generated.intersection(raw.columns):
        raise ValueError('input contains reserved Pgen output columns')
    normalized, ingestion = normalize_query(raw, valid_aa=False, source='airr', return_report=True)
    fields = _resolve(raw)
    if 'locus' in fields:
        locus = pl.col(fields['locus']).cast(pl.String).str.strip_chars().str.to_uppercase()
    elif 'v' in fields:
        v_prefix = pl.col(fields['v']).cast(pl.String).str.slice(0, 3).str.to_uppercase()
        locus = pl.when(v_prefix.is_in(['TRA', 'TRB'])).then(v_prefix).otherwise(None)
    else:
        locus = pl.lit(None, dtype=pl.String)
    raw = raw.with_columns(locus.alias('_pgen_locus'))
    rows = raw.join(normalized.select('query_id', pl.col('cdr3').alias('junction_aa_normalized')),
                    on='query_id', how='left', validate='1:1', maintain_order='left')
    species_ok = (pl.col('species').str.strip_chars().str.to_lowercase().is_in(
        ['human', 'homosapiens', 'homo sapiens'] if organism == 'human' else
        ['mouse', 'musmusculus', 'mus musculus']) if 'species' in rows.columns else pl.lit(True))
    # --locus supplies absent locus metadata, without requiring a V call.
    locus_ok = pl.col('_pgen_locus').fill_null(a.locus) == a.locus
    rows = rows.with_columns(
        pl.when(~species_ok.fill_null(False)).then(pl.lit('unselected_species'))
        .when(~locus_ok).then(pl.lit('unselected_locus'))
        .when(pl.col('junction_aa_normalized').is_null()).then(pl.lit('missing_junction'))
        .when(~pl.col('junction_aa_normalized').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$'))
        .then(pl.lit('invalid_junction')).otherwise(pl.lit('scored')).alias('status'))
    selected = rows.filter(pl.col('status') == 'scored')
    sequences = selected.select('junction_aa_normalized')
    control = None
    if a.control:
        from vdjtools.io import read_vdjtools
        control_rows = read_vdjtools(a.control)
        # Pooled AA gene calls are descriptive; these keys are scored under
        # the explicitly selected model, without conditioning on V/J.
        j_locus, v_locus = pl.col('j_call').str.slice(0, 3), pl.col('locus')
        other_chain = control_rows.filter(v_locus.is_in(['TRA', 'TRB']) &
            (v_locus == j_locus) & (j_locus != a.locus)).height
        conflicts = control_rows.filter(v_locus.is_in(['TRA', 'TRB']) &
            j_locus.is_in(['TRA', 'TRB']) & (v_locus != j_locus)).height
        control_rows = control_rows.select(
            pl.col('junction_aa').str.strip_chars().str.to_uppercase().alias('junction_aa_normalized'))
        canonical = control_rows.filter(pl.col('junction_aa_normalized').str.contains(
            r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False))
        control = dict(path=str(Path(a.control).resolve()), sha256=sha256(Path(a.control)),
            source_rows=control_rows.height, canonical_rows=canonical.height,
            excluded_rows=control_rows.height-canonical.height,
            conflicting_gene_rows=conflicts,
            other_chain_first_call_rows=other_chain,
            reader_sha256=sha256(Path(read_vdjtools.__code__.co_filename)),
            unit='AA junction keys only; control read abundance and AA incidence are not probability weights')
        sequences = pl.concat([sequences, canonical])
        del control_rows, canonical
    sequences = sequences.unique().sort('junction_aa_normalized')
    print(f'pgen: {sequences.height} distinct {organism} {a.locus} junctions, Hamming <= {a.mismatches}', flush=True)
    native_start = time.perf_counter()
    probabilities = np.asarray(native.pgen_aa_batch(model, sequences['junction_aa_normalized'].to_list(),
                               mismatches=a.mismatches, threads=a.threads) if sequences.height else [], dtype=np.float64)
    native_wall = time.perf_counter() - native_start
    if probabilities.shape != (sequences.height,) or not np.all(
            np.isfinite(probabilities) & (probabilities >= 0) & (probabilities <= 1)):
        raise RuntimeError('native Pgen returned an invalid probability; no scores written')
    scored = selected.select('query_id', 'junction_aa_normalized').join(
        sequences.with_columns(pl.Series('pgen', probabilities)), on='junction_aa_normalized',
        how='left', validate='m:1').select('query_id', 'pgen')
    rows = rows.join(scored, on='query_id', how='left', validate='1:1', maintain_order='left').with_columns(
        pl.when(pl.col('pgen') == 0).then(pl.lit('zero_mass')).otherwise(pl.col('status')).alias('status'),
        pl.when(pl.col('pgen') > 0).then(pl.col('pgen').log10()).otherwise(None).alias('log10_pgen'),
        pl.col('pgen').is_not_null().alias('available')).drop('_pgen_locus')
    prefix.parent.mkdir(parents=True, exist_ok=True)
    from .__main__ import _flat_table
    _flat_table(rows).write_csv(str(prefix) + '.scores.tsv', separator='\t')
    features = None
    if a.control:
        feature_path = Path(str(prefix) + '.features.parquet')
        sequences.rename({'junction_aa_normalized': 'junction_aa'}).with_columns(
            pl.Series('pgen', probabilities), pl.lit(organism).alias('species'),
            pl.lit(a.locus).alias('locus')).write_parquet(feature_path)
        features = dict(path=feature_path.name, sha256=sha256(feature_path))
    manifest = {
        'estimator': 'marginal-full-junction-pgen-v1',
        'probability': f'P(Hamming(junction_aa, generated_junction_aa) <= {a.mismatches}); '
                       'same length, no indels; summed over all model V/J and recombination scenarios',
        'sequence_contract': 'full anchor-inclusive junction; AIRR junction_aa or normalized VDJdb cdr3; no anchor fabrication',
        'input_format': 'vdjdb' if a.vdjdb else 'airr',
        'reference': reference,
        'control': control, 'features': features,
        'sample': {'sha256': sha256(Path(a.sample))} if Path(a.sample).is_file() else reference,
        'model': {'sha256': model_hash, 'files': hashes, 'metadata': json.loads(model.manifest.to_json()),
                  'collapse_alleles': False, 'anchor_repair': 'vdjtools load_model'},
        'software': {'vdjmatch': version('vdjmatch'), 'vdjtools': version('vdjtools'),
                     'native_sha256': sha256(Path(_core.__file__)), 'cli_sha256': sha256(Path(__file__)),
                     'model_io_sha256': sha256(Path(model_io.__file__)), 'native_bridge_sha256': sha256(Path(native.__file__)),
                     'formatter_sha256': sha256(Path(_flat_table.__code__.co_filename))},
        'options': {'species': organism, 'locus': a.locus, 'threads': a.threads,
                    'mismatches': a.mismatches,
                    'model_path': str(model_path.resolve()), 'sample': str(Path(a.sample).resolve()), 'vdjdb': a.vdjdb},
        'counts': {'input_rows': raw.height, 'selected_rows': selected.height, 'native_rows': sequences.height,
                   'dispositions': dict(rows.group_by('status').len().iter_rows())},
        'ingestion': ingestion, 'timing_seconds': {'native_batch': native_wall, 'total': time.perf_counter() - start},
    }
    Path(str(prefix) + '.manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return 0
