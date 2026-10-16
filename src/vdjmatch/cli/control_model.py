"""Prepare and validate explicitly stored observed-control diagnostic models."""
from pathlib import Path
import json


def register(subparsers):
    p = subparsers.add_parser('prepare-control', help=__doc__)
    p.add_argument('--control', required=True)
    p.add_argument('--species', default='human')
    p.add_argument('--locus', choices=['TRA', 'TRB'], required=True)
    p.add_argument('--output-dir', required=True, help='new directory; existing paths are refused')
    modes = ['exact3', 'imgt3', 'wild4', 'imgt_wild4']
    p.add_argument('--motif-mode', nargs='+', choices=modes, default=modes)
    p.add_argument('--threads', type=int, default=1)
    p.set_defaults(func=main)


def main(a):
    import numpy as np
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from importlib.metadata import version
    from .combined import _raw_control
    from ..evalue.mock_reference import control_events
    from ..aggregate.structure import control_gene_usage
    from ..aggregate.motif_information import motif_population
    from ..match.motifs import masked_motif_presence_counts
    from ..evalue.epitope_rules import _tokens

    directory = Path(a.output_dir).resolve()
    modes = list(a.motif_mode)
    if directory.exists() or a.threads < 1 or len(set(modes)) != len(modes):
        raise ValueError('prepare-control requires a new output directory, positive threads and distinct modes')
    organism = _organism(a.species)
    source_functions = (_raw_control, control_events, control_gene_usage, motif_population, masked_motif_presence_counts, _tokens)
    sources_before = {fn.__name__: sha256(Path(fn.__code__.co_filename)) for fn in source_functions}
    producer_before, raw_before = sha256(Path(__file__)), sha256(Path(a.control))
    directory.mkdir(parents=True, exist_ok=False)
    raw, info = _raw_control(a.control)
    events = control_events(raw, species=organism, locus=a.locus)
    del raw
    files = {}
    def write(name, frame):
        path = directory/name
        frame.write_parquet(path)
        files[name] = dict(sha256=sha256(path), rows=frame.height,
            schema={c: str(t) for c, t in frame.schema.items()})
    write('events.parquet', events)
    write('genes.parquet', control_gene_usage(events, locus=a.locus))
    for mode in modes:
        population = motif_population(events.rename({'junction_aa': 'cdr3'}), locus=a.locus, mode=mode)
        write(mode+'.vocabulary.parquet', population['vocabulary'])
        arrays = {}
        for scope in ('exact_length', 'all_lengths'):
            strata, arrays[scope] = population[scope]
            write(mode+'.'+scope+'.parquet', strata)
        path = directory/(mode+'.counts.npz')
        np.savez_compressed(path, **arrays)
        files[path.name] = dict(sha256=sha256(path), shapes={k: list(v.shape) for k, v in arrays.items()})
        del population, arrays
    raw_after = sha256(Path(a.control))
    sources_after = {fn.__name__: sha256(Path(fn.__code__.co_filename)) for fn in source_functions}
    producer_after = sha256(Path(__file__))
    if raw_before != raw_after or sources_before != sources_after or producer_before != producer_after:
        raise ValueError('control input or producing source changed during preparation; no complete manifest written')
    info.update(path=str(Path(a.control).resolve()), sha256=raw_before,
        event_keys=events.height, total_events=int(events['event_weight'].sum()),
        reader_sha256=sha256(Path(_raw_control.__code__.co_filename)),
        event_normalizer_sha256=sha256(Path(control_events.__code__.co_filename)))
    manifest = dict(format='vdjmatch-control-model-v1', species=organism, locus=a.locus, modes=modes,
        raw=dict(path=info['path'], sha256=raw_before, sha256_before=raw_before, sha256_after=raw_after), control_info=info, files=files,
        producer=dict(sha256=producer_before, sha256_after=producer_after,
            sources_before=sources_before, sources_after=sources_after,
            software={name: version(name) for name in ('vdjmatch', 'vdjtools', 'polars', 'numpy')}))
    (directory/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0


def load_control_model(directory, raw_path, *, species, locus, modes):
    """Load only requested stored views; ``modes=()`` loads events alone.

    Return events, gene usage (None for events-only), motif populations and the
    original control report plus this model's identity. Never rebuild a view.
    """
    from zipfile import ZipFile
    import numpy as np
    import polars as pl
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..evalue.epitope_rules import MOTIFS
    directory = Path(directory).resolve()
    manifest_path = directory/'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if (not isinstance(manifest, dict) or manifest.get('format') != 'vdjmatch-control-model-v1' or
            not isinstance(manifest.get('raw'), dict) or not isinstance(manifest.get('modes'), list) or
            not isinstance(manifest.get('producer'), dict) or
            manifest.get('species') != _organism(species) or manifest.get('locus') != locus or
            not set(modes) <= set(MOTIFS) or not set(modes) <= set(manifest.get('modes', [])) or
            not isinstance(manifest.get('files'), dict) or not isinstance(manifest.get('control_info'), dict) or
            manifest['raw'].get('sha256') != sha256(Path(raw_path)) or
            manifest['raw'].get('sha256_before') != manifest['raw'].get('sha256_after') or
            manifest['raw'].get('sha256') != manifest['control_info'].get('sha256') or
            not isinstance(manifest['producer'].get('sources_before'), dict) or
            not manifest['producer']['sources_before'] or
            manifest['producer']['sources_before'] != manifest['producer'].get('sources_after') or
            manifest['producer'].get('sha256') != manifest['producer'].get('sha256_after')):
        raise ValueError('control model manifest, raw identity, species/locus or requested modes do not match')
    files = manifest['files']
    def checked(name):
        path = (directory/name).resolve()
        entry = files.get(name)
        if path.parent != directory or not isinstance(entry, dict) or entry.get('sha256') != sha256(path):
            raise ValueError('control model file identity mismatch: '+name)
        return path
    def table(name, columns):
        frame = pl.read_parquet(checked(name))
        entry = files[name]
        if (set(frame.columns) != set(columns) or
                (isinstance(columns, dict) and dict(frame.schema) != columns) or frame.height != entry.get('rows') or
                {c: str(t) for c, t in frame.schema.items()} != entry.get('schema')):
            raise ValueError('control model table schema/count mismatch: '+name)
        return frame
    events = table('events.parquet', dict(junction_aa=pl.String, v=pl.String, j=pl.String,
        event_weight=pl.Int64, species=pl.String, locus=pl.String))
    if (not events.height or not events.schema['event_weight'].is_integer() or
            any(events[c].null_count() for c in events.columns) or
            events.filter((pl.col('event_weight') <= 0) | (pl.col('species') != _organism(species)) |
                (pl.col('locus') != locus) | ~pl.col('junction_aa').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$')).height or
            events.select('junction_aa', 'v', 'j').n_unique() != events.height or
            events.height != manifest['control_info'].get('event_keys') or
            int(events['event_weight'].sum()) != manifest['control_info'].get('total_events')):
        raise ValueError('control model events/counts do not match the declared population')
    info = dict(manifest['control_info'], path=str(Path(raw_path).resolve()),
        control_model=dict(path=str(directory), manifest_sha256=sha256(manifest_path)))
    if not modes:
        return events, None, {}, info
    usage = table('genes.parquet', ('segment', 'length_scope', 'junction_length', 'gene',
        'control_junction_gene_keys', 'total_control_junction_gene_keys', 'control_events',
        'total_control_events', 'control_event_fraction', 'control_key_fraction'))
    populations = {}
    for mode in modes:
        vocabulary = table(mode+'.vocabulary.parquet', dict(_token_id=pl.UInt32, t0=pl.String))
        width = {'exact3': 8000, 'imgt3': 1331, 'wild4': 32000, 'imgt_wild4': 5324}[mode]
        if (vocabulary.height != width or vocabulary['t0'].null_count() or vocabulary['t0'].n_unique() != width or
                not np.array_equal(vocabulary['_token_id'].to_numpy(), np.arange(width))):
            raise ValueError('control model vocabulary shape/identity mismatch: '+mode)
        path = checked(mode+'.counts.npz')
        with ZipFile(path) as archive:
            if (set(archive.namelist()) != {'exact_length.npy', 'all_lengths.npy'} or
                    any(item.file_size > 1024**3+4096 for item in archive.infolist())):
                raise ValueError('control model count archive exceeds the bounded array schema')
        population = dict(vocabulary=vocabulary, mode=mode)
        with np.load(path, allow_pickle=False) as arrays:
            for scope in ('exact_length', 'all_lengths'):
                strata = table(mode+'.'+scope+'.parquet', dict(_stratum=pl.UInt32, v_gene=pl.String,
                    junction_length=pl.UInt32, keys=pl.UInt32))
                counts = arrays[scope]
                shape = (strata.height, width)
                if (counts.dtype != np.int64 or counts.shape != shape or counts.nbytes > 1024**3 or
                        list(shape) != files[path.name].get('shapes', {}).get(scope) or counts.min(initial=0) < 0 or
                        not np.array_equal(strata['_stratum'].to_numpy(), np.arange(strata.height)) or
                        strata['keys'].null_count() or (strata['keys'] <= 0).any() or
                        np.any(counts.max(axis=1) > strata['keys'].to_numpy())):
                    raise ValueError('control model count array/stratum mismatch: '+mode+'/'+scope)
                population[scope] = strata, counts
        populations[mode] = population
    return events, usage, populations, info
