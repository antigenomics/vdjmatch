"""Prepared reference assets: normalized observations, native indexes and profiles."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import polars as pl
from seqtree import Index

from .cache import locked, publish, sha256, staging, valid


def _identity(source, *, asset, species, filters):
    from importlib import import_module
    from importlib.metadata import version
    import seqtree._core as native
    from ..aggregate import fast_profiles, structure, profile_noise
    from ..cli import epitope
    from ..match import engine
    from ..evalue import epitope_rules
    from . import schema, vdjdb
    overlap = import_module('vdjmatch.cluster.overlap')

    if not all(callable(getattr(Index, name, None)) for name in ('build', 'save', 'load')):
        raise RuntimeError('prepared references require native seqtree index save/load support')
    validation = [source.parent/'release.json']
    if source.is_file():
        validation += [source.with_name(source.name+'.sha256'), source.parent/'SHA256SUMS']
    return dict(format='vdjmatch-reference-index-v1', source=vdjdb.provenance(source),
        validation_files={p.name: sha256(p) for p in validation if p.is_file()},
        asset=asset, species=species, filters=json.loads(json.dumps(
            {k: v for k, v in filters.items() if v is not None}, default=os.fspath)),
        profile_recipe=fast_profiles.PROFILE_RECIPE,
        source_sha256={name: sha256(Path(module.__file__)) for name, module in
            [('normalizer', schema), ('loader', vdjdb), ('engine', engine), ('profiles', fast_profiles),
             ('profile_genes', structure), ('profile_populations', epitope), ('profile_graph', overlap),
             ('profile_noise', profile_noise), ('profile_tokens', epitope_rules)]} |
            {'bundle': sha256(Path(__file__))},
        numerics={name: version(name) for name in ['numpy', 'scipy', 'polars']},
        native=dict(version=version('seqtree'), sha256=sha256(Path(native.__file__)),
            capabilities=['Index.build', 'Index.save', 'Index.load']))


def _restore(directory, identity):
    from ..match.engine import VdjdbIndex

    manifest_path = directory/'manifest.json'
    try:
        if not valid(manifest_path):
            return None
        manifest = json.loads(manifest_path.read_text())
        genes = manifest['genes']
        if (manifest['identity'] != identity or genes != sorted(set(genes)) or
                not set(genes) <= {'TRA', 'TRB'} or manifest['artifact_path'] != str(directory)):
            return None
        expected = {'normalized-reference.parquet', 'profiles.parquet', 'profiles.tsv', 'joint-profiles.tsv'} | {
            gene+suffix for gene in genes for suffix in ('.seqtree', '.keys.parquet')}
        if set(manifest['files']) != expected or not all(
                valid(directory/name, manifest['files'][name]['sha256']) for name in expected):
            return None
        reference = pl.read_parquet(directory/'normalized-reference.parquet')
        profiles = pl.read_parquet(directory/'profiles.parquet')
        records = reference.filter(pl.col('reference_valid'))
        if 'record_id' not in records.columns:
            records = records.with_row_index('record_id').with_columns(pl.col('record_id').cast(pl.String))
        if sorted(records['gene'].drop_nulls().unique().to_list()) != genes:
            return None
        by_gene = {}
        for gene in genes:
            part = records.filter(pl.col('gene') == gene)
            keys = pl.read_parquet(directory/(gene+'.keys.parquet'))
            expected_keys = part.select('cdr3').unique(maintain_order=True).with_row_index('ref_id')
            native = Index.load(str(directory/(gene+'.seqtree')))
            if not keys.equals(expected_keys) or len(native) != keys.height:
                return None
            by_gene[gene] = native, keys, part
        if reference.height != manifest['reference_rows'] or profiles.height != manifest['profile_rows']:
            return None
        return VdjdbIndex(by_gene, reference), profiles, manifest
    except Exception:
        # Native format and parquet readers use different exception classes for
        # incomplete/corrupt assets; all failed loads rebuild from the raw source.
        return None


def load_or_build_index(source, *, asset='default', species=None, threads=0, **filters):
    """Return ``(VdjdbIndex, flat_profiles, manifest)`` from a verified prepared bundle.

    Bundles live beside the source, including directory sources, so a bundle
    never enters its own source identity. Every input file is content-hashed;
    unchanged assets bypass ingestion, native building and classification.
    Observation identities, pairing and conflicting labels remain intact.
    ``filters`` are forwarded unchanged to :func:`vdjmatch.db.load`.
    """
    from .. import db
    from ..aggregate.fast_profiles import reference_profiles, joint_reference_profiles
    from ..match.engine import VdjdbIndex

    from vdjtools.cores import kernel_threads
    threads = kernel_threads(threads)
    source = Path(source).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(source)
    identity = _identity(source, asset=asset, species=species, filters=filters)
    encoded = json.dumps(identity, sort_keys=True, separators=(',', ':'))
    digest = hashlib.sha256(encoded.encode()).hexdigest()
    root = source.with_name(source.name+'.vdjmatch')
    root.mkdir(parents=True, exist_ok=True)
    directory = root/digest
    with locked(directory):
        if _identity(source, asset=asset, species=species, filters=filters) != identity:
            raise ValueError('reference source or implementation changed while waiting for its lock')
        restored = _restore(directory, identity)
        if restored is not None:
            return restored
        directory.mkdir(exist_ok=True)
        (directory/'manifest.json').unlink(missing_ok=True)
        reference = db.load(source, asset=asset, species=species, **filters)
        index = VdjdbIndex.build(reference)
        reference = index.reference
        profiles = reference_profiles(reference, threads=threads)
        files = {}
        def write(name, writer):
            target = directory/name
            with staging(target) as temporary:
                writer(temporary)
                publish(temporary, target)
            files[name] = dict(sha256=target.with_name(target.name+'.sha256').read_text().strip())
        write('normalized-reference.parquet', reference.write_parquet)
        write('profiles.parquet', profiles.write_parquet)
        write('profiles.tsv', lambda path: profiles.write_csv(path, separator='\t'))
        write('joint-profiles.tsv', lambda path: joint_reference_profiles(profiles).write_csv(path, separator='\t'))
        for gene, (native, keys, _) in index._by_gene.items():
            write(gene+'.seqtree', lambda path, native=native: native.save(str(path)))
            write(gene+'.keys.parquet', keys.write_parquet)
        if _identity(source, asset=asset, species=species, filters=filters) != identity:
            raise ValueError('reference source or implementation changed during preparation; bundle not committed')
        manifest = dict(identity=identity, artifact_path=str(directory), files=files,
            genes=index.genes, reference_rows=reference.height, profile_rows=profiles.height,
            threads=threads)
        with staging(directory/'manifest.json') as temporary:
            temporary.write_text(json.dumps(manifest, indent=2)+'\n')
            publish(temporary, directory/'manifest.json')
        return index, profiles, manifest
