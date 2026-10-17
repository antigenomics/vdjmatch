"""Rebuild reference classification tables freshly from raw reference observations."""
from pathlib import Path


def register(subparsers):
    parser = subparsers.add_parser('reference-profiles', help=__doc__)
    parser.add_argument('--vdjdb', required=True)
    parser.add_argument('--species', choices=['human', 'mouse', 'HomoSapiens', 'MusMusculus'])
    parser.add_argument('--reference-population', choices=['curated-benchmark', 'all-valid'], default='curated-benchmark')
    parser.add_argument('--shuffles', type=int, default=999, help='label-shuffle null replicates (default: 999)')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--output-prefix', required=True)
    parser.set_defaults(func=main)


def main(a):
    import json
    import time
    from importlib.metadata import version
    import polars as pl
    import seqtree._core as native
    from .. import db
    from ..aggregate.fast_profiles import reference_profiles, joint_reference_profiles, PROFILE_RECIPE
    from ..db.cache import sha256
    from ..runner.annotation import _validate_output_paths

    started = time.monotonic()
    if a.threads < 1:
        raise ValueError('positive threads required')
    if a.shuffles < 200:
        raise ValueError('at least 200 label shuffles required for scientific profiles')
    prefix = _validate_output_paths(a.output_prefix, [a.vdjdb])
    def freeze():
        paths = sorted(Path(__file__).resolve().parents[1].rglob('*.py'))+[Path(native.__file__)]
        return {str(path): sha256(path) for path in paths}
    sources = freeze()
    identity = db.provenance(a.vdjdb)
    species = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}.get(a.species, a.species)
    reference = db.load(a.vdjdb, species=species)
    if a.reference_population == 'all-valid':
        reference = reference.with_columns(pl.col('cdr3').str.contains(
            r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).alias('reference_valid'))
    profiles = reference_profiles(reference, threads=a.threads, shuffles=a.shuffles)
    joint = joint_reference_profiles(profiles)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for suffix, frame in [('profiles', profiles), ('joint-profiles', joint)]:
        path = Path(str(prefix)+'.'+suffix+'.tsv')
        frame.write_csv(path, separator='\t')
        outputs[path.name] = sha256(path)
    after = freeze()
    if identity != db.provenance(a.vdjdb) or sources != after:
        raise ValueError('reference or producing source/native changed during fresh profiling')
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(dict(reference=identity,
        reference_population=a.reference_population, species=species, profile_recipe=PROFILE_RECIPE.replace('999-label-shuffles', f'{a.shuffles}-label-shuffles'), shuffles=a.shuffles,
        source_before=sources, source_after=after, source_frozen=True, reference_rows=reference.height,
        profile_rows=profiles.height, joint_rows=joint.height, threads=a.threads,
        calculation='fresh source observations, reference graphs and label shuffles; no saved diagnostic rows',
        software={name: version(name) for name in ('vdjmatch', 'seqtree', 'polars', 'numpy', 'scipy')},
        tables=outputs, seconds=time.monotonic()-started), indent=2)+'\n')
    return 0
