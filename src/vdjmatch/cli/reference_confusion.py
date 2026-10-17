"""Summarize fresh benchmark winning labels against annotated reference truth."""
from pathlib import Path


def register(subparsers):
    parser = subparsers.add_parser('score-confusion', help=__doc__)
    parser.add_argument('scores', help='long TSV: species,locus,query_id,epitope,scorer,score,supported')
    parser.add_argument('--truth', required=True, help='TSV of all known query/epitope memberships')
    parser.add_argument('--profiles', required=True, help='pooled reference profiles with source-record eligibility')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--output-prefix', required=True)
    parser.set_defaults(func=main)


def main(a):
    import json
    import time
    import polars as pl
    from ..aggregate.reference_confusion import reference_confusion
    from ..db.cache import sha256
    from ..runner.annotation import _validate_output_paths

    started = time.monotonic()
    paths = [Path(path) for path in (a.scores, a.truth, a.profiles)]
    prefix = _validate_output_paths(a.output_prefix, paths)
    identities = {name: sha256(path) for name, path in zip(('scores', 'truth', 'profiles'), paths)}
    frames = [pl.read_csv(path, separator='\t', infer_schema_length=None,
        schema_overrides={'species': pl.String, 'locus': pl.String, 'epitope': pl.String,
            **({'query_id': pl.String} if i < 2 else {})}) for i, path in enumerate(paths)]
    metrics, calls = reference_confusion(*frames)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for suffix, table in [('confusion', metrics), ('calls', calls)]:
        path = Path(str(prefix)+'.'+suffix+'.tsv')
        table.write_csv(path, separator='\t')
        outputs[path.name] = sha256(path)
    if identities != {name: sha256(path) for name, path in zip(('scores', 'truth', 'profiles'), paths)}:
        raise ValueError('confusion inputs changed during reduction')
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(dict(inputs=identities, tables=outputs,
        source_sha256=sha256(Path(__file__)), reducer_sha256=sha256(Path(reference_confusion.__code__.co_filename)),
        population='complete explicit scorer query dispositions; eligible targets have >=30 valid source records',
        decision='all supported finite maximum-score labels win; exact ties retained; explicit no-hit rows retained',
        interpretation='one-vs-other annotation retrieval; other-class labels do not demonstrate experimental nonbinding',
        exact_policy='inherited from producing scores; benchmark producer must filter exact hits only and keep models fixed',
        seconds=time.monotonic()-started), indent=2)+'\n')
    return 0
