"""Finite-population discovery of observed V-conditioned core three-residue motifs."""


def register(subparsers):
    parser = subparsers.add_parser('motif-depth', help=__doc__)
    parser.add_argument('--vdjdb', required=True)
    parser.add_argument('--depth', type=int, required=True, action='append')
    parser.add_argument('--output-prefix', required=True)
    parser.set_defaults(func=main)


def catalog_depth(reference, depths):
    """Expected observed catalog size in a uniform subsample of junction/V keys."""
    from bisect import bisect_left
    import numpy as np
    import polars as pl
    from scipy.stats import hypergeom
    from ..aggregate.structure import _gene_expression
    from ..evalue.epitope_rules import _tokens

    if not depths or any(type(n) is not int or n < 1 for n in depths):
        raise ValueError('positive integer sample depths required')
    groups = ['species', 'gene', 'epitope']
    schema = {**dict.fromkeys(groups, pl.String),
        **dict.fromkeys(['source_records', 'junction_v_keys', 'unresolved_v_keys',
                         'observed_motifs', 'guaranteed_catalog_depth', 'expected_90pct_depth', 'sample_keys'], pl.Int64),
        'expected_observed_motifs': pl.Float64}
    rows = []
    valid = reference.filter(pl.col('reference_valid') &
        pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$') &
        pl.col('gene').is_in(['TRA', 'TRB']))
    for (species, locus, epitope), group in valid.partition_by(groups, as_dict=True).items():
        keys = group.select(pl.col('cdr3').alias('seq0'),
            _gene_expression(locus, 'V').alias('v0')).unique().with_row_index('_id')
        # Unresolved V keys remain sampling opportunities but do not identify a V-specific motif.
        tokens = _tokens(keys.filter(pl.col('v0').is_not_null()), [0], 'exact3')
        carriers = tokens.group_by('v0', 't0').len()['len'].to_numpy()
        sizes, frequencies = np.unique(carriers, return_counts=True)
        # A sample can omit the rarest motif while it contains only noncarriers.
        guaranteed_depth = keys.height - int(sizes[0]) + 1 if sizes.size else 0
        def expected_count(depth):
            return float(np.dot(frequencies, hypergeom.sf(0, keys.height, sizes, depth)))

        depth90 = bisect_left(range(guaranteed_depth + 1), 0.9 * len(carriers),
                              key=expected_count) if sizes.size else 0
        for depth in sorted(set(depths)):
            expected = None
            if depth <= keys.height:
                expected = expected_count(depth)
            rows.append((species, locus, epitope, group.height, keys.height,
                keys['v0'].null_count(), len(carriers), guaranteed_depth, depth90, depth, expected))
    return pl.DataFrame(rows, schema=schema, orient='row').sort(*groups, 'sample_keys')


def main(args):
    import json
    from pathlib import Path
    from .. import db
    from ..db.cache import sha256
    from ..runner.annotation import _validate_output_paths
    prefix = _validate_output_paths(args.output_prefix, [args.vdjdb])
    output = Path(str(prefix) + '.motif-depth.tsv')
    reference = db.provenance(args.vdjdb)
    table = catalog_depth(db.load(args.vdjdb), args.depth)
    if db.provenance(args.vdjdb) != reference:
        raise ValueError('reference changed during motif depth calculation')
    output.parent.mkdir(parents=True, exist_ok=True)
    table.write_csv(output, separator='\t')
    Path(str(prefix) + '.manifest.json').write_text(json.dumps(dict(reference=reference,
        reference_population='curated-benchmark', depths=sorted(set(args.depth)),
        catalog='V gene and core exact3; one presence per distinct junction/V key',
        output_sha256=sha256(output)), indent=2) + '\n')
    return 0
