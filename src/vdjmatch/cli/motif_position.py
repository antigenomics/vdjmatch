"""Contrast junction positions against prepared controls within one V/J/length stratum."""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('motif-position', help=__doc__)
    p.add_argument('--vdjdb', required=True)
    p.add_argument('--control', required=True, help='original immutable raw control, for model identity validation')
    p.add_argument('--control-model', required=True, help='prepared control directory; no model is rebuilt')
    p.add_argument('--species', required=True)
    p.add_argument('--locus', required=True, choices=['TRA', 'TRB'])
    p.add_argument('--epitope', required=True)
    p.add_argument('--v-gene', required=True)
    p.add_argument('--j-gene', required=True)
    p.add_argument('--junction-length', required=True, type=int)
    p.add_argument('--motif', help='optional exact amino-acid motif, including its start-position counts')
    p.add_argument('--reference-population', choices=['curated-benchmark', 'all-valid'], default='curated-benchmark')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--output-prefix', required=True)
    p.set_defaults(func=main)


def main(a):
    import json
    from importlib.metadata import version
    import polars as pl
    import seqtree._core as native
    from .. import db
    from ..aggregate.motif_position import position_contrast
    from ..aggregate.structure import _gene_expression
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..runner.annotation import _validate_output_paths
    from .control_model import load_control_model

    prefix = _validate_output_paths(a.output_prefix, [a.vdjdb, a.control, a.control_model])
    organism = _organism(a.species)
    species = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[organism]
    raw = db.load(a.vdjdb, species=species, gene=a.locus, epitope=a.epitope)
    selected = raw.filter(pl.col('reference_valid')) if a.reference_population == 'curated-benchmark' else raw
    events, _, _, control_info = load_control_model(a.control_model, a.control,
        species=organism, locus=a.locus, modes=())
    positions, motifs = position_contrast(selected, events, locus=a.locus, v_gene=a.v_gene,
        j_gene=a.j_gene, junction_length=a.junction_length, motif=a.motif)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    tables = {}
    for suffix, table in [('positions', positions), ('motifs', motifs)]:
        if table is None:
            continue
        table = table.with_columns(pl.lit(species).alias('species'), pl.lit(a.locus).alias('locus'),
            pl.lit(a.epitope).alias('epitope'))
        path = Path(str(prefix)+'.'+suffix+'.tsv')
        table.write_csv(path, separator='\t')
        tables[suffix] = dict(path=str(path), sha256=sha256(path), rows=table.height)
    manifest = dict(command='motif-position', options=vars(a) | {'func': 'motif-position'},
        reference=db.provenance(a.vdjdb), control=control_info, tables=tables,
        reference_records_before_stratum=selected.height, reference_records_loaded=raw.height,
        population_unit='unique full junction/allele-pooled V gene/allele-pooled J gene keys; event weights ignored',
        positions='1-based full junction including conserved anchors; motif position 0 is presence anywhere',
        comparison='reference fraction minus prepared-control fraction conditional on both V/J genes and junction length',
        interpretation='descriptive background contrast; not calibrated selection or binding evidence',
        native_batches=0, software={name: version(name) for name in ('vdjmatch', 'seqtree', 'vdjtools', 'polars')},
        native_sha256=sha256(Path(native.__file__)),
        source_sha256={fn.__name__: sha256(Path(fn.__code__.co_filename))
            for fn in (main, position_contrast, _gene_expression, db.load, load_control_model)})
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0
