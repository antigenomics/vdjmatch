"""One fresh execution-owned preparation for motif and community diagnostics."""
from argparse import Namespace


def register(subparsers):
    from .structure import register as register_structure
    from .coverage import add_graph_options
    parser = register_structure(subparsers, 'reference-diagnostics')
    add_graph_options(parser, combined=True)
    parser.set_defaults(func=main, motif_mode='all', mock_replicates=50)


def main(a):
    from . import structure, coverage
    if not a.control or a.mock_replicates < 2:
        raise ValueError('reference-diagnostics requires --control and at least two mock replicates')
    inputs = structure.main(a, return_inputs=True)
    for geometry in (('levenshtein', 'sw') if a.geometry == 'both' else (a.geometry,)):
        options = vars(a) | dict(geometry=geometry, output_prefix=a.output_prefix+'.'+geometry,
            scope=a.scope if geometry == 'levenshtein' else '1,1,1,1',
            sw_radius=a.sw_radius if geometry == 'sw' else None,
            sw_radius_knee=a.sw_radius_knee if geometry == 'sw' else False)
        coverage.main(Namespace(**options), inputs=inputs)
    return 0
