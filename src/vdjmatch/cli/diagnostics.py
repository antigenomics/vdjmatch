"""One fresh execution-owned preparation for motif and community diagnostics."""
from argparse import Namespace


def register(subparsers):
    from .structure import register as register_structure
    from .coverage import add_graph_options
    parser = register_structure(subparsers, 'reference-diagnostics')
    add_graph_options(parser, combined=True)
    parser.set_defaults(func=main, motif_mode='all', mock_replicates=200)


def main(a):
    from . import structure, coverage
    if not a.control or a.mock_replicates < 2:
        raise ValueError('reference-diagnostics requires --control and at least two mock replicates')
    if a.geometry == 'levenshtein' and (a.sw_radius is not None or a.sw_radius_knee):
        raise ValueError('SW radius options require SW geometry')
    passes = []
    for geometry in (('levenshtein', 'sw') if a.geometry == 'both' else (a.geometry,)):
        options = vars(a) | dict(geometry=geometry, output_prefix=a.output_prefix+'.'+geometry,
            scope='1,1,1,1' if a.geometry == 'both' and geometry == 'sw' else a.scope,
            sw_radius=a.sw_radius if geometry == 'sw' else None,
            sw_radius_knee=a.sw_radius_knee if geometry == 'sw' else False)
        options = Namespace(**options)
        coverage.validate_options(options)
        passes.append(options)
    inputs = structure.main(a, return_inputs=True)
    for options in passes:
        coverage.main(options, inputs=inputs)
    return 0
