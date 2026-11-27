"""Record publicity from one exact batched Hamming<=1 control search."""
import numpy as np
import polars as pl


def junction_publicity(junctions, control, *, threads=1):
    """Return distinct control-junction matches, retaining exact matches.

    Inputs are valid full anchor-inclusive junctions from one species/locus.
    Control repeats collapse to unique junctions; no donor/publicity prevalence
    is inferred from a pooled repertoire. Output follows query order.
    """
    from seqtree import Index
    from ..match.scope import search_params

    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError('positive integer threads required')
    queries, background = (
        seqs if isinstance(seqs, pl.Series) else pl.Series(list(seqs))
        for seqs in (junctions, control)
    )
    for series in (queries, background):
        if len(series) and (series.dtype != pl.String or
                not series.str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).all()):
            raise ValueError('publicity requires valid full amino-acid junctions')
    if not len(background):
        raise ValueError('publicity requires nonempty control junctions')
    background = background.unique().sort()
    index = Index.build(background.to_list(), 'aa')
    counts = np.asarray(index.edit_histogram_batch(queries.to_list(), search_params('1,0,0,1'), threads, False), dtype=np.int64).reshape(-1, 2)
    return pl.DataFrame({'cdr3': queries.cast(pl.String),
        'control_exact_junctions': counts[:, 0],
        'control_hamming1_junctions': counts.sum(axis=1)}).with_columns(
            pl.lit(len(background), dtype=pl.Int64).alias('control_junctions'))
