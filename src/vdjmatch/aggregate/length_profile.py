"""Empirical junction spectratypes on allele-pooled joint V/J keys.

CDF distances are descriptive effects against a supplied control, not Gaussian
normality tests or binding confidence. Junction lengths include both anchors.
"""
import numpy as np
import polars as pl

from .structure import _gene_expression


def length_histogram(records, *, locus, sequence_column='cdr3'):
    """Return exact length counts and retained/dropped counts; keep missing genes.

    Each distinct junction/V-gene/J-gene key weighs one. Unresolved genes form
    explicit null categories; alleles are pooled through the existing gene parser.
    Neither source-record multiplicity nor control-event weights supply weights.
    """
    if locus not in {'TRA', 'TRB'} or not {sequence_column, 'v', 'j'}.issubset(records.columns):
        raise ValueError('spectratype requires TRA/TRB and junction/V/J columns')
    selected = records.select(pl.col(sequence_column).cast(pl.String).alias('cdr3'),
        _gene_expression(locus, 'V').alias('v_gene'),
        _gene_expression(locus, 'J').alias('j_gene')).filter(
            pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False))
    keys = selected.unique().with_columns(
        pl.col('cdr3').str.len_chars().alias('junction_length'),
        (pl.col('v_gene').is_null() | pl.col('j_gene').is_null()).alias('_unresolved'))
    histogram = keys.group_by('junction_length').agg(pl.len().alias('keys'),
        pl.col('_unresolved').sum().alias('unresolved_gene_keys')).sort('junction_length')
    return histogram, dict(source_rows=records.height, selected_rows=selected.height,
        dropped_rows=records.height-selected.height, selected_keys=keys.height,
        unresolved_gene_keys=int(keys['_unresolved'].sum()))


def _distribution(histogram):
    if not {'junction_length', 'keys'}.issubset(histogram.columns):
        raise ValueError('length histogram requires junction_length and keys')
    ordered = histogram.sort('junction_length')
    if ordered.is_empty():
        return None
    if (ordered['junction_length'].null_count() or ordered['keys'].null_count() or
            not ordered['junction_length'].dtype.is_integer() or not ordered['keys'].dtype.is_integer() or
            (ordered['junction_length'] <= 0).any() or (ordered['keys'] <= 0).any() or
            ordered['junction_length'].n_unique() != ordered.height):
        raise ValueError('length histogram requires unique positive integer lengths and counts')
    lengths, counts = ordered['junction_length'].to_numpy(), ordered['keys'].to_numpy()
    cumulative = counts.cumsum()
    total = int(cumulative[-1])
    median = (float(lengths[np.searchsorted(cumulative, (total+1)//2)]) +
              float(lengths[np.searchsorted(cumulative, (total+2)//2)]))/2
    return dict(lengths=lengths, counts=counts, cumulative=cumulative,
        cdf=cumulative/total, total=total, mean=float(np.dot(lengths, counts/total)), median=median)


def length_control(histogram):
    """Prepare one shared control CDF and its inclusive discrete central interval."""
    result = _distribution(histogram)
    if result is None:
        raise ValueError('spectratype requires nonempty control keys')
    n, lengths, cumulative = result['total'], result['lengths'], result['cumulative']
    result['lower'] = int(lengths[np.searchsorted(cumulative, (n+39)//40)])
    result['upper'] = int(lengths[np.searchsorted(cumulative, (39*n+39)//40)])
    result['tail_keys'] = int(result['counts'][(lengths < result['lower']) | (lengths > result['upper'])].sum())
    return result


def length_summary(histogram, control):
    """Compare empirical CDFs without bin matching, smoothing or a p-value.

    Wasserstein-1 integrates absolute CDF differences over integer length gaps,
    including lengths absent in either population. Median is the ordinary sample
    median (mean of the middle two observations for even key counts).
    """
    reference = _distribution(histogram)
    row = dict(reference_keys=reference['total'] if reference else 0, control_keys=control['total'],
        control_central95_lower_length=control['lower'], control_central95_upper_length=control['upper'],
        control_mean_length=control['mean'], control_median_length=control['median'],
        control_tail_keys=control['tail_keys'], control_tail_fraction=control['tail_keys']/control['total'])
    if reference is None:
        return row | dict(reference_mean_length=None, reference_median_length=None,
            cdf_max_distance=None, wasserstein1_residues=None, mean_length_shift_residues=None,
            median_length_shift_residues=None, reference_tail_keys=0, reference_tail_fraction=None,
            spectratype_status='empty_reference')
    grid = np.union1d(reference['lengths'], control['lengths'])
    def cdf_at(population):
        index = np.searchsorted(population['lengths'], grid, side='right')
        return np.r_[0., population['cdf']][index]
    gap = np.abs(cdf_at(reference)-cdf_at(control))
    tails = int(reference['counts'][(reference['lengths'] < control['lower']) |
                                   (reference['lengths'] > control['upper'])].sum())
    return row | dict(reference_mean_length=reference['mean'], reference_median_length=reference['median'],
        cdf_max_distance=float(gap.max()), wasserstein1_residues=float(np.dot(gap[:-1], np.diff(grid))),
        mean_length_shift_residues=reference['mean']-control['mean'],
        median_length_shift_residues=reference['median']-control['median'],
        reference_tail_keys=tails, reference_tail_fraction=tails/reference['total'], spectratype_status='available')
