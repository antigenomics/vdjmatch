"""Database-wide equal-tail selection over supported candidate pMHCs and receptor axes."""
import numpy as np
import polars as pl

from ..aggregate.candidates import PMHC
from .equal_tail import equal_tail_scores


def equal_tail_selection(queries, reference, *, locus, species, threads, exclude_exact=False,
                         distance_geometry='native', background_scope='all-hla', withheld_column=None):
    """Return one selected any-call score per query, retaining tied and unsupported cases.

    For a paired query, take the maximum over linked, TRA, and TRB selections.
    Withholding omits the annotated peptide at every restriction only from selection.
    """
    if queries['query_id'].null_count() or queries['query_id'].n_unique() != len(queries):
        raise ValueError('database-wide selection requires unique nonnull query IDs')
    held = queries[withheld_column].to_numpy() if withheld_column is not None else None
    n = len(queries)
    numerator = np.ones(n, dtype=np.uint64)
    denominator = np.ones(n, dtype=np.uint64)
    winner = np.full(n, -1, dtype=np.int64)
    ambiguous = np.zeros(n, dtype=bool)
    catalog = {}
    axes = ('paired', 'TRA', 'TRB') if locus == 'paired' else (locus,)

    def consume(task, scores):
        candidate = tuple(task[c] for c in PMHC)
        label = catalog.setdefault(candidate, len(catalog))
        supported = (scores['has_evidence'] & scores['calibration_available']).to_numpy()
        if held is not None:
            supported = supported & (held != task['epitope'])
        ix = np.flatnonzero(supported)
        if not len(ix):
            return
        nums = scores['product_numerator'].to_numpy()[ix]
        dens = scores['tail_denominator'].to_numpy()[ix]
        # Cross-products can exceed uint64. Python integers keep exact ties across
        # candidate-specific background denominators, including ties across axes.
        left = nums.astype(object) * denominator[ix].astype(object)**2
        right = numerator[ix].astype(object) * dens.astype(object)**2
        prior = winner[ix]
        tied = ix[(prior >= 0) & (prior != label) & (left == right)]
        ambiguous[tied] = True
        improved = (prior < 0) | (left < right)
        replace = ix[improved]
        winner[replace] = label
        ambiguous[replace] = False
        numerator[replace], denominator[replace] = nums[improved], dens[improved]

    for axis in axes:
        equal_tail_scores(queries, reference, locus=axis, species=species, threads=threads,
                          exclude_exact=exclude_exact, distance_geometry=distance_geometry,
                          background_scope=background_scope, _consume=consume)
    labels = list(catalog)
    any_hit = winner >= 0
    score = -np.log(np.asarray(numerator.astype(object) / denominator.astype(object)**2, dtype=float))
    base = queries.select('query_id', *[c for c in ('pair_id','sequence_id') if c in queries])
    return base.with_columns(
        pl.Series('score', score).set(pl.Series(~any_hit), None),
        pl.Series('product_numerator', numerator).set(pl.Series(~any_hit), None),
        pl.Series('tail_denominator', denominator).set(pl.Series(~any_hit), None),
        pl.Series('selection_status', np.where(~any_hit, 'no_support', np.where(ambiguous, 'ambiguous', 'unique'))),
        *[pl.Series('selected_'+c, [labels[i][j] if i >= 0 and not tie else None
                                   for i,tie in zip(winner, ambiguous)], dtype=pl.String)
          for j,c in enumerate(PMHC)],
    )
