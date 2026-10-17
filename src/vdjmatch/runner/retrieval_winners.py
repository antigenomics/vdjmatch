"""Streaming all-label maxima with bounded, bit-packed exact ties."""

from numbers import Integral

import numpy as np
import polars as pl


class WinnerAccumulator:
    """Retain every finite supported maximum without a query-by-label score grid.

    Repeated coordinates, including different scores, reduce by maximum. Support
    is candidate-local: an unsupported candidate never erases an earlier winner.
    ``max_bytes`` bounds persistent numeric state, not caller batches or total RSS.
    """

    def __init__(self, n_queries, targets: list[str], scorers: list[str],
                 max_bytes=12 * 1024**3):
        if (not isinstance(n_queries, Integral) or isinstance(n_queries, bool)
                or not 0 <= n_queries <= 2**32):
            raise ValueError('n_queries must be an integer in [0, 2**32]')
        if (not isinstance(max_bytes, Integral) or isinstance(max_bytes, bool)
                or max_bytes <= 0):
            raise ValueError('max_bytes must be a positive integer')
        for name, values in [('targets', targets), ('scorers', scorers)]:
            if (not values or any(not isinstance(x, str) or not x for x in values)
                    or len(set(values)) != len(values)):
                raise ValueError(f'{name} must contain distinct nonempty strings')
        if len(targets) > 65536:
            raise ValueError('at most 65536 targets are supported')
        n_queries = int(n_queries)
        n_bytes = (len(targets) + 7) // 8
        required = len(scorers) * n_queries * (8 + n_bytes)
        if required > max_bytes:
            raise ValueError(f'winner array budget exceeded: {required} bytes required')
        self.n_queries = n_queries
        self.targets = tuple(targets)
        self.scorers = tuple(scorers)
        self.maxima = np.full((len(scorers), n_queries), -np.inf, dtype=np.float64)
        self.ties = np.zeros((len(scorers), n_queries, n_bytes), dtype=np.uint8)

    def consume(self, scorer, query_indices, target_indices, scores, supported=None):
        """Reduce unordered one-dimensional candidate arrays, retaining exact ties."""
        try:
            scorer_index = self.scorers.index(scorer)
        except ValueError:
            raise ValueError(f'unknown scorer: {scorer!r}') from None
        q, t = np.asarray(query_indices), np.asarray(target_indices)
        raw = np.asarray(scores)
        if q.ndim != 1 or t.shape != q.shape or raw.shape != q.shape:
            raise ValueError('candidate arrays must be one-dimensional with equal lengths')
        for name, coords, size in [('query', q, self.n_queries),
                                   ('target', t, len(self.targets))]:
            if (coords.dtype.kind not in 'iu' or np.any(coords < 0)
                    or np.any(coords >= size)):
                raise ValueError(f'{name} indices must be integers within the domain')
        if raw.dtype.kind not in 'fiu':
            raise ValueError('scores must be numeric')
        values = raw.astype(np.float64, copy=False)
        valid = np.isfinite(values)
        if supported is not None:
            support = np.asarray(supported)
            if support.shape != q.shape or support.dtype.kind != 'b':
                raise ValueError('supported must be a Boolean array matching scores')
            valid &= support
        if not valid.any():
            return
        q, t, values = q[valid], t[valid], values[valid]
        rows, inverse = np.unique(q, return_inverse=True)
        batch_maxima = np.full(len(rows), -np.inf)
        np.maximum.at(batch_maxima, inverse, values)
        maxima, ties = self.maxima[scorer_index], self.ties[scorer_index]
        improved = rows[batch_maxima > maxima[rows]]
        ties[improved] = 0
        maxima[rows] = np.maximum(maxima[rows], batch_maxima)
        winners = values == maxima[q]
        q, t = q[winners], t[winners]
        np.bitwise_or.at(ties, (q, t // 8), np.left_shift(np.uint8(1), (t % 8).astype(np.uint8)))

    def frames(self, query_ids=None, max_rows=65536):
        """Yield bounded long tables, including one explicit row for each no-hit.

        ``query_ids`` supplies one unique, non-null stable identity per query in
        positional order. Its Polars dtype is preserved. Defaults are UInt32.
        ``max_rows`` must fit every target tie of at least one query.
        """
        n_targets = len(self.targets)
        if (not isinstance(max_rows, Integral) or isinstance(max_rows, bool)
                or max_rows < n_targets):
            raise ValueError('max_rows must be an integer at least the target count')
        ids = None if query_ids is None else pl.Series('query_id', query_ids)
        if ids is not None and (len(ids) != self.n_queries or ids.null_count()
                                or ids.n_unique() != self.n_queries):
            raise ValueError('query_ids must contain one unique non-null identity per query')
        targets = pl.Series('epitope', self.targets, dtype=pl.String)
        for s, scorer in enumerate(self.scorers):
            for start in range(0, self.n_queries, int(max_rows) // n_targets):
                end = min(start + int(max_rows) // n_targets, self.n_queries)
                bits = np.unpackbits(self.ties[s, start:end], axis=1,
                                     count=n_targets, bitorder='little')
                nohit = ~np.isfinite(self.maxima[s, start:end])
                bits[nohit, 0] = 1
                rows, labels = np.nonzero(bits)
                hit = ~nohit[rows]
                yield pl.DataFrame({
                    'query_id': (pl.Series('query_id', (rows + start).astype(np.uint32))
                                 if ids is None else ids.gather(rows + start)),
                    'epitope': targets.gather(labels).set(pl.Series(~hit), None),
                    'scorer': pl.Series('scorer', [scorer] * len(rows), dtype=pl.String),
                    'score': pl.Series(self.maxima[s, rows + start]).set(pl.Series(~hit), None),
                    'supported': hit,
                })
