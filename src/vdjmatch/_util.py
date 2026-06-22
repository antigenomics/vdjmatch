"""Small shared utilities (progress / verbosity)."""
from __future__ import annotations

import sys
from collections.abc import Iterable


def progress(iterable: Iterable, *, total: int | None = None, desc: str = "",
             enable: bool = False) -> Iterable:
    """Wrap an iterable in a ``tqdm`` progress bar when ``enable`` and tqdm is importable; otherwise
    pass through (a terse stderr counter if tqdm is missing). Keeps progress optional and dependency-
    light — the core never hard-requires tqdm."""
    if not enable:
        return iterable
    try:
        from tqdm import tqdm
        return tqdm(iterable, total=total, desc=desc, leave=False)
    except ImportError:
        return _counter(iterable, total, desc)


def _counter(iterable, total, desc):
    n = 0
    for x in iterable:
        n += 1
        print(f"\r{desc} {n}/{total or '?'}", end="", file=sys.stderr, flush=True)
        yield x
    print("", file=sys.stderr)


def chunked(seq, size: int):
    """Yield ``seq`` in lists of at most ``size`` (size<=0 -> one chunk)."""
    seq = list(seq)
    if size <= 0:
        yield seq
        return
    for i in range(0, len(seq), size):
        yield seq[i:i + size]
