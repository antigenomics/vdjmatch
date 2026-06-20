"""Cache directory resolution for downloaded VDJdb releases."""
from __future__ import annotations

import os
from pathlib import Path


def cache_dir(override: str | os.PathLike | None = None) -> Path:
    """Return (and create) the vdjmatch cache directory.

    Order: explicit ``override`` → ``$VDJMATCH_CACHE`` → ``~/.cache/vdjmatch``.
    """
    base = override or os.environ.get("VDJMATCH_CACHE") or (Path.home() / ".cache" / "vdjmatch")
    p = Path(base)
    p.mkdir(parents=True, exist_ok=True)
    return p
