"""Immutable reference cache with bounded, process-coordinated atomic publication."""

from __future__ import annotations

import hashlib
import os
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path


def cache_dir(override: str | os.PathLike | None = None) -> Path:
    """Explicit override → VDJMATCH_CACHE → ~/.cache/vdjmatch."""
    p = Path(
        override
        or os.environ.get("VDJMATCH_CACHE")
        or Path.home() / ".cache" / "vdjmatch"
    )
    p.mkdir(parents=True, exist_ok=True)
    return p


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()


@contextmanager
def locked(path: Path, timeout: float = 120):
    """OS-owned file lock; a killed owner releases it without stale-lock deletion."""
    with path.with_name(path.name + ".lock").open("a+b") as fh:
        deadline = time.monotonic() + timeout
        if os.name == "nt":
            import msvcrt

            fh.seek(0)
            fh.write(b"\0")
            fh.flush()

            def lock():
                msvcrt.locking(fh.fileno(), msvcrt.LK_NBLCK, 1)

            def unlock():
                msvcrt.locking(fh.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl

            def lock():
                fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)

            def unlock():
                fcntl.flock(fh, fcntl.LOCK_UN)

        while True:
            try:
                fh.seek(0)
                lock()
                break
            except (BlockingIOError, OSError):
                if time.monotonic() >= deadline:
                    raise TimeoutError(
                        f"timed out waiting for reference cache lock: {path}"
                    )
                time.sleep(0.05)
        try:
            yield
        finally:
            fh.seek(0)
            unlock()


def valid(path: Path, expected: str | None = None) -> bool:
    if not path.is_file():
        return False
    digest = path.with_name(path.name + ".sha256")
    if expected is None:
        if not digest.is_file():
            return False
        expected = digest.read_text().strip()
    return sha256(path) == expected


@contextmanager
def staging(path: Path):
    """Unique sibling file, cleaned only by its owner."""
    fd, name = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    os.close(fd)
    stage = Path(name)
    try:
        yield stage
    finally:
        stage.unlink(missing_ok=True)


def publish(stage: Path, target: Path):
    digest = sha256(stage)
    os.replace(stage, target)
    checksum = target.with_name(target.name + ".sha256")
    with staging(checksum) as sidecar:
        sidecar.write_text(digest + "\n")
        os.replace(sidecar, checksum)
