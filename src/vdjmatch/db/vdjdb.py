"""Fetch and parse VDJdb releases (github.com/antigenomics/vdjdb-db).

The latest release ships a single ZIP (e.g. ``vdjdb-2026-06-03.zip``) that contains the
``vdjdb.txt`` / ``vdjdb.slim.txt`` / ``vdjdb_full.txt`` tables. ``fetch_latest`` downloads and
caches the ZIP by release tag and extracts the requested table; ``load`` parses any VDJdb table
(downloaded or a local snapshot) into a normalized polars frame.
"""
from __future__ import annotations

import io
import json
import os
import urllib.request
import zipfile
from pathlib import Path

import polars as pl

from . import schema
from .cache import cache_dir

_REPO = "antigenomics/vdjdb-db"
_UA = {"User-Agent": "vdjmatch", "Accept": "application/vnd.github+json"}
# requested table -> substring matched against ZIP member basenames (most specific first)
_MEMBER = {"full": "vdjdb_full.txt", "slim": "vdjdb.slim.txt", "default": "vdjdb.txt"}


def _release_json(pin: str | None) -> dict:
    url = (f"https://api.github.com/repos/{_REPO}/releases/tags/{pin}" if pin
           else f"https://api.github.com/repos/{_REPO}/releases/latest")
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=60) as r:  # noqa: S310 (trusted host)
        return json.load(r)


def _zip_asset(rel: dict) -> str:
    for a in rel.get("assets", []):
        if a["name"].endswith(".zip"):
            return a["browser_download_url"]
    raise RuntimeError(f"no .zip asset in VDJdb release {rel.get('tag_name')!r}")


def fetch_latest(asset: str = "slim", cache: str | os.PathLike | None = None,
                 pin: str | None = None, force: bool = False) -> Path:
    """Download (and cache) a VDJdb table from the latest release (or ``pin`` tag).

    Args:
        asset: which table to extract — ``"slim"`` (default), ``"full"``, or ``"default"``.
        cache: cache directory (default ``~/.cache/vdjmatch`` or ``$VDJMATCH_CACHE``).
        pin: pin a specific release tag for reproducibility; ``None`` = latest.
        force: re-download even if cached.

    Returns:
        Path to the extracted table file, named ``vdjdb-<tag>.<asset>.txt`` in the cache.
    """
    cdir = cache_dir(cache)
    rel = _release_json(pin)
    tag = rel["tag_name"]
    out = cdir / f"vdjdb-{tag}.{asset}.txt"
    if out.exists() and not force:
        return out
    zip_path = cdir / f"vdjdb-{tag}.zip"
    if not zip_path.exists() or force:
        req = urllib.request.Request(_zip_asset(rel), headers=_UA)
        with urllib.request.urlopen(req, timeout=300) as r:  # noqa: S310
            data = r.read()
        zip_path.write_bytes(data)
    want = _MEMBER.get(asset, _MEMBER["default"])
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        member = next((n for n in names if Path(n).name == want), None)
        if member is None:  # fall back to substring match
            key = want.replace("vdjdb", "").strip("._") or "vdjdb.txt"
            member = next((n for n in names if key in Path(n).name), None)
        if member is None:
            raise RuntimeError(f"{want!r} not found in {zip_path.name}; members={names}")
        out.write_bytes(zf.read(member))
    return out


def load(source: str | os.PathLike | None = None, *, asset: str = "slim",
         species: str | None = None, gene: str | None = None, mhc_class: str | None = None,
         min_score: int = 0, paired_only: bool = False, pin: str | None = None) -> pl.DataFrame:
    """Parse a VDJdb table into a normalized polars frame (see ``schema.CANONICAL``).

    ``source`` is a path to an existing table (e.g. a local snapshot); if ``None`` the latest
    release is fetched (``asset``/``pin``). Optional filters: ``species`` (e.g. ``"HomoSapiens"``),
    ``gene`` (``"TRA"``/``"TRB"``), ``mhc_class`` (``"MHCI"``/``"MHCII"``), ``min_score`` (VDJdb
    confidence), ``paired_only`` (keep only rows with a non-zero ``complex_id``).
    """
    path = Path(source) if source is not None else fetch_latest(asset=asset, pin=pin)
    df = pl.read_csv(path, separator="\t", quote_char=None, infer_schema_length=0)  # all str
    df = schema.normalize(df)
    df = df.with_columns(
        pl.col("vdjdb_score").cast(pl.Int64, strict=False).fill_null(0),
        pl.col("complex_id").cast(pl.Int64, strict=False).fill_null(0),
    )
    if species is not None:
        df = df.filter(pl.col("species") == species)
    if gene is not None:
        df = df.filter(pl.col("gene") == gene)
    if mhc_class is not None:
        df = df.filter(pl.col("mhc_class") == mhc_class)
    if min_score > 0:
        df = df.filter(pl.col("vdjdb_score") >= min_score)
    if paired_only:
        df = df.filter(pl.col("complex_id") != 0)
    return df
