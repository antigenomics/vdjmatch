"""Resolve legacy, primary and AIRR VDJdb releases through one offline-capable loader."""

from __future__ import annotations

import gzip
import json
import os
import shutil
import tempfile
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath

import polars as pl

from . import schema
from .cache import cache_dir, locked, publish, sha256, staging, valid

_REPO = "antigenomics/vdjdb-db"
_UA = {"User-Agent": "vdjmatch", "Accept": "application/vnd.github+json"}
_MEMBER = {
    "full": "vdjdb_full.txt",
    "slim": "vdjdb.slim.txt",
    "default": "vdjdb.txt",
    "legacy": "vdjdb.txt",
}
_HF_REPO = "isalgo/airr_benchmark"
_HF_TAG = "2026-06-11-ZENODO"


def _json_url(url: str) -> dict:
    with urllib.request.urlopen(
        urllib.request.Request(url, headers=_UA), timeout=60
    ) as r:
        return json.load(r)


def _release_json(pin: str | None) -> dict:
    suffix = "tags/" + urllib.parse.quote(pin, safe="") if pin else "latest"
    return _json_url(f"https://api.github.com/repos/{_REPO}/releases/{suffix}")


def _select_asset(rel: dict, asset: str) -> tuple[dict, dict | None]:
    if asset not in {"default", "primary", "legacy", "slim", "full", "airr"}:
        raise ValueError(f"unsupported VDJdb asset: {asset!r}")
    assets = rel.get("assets", [])
    manifests = [a for a in assets if a["name"] == "manifest.json"]
    if manifests:
        if len(manifests) != 1:
            raise ValueError("release has multiple manifests")
        manifest = _json_url(manifests[0]["browser_download_url"])
        if manifest.get("tag") != rel["tag_name"]:
            raise ValueError("manifest tag does not match GitHub release")
        role = (
            "primary"
            if asset in {"default", "primary"}
            else "airr"
            if asset == "airr"
            else "legacy"
        )
        entries = [b for b in manifest.get("bundles", []) if b.get("role") == role]
        if len(entries) != 1:
            raise ValueError(f"manifest must declare one {role!r} bundle")
        entry = entries[0]
        matches = [a for a in assets if a["name"] == entry.get("file")]
        if len(matches) != 1 or not matches[0]["name"].endswith(".zip"):
            raise ValueError(
                f"manifest bundle {entry.get('file')!r} missing from release"
            )
        if not isinstance(entry.get("sha256"), str) or len(entry["sha256"]) != 64:
            raise ValueError("manifest bundle requires SHA256 digest")
        return matches[0], entry
    zips = [a for a in assets if a["name"].endswith(".zip")]
    if len(zips) != 1:
        raise ValueError(
            "release requires manifest.json when it does not have one historical ZIP"
        )
    if asset in {"primary", "airr"}:
        raise ValueError(f"historical release does not declare a {asset!r} bundle")
    return zips[0], None


def _release_checksum(rel: dict, remote: dict) -> str | None:
    checksums = [a for a in rel.get("assets", []) if a["name"] == "SHA256SUMS"]
    if len(checksums) > 1:
        raise ValueError("release has multiple SHA256SUMS assets")
    if not checksums:
        return None
    req = urllib.request.Request(checksums[0]["browser_download_url"], headers=_UA)
    with urllib.request.urlopen(req, timeout=60) as src:
        lines = src.read().decode("utf-8").splitlines()
    matches = [
        line.split()[0]
        for line in lines
        if len(line.split()) == 2 and line.split()[1].lstrip("*") == remote["name"]
    ]
    if len(matches) != 1 or len(matches[0]) != 64:
        raise ValueError("SHA256SUMS does not identify the selected release asset")
    return matches[0]


def _archive_members(zf: zipfile.ZipFile) -> dict[str, str]:
    members = {}
    expanded = 0
    for info in zf.infolist():
        expanded += info.file_size
        if info.file_size > 2 * 1024**3 or expanded > 8 * 1024**3:
            raise ValueError(
                "reference archive exceeds supported expansion budget (2 GiB member, 8 GiB total)"
            )
        p = PurePosixPath(info.filename)
        if (
            p.is_absolute()
            or ".." in p.parts
            or "\\" in info.filename
            or (info.external_attr >> 16) & 0o170000 == 0o120000
        ):
            raise ValueError(f"unsafe archive member: {info.filename!r}")
        if info.is_dir():
            continue
        if p.name in members:
            raise ValueError(f"ambiguous archive basename: {p.name}")
        members[p.name] = info.filename
    return members


def _validate_zip(path: Path, entry: dict | None = None):
    with zipfile.ZipFile(path) as zf:
        names = _archive_members(zf)
        if not names:
            raise ValueError("empty VDJdb archive")
        if entry is not None:
            if set(entry.get("members", [])) != set(zf.namelist()):
                raise ValueError("archive members differ from manifest")
            if entry.get("bytes") is not None and path.stat().st_size != entry["bytes"]:
                raise ValueError("archive size differs from manifest")
        if bad := zf.testzip():
            raise ValueError(f"corrupt archive member: {bad}")
        if not any(
            n in names
            for n in (
                "vdjdb.txt",
                "vdjdb.slim.txt",
                "vdjdb_full.txt",
                "vdjdb.parquet",
                "records.tsv",
                "records.parquet",
                "vdjdb.rearrangement.tsv",
            )
        ):
            raise ValueError("archive contains no supported VDJdb reference")


def fetch_latest(
    asset: str = "default",
    cache: str | os.PathLike | None = None,
    pin: str | None = None,
    force: bool = False,
) -> Path:
    """Fetch latest/pinned release; default selects primary, or historical long table.

    ``primary``/``airr`` return a ZIP suitable for :func:`load`; historical and explicit
    ``legacy``/``slim``/``full`` requests return an extracted table. Downloads stream into
    unique staging files and are checksummed and published under an OS process lock.
    """
    rel = _release_json(pin)
    remote, entry = _select_asset(rel, asset)
    # URL identity separates assets/tags and avoids untrusted names becoming filesystem paths.
    import hashlib

    identity = hashlib.sha256(remote["browser_download_url"].encode()).hexdigest()[:24]
    folder = cache_dir(cache) / "github" / identity
    folder.mkdir(parents=True, exist_ok=True)
    archive = folder / "release.zip"
    expected = entry["sha256"] if entry else remote.get("digest", "")
    if expected and expected.startswith("sha256:"):
        expected = expected[7:]
    sums = _release_checksum(rel, remote)
    if sums and expected and sums != expected:
        raise ValueError("release checksum declarations disagree")
    expected = expected or sums or None
    with locked(archive):
        if force or not valid(archive, expected):
            with staging(archive) as stage:
                req = urllib.request.Request(
                    remote["browser_download_url"], headers=_UA
                )
                with (
                    urllib.request.urlopen(req, timeout=300) as src,
                    stage.open("wb") as dst,
                ):
                    shutil.copyfileobj(src, dst, length=1024 * 1024)
                if expected and sha256(stage) != expected:
                    raise ValueError("downloaded VDJdb SHA256 does not match release")
                _validate_zip(stage, entry)
                layout = entry["role"] if entry and asset == "default" else asset
                schema.normalize(_load_archive(stage, layout))
                publish(stage, archive)
        _validate_zip(archive, entry)
        metadata = {
            "source": "github",
            "tag": rel["tag_name"],
            "url": remote["browser_download_url"],
            "bundle": entry["role"] if entry else "historical",
            "sha256": sha256(archive),
        }
        with staging(folder / "release.json") as stage:
            stage.write_text(json.dumps(metadata, indent=2) + "\n")
            os.replace(stage, folder / "release.json")
        if entry and entry["role"] in {"primary", "airr"}:
            return archive
        member = _MEMBER.get(asset, "vdjdb.txt")
        out = folder / member
        if force or not valid(out):
            with zipfile.ZipFile(archive) as zf:
                names = _archive_members(zf)
                if member not in names:
                    raise ValueError(
                        f"required table {member!r} missing from release archive"
                    )
                with staging(out) as stage:
                    with zf.open(names[member]) as src, stage.open("wb") as dst:
                        shutil.copyfileobj(src, dst)
                    schema.normalize(_read_table(stage, parquet=False, legacy=True))
                    publish(stage, out)
        return out


def fetch_hf(
    tag: str = _HF_TAG,
    asset: str = "default",
    cache: str | os.PathLike | None = None,
    force: bool = False,
    repo: str = _HF_REPO,
) -> Path:
    """Fetch a pinned gzipped legacy table, with cache identity scoped to HF repo/tag/table."""
    if asset not in _MEMBER:
        raise ValueError(f"unsupported HF legacy asset: {asset!r}")
    import hashlib

    member = _MEMBER[asset]
    filename = f"vdjdb/vdjdb-{tag}/{member}.gz"
    identity = hashlib.sha256(f"{repo}/{filename}".encode()).hexdigest()[:24]
    folder = cache_dir(cache) / "huggingface" / identity
    folder.mkdir(parents=True, exist_ok=True)
    out = folder / member
    with locked(out):
        if force or not valid(out):
            from huggingface_hub import hf_hub_download

            gz = hf_hub_download(
                repo_id=repo,
                repo_type="dataset",
                filename=filename,
                cache_dir=str(cache_dir(cache) / "hf_hub"),
                force_download=force,
            )
            with staging(out) as stage:
                with gzip.open(gz, "rb") as src, stage.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
                schema.normalize(_read_table(stage, parquet=False, legacy=True))
                publish(stage, out)
        metadata = {
            "source": "huggingface",
            "repo": repo,
            "tag": tag,
            "bundle": "legacy",
            "file": filename,
            "sha256": sha256(out),
        }
        with staging(folder / "release.json") as stage:
            stage.write_text(json.dumps(metadata, indent=2) + "\n")
            os.replace(stage, folder / "release.json")
    return out


def provenance(source: str | os.PathLike) -> dict:
    """Describe a resolved local reference, including fetched release identity when available."""
    path = Path(source)
    metadata = path.parent / "release.json"
    out = (
        json.loads(metadata.read_text()) if metadata.is_file() else {"source": "local"}
    )
    out = {**out, "path": str(path.resolve())}
    if path.is_file():
        out["input_sha256"] = sha256(path)
        if path.suffix == ".zip":
            import hashlib

            with zipfile.ZipFile(path) as zf:
                member = _archive_members(zf).get("vdjdb.schema.json")
                if member:
                    out["schema_sha256"] = hashlib.sha256(zf.read(member)).hexdigest()
    else:
        out["files"] = {
            str(p.relative_to(path)): sha256(p)
            for p in sorted(path.rglob("*"))
            if p.is_file()
        }
    return out


def _read_table(
    path: Path, *, parquet: bool | None = None, legacy: bool | None = None
) -> pl.DataFrame:
    is_parquet = parquet if parquet is not None else path.suffix == ".parquet"
    if is_parquet:
        return pl.read_parquet(path)
    if legacy is None:
        legacy = path.suffix == ".txt"
    return pl.read_csv(
        path, separator="\t", quote_char=None if legacy else '"', infer_schema_length=0
    )


def _join(
    left: pl.DataFrame,
    right: pl.DataFrame,
    keys: list[str],
    name: str,
    *,
    unique: bool = True,
) -> pl.DataFrame:
    for frame in (left, right):
        if not set(keys).issubset(frame.columns):
            raise ValueError(f"{name} requires join keys {keys}")
        if frame.select(
            pl.any_horizontal(
                *(pl.col(k).is_null() | (pl.col(k).cast(pl.String) == "") for k in keys)
            ).any()
        ).item():
            raise ValueError(f"{name} has missing join keys")
    if unique and right.select(keys).n_unique() != right.height:
        raise ValueError(f"{name} has duplicate join keys {keys}")
    if left.join(right.select(keys).unique(), on=keys, how="anti").height:
        raise ValueError(f"{name} has orphan reference rows")
    return left.join(
        right, on=keys, how="left", suffix="_" + name, maintain_order="left"
    )


def _require(frame: pl.DataFrame, columns: set[str], name: str):
    if missing := columns - set(frame.columns):
        raise ValueError(f"{name} requires columns {sorted(missing)}")


def _load_directory(path: Path, asset: str) -> pl.DataFrame:
    files: dict[str, Path] = {}
    for file in path.rglob("*"):
        if file.is_file() and file.name in {
            "vdjdb.parquet",
            "records.parquet",
            "records.tsv",
            "chains.parquet",
            "chains.tsv",
            "evidence.parquet",
            "evidence.tsv",
            "vdjdb.schema.json",
            "vdjdb.rearrangement.tsv",
            "vdjdb.reactivity.tsv",
            *_MEMBER.values(),
        }:
            if file.name in files:
                raise ValueError(f"ambiguous release member: {file.name}")
            files[file.name] = file
    if asset not in {"default", "primary", "legacy", "slim", "full", "airr"}:
        raise ValueError(f"unsupported VDJdb asset: {asset!r}")
    if asset in {"legacy", "slim", "full"}:
        member = _MEMBER[asset]
        if member not in files:
            raise ValueError(f"release directory missing {member}")
        return _read_table(files[member])
    registry = None
    if "vdjdb.schema.json" in files:
        registry = json.loads(files["vdjdb.schema.json"].read_text())
        if not isinstance(registry.get("fields"), list) or not isinstance(
            registry.get("tables"), dict
        ):
            raise ValueError("invalid VDJdb schema registry")
        if any(
            not isinstance(cols, list) or any(not isinstance(c, str) for c in cols)
            for cols in registry["tables"].values()
        ):
            raise ValueError("invalid VDJdb schema table columns")
        for field in registry["fields"]:
            if (
                not isinstance(field, dict)
                or not isinstance(field.get("name"), str)
                or not isinstance(field.get("ships_as"), str)
            ):
                raise ValueError("invalid VDJdb schema field names")
        if len({f["name"] for f in registry["fields"]}) != len(registry["fields"]):
            raise ValueError("duplicate VDJdb schema field names")

    def read(name: str):
        file = next(
            (files[k] for k in (name + ".parquet", name + ".tsv") if k in files), None
        )
        if file is None:
            raise ValueError(f"incomplete rich release: missing {name} table")
        frame = _read_table(file)
        if registry and name in registry["tables"]:
            names = {f["name"]: f["ships_as"] for f in registry["fields"]}
            required = {names.get(c, c) for c in registry["tables"][name]}
            if missing := required - set(frame.columns):
                raise ValueError(
                    f"{name} does not match schema; missing {sorted(missing)}"
                )
        return frame

    def evidence(out: pl.DataFrame, records: pl.DataFrame):
        if not any(k.startswith("evidence.") for k in files):
            return out
        ev = read("evidence")
        if not {"record_id", "evidence_id"}.issubset(ev.columns):
            raise ValueError("evidence requires record_id and evidence_id")
        if ev.select("record_id", "evidence_id").n_unique() != ev.height:
            raise ValueError("evidence has duplicate observation IDs")
        if ev.join(
            records.select("record_id").unique(), on="record_id", how="anti"
        ).height:
            raise ValueError("evidence has orphan record_id")
        columns = [c for c in ev.columns if c != "record_id"]
        grouped = ev.group_by("record_id").agg(pl.struct(columns).alias("evidence"))
        return out.join(grouped, on="record_id", how="left", maintain_order="left")

    if asset != "airr" and "vdjdb.parquet" in files:
        out = _read_table(files["vdjdb.parquet"])
        if not {"record_id", "gene"}.issubset(out.columns):
            raise ValueError("joined primary table requires record_id and gene")
        if out.select("record_id", "gene").n_unique() != out.height:
            raise ValueError("joined primary table has duplicate (record_id, gene)")
        return evidence(out, out)
    if asset != "airr" and any(k.startswith("records.") for k in files):
        chains, records = read("chains"), read("records")
        _require(chains, {"record_id", "gene", "cdr3"}, "chains")
        _require(records, {"record_id"}, "records")
        if chains.select("record_id", "gene").n_unique() != chains.height:
            raise ValueError("chains has duplicate (record_id, gene)")
        if records.join(
            chains.select("record_id").unique(), on="record_id", how="anti"
        ).height:
            raise ValueError("records has rows without chains")
        out = _join(chains, records, ["record_id"], "records")
        return evidence(out, records)
    if "vdjdb.rearrangement.tsv" in files or asset == "airr":
        required = {"vdjdb.rearrangement.tsv", "vdjdb.reactivity.tsv"}
        if missing := required - files.keys():
            raise ValueError(f"incomplete AIRR release: missing {sorted(missing)}")
        chains = _read_table(files["vdjdb.rearrangement.tsv"])
        activity = _read_table(files["vdjdb.reactivity.tsv"])
        if "junction_aa" not in chains.columns:
            raise ValueError("AIRR reference requires anchored junction_aa")
        if not {"reactivity_readout", "reactivity_unit", "reactivity_value"}.issubset(
            activity.columns
        ):
            raise ValueError("AIRR Reactivity missing confidence fields")
        _require(chains, {"sequence_id", "cell_id", "locus"}, "AIRR Rearrangement")
        _require(activity, {"cell_id"}, "AIRR Reactivity")
        if (
            chains["sequence_id"].n_unique() != chains.height
            or chains.select("cell_id", "locus").n_unique() != chains.height
        ):
            raise ValueError(
                "AIRR reference has duplicate sequence/cell-locus identities"
            )
        if activity.join(
            chains.select("cell_id").unique(), on="cell_id", how="anti"
        ).height:
            raise ValueError("AIRR Reactivity contains cells without rearrangements")
        return _join(chains, activity, ["cell_id"], "reactivity")
    if asset == "primary":
        raise ValueError("incomplete primary release: missing rich tables")
    if "vdjdb.txt" in files:
        return _read_table(files["vdjdb.txt"])
    raise ValueError("directory contains no supported reference tables")


def _load_archive(path: Path, asset: str) -> pl.DataFrame:
    # Extraction only into an owned temporary directory: read-only inputs work.
    with tempfile.TemporaryDirectory(prefix="vdjmatch-release-") as d:
        with zipfile.ZipFile(path) as zf:
            zf.extractall(d)
        return _load_directory(Path(d), asset)


def load(
    source: str | os.PathLike | None = None,
    *,
    asset: str = "default",
    species: str | None = None,
    gene: str | None = None,
    mhc_class: str | None = None,
    epitope: str | list[str] | None = None,
    mhc_a: str | list[str] | None = None,
    mhc_b: str | list[str] | None = None,
    mhc_match: str = "exact",
    reference_id: str | list[str] | None = None,
    exclude_reference_ids: str | list[str] | None = None,
    evidence_type: str | list[str] | None = None,
    min_score: int = 0,
    paired_only: bool = False,
    pin: str | None = None,
    cache: str | os.PathLike | None = None,
) -> pl.DataFrame:
    """Load a table, release ZIP or directory; local inputs never access the network.

    Metadata and independent observations survive normalization. Primary tables use stable
    record IDs for pairing. AIRR requires Rearrangement + Reactivity linked by cell_id and
    junction_aa; species filtering errors when the export has no receptor species.
    Peptide, MHC and study selectors accept a string or list of strings (OR within
    each selector, AND across selectors). Evidence selects any matching observation
    attached to the record and retains all of that record's chain/evidence metadata.
    A requested evidence predicate requires a supplied evidence_type field.
    ``mhc_match="compatible"`` matches family/allele restrictions with intersecting
    declared resolution, retaining original restriction labels. Exact matching is
    the default for backwards compatibility.
    """
    if mhc_match not in {"exact", "compatible"}:
        raise ValueError("mhc_match must be exact or compatible")
    path = (
        Path(source)
        if source is not None
        else fetch_latest(asset=asset, pin=pin, cache=cache)
    )
    if not path.exists():
        raise FileNotFoundError(f"reference does not exist: {path}")
    if (
        path.is_file()
        and path.with_name(path.name + ".sha256").is_file()
        and not valid(path)
    ):
        raise ValueError(f"reference local checksum mismatch: {path}")
    checksums = path.parent / "SHA256SUMS"
    if path.is_file() and checksums.is_file():
        matches = [
            line.split()[0]
            for line in checksums.read_text().splitlines()
            if len(line.split()) == 2 and line.split()[1].lstrip("*") == path.name
        ]
        if matches and (len(matches) != 1 or sha256(path) != matches[0]):
            raise ValueError(f"reference SHA256SUMS mismatch: {path}")
    if path.is_dir():
        raw = _load_directory(path, asset)
    elif path.suffix == ".zip":
        _validate_zip(path)
        raw = _load_archive(path, asset)
    elif path.name in {
        "records.tsv",
        "records.parquet",
        "chains.tsv",
        "chains.parquet",
        "vdjdb.rearrangement.tsv",
        "vdjdb.reactivity.tsv",
    }:
        raw = _load_directory(path.parent, asset)
    else:
        raw = _read_table(path)
    df = schema.normalize(raw)
    selectors = [
        ("epitope", epitope),
        ("mhc_a", mhc_a),
        ("mhc_b", mhc_b),
        ("reference_id", reference_id),
        ("exclude_reference_ids", exclude_reference_ids),
        ("evidence_type", evidence_type),
    ]
    for name, selection in selectors:
        if selection is None:
            continue
        values = [selection] if isinstance(selection, str) else selection
        if not isinstance(values, list) or any(
            not isinstance(value, str) for value in values
        ):
            raise ValueError(f"{name} must be a string or list of strings")
        if name in {"mhc_a", "mhc_b"} and mhc_match == "compatible":
            df = df.filter(pl.any_horizontal(*(schema.mhc_compatible(name, v) for v in values)) if values else pl.lit(False))
        elif name == "exclude_reference_ids":
            df = df.filter(~pl.col("reference_id").is_in(values).fill_null(False))
        elif name == "evidence_type":
            if "evidence_type" in df.columns:
                predicate = pl.col("evidence_type").is_in(values)
            elif (
                "evidence" in df.columns
                and isinstance(df.schema["evidence"], pl.List)
                and isinstance(df.schema["evidence"].inner, pl.Struct)
                and "evidence_type"
                in [field.name for field in df.schema["evidence"].inner.fields]
            ):
                predicate = (
                    pl.col("evidence")
                    .list.eval(pl.element().struct.field("evidence_type").is_in(values))
                    .list.any()
                )
            else:
                raise ValueError(
                    "reference has no evidence_type metadata; cannot apply evidence predicate"
                )
            df = df.filter(predicate.fill_null(False))
        else:
            df = df.filter(pl.col(name).is_in(values))
    if species is not None:
        if df.height and df["species"].null_count() == df.height:
            raise ValueError(
                "reference has no receptor species; omit species filter or supply a species-annotated reference"
            )
        df = df.filter(pl.col("species") == species)
    if gene is not None:
        df = df.filter(pl.col("gene") == gene)
    if mhc_class is not None:
        df = df.filter(pl.col("mhc_class") == mhc_class)
    if min_score > 0:
        df = df.filter(pl.col("vdjdb_score") >= min_score)
    if paired_only:
        df = df.filter(pl.col("complex_id").cast(pl.String) != "0")
    return df


def replicated(df: pl.DataFrame, min_refs: int = 2) -> pl.DataFrame:
    """Unique receptor/epitope keys supported by at least min_refs distinct references."""
    key = ["gene", "cdr3", "v", "j", "epitope"]
    ref = pl.col("reference_id")
    keep = (
        df.group_by(key)
        .agg(ref.filter(ref.is_not_null() & (ref != "")).n_unique().alias("n_refs"))
        .filter(pl.col("n_refs") >= min_refs)
    )
    return (
        df.join(keep, on=key, how="inner").select(*key, "mhc_class", "n_refs").unique()
    )
