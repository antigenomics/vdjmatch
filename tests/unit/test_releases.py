"""Tiny release contracts: offline formats, identity, joins and atomic bootstrap."""

import hashlib
import io
import json
import multiprocessing
import sys
import zipfile
from pathlib import Path

import polars as pl
import pytest

from vdjmatch import db
from vdjmatch.db import vdjdb
from vdjmatch.db.cache import locked, publish, staging, valid


def rich_tables(directory, parquet=False):
    records = pl.DataFrame(
        {
            "record_id": ["R1", "R2"],
            "antigen_epitope": ["PEP", "PEP"],
            "mhc_a": ["HLA-A*01:01", "HLA-A*02:01"],
            "mhc_class": ["MHCI"] * 2,
            "vdjdb_score": [2, 1],
            "reference_id": ["PMID:1", "PMID:2"],
            "species": ["HomoSapiens"] * 2,
            "method_verification": ["positive", "unknown"],
        }
    )
    chains = pl.DataFrame(
        {
            "record_id": ["R1", "R1", "R2"],
            "gene": ["TRA", "TRB", "TRB"],
            "cdr3": ["CAVF", "CASSF", "CASSF"],
            "v_segm": ["TRAV1", "TRBV1", "TRBV2"],
            "j_segm": ["TRAJ1", "TRBJ1", "TRBJ2"],
            "clonotype_id": ["CT1", "CT2", "CT3"],
            "clone_id": ["CX1", "CX1", ""],
            "cdr3nt_pgen": [0.1, 0.2, 0.3],
        }
    )
    directory.mkdir(parents=True, exist_ok=True)
    for name, frame in {"records": records, "chains": chains}.items():
        if parquet:
            frame.write_parquet(directory / (name + ".parquet"))
        else:
            frame.write_csv(directory / (name + ".tsv"), separator="\t")
    return chains.join(records, on="record_id", how="left")


@pytest.mark.parametrize("parquet", [True, False])
def test_rich_join_preserves_observations_and_metadata(tmp_path, parquet):
    rich_tables(tmp_path, parquet)
    pl.DataFrame(
        {
            "record_id": ["R1", "R1"],
            "evidence_id": ["E1", "E2"],
            "gene": ["TRA", "TRB"],
            "evidence_type": ["structure_native", "independent_study"],
        }
    ).write_csv(tmp_path / "evidence.tsv", separator="\t")
    out = db.load(tmp_path)
    assert out.height == 3
    assert out["complex_id"].to_list() == ["R1", "R1", "0"]
    assert out["v"].to_list() == ["TRAV1", "TRBV1", "TRBV2"]
    assert out["clonotype_id"].to_list() == ["CT1", "CT2", "CT3"]
    assert out["method_verification"].to_list() == ["positive", "positive", "unknown"]
    assert len(out["evidence"][0]) == 2
    assert db.load(tmp_path, paired_only=True).height == 2
    assert db.load(tmp_path, min_score=2, gene="TRB").height == 1


def test_joined_and_split_equivalent(tmp_path):
    joined = rich_tables(tmp_path)
    expected = db.load(tmp_path)
    joined.write_parquet(tmp_path / "vdjdb.parquet")
    observed = db.load(tmp_path)
    assert observed.select(
        "record_id", "gene", "cdr3", "v", "epitope", "complex_id"
    ).equals(expected.select("record_id", "gene", "cdr3", "v", "epitope", "complex_id"))
    assert db.load(tmp_path / "records.tsv").height == 3


def test_schema_registry_and_quoted_tsv(tmp_path):
    rich_tables(tmp_path)
    records = pl.read_csv(tmp_path / "records.tsv", separator="\t").with_columns(
        pl.lit('a\tb\n"c"').alias("comment")
    )
    records.write_csv(tmp_path / "records.tsv", separator="\t")
    registry = {
        "tables": {
            "records": ["record_id", "antigen.epitope"],
            "chains": ["gene", "v.segm"],
        },
        "fields": [
            {"name": "record_id", "ships_as": "record_id"},
            {"name": "antigen.epitope", "ships_as": "antigen_epitope"},
            {"name": "gene", "ships_as": "gene"},
            {"name": "v.segm", "ships_as": "v_segm"},
        ],
    }
    (tmp_path / "vdjdb.schema.json").write_text(json.dumps(registry))
    assert db.load(tmp_path)["comment"][0] == 'a\tb\n"c"'
    registry["tables"]["records"].append("missing")
    (tmp_path / "vdjdb.schema.json").write_text(json.dumps(registry))
    with pytest.raises(ValueError, match="match schema"):
        db.load(tmp_path)


@pytest.mark.parametrize("defect", ["orphan", "duplicate", "missing"])
def test_rich_join_validation(tmp_path, defect):
    rich_tables(tmp_path)
    if defect == "missing":
        (tmp_path / "records.tsv").unlink()
    else:
        table = pl.read_csv(tmp_path / "records.tsv", separator="\t")
        if defect == "duplicate":
            table = pl.concat([table, table.head(1)])
        else:
            table = table.filter(pl.col("record_id") != "R1")
        table.write_csv(tmp_path / "records.tsv", separator="\t")
    with pytest.raises(ValueError):
        db.load(tmp_path)


def airr_tables(directory):
    pl.DataFrame(
        {
            "sequence_id": ["S1", "S2"],
            "cell_id": ["R1", "R1"],
            "locus": ["TRA", "TRB"],
            "junction_aa": ["CAVF", "CASSF"],
            "cdr3_aa": ["AV", "ASS"],
            "v_call": ["TRAV1", "TRBV1"],
            "j_call": ["TRAJ1", "TRBJ1"],
        }
    ).write_csv(directory / "vdjdb.rearrangement.tsv", separator="\t")
    pl.DataFrame(
        {
            "cell_id": ["R1"],
            "reactivity_id": ["R1"],
            "peptide_sequence_aa": ["PEP"],
            "mhc_class": ["MHC-I"],
            "mhc_allele_1": ["HLA-A*02:01"],
            "reactivity_value": [3],
            "reactivity_readout": ["confidence"],
            "reactivity_unit": ["vdjdb.score"],
            "reactivity_refs": ["PMID:1"],
        }
    ).write_csv(directory / "vdjdb.reactivity.tsv", separator="\t")


def test_airr_uses_junction_and_cell_identity(tmp_path):
    airr_tables(tmp_path)
    out = db.load(tmp_path, asset="airr")
    assert out["cdr3"].to_list() == ["CAVF", "CASSF"]
    assert out["cdr3_aa"].to_list() == ["AV", "ASS"]
    assert out["complex_id"].to_list() == ["R1", "R1"]
    assert out["mhc_class"].to_list() == ["MHCI", "MHCI"]
    assert out["vdjdb_score"].to_list() == [3, 3]
    with pytest.raises(ValueError, match="no receptor species"):
        db.load(tmp_path, species="HomoSapiens")
    chains = pl.read_csv(tmp_path / "vdjdb.rearrangement.tsv", separator="\t").drop(
        "junction_aa"
    )
    chains.write_csv(tmp_path / "vdjdb.rearrangement.tsv", separator="\t")
    with pytest.raises(ValueError, match="junction_aa"):
        db.load(tmp_path)


def test_local_legacy_zip_and_read_only_directory(tmp_path, monkeypatch):
    monkeypatch.setattr(
        vdjdb, "_release_json", lambda pin: pytest.fail("offline network access")
    )
    archive = tmp_path / "legacy.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr(
            "release/vdjdb.txt",
            "gene\tcdr3\tantigen.epitope\tmethod\nTRB\tCASSF\tPEP\tcustom\n",
        )
    assert db.load(archive)["method"].to_list() == ["custom"]
    directory = tmp_path / "readonly"
    rich_tables(directory)
    for p in directory.iterdir():
        p.chmod(0o444)
    directory.chmod(0o555)
    try:
        assert db.load(directory).height == 3
    finally:
        directory.chmod(0o755)
        for p in directory.iterdir():
            p.chmod(0o644)


@pytest.mark.parametrize("member", ["../vdjdb.txt", "/vdjdb.txt", "dir\\vdjdb.txt"])
def test_unsafe_archive(tmp_path, member):
    archive = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr(member, "bad")
    with pytest.raises(ValueError, match="unsafe"):
        db.load(archive)


def test_empty_and_bare_cdr3(tmp_path):
    pl.DataFrame(
        schema={"gene": pl.String, "cdr3": pl.String, "antigen.epitope": pl.String}
    ).write_csv(tmp_path / "empty.tsv", separator="\t")
    assert db.load(tmp_path / "empty.tsv").height == 0
    pl.DataFrame({"locus": ["TRB"], "cdr3_aa": ["ASS"], "epitope": ["PEP"]}).write_csv(
        tmp_path / "bare.tsv", separator="\t"
    )
    with pytest.raises(ValueError, match="bare cdr3_aa"):
        db.load(tmp_path / "bare.tsv")


def mock_release(monkeypatch, payload, manifest=None):
    rel = {
        "tag_name": "test-tag",
        "assets": [
            {
                "name": "test.zip",
                "browser_download_url": "https://example.test/test.zip",
            }
        ],
    }
    if manifest:
        rel["assets"].append(
            {
                "name": "manifest.json",
                "browser_download_url": "https://example.test/manifest.json",
            }
        )
    monkeypatch.setattr(vdjdb, "_release_json", lambda pin: rel)
    monkeypatch.setattr(vdjdb, "_json_url", lambda url: manifest)
    calls = []

    def open_url(req, timeout):
        calls.append(req.full_url)
        return io.BytesIO(payload)

    monkeypatch.setattr(vdjdb.urllib.request, "urlopen", open_url)
    return calls


def legacy_zip():
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr(
            "release/vdjdb.txt", "gene\tcdr3\tantigen.epitope\nTRB\tCASSF\tPEP\n"
        )
    return buf.getvalue()


def test_atomic_fetch_cache_corruption_and_interruption(tmp_path, monkeypatch):
    calls = mock_release(monkeypatch, legacy_zip())
    path = db.fetch_latest(cache=tmp_path)
    assert db.load(path).height == 1
    assert db.fetch_latest(cache=tmp_path) == path and len(calls) == 1
    path.write_text("corrupted table")
    assert db.load(db.fetch_latest(cache=tmp_path)).height == 1
    archive = path.parent / "release.zip"
    archive.write_bytes(b"truncated zip")
    assert db.fetch_latest(cache=tmp_path) == path and len(calls) == 2

    def interrupted(req, timeout):
        raise OSError("interrupted")

    monkeypatch.setattr(vdjdb.urllib.request, "urlopen", interrupted)
    with pytest.raises(OSError, match="interrupted"):
        db.fetch_latest(cache=tmp_path, force=True)
    assert valid(archive) and db.load(path).height == 1
    assert not list(path.parent.glob(".release.zip.*"))
    assert db.provenance(path)["tag"] == "test-tag"


def test_manifest_role_and_checksum(tmp_path, monkeypatch):
    rich = tmp_path / "rich"
    rich_tables(rich)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for p in rich.iterdir():
            zf.write(p, "release/" + p.name)
    payload = buf.getvalue()
    manifest = {
        "tag": "test-tag",
        "bundles": [
            {
                "role": "primary",
                "file": "test.zip",
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
                "members": ["release/chains.tsv", "release/records.tsv"],
            }
        ],
    }
    calls = mock_release(monkeypatch, payload, manifest)
    path = db.fetch_latest(cache=tmp_path)
    assert path.suffix == ".zip" and len(calls) == 1
    manifest["bundles"][0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="SHA256"):
        db.fetch_latest(cache=tmp_path, force=True)
    assert db.load(path).height == 3
    manifest["tag"] = "wrong-tag"
    with pytest.raises(ValueError, match="tag"):
        db.fetch_latest(cache=tmp_path)


def test_multiple_assets_require_manifest(monkeypatch):
    monkeypatch.setattr(
        vdjdb,
        "_release_json",
        lambda pin: {
            "tag_name": "test",
            "assets": [{"name": "one.zip"}, {"name": "two.zip"}],
        },
    )
    with pytest.raises(ValueError, match="manifest"):
        db.fetch_latest()


def _bootstrap_worker(path, queue):
    path = Path(path)
    with locked(path, timeout=10):
        if not valid(path):
            with staging(path) as stage:
                stage.write_bytes(b"complete resource")
                publish(stage, path)
            queue.put("writer")
        else:
            queue.put("reader")


def test_process_coordinated_publication(tmp_path):
    context = multiprocessing.get_context("spawn")
    queue = context.Queue()
    workers = [
        context.Process(
            target=_bootstrap_worker, args=(str(tmp_path / "resource"), queue)
        )
        for _ in range(3)
    ]
    for worker in workers:
        worker.start()
    for worker in workers:
        worker.join(15)
        if worker.is_alive():
            worker.terminate()
            worker.join(2)
        assert worker.exitcode == 0
    outcomes = [queue.get(timeout=2) for _ in workers]
    assert outcomes.count("writer") == 1
    assert valid(tmp_path / "resource")


def test_hf_source_isolation_and_local_digest(tmp_path, monkeypatch):
    import gzip
    from types import SimpleNamespace

    source = tmp_path / "source.gz"
    with gzip.open(source, "wb") as fh:
        fh.write(b"gene\tcdr3\tantigen.epitope\nTRB\tCASSF\tPEP\n")
    calls = []

    def download(**kwargs):
        calls.append(kwargs)
        return str(source)

    monkeypatch.setitem(
        sys.modules, "huggingface_hub", SimpleNamespace(hf_hub_download=download)
    )
    first = db.fetch_hf(tag="tag", repo="owner/one", cache=tmp_path)
    second = db.fetch_hf(tag="tag", repo="owner/two", cache=tmp_path)
    assert first != second
    assert first.parent.parent.name == "huggingface"
    assert db.fetch_hf(tag="tag", repo="owner/one", cache=tmp_path) == first
    assert len(calls) == 2
    first.write_text("broken")
    with pytest.raises(ValueError, match="checksum mismatch"):
        db.load(first)
    assert db.load(db.fetch_hf(tag="tag", repo="owner/one", cache=tmp_path)).height == 1
    assert len(calls) == 3
    assert calls[0]["cache_dir"] == str(tmp_path / "hf_hub")


def test_concurrent_fetch_downloads_once(tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor

    calls = mock_release(monkeypatch, legacy_zip())
    with ThreadPoolExecutor(max_workers=3) as pool:
        paths = list(pool.map(lambda _: db.fetch_latest(cache=tmp_path), range(3)))
    assert len(set(paths)) == 1 and len(calls) == 1
    assert db.load(paths[0]).height == 1


def test_malformed_new_download_not_published(tmp_path, monkeypatch):
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w") as zf:
        zf.writestr("release/vdjdb.txt", "unexpected\theader\nno\tjunction\n")
    mock_release(monkeypatch, payload.getvalue())
    with pytest.raises(ValueError, match="junction"):
        db.fetch_latest(cache=tmp_path)
    assert not list(tmp_path.rglob("release.zip"))
    assert not list(tmp_path.rglob(".release.zip.*"))


def test_historical_sha256sums_checked(tmp_path, monkeypatch):
    payload = legacy_zip()
    release = {
        "tag_name": "test",
        "assets": [
            {
                "name": "test.zip",
                "browser_download_url": "https://example.test/test.zip",
            },
            {
                "name": "SHA256SUMS",
                "browser_download_url": "https://example.test/SHA256SUMS",
            },
        ],
    }
    monkeypatch.setattr(vdjdb, "_release_json", lambda pin: release)
    digest = hashlib.sha256(payload).hexdigest()

    def response(req, timeout):
        data = (
            f"{digest}  test.zip\n".encode()
            if req.full_url.endswith("SHA256SUMS")
            else payload
        )
        return io.BytesIO(data)

    monkeypatch.setattr(vdjdb.urllib.request, "urlopen", response)
    assert db.load(db.fetch_latest(cache=tmp_path)).height == 1
    digest = "0" * 64
    with pytest.raises(ValueError, match="SHA256"):
        db.fetch_latest(cache=tmp_path, force=True)


def test_reference_selectors_preserve_metadata_and_observation_order(tmp_path):
    rich_tables(tmp_path)
    selected = db.load(tmp_path, epitope="PEP", mhc_a=["HLA-A*02:01"], reference_id="PMID:2")
    assert selected["record_id"].to_list() == ["R2"]
    assert selected["method_verification"].to_list() == ["unknown"]
    assert selected["clonotype_id"].to_list() == ["CT3"]
    assert db.load(tmp_path, exclude_reference_ids=["PMID:2"])["record_id"].to_list() == ["R1", "R1"]
    assert db.load(tmp_path, reference_id=["PMID:1", "PMID:2"]).height == 3
    assert db.load(tmp_path, epitope=[]).height == 0
    with pytest.raises(ValueError, match="list of strings"):
        db.load(tmp_path, mhc_a=[2])


def test_reference_mhc_b_selection(tmp_path):
    joined = rich_tables(tmp_path).with_columns(pl.lit("B2M").alias("mhc_b"))
    path = tmp_path / "joined.tsv"
    joined.write_csv(path, separator="\t")
    assert db.load(path, mhc_b="B2M").height == 3
    assert db.load(path, mhc_b=["OTHER"]).height == 0


def test_evidence_selector_matches_record_without_discarding_other_evidence(tmp_path):
    rich_tables(tmp_path)
    pl.DataFrame({"record_id": ["R1", "R1", "R2"],
                  "evidence_id": ["E1", "E2", "E3"],
                  "evidence_type": ["structure_native", "independent_study", "other"]})\
        .write_csv(tmp_path / "evidence.tsv", separator="\t")
    selected = db.load(tmp_path, evidence_type="structure_native")
    assert selected["record_id"].to_list() == ["R1", "R1"]
    assert len(selected["evidence"][0]) == 2
    assert db.load(tmp_path, evidence_type=["structure_native", "other"]).height == 3
    assert db.load(tmp_path, evidence_type="absent").height == 0


def test_evidence_predicate_requires_metadata_even_after_other_filters_empty(tmp_path):
    rich_tables(tmp_path)
    with pytest.raises(ValueError, match="no evidence_type"):
        db.load(tmp_path, epitope="absent", evidence_type="structure_native")


def test_flat_evidence_type_selector(tmp_path):
    joined = rich_tables(tmp_path).with_columns(
        pl.when(pl.col("record_id") == "R1").then(pl.lit("positive"))
        .otherwise(pl.lit("unknown")).alias("evidence_type"))
    path = tmp_path / "joined.tsv"
    joined.write_csv(path, separator="\t")
    assert db.load(path, evidence_type="positive").height == 2
