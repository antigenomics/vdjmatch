"""Unit tests for vdjmatch.config.Params and the CLI param-resolution / reproducibility wiring."""
import argparse
import json

import polars as pl
import pytest

from vdjmatch.config import Params
from vdjmatch.cli import __main__ as cli


def test_cli_fresh_raw_control_and_exact_exclusion(tmp_path):
    reference, query, control = (tmp_path / name for name in ("reference.tsv", "query.tsv", "control.tsv"))
    pl.DataFrame({"cdr3": ["CASSF", "CATSF"], "gene": ["TRB"] * 2,
                  "species": ["HomoSapiens"] * 2,
                  "antigen.epitope": ["PEP"] * 2}).write_csv(reference, separator="\t")
    pl.DataFrame({"sequence_id": ["input"], "junction_aa": ["CASSF"],
                  "locus": ["TRB"]}).write_csv(query, separator="\t")
    pl.DataFrame({"cdr3aa": ["CASSF", "CASRF"]}).write_csv(control, separator="\t")
    prefix = tmp_path / "out"
    arguments = ["match", str(query), "--vdjdb", str(reference), "--control", f"TRB={control}",
                 "--fresh-control", "--exclude-exact", "--no-align", "-o", str(prefix)]
    assert cli.main(arguments) == 0
    candidates = pl.read_csv(tmp_path / "out.query.candidates.txt", separator="\t")
    assert candidates["n_target"].to_list() == [1]
    assert candidates["n_control"].to_list() == [1]
    manifest = json.loads((tmp_path / "out.query.manifest.json").read_text())
    assert manifest["exclude_exact"] is True
    assert manifest["controls"]["TRB"]["source_rows"] == 2
    with pytest.raises(SystemExit):
        cli.main(arguments + ["--control", f"TRB={control}"])
    for flags in (["--search-mode", "ranked"], ["--top-k", "10"]):
        with pytest.raises(SystemExit):
            cli.main(arguments + flags)
    assert cli.main(arguments + ["--search-mode", "ball", "-o", str(tmp_path / "ball")]) == 0
    manifest = json.loads((tmp_path / "ball.query.manifest.json").read_text())
    assert manifest["radii"] == [1, 2, 3, 4, 5]
    assert "ranked_model" not in manifest and "top_k" not in manifest



# --- Params dataclass round-trip ---
def test_params_defaults():
    p = Params.defaults()
    assert p.scope == "1,0,0,1" and p.matrix == "vdjam" and p.min_score == 0
    assert p.evalue is True and p.align is True and p.species == "HomoSapiens" and p.threads == 1


def test_params_json_roundtrip(tmp_path):
    p = Params(scope="2,1,1,2", matrix="none", min_score=1, match_v=True, match_j=True,
               evalue=False, align=False, species="MusMusculus", threads=8)
    path = tmp_path / "params.json"
    p.to_json(path)
    assert path.exists()
    loaded = Params.from_json(path)
    assert loaded == p


def test_params_from_json_ignores_unknown_keys(tmp_path):
    path = tmp_path / "extra.json"
    path.write_text(json.dumps({"scope": "3,0,0,3", "bogus": 42}))
    p = Params.from_json(path)
    assert p.scope == "3,0,0,3" and p.matrix == "vdjam"   # unknown key dropped, rest defaulted


# --- CLI override-beats-JSON ---
def _ns(**kw):
    base = {"config": None}
    base.update(kw)
    return argparse.Namespace(**base)


@pytest.mark.parametrize("mode", ["ball"])
def test_nonfixed_mode_rejects_custom_scope(mode):
    with pytest.raises(ValueError, match="scope applies to fixed"):
        cli._cmd_match(_ns(search_mode=mode, scope="2,0,0,2"))


def test_resolve_params_no_config_no_flags_is_defaults():
    assert cli._resolve_params(_ns()) == Params.defaults()


def test_resolve_params_loads_config(tmp_path):
    cfg = tmp_path / "c.json"
    Params(scope="2,1,1,2", threads=4).to_json(cfg)
    p = cli._resolve_params(_ns(config=str(cfg)))
    assert p.scope == "2,1,1,2" and p.threads == 4


def test_resolve_params_cli_overrides_json(tmp_path):
    cfg = tmp_path / "c.json"
    Params(scope="2,1,1,2", matrix="none", threads=4, evalue=False).to_json(cfg)
    # explicitly-passed flags (present on the namespace) win over the JSON values
    p = cli._resolve_params(_ns(config=str(cfg), scope="1,0,0,1", threads=16))
    assert p.scope == "1,0,0,1"          # CLI override
    assert p.threads == 16               # CLI override
    assert p.matrix == "none"            # untouched -> JSON value kept
    assert p.evalue is False             # untouched -> JSON value kept


def test_resolve_params_store_const_flags():
    # --no-evalue / --match-v map to const False/True on the namespace
    p = cli._resolve_params(_ns(evalue=False, match_v=True))
    assert p.evalue is False and p.match_v is True
    assert p.align is True                # not passed -> default kept


# --- <prefix>.params.json is written by a run ---
def test_cmd_match_writes_params_json(tmp_path, monkeypatch):
    # stub the heavy/network bits so the test is fast, deterministic, and offline
    monkeypatch.setattr(cli.db, "fetch_latest", lambda **k: tmp_path / "s.tsv")
    monkeypatch.setattr(cli.db, "provenance", lambda path: {"source":"test"})
    monkeypatch.setattr(cli.db, "load", lambda *a, **k: pl.DataFrame({"gene": ["TRB"]}))

    class _Idx:
        genes = ["TRB"]
    monkeypatch.setattr(cli.match.VdjdbIndex, "build", classmethod(lambda cls, *a, **k: _Idx()))
    monkeypatch.setattr(cli.match, "load_vdjam", lambda: "M")
    empty = pl.DataFrame()
    monkeypatch.setattr(cli, "annotate_sample",
                        lambda *a, **k: {"hits": empty, "summary": empty, "calls": empty})

    sample = tmp_path / "s.tsv"
    sample.write_text("junction_aa\tlocus\nCASSF\tTRB\n")
    prefix = tmp_path / "run" / "out"
    cfg = tmp_path / "c.json"
    Params(scope="2,1,1,2", min_score=1).to_json(cfg)

    rc = cli._cmd_match(_ns(config=str(cfg), scope="1,0,0,1", verbose=False,
                            samples=[str(sample)], output_prefix=str(prefix),
                            vdjdb=None, asset="full", pin=None))
    assert rc == 0
    written = Params.from_json(f"{prefix}.params.json")
    assert written.scope == "1,0,0,1"    # CLI override recorded
    assert written.min_score == 1        # from JSON config


def test_raw_airr_control_rejects_bare_cdr3_with_only_j_call(tmp_path):
    from vdjmatch.evalue.control import raw_background

    path = tmp_path / "control.tsv"
    pl.DataFrame({"cdr3_aa": ["ASS"], "j_call": ["TRBJ1-1"]}).write_csv(
        path, separator="\t"
    )
    with pytest.raises(ValueError, match="requires junction_aa"):
        raw_background("TRB", "human", path)
