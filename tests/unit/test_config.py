"""Unit tests for vdjmatch.config.Params and the CLI param-resolution / reproducibility wiring."""
import argparse
import json

import polars as pl

from vdjmatch.config import Params
from vdjmatch.cli import __main__ as cli


# --- Params dataclass round-trip ---
def test_params_defaults():
    p = Params.defaults()
    assert p.scope == "1,0,0,1" and p.matrix == "vdjam" and p.min_score == 0
    assert p.evalue is True and p.align is True and p.species == "HomoSapiens" and p.threads == 0


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
