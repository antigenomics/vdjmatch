"""Ranked CLI runner retains the shared validation and identity contracts."""

import polars as pl
import pytest

from vdjmatch import io
from vdjmatch.match.engine import VdjdbIndex
from vdjmatch.runner.multisample import annotate_sample


def test_automatic_ranked_control_rejects_species_mismatch(monkeypatch):
    index = VdjdbIndex.build(pl.DataFrame({"gene": ["TRB"], "cdr3": ["CASSF"],
                                         "epitope": ["E"], "species": ["MusMusculus"]}))
    q = pl.DataFrame({"query_id": [0], "cdr3": ["CASSF"], "locus": ["TRB"]})
    monkeypatch.setattr(io, "read_rearrangement", lambda *args, **kw: (q, {}))
    with pytest.raises(ValueError, match="species does not match"):
        annotate_sample(index, "unused", search_mode="ranked", species="human")


def test_paired_invalid_beta_remains_invalid_query(monkeypatch):
    index = VdjdbIndex.build(pl.DataFrame({"complex_id": [1, 1], "gene": ["TRA", "TRB"],
                                         "cdr3": ["CAVF", "CASSF"], "epitope": ["E", "E"]}))
    q = pl.DataFrame({"query_id": [0], "cdr3a": ["CAVF"], "cdr3b": [None]},
                     schema_overrides={"cdr3b": pl.String})
    monkeypatch.setattr(io, "read_cell", lambda *args, **kw: (q, {}))
    tables = annotate_sample(index, "unused", search_mode="ranked", paired=True, with_evalue=False)
    assert tables["calls"]["vdjmatch_status"].to_list() == ["invalid_query"]
