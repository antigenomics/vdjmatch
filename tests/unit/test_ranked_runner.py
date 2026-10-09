"""Neighbourhood runner rejects retired selection and preserves paired identities."""

import polars as pl
import pytest

from vdjmatch import io
from vdjmatch.match.engine import VdjdbIndex
from vdjmatch.runner.multisample import annotate_sample


def test_retired_ranked_runner_rejected():
    with pytest.raises(ValueError, match="search_mode must be fixed or ball"):
        annotate_sample(None, "unused", search_mode="ranked")


def test_paired_invalid_beta_remains_invalid_query(monkeypatch):
    index = VdjdbIndex.build(pl.DataFrame({"complex_id": [1, 1], "gene": ["TRA", "TRB"],
                                         "cdr3": ["CAVF", "CASSF"], "epitope": ["E", "E"]}))
    q = pl.DataFrame({"query_id": [0], "cdr3a": ["CAVF"], "cdr3b": [None]},
                     schema_overrides={"cdr3b": pl.String})
    monkeypatch.setattr(io, "read_cell", lambda *args, **kw: (q, {}))
    tables = annotate_sample(index, "unused", search_mode="fixed", paired=True, with_evalue=False)
    assert tables["calls"]["vdjmatch_status"].to_list() == ["invalid_query"]
