"""Declared MHC precision selects compatible records without changing source labels."""
import polars as pl

from vdjmatch import db
from vdjmatch.db.schema import mhc_compatible


def test_family_compatibility_preserves_restriction_and_rejects_siblings(tmp_path):
    path = tmp_path / "ref.tsv"
    pl.DataFrame({
        "gene": ["TRB"] * 5, "cdr3": ["CASSLGQAYEQYF"] * 5,
        "epitope": ["E"] * 5,
        "mhc_a": ["HLA-A*02", "HLA-A*02:01", "HLA-A*02:01:48", "HLA-A*02:06", "HLA-A*03:01"],
    }).write_csv(path, separator="\t")
    family = db.load(path, mhc_a="A*02", mhc_match="compatible")
    assert family["mhc_a"].to_list() == ["HLA-A*02", "HLA-A*02:01", "HLA-A*02:01:48", "HLA-A*02:06"]
    allele = db.load(path, mhc_a="A*02:01", mhc_match="compatible")
    assert allele["mhc_a"].to_list() == ["HLA-A*02", "HLA-A*02:01", "HLA-A*02:01:48"]
    assert db.load(path, mhc_a="HLA-A*02:01")["mhc_a"].to_list() == ["HLA-A*02:01"]
    assert db.load(path, mhc_a=[], mhc_match="compatible").is_empty()


def test_compatible_class_ii_requires_both_chains_and_missing_is_not_a_match():
    frame = pl.DataFrame({
        "mhc_a": ["HLA-DQA1*05:01", "DQA1*05", "DQA1*05:01", None],
        "mhc_b": ["HLA-DQB1*02:01", "DQB1*02", "DQB1*03:01", "DQB1*02:01"],
    })
    selected = frame.filter(mhc_compatible("mhc_a", "DQA1*05") & mhc_compatible("mhc_b", "DQB1*02"))
    assert selected.height == 2
    mouse = pl.DataFrame({"mhc_a": ["H-2Kb", "H-2Kd", None]})
    assert mouse.filter(mhc_compatible("mhc_a", "H-2Kb"))["mhc_a"].to_list() == ["H-2Kb"]
