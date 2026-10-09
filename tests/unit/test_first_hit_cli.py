"""Small offline parity/identity checks for the published global first-hit CLI."""

import argparse
import json

import polars as pl
import pytest
from seqtree import Index

from vdjmatch.cli import first_hit as cli
from vdjmatch.evalue import first_hit
from vdjmatch.match.scope import search_params


def files(tmp_path, samples=None):
    reference = tmp_path / "reference.tsv"
    pl.DataFrame(
        {
            "cdr3": ["CASSF"] * 3 + ["CATSF", "CATSF", "CAVF"],
            "gene": ["TRB"] * 6,
            "species": ["HomoSapiens"] * 6,
            "v.segm": ["TRBV1"] * 6,
            "j.segm": ["TRBJ1"] * 6,
            "antigen.epitope": ["PEP"] * 6,
            "reference.id": ["P1", "P1", "P2", "P1", "P2", "P1"],
        }
    ).write_csv(reference, separator="\t")
    sample = tmp_path / "query.tsv"
    if samples is None:
        samples = {
            "query_id": ["same", "same", "far", "invalid", "missing", "other"],
            "sequence_id": ["s1", "s2", "s3", "s4", "s5", "s6"],
            "junction_aa": ["CASSF", "CASSF", "WWWWWW", "CA*F", None, "CASSF"],
            "locus": ["TRB"] * 5 + ["TRA"],
            "count": [1, 2, 3, 4, 5, 6],
        }
    pl.DataFrame(samples).write_csv(sample, separator="\t")
    control = tmp_path / "control.tsv"
    pl.DataFrame({"cdr3aa": ["CASSF", "CASRF", "CASRF", "CA*F", None]}).write_csv(
        control, separator="\t"
    )
    return sample, reference, control


def args(tmp_path, samples=None, *extra):
    sample, reference, control = files(tmp_path, samples)
    parser = argparse.ArgumentParser()
    cli.register(parser.add_subparsers())
    return parser.parse_args(
        [
            "first-hit",
            str(sample),
            "--vdjdb",
            str(reference),
            "--control",
            str(control),
            "--scope",
            "1,0,0,1",
            "--output-prefix",
            str(tmp_path / "out"),
            *extra,
        ]
    )


def test_oracle_identity_shortlist_and_statuses(tmp_path):
    options = args(tmp_path)
    assert cli.main(options) == 0
    out = pl.read_csv(tmp_path / "out.evidence.tsv", separator="\t")
    assert out["query_id"].to_list() == [
        "same",
        "same",
        "far",
        "invalid",
        "missing",
        "other",
    ]
    assert out["sequence_id"].to_list() == ["s1", "s2", "s3", "s4", "s5", "s6"]
    assert out["count"].to_list() == [1, 2, 3, 4, 5, 6]
    assert out["status"].to_list() == [
        "matched",
        "matched",
        "no_hit",
        "invalid_sequence",
        "missing_sequence",
        "locus_mismatch",
    ]
    assert out["p_any"][2] is None  # No selected radius: the public oracle omits p_any.
    target = Index.build(["CASSF", "CATSF"], "aa")
    control = Index.build(["CASRF", "CASSF"], "aa")
    th, cc = first_hit.scan(
        target,
        [""] * 2,
        control,
        ["CASSF", "CASSF", "WWWWWW"],
        params=search_params("1,0,0,1"),
        threads=1,
    )
    for i in range(3):
        expected = first_hit.pvalue(th[i], cc[i], 2, 2)
        for column in ("radius", "n_target", "n_control", "E", "p_enrichment"):
            assert out[column][i] == expected[column]
    manifest = json.loads((tmp_path / "out.manifest.json").read_text())
    assert manifest["reference"] == {
        "source_records": 6,
        "selected_records": 5,
        "valid_records": 5,
        "unique_junctions": 2,
    }
    assert manifest["control"]["excluded_rows"] == 2
    assert manifest["control"]["unique_junctions"] == 2
    assert manifest["queries"]["searched_rows"] == 3


def test_exact_exclusion_and_zero_control_rule(tmp_path):
    options = args(
        tmp_path, {"junction_aa": ["CASSF"], "locus": ["TRB"]}, "--exclude-exact"
    )
    # Both query/control exact hits are excluded; the nearest target is distance1.
    pl.DataFrame({"cdr3aa": ["CASSF"]}).write_csv(options.control, separator="\t")
    cli.main(options)
    row = pl.read_csv(tmp_path / "out.evidence.tsv", separator="\t").row(0, named=True)
    assert row["radius"] == 1
    assert row["n_target"] == 1
    assert row["n_control"] == 0
    assert row["rule_of_three"] is True
    assert row["E"] == 6.0


def test_single_batch_for_thousand_original_rows(tmp_path, monkeypatch):
    options = args(tmp_path, {"junction_aa": ["CASSF"] * 1000, "locus": ["TRB"] * 1000})
    original = first_hit.scan
    calls = []

    def scan(*a, **kw):
        calls.append((len(a[3]), kw["chunk"]))
        return original(*a, **kw)

    monkeypatch.setattr(first_hit, "scan", scan)
    cli.main(options)
    assert calls == [(1000, 1000)]
    assert pl.read_csv(tmp_path / "out.evidence.tsv", separator="\t").height == 1000


@pytest.mark.parametrize("value", ["0", "-1"])
def test_invalid_threads_fail_before_input_read(tmp_path, value):
    options = args(tmp_path, None, "--threads", value)
    with pytest.raises(ValueError, match="threads"):
        cli.main(options)


def test_nonbundled_control_requires_explicit_raw_input(tmp_path):
    options = args(tmp_path)
    options.control = None
    options.locus = "TRA"
    with pytest.raises(ValueError, match="supply --control"):
        cli._control(options)


def test_thread_equality_and_all_missing_rows(tmp_path):
    options = args(tmp_path)
    cli.main(options)
    serial = pl.read_csv(tmp_path / "out.evidence.tsv", separator="\t")
    options.threads = 4
    cli.main(options)
    assert pl.read_csv(tmp_path / "out.evidence.tsv", separator="\t").equals(serial)
    pl.DataFrame(
        {"sequence_id": ["empty"], "junction_aa": [None], "count": [7]},
        schema_overrides={"junction_aa": pl.String},
    ).write_csv(options.sample, separator="\t")
    cli.main(options)
    out = pl.read_csv(tmp_path / "out.evidence.tsv", separator="\t")
    assert out["status"].to_list() == ["missing_sequence"]
    assert out["count"].to_list() == [7]
    assert out["query_id"].to_list() == [0]


def test_invalid_counts_fail_instead_of_disappearing(tmp_path):
    options = args(tmp_path, {"junction_aa": ["CA*F"], "count": [-1]})
    with pytest.raises(ValueError, match="non-negative integer"):
        cli.main(options)


def test_native_vdjtools_counts_are_preserved(tmp_path):
    options = args(tmp_path)
    pl.DataFrame(
        {
            "#count": [7],
            "freq": [1.0],
            "cdr3nt": ["TGT"],
            "cdr3aa": ["CASSF"],
            "v": ["TRBV1"],
            "d": ["TRBD1"],
            "j": ["TRBJ1"],
        }
    ).write_csv(options.sample, separator="\t")
    cli.main(options)
    out = pl.read_csv(tmp_path / "out.evidence.tsv", separator="\t")
    assert out["count"].to_list() == [7]
