"""Offline native-oracle tests for generic batch-search row/count contracts."""

import argparse
import json

import polars as pl
import pytest
from seqtree import Index

from vdjmatch.cli import search
from vdjmatch.match import search_params


def options(tmp_path, queries=None, references=None, *extra):
    query = tmp_path / "query.tsv"
    reference = tmp_path / "reference.tsv"
    pl.DataFrame(
        {
            "aa": pl.Series(
                queries if queries is not None else ["ASS", "WWW"], dtype=pl.String
            )
        }
    ).write_csv(query, separator="\t")
    pl.DataFrame(
        {
            "reference_aa": pl.Series(
                references if references is not None else ["ASS", "ASS", "ATS"],
                dtype=pl.String,
            )
        }
    ).write_csv(reference, separator="\t")
    parser = argparse.ArgumentParser()
    search.register(parser.add_subparsers())
    return parser.parse_args(
        [
            "search",
            str(query),
            "--reference",
            str(reference),
            "--sequence-col",
            "aa",
            "--reference-sequence-col",
            "reference_aa",
            "--output-prefix",
            str(tmp_path / "out"),
            *extra,
        ]
    )


def test_native_oracle_preserves_reference_duplicates_and_unanchored_sequences(
    tmp_path,
):
    args = options(tmp_path)
    search.main(args)
    pairs = pl.read_csv(tmp_path / "out.pairs.tsv", separator="\t")
    native = Index.build(["ASS", "ASS", "ATS"], "aa").search_batch(
        ["ASS", "WWW"], search_params("1"), 1
    )
    expected = sorted(
        (row, h.ref_id, h.score, h.n_subs, h.n_ins, h.n_dels)
        for row, hs in enumerate(native)
        for h in hs
    )
    assert sorted(pairs.rows()) == expected
    assert set(pairs["reference_row"]) == {0, 1, 2}
    counts = pl.read_csv(tmp_path / "out.counts.tsv", separator="\t")
    assert counts.rows() == [(0, 3), (1, 0)]
    manifest = json.loads((tmp_path / "out.manifest.json").read_text())
    assert (
        manifest["reference_rows"] == 3 and manifest["reference_unique_sequences"] == 2
    )
    assert manifest["sequence_conventions"]["query"] == "generic amino-acid sequence"
    assert manifest["resources"]["native_batches"] == 1


def test_same_key_filter_preserves_row_identity_and_counts(tmp_path):
    args = options(tmp_path, None, None, "--same-key-col", "v")
    pl.DataFrame({"aa": ["ASS", "WWW"], "v": ["V1", "V2"]}).write_csv(
        args.query, separator="\t"
    )
    pl.DataFrame({"reference_aa": ["ASS", "ASS", "ATS"], "v": ["V2", "V1", "V1"]}).write_csv(
        args.reference, separator="\t"
    )
    search.main(args)
    pairs = pl.read_csv(tmp_path / "out.pairs.tsv", separator="\t")
    assert pairs["reference_row"].to_list() == [1, 2]
    assert pl.read_csv(tmp_path / "out.counts.tsv", separator="\t").rows() == [(0, 2), (1, 0)]
    pl.DataFrame({"aa": ["ASS"], "v": [None]}).write_csv(args.query, separator="\t")
    with pytest.raises(ValueError, match="non-null"):
        search.main(args)


def test_counts_only_exact_exclusion_and_one_batch(tmp_path, monkeypatch):
    import seqtree

    args = options(tmp_path, None, None, "--counts-only", "--exclude-exact")
    calls = []
    native_class = seqtree.Index

    class CountedIndex:
        @staticmethod
        def build(*a, **kw):
            native = native_class.build(*a, **kw)
            return CountedIndex(native)

        def __init__(self, native):
            self.native = native

        def __len__(self):
            return len(self.native)

        def search_batch(self, queries, params, threads):
            calls.append((len(queries), threads))
            return self.native.search_batch(queries, params, threads)

    monkeypatch.setattr(seqtree, "Index", CountedIndex)
    search.main(args)
    assert calls == [(2, 1)]
    assert pl.read_csv(tmp_path / "out.counts.tsv", separator="\t").rows() == [
        (0, 1),
        (1, 0),
    ]
    assert not (tmp_path / "out.pairs.tsv").exists()


def test_empty_reference_and_empty_query_have_typed_outputs(tmp_path):
    args = options(tmp_path, ["ASS"], [])
    search.main(args)
    assert pl.read_csv(tmp_path / "out.counts.tsv", separator="\t").rows() == [(0, 0)]
    assert pl.read_csv(tmp_path / "out.pairs.tsv", separator="\t").columns == [
        "query_row",
        "reference_row",
        "score",
        "n_subs",
        "n_ins",
        "n_dels",
    ]
    args = options(tmp_path, [], ["ASS"])
    search.main(args)
    assert pl.read_csv(tmp_path / "out.counts.tsv", separator="\t").height == 0
    assert pl.read_csv(tmp_path / "out.pairs.tsv", separator="\t").height == 0


@pytest.mark.parametrize("role", ["query", "reference"])
def test_invalid_rows_fail_with_counts(tmp_path, role):
    args = options(
        tmp_path,
        [None, "AX", ""] if role == "query" else ["ASS"],
        [None, "AX", ""] if role == "reference" else ["ASS"],
    )
    with pytest.raises(
        ValueError,
        match=role + " has 3 invalid amino-acid rows: 1 null, 2 empty/nonstandard",
    ):
        search.main(args)
    assert not (tmp_path / "out.counts.tsv").exists()


def test_threads_engines_and_supported_matrix_are_explicit(tmp_path):
    args = options(tmp_path)
    search.main(args)
    serial = pl.read_csv(tmp_path / "out.pairs.tsv", separator="\t").sort(
        "query_row", "reference_row"
    )
    args.threads = 4
    search.main(args)
    assert serial.equals(
        pl.read_csv(tmp_path / "out.pairs.tsv", separator="\t").sort(
            "query_row", "reference_row"
        )
    )
    args.threads = 0
    with pytest.raises(ValueError, match="positive integer"):
        search.main(args)
    args.threads, args.engine, args.matrix = 1, "seqtrie", "vdjam"
    with pytest.raises(ValueError, match="requires --engine seqtm"):
        search.main(args)
    args.engine = "seqtm"
    search.main(args)
    from vdjmatch.match.scoring import DEFAULT_SCALE, load_vdjam

    native = Index.build(["ASS", "ASS", "ATS"], "aa").search_batch(
        ["ASS", "WWW"],
        search_params(
            "1,0,0,1",
            matrix=load_vdjam(),
            gap_open=DEFAULT_SCALE,
            gap_extend=DEFAULT_SCALE,
        ),
        1,
    )
    expected = sorted(
        (row, h.ref_id, h.score, h.n_subs, h.n_ins, h.n_dels)
        for row, hs in enumerate(native)
        for h in hs
    )
    assert (
        sorted(pl.read_csv(tmp_path / "out.pairs.tsv", separator="\t").rows())
        == expected
    )


def test_seqtrie_total_radius_null_edits_and_exact_exclusion(tmp_path):
    args = options(
        tmp_path,
        ["ASS"],
        ["ASS", "ATS", "ASSS", "AS"],
        "--engine",
        "seqtrie",
        "--exclude-exact",
    )
    search.main(args)
    pairs = pl.read_csv(tmp_path / "out.pairs.tsv", separator="\t")
    assert sorted(pairs.select("query_row", "reference_row", "score").rows()) == [
        (0, 1, 1),
        (0, 2, 1),
        (0, 3, 1),
    ]
    assert all(pairs[c].null_count() == 3 for c in ("n_subs", "n_ins", "n_dels"))
    assert pl.read_csv(tmp_path / "out.counts.tsv", separator="\t").rows() == [(0, 3)]
    manifest = json.loads((tmp_path / "out.manifest.json").read_text())
    assert "per-edit caps ignored" in manifest["native_scope_contract"]
    assert manifest["edit_decomposition_available"] is False
    args.threads = 4
    search.main(args)
    assert pairs.equals(pl.read_csv(tmp_path / "out.pairs.tsv", separator="\t"))


def test_input_table_cannot_be_overwritten_by_output_prefix(tmp_path):
    args = options(tmp_path)
    query = tmp_path / "out.counts.tsv"
    query.write_bytes((tmp_path / "query.tsv").read_bytes())
    args.query = str(query)
    before = query.read_bytes()
    with pytest.raises(ValueError, match="overwrite an input"):
        search.main(args)
    assert query.read_bytes() == before


def test_duplicate_query_rows_default_reference_column_and_junction_label(tmp_path):
    args = options(tmp_path, ["CASSF", "CASSF"], ["CASSF"])
    pl.DataFrame({"junction_aa": ["CASSF", "CASSF"]}).write_csv(
        args.query, separator="\t"
    )
    pl.DataFrame({"junction_aa": ["CASSF"]}).write_csv(args.reference, separator="\t")
    args.sequence_col, args.reference_sequence_col = "junction_aa", None
    search.main(args)
    assert pl.read_csv(tmp_path / "out.counts.tsv", separator="\t").rows() == [
        (0, 1),
        (1, 1),
    ]
    assert pl.read_csv(tmp_path / "out.pairs.tsv", separator="\t").select(
        "query_row", "reference_row"
    ).rows() == [(0, 0), (1, 0)]
    manifest = json.loads((tmp_path / "out.manifest.json").read_text())
    assert manifest["sequence_conventions"] == {
        "query": "declared junction_aa",
        "reference": "declared junction_aa",
    }


@pytest.mark.parametrize("matrix_name", ["blosum62", "vdjam"])
@pytest.mark.parametrize("exclude_exact", [False, True])
def test_position_significance_exact_native_pssm_agreement(
    tmp_path, matrix_name, exclude_exact
):
    from seqtree import PositionalMatrix, SubstitutionMatrix
    from vdjmatch.match.regions import load_significance, significance_weights
    from vdjmatch.match.scoring import load_vdjam

    queries = ["ASS", "ATSS", "ASS", "WWW"]
    references = ["ASS", "ATS", "AAS", "ATS", "ATSS", "ATTS", "ASSS", "WWA", "AS"]
    args = options(
        tmp_path,
        queries,
        references,
        "--scope",
        "2,0,0,2",
        "--matrix",
        matrix_name,
        "--position-significance",
    )
    args.exclude_exact = exclude_exact
    search.main(args)
    matrix = load_vdjam() if matrix_name == "vdjam" else SubstitutionMatrix.blosum62()
    expected = []
    # Native oracle groups by length solely because a PSSM has fixed width.
    sig = load_significance()
    for length in sorted({len(q) for q in queries}):
        weights = [max(1, round(100 * w)) for w in significance_weights(length, sig)]
        pssm = PositionalMatrix.from_weights(matrix, weights)
        qrows = [i for i, q in enumerate(queries) if len(q) == length]
        params = search_params("2,0,0,2", pos_matrix=pssm)
        native = Index.build(references, "aa").search_batch(
            [queries[i] for i in qrows], params, 1
        )
        expected.extend(
            (qrows[row], h.ref_id, h.score, h.n_subs, h.n_ins, h.n_dels)
            for row, hs in enumerate(native)
            for h in hs
            if not exclude_exact or h.n_subs > 0
        )
    pairs = pl.read_csv(tmp_path / "out.pairs.tsv", separator="\t")
    assert sorted(pairs.rows()) == sorted(expected)
    counts = pl.read_csv(tmp_path / "out.counts.tsv", separator="\t")
    assert counts["n_hits"].to_list() == [
        sum(row[0] == i for row in expected) for i in range(len(queries))
    ]


def test_position_significance_rejects_unsupported_options_and_counts_skip_scoring(
    tmp_path, monkeypatch
):
    args = options(tmp_path, None, None, "--position-significance")
    for matrix, scope in (("none", "1,0,0,1"), ("blosum62", "1,1,0,1")):
        args.matrix, args.scope = matrix, scope
        with pytest.raises(ValueError, match="position significance requires"):
            search.main(args)
    args.matrix, args.scope, args.counts_only = "blosum62", "1,0,0,1", True

    def unexpected(*a):
        pytest.fail("counts-only must not construct positional scores")

    monkeypatch.setattr(search, "_positional_scores", unexpected)
    search.main(args)
    assert pl.read_csv(tmp_path / "out.counts.tsv", separator="\t").rows() == [
        (0, 3),
        (1, 0),
    ]
    manifest = json.loads((tmp_path / "out.manifest.json").read_text())
    assert manifest["scores_reported"] is False
