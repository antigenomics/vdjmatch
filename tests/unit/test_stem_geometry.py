"""Joint V-J hairpin stem term: semantics, units and TCRdist3 reproduction.

The stem is the one geometry addition the goal scoreboard selects (``GOAL.md``). These
tests pin the properties the selection rests on, so a later change that breaks them
fails here rather than silently producing a different scorer under the same name.
"""

import numpy as np
import pytest

from vdjmatch.match.tcrdist import (STEM_DEPTH_DEFAULT, _substitution_costs,
                                    distance_matrix, paired_distance_matrix,
                                    stem_penalty)

A = ["CASSIVGGNEQFF"]
B = ["CASSIVGGNEQFF", "CASSMRSTGELFF", "WASSIVGGNEQFF", "CASSIVGGNEQFW"]
V = ["TRBV19"]


def test_identity_costs_nothing():
    assert (stem_penalty(A, A, depth=4) == 0).all()


def test_depth_zero_disables_the_term():
    assert (stem_penalty(A, B, depth=0) == 0).all()
    assert (stem_penalty(A, B, depth=4, weight=0) == 0).all()


def test_default_depth_is_four():
    assert STEM_DEPTH_DEFAULT == 4


def test_shape_and_integer_units():
    out = stem_penalty(B, A, depth=4)
    assert out.shape == (len(B), len(A))
    assert out.dtype == np.int32


def test_symmetric_in_its_arguments():
    assert (stem_penalty(A, B, depth=4).T == stem_penalty(B, A, depth=4)).all()


def test_a_level_breaks_when_either_end_differs():
    """Both the N-side and the C-side residue of level 0 gate the same level."""
    nterm = stem_penalty(A, ["WASSIVGGNEQFF"], depth=1)[0, 0]
    cterm = stem_penalty(A, ["CASSIVGGNEQFW"], depth=1)[0, 0]
    assert nterm > 0 and cterm > 0


def test_combines_as_maximum_not_sum():
    """Breaking both residues of one level costs no more than the worse of the two."""
    costs = _substitution_costs()
    both = stem_penalty(A, ["WASSIVGGNEQFW"], depth=1, weight=1.0)[0, 0]
    nterm = stem_penalty(A, ["WASSIVGGNEQFF"], depth=1, weight=1.0)[0, 0]
    cterm = stem_penalty(A, ["CASSIVGGNEQFW"], depth=1, weight=1.0)[0, 0]
    assert both == max(nterm, cterm)
    assert both < nterm + cterm


def test_bounded_by_the_shared_matrix_scale():
    """weight=1 spans the same range as one aligned residue, never Pmax times it."""
    costs = _substitution_costs()
    out = stem_penalty(B, B, depth=STEM_DEPTH_DEFAULT, weight=1.0)
    assert 0 <= out.min() and out.max() <= int(costs.max())


def test_distance_matrix_default_reproduces_tcrdist3():
    """stem_depth defaults to 0, so the shipped distance is unchanged."""
    plain = distance_matrix(A, B, V, V * len(B))
    assert (distance_matrix(A, B, V, V * len(B), stem_depth=0) == plain).all()


def test_distance_matrix_adds_exactly_the_stem():
    plain = distance_matrix(A, B, V, V * len(B))
    withstem = distance_matrix(A, B, V, V * len(B), stem_depth=4)
    assert (withstem - plain == stem_penalty(A, B, depth=4)).all()


def test_stem_never_shortens_a_distance():
    """Non-negativity is what makes re-ranking a retrieved ball exact."""
    plain = distance_matrix(A, B, V, V * len(B))
    withstem = distance_matrix(A, B, V, V * len(B), stem_depth=4)
    assert (withstem >= plain).all()


def test_rejects_junctions_shorter_than_the_stem():
    with pytest.raises(ValueError, match="length"):
        stem_penalty(["CASSIVGG"], ["CASSIVGG"], depth=8)


@pytest.mark.parametrize("bad", [-1, 2.5, True])
def test_rejects_invalid_depth(bad):
    with pytest.raises(ValueError, match="depth"):
        stem_penalty(A, B, depth=bad)


def test_rejects_negative_weight():
    with pytest.raises(ValueError, match="weight"):
        stem_penalty(A, B, depth=4, weight=-1)


def test_empty_inputs_are_empty_not_an_error():
    assert stem_penalty([], B, depth=4).shape == (0, len(B))
    assert stem_penalty(A, [], depth=4).shape == (len(A), 0)


def test_cli_exposes_the_stem_options():
    """The geometry GOAL.md selects must be reachable without editing code."""
    import argparse

    from vdjmatch.cli.tcrdist import register

    parser = argparse.ArgumentParser()
    register(parser.add_subparsers(dest="cmd"))
    args = parser.parse_args(
        ["tcrdist-neighbours", "s.tsv", "--vdjdb", "db.tsv", "--locus", "TRB", "--output-prefix", "o",
         "--stem-depth", "4"])
    assert args.stem_depth == 4
    assert args.stem_weight == 1.0


def test_cli_stem_defaults_to_tcrdist3():
    import argparse

    from vdjmatch.cli.tcrdist import register

    parser = argparse.ArgumentParser()
    register(parser.add_subparsers(dest="cmd"))
    args = parser.parse_args(
        ["tcrdist-neighbours", "s.tsv", "--vdjdb", "db.tsv", "--locus", "TRB", "--output-prefix", "o"])
    assert args.stem_depth == 0


# --- paired: the priority view, and where stem4's largest margin is measured ---

QA, QB = ["CASSIVGGNEQFF"], ["CASSMRSTGELFF"]
RA = ["CASSIVGGNEQFF", "CASSLRSTGELFF"]
RB = ["CASSMRSTGELFF", "CASSIVGGNEQFF"]
PV = dict(query_v_alpha=["TRAV12-1"], query_v_beta=["TRBV19"],
          reference_v_alpha=["TRAV12-1"] * 2, reference_v_beta=["TRBV19"] * 2)


def test_paired_default_reproduces_tcrdist3():
    plain = paired_distance_matrix(QA, QB, RA, RB, **PV)
    assert (paired_distance_matrix(QA, QB, RA, RB, **PV, stem_depth=0) == plain).all()


def test_paired_carries_one_stem_per_chain():
    """A linked pair sums two independent hairpins, not one term for the pair."""
    plain = paired_distance_matrix(QA, QB, RA, RB, **PV)
    withstem = paired_distance_matrix(QA, QB, RA, RB, **PV, stem_depth=4)
    expected = stem_penalty(QA, RA, depth=4) + stem_penalty(QB, RB, depth=4)
    assert (withstem - plain == expected).all()


def test_paired_stem_never_shortens_a_distance():
    plain = paired_distance_matrix(QA, QB, RA, RB, **PV)
    withstem = paired_distance_matrix(QA, QB, RA, RB, **PV, stem_depth=4)
    assert (withstem >= plain).all()


def test_paired_identity_pair_costs_nothing():
    out = paired_distance_matrix(QA, QB, QA, QB, query_v_alpha=["TRAV12-1"],
                                 query_v_beta=["TRBV19"],
                                 reference_v_alpha=["TRAV12-1"],
                                 reference_v_beta=["TRBV19"], stem_depth=4)
    assert out[0, 0] == 0


# --- CLI contract ---

def _args(*extra):
    import argparse

    from vdjmatch.cli.tcrdist import register

    parser = argparse.ArgumentParser()
    register(parser.add_subparsers(dest="cmd"))
    return parser.parse_args(
        ["tcrdist-neighbours", "s.tsv", "--vdjdb", "db.tsv", "--locus", "TRB",
         "--output-prefix", "o", *extra])


@pytest.mark.parametrize("flag,value", [("--stem-depth", "-1"), ("--stem-weight", "-0.5")])
def test_cli_rejects_negative_stem_parameters(flag, value, tmp_path):
    """Caught at the CLI boundary, not deep inside a distance batch."""
    from vdjmatch.cli.tcrdist import main as run

    with pytest.raises(SystemExit, match="nonnegative"):
        run(_args(flag, value))


def test_cli_refuses_stem_with_native_background():
    """A geometry silently dropped by the native counter must not look like success."""
    from vdjmatch.cli.tcrdist import main as run

    with pytest.raises(SystemExit, match="stem-depth"):
        run(_args("--stem-depth", "4", "--background", "real"))


def test_cli_manifest_declares_the_stem_geometry(tmp_path):
    """A stem run must not describe itself as the TCRdist3 default.

    The raw flags always reached manifest['parameters'], but the interpreted
    geometry did not, so a stem run and a TCRdist3 run were indistinguishable to
    anyone reading the geometry description rather than the argument dump.
    """
    import json

    import polars as pl

    from vdjmatch.cli.__main__ import main

    q, r = tmp_path / "query.tsv", tmp_path / "reference.tsv"
    meta = dict(species="HomoSapiens", epitope="E", mhc_a="HLA-A*02:01",
                mhc_b="B2M", mhc_class="MHCI")
    pl.DataFrame([dict(junction_aa="CASSLGQAYEQYF", v_call="TRBV19*01",
                       locus="TRB", **meta)]).write_csv(q, separator="\t")
    pl.DataFrame([dict(cdr3="CASSAGQAYEQYF", v="TRBV7-9*01", gene="TRB", **meta)]
                 ).write_csv(r, separator="\t")

    seen = {}
    for depth in (0, 4):
        prefix = tmp_path / f"d{depth}"
        assert main(["tcrdist-neighbours", str(q), "--vdjdb", str(r), "--locus", "TRB",
                     "--targets-from-sample", "--radius", "12", "--unbounded-nearest",
                     "--stem-depth", str(depth),
                     "--output-prefix", str(prefix)]) == 0
        seen[depth] = json.loads((tmp_path / f"d{depth}.manifest.json").read_text())

    for depth, manifest in seen.items():
        assert manifest["parameters"]["stem_depth"] == depth
        assert manifest["distance_scale"]["stem_depth"] == depth
        assert ("stem" in manifest["distance_scale"]) == bool(depth)
        assert ("joint stem" in manifest["distance"]) == bool(depth)
    assert "TCRdist3 default" in seen[0]["distance"]
    assert seen[0]["distance"] != seen[4]["distance"]


def test_reproduces_the_independently_implemented_reference_arm():
    """Benchmark-grade gate: match the measured `stem4` arm's own arithmetic.

    This is the check that caught a real defect. An earlier version scaled the term
    by Pmax, which is what the J and V-gene terms of the sweep do; the stem arm does
    not, and the extra factor made the term four times too strong -- a geometry that
    would have been reported under the name of one that was never measured.

    The reference below is transcribed from the sweep's own inner loop rather than
    from this module, and the sequences are deterministic, so the test runs anywhere
    and needs no reference database.
    """
    import random

    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    rng = random.Random(20261013)
    seqs = ["C" + "".join(rng.choice(alphabet) for _ in range(rng.randint(8, 18))) + "F"
            for _ in range(120)]
    queries, references = seqs[:60], seqs[60:]

    costs = _substitution_costs("full")
    index = {c: i for i, c in enumerate(__import__("seqtree").amino_acids())}
    depth = STEM_DEPTH_DEFAULT

    def ends(batch):
        head = np.array([[index[s[i]] for i in range(depth)] for s in batch])
        tail = np.array([[index[s[len(s) - 1 - i]] for i in range(depth)] for s in batch])
        return head, tail

    qh, qt = ends(queries)
    rh, rt = ends(references)
    stem = np.zeros((len(queries), len(references)), dtype=np.int64)
    for level in range(depth):
        stem += np.maximum(costs[qh[:, level][:, None], rh[None, :, level]],
                           costs[qt[:, level][:, None], rt[None, :, level]])
    expected = np.rint(1.0 * stem / depth).astype(np.int64)

    got = stem_penalty(queries, references, depth=depth, weight=1.0)
    assert got.shape == expected.shape
    assert np.array_equal(got, expected), "stem no longer reproduces the measured arm"
    # The Pmax regression specifically: the term must span one residue's cost range.
    assert got.max() <= int(costs.max())
