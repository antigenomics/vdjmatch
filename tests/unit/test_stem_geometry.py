"""Joint V-J hairpin stem term: semantics, units and clipped BLOSUM62 geometry reproduction.

The stem is the one geometry addition the goal scoreboard selects (``GOAL.md``). These
tests pin the properties the selection rests on, so a later change that breaks them
fails here rather than silently producing a different scorer under the same name.
"""

import numpy as np
import pytest

from vdjmatch.match.geometry import (STEM_DEPTH_DEFAULT, _substitution_costs,
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


def test_distance_matrix_default_reproduces_geometry3():
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

    from vdjmatch.cli.neighbours import register

    parser = argparse.ArgumentParser()
    register(parser.add_subparsers(dest="cmd"))
    args = parser.parse_args(
        ["neighbours", "s.tsv", "--vdjdb", "db.tsv", "--locus", "TRB", "--output-prefix", "o",
         "--stem-depth", "4"])
    assert args.stem_depth == 4
    assert args.stem_weight == 1.0


def test_cli_stem_defaults_to_geometry3():
    import argparse

    from vdjmatch.cli.neighbours import register

    parser = argparse.ArgumentParser()
    register(parser.add_subparsers(dest="cmd"))
    args = parser.parse_args(
        ["neighbours", "s.tsv", "--vdjdb", "db.tsv", "--locus", "TRB", "--output-prefix", "o"])
    assert args.stem_depth == 0


# --- paired: the priority view, and where stem4's largest margin is measured ---

QA, QB = ["CASSIVGGNEQFF"], ["CASSMRSTGELFF"]
RA = ["CASSIVGGNEQFF", "CASSLRSTGELFF"]
RB = ["CASSMRSTGELFF", "CASSIVGGNEQFF"]
PV = dict(query_v_alpha=["TRAV12-1"], query_v_beta=["TRBV19"],
          reference_v_alpha=["TRAV12-1"] * 2, reference_v_beta=["TRBV19"] * 2)


def test_paired_default_reproduces_geometry3():
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

    from vdjmatch.cli.neighbours import register

    parser = argparse.ArgumentParser()
    register(parser.add_subparsers(dest="cmd"))
    return parser.parse_args(
        ["neighbours", "s.tsv", "--vdjdb", "db.tsv", "--locus", "TRB",
         "--output-prefix", "o", *extra])


@pytest.mark.parametrize("flag,value", [("--stem-depth", "-1"), ("--stem-weight", "-0.5")])
def test_cli_rejects_negative_stem_parameters(flag, value, tmp_path):
    """Caught at the CLI boundary, not deep inside a distance batch."""
    from vdjmatch.cli.neighbours import main as run

    with pytest.raises(SystemExit, match="nonnegative"):
        run(_args(flag, value))


def test_cli_manifest_declares_the_stem_geometry(tmp_path):
    """A stem run must not describe itself as the clipped BLOSUM62 geometry default.

    The raw flags always reached manifest['parameters'], but the interpreted
    geometry did not, so a stem run and a clipped BLOSUM62 geometry run were indistinguishable to
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
        assert main(["neighbours", str(q), "--vdjdb", str(r), "--locus", "TRB",
                     "--targets-from-sample", "--radius", "12", "--unbounded-nearest",
                     "--stem-depth", str(depth),
                     "--output-prefix", str(prefix)]) == 0
        seen[depth] = json.loads((tmp_path / f"d{depth}.manifest.json").read_text())

    for depth, manifest in seen.items():
        assert manifest["parameters"]["stem_depth"] == depth
        assert manifest["distance_scale"]["stem_depth"] == depth
        assert ("stem" in manifest["distance_scale"]) == bool(depth)
        assert ("joint stem" in manifest["distance"]) == bool(depth)
    assert "clipped BLOSUM62 geometry default" in seen[0]["distance"]
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


@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('weight', [0.0, 1.0, 1.5])
def test_native_stem_counts_and_mass_equal_dense(paired, weight):
    from vdjmatch.match.geometry import total_distance_count_batch, paired_total_distance_mass_batch
    q = ['CASSLGQAYEQYF', 'WASSLGQAYEQYF']
    r = q + ['CASSAGQAYEQWF', q[0], 'CASSLGQAPAYEQYF']
    qv, rv = ['TRBV19'] * 2, ['TRBV19', 'TRBV7-9', 'TRBV19', 'TRBV19', 'TRBV19']
    opts = dict(stem_depth=4, stem_weight=weight, ctrim=3)
    if paired:
        qa, ra = q[::-1], r[::-1]
        qva, rva = ['TRAV1-1'] * 2, ['TRAV1-1'] * 5
        d = 200 * paired_distance_matrix(qa, q, ra, r, qva, qv, rva, rv, **opts)
        eligible = np.array([[not (a == c and b == e) for c, e in zip(ra, r)]
                             for a, b in zip(qa, q)])
    else:
        d = 200 * distance_matrix(q, r, qv, rv, **opts)
        eligible = np.array([[a != b for b in r] for a in q])
    cuts = [[-1, 0, 200, 18000, *row.tolist(), 18000] for row in d]
    for threads in (1, 4):
        if paired:
            got, mass = paired_total_distance_mass_batch(qa, q, ra, r, qva, qv, rva, rv,
                                                        cuts, threads=threads, exclude_exact=True, **opts)
        else:
            got, mass = total_distance_count_batch(q, r, qv, rv, cuts, threads=threads,
                                                  exclude_exact=True, linear_mass=True, **opts)
        assert got == [[int(np.sum(eligible[i] & (d[i] <= t))) for t in row]
                       for i, row in enumerate(cuts)]
        assert mass == [[int(np.maximum(0, t - d[i][eligible[i]]).sum()) for t in row]
                        for i, row in enumerate(cuts)]
        shared_cuts=[-1,0,200,18000,200,0]
        if paired:
            shared=paired_total_distance_mass_batch(qa,q,ra,r,qva,qv,rva,rv,
                shared_cuts,threads=threads,exclude_exact=True,shared_thresholds=True,**opts)
            expanded=paired_total_distance_mass_batch(qa,q,ra,r,qva,qv,rva,rv,
                [shared_cuts]*len(q),threads=threads,exclude_exact=True,**opts)
        else:
            shared=total_distance_count_batch(q,r,qv,rv,shared_cuts,threads=threads,
                exclude_exact=True,linear_mass=True,shared_thresholds=True,**opts)
            expanded=total_distance_count_batch(q,r,qv,rv,[shared_cuts]*len(q),threads=threads,
                exclude_exact=True,linear_mass=True,**opts)
        assert shared==expanded


@pytest.mark.parametrize('weight', [float('nan'), float('inf')])
def test_stem_rejects_nonfinite_weight(weight):
    with pytest.raises(ValueError, match='weight'):
        stem_penalty(A, B, depth=4, weight=weight)


def test_unclipped_matrix_native_stem_matches_dense():
    from seqtree import SubstitutionMatrix, amino_acids
    from vdjmatch.match.geometry import total_distance_count_batch
    b = SubstitutionMatrix.blosum62()
    costs = np.array([[b.penalty(a, c) for c in amino_acids()] for a in amino_acids()])
    q, r = A, B
    opts = dict(stem_depth=4, ctrim=3, matrix='blosum62')
    d = distance_matrix(q, r, V, V * len(r), **opts)
    clipped = distance_matrix(q, r, V, V * len(r), stem_depth=4, ctrim=3)
    assert costs.max() == 26 and not np.array_equal(d, clipped)
    cuts = [[int(t) * 200 for t in d[0]]]
    assert total_distance_count_batch(q, r, V, V * len(r), cuts, **opts) == [
        [int((d[0] <= t // 200).sum()) for t in cuts[0]]]


def test_same_gene_counts_retain_allele_loop_cost_in_one_batch():
    from vdjmatch.match.geometry import total_distance_count_batch,distance_matrix
    q=['CASSLGQAYEQYF','CASSAGQAYEQYF']
    r=['CASSLGQAYEQYF','CASSLGQAYEQWF','CASSAGQAYEQYF']
    qv=['TRBV19*01','TRBV7-9*01']
    rv=['TRBV19*02','TRBV7-9*01','TRBV19*01']
    d=distance_matrix(q,r,qv,rv,ctrim=3,matrix='blosum62',stem_depth=4)
    t=[0,200*26,200*585,200*2600]
    same=np.array([[a.split('*')[0]==b.split('*')[0] for b in rv] for a in qv])
    expected=np.array([np.sum((200*d<=cut)&same,axis=1) for cut in t]).T
    for threads in (1,4):
        actual=total_distance_count_batch(q,r,qv,rv,[t]*len(q),ctrim=3,
            matrix='blosum62',stem_depth=4,same_v_gene=True,threads=threads)
        np.testing.assert_array_equal(actual,expected)
