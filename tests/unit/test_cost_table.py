"""Explicit substitution-cost tables: lossless round trip and rejected malformed input."""
import numpy as np
import pytest
from seqtree import amino_acids

from vdjmatch.match import geometry
from vdjmatch.match.scoring import load_cost_table

AA = "ACDEFGHIKLMNPQRSTVWY"


def write_table(path, costs, symbols):
    position = {a: i for i, a in enumerate(symbols)}
    with open(path, "w") as fh:
        fh.write("aa.1\taa.2\tcost\n")
        for i, a in enumerate(AA):
            for b in AA[i + 1:]:
                fh.write(f"{a}\t{b}\t{int(costs[position[a], position[b]])}\n")
    return path


def test_round_trip_reproduces_the_shipped_clipped_geometry(tmp_path):
    """A cost table written from the default geometry must read back bit-identical.

    This is the contract a derived matrix relies on: costs are consumed directly, never through
    the lossy Gram transform, so what was fitted is what the engine scores.
    """
    symbols = amino_acids()
    base = geometry._substitution_costs()
    path = write_table(tmp_path / "costs.tsv", base, symbols)
    back = geometry._substitution_costs("full", str(path))
    index = [symbols.index(c) for c in AA]
    np.testing.assert_array_equal(base[np.ix_(index, index)], back[np.ix_(index, index)])
    assert back.dtype == np.int32
    np.testing.assert_array_equal(back, back.T)
    assert not np.diag(back).any()
    assert geometry.distance_units(back) == geometry.distance_units(base)


def test_blosum62_cost_table_round_trips_over_the_whole_alphabet(tmp_path):
    """A table written from BLOSUM62 must read back as BLOSUM62 across all 24 symbols.

    The padding rows matter: V-loop templates carry characters outside the standard twenty, so a
    table that guessed those rows would change the V-loop distances feeding the kernel and make two
    matrices incomparable for reasons that have nothing to do with the 20x20 under study.
    """
    symbols = amino_acids()
    base = geometry._substitution_costs("full", "blosum62")
    back = geometry._substitution_costs(
        "full", str(write_table(tmp_path / "b62.tsv", base, symbols)))
    np.testing.assert_array_equal(base, back)


def test_padding_rows_are_shared_between_different_tables(tmp_path):
    """Two different 20x20 tables differ only on the standard residues, never on padding."""
    symbols = amino_acids()
    one = geometry._substitution_costs("full", "blosum62")
    other = geometry._substitution_costs("full", "pam250")
    a = geometry._substitution_costs("full", str(write_table(tmp_path / "a.tsv", one, symbols)))
    b = geometry._substitution_costs("full", str(write_table(tmp_path / "b.tsv", other, symbols)))
    padding = [symbols.index(s) for s in "BZX*"]
    standard = [symbols.index(s) for s in AA]
    scale = b[np.ix_(standard, standard)].max() / a[np.ix_(standard, standard)].max()
    np.testing.assert_allclose(b[np.ix_(padding, standard)],
                               np.rint(a[np.ix_(padding, standard)] * scale), atol=1)


@pytest.mark.parametrize("body,message", [
    ("aa.1\taa.2\tcost\nA\tC\t-1\n", "non-negative integers"),
    ("aa.1\taa.2\tcost\nA\tC\t1.5\n", "non-negative integers"),
    ("aa.1\taa.2\tcost\nA\tA\t3\n", "zero on the diagonal"),
    ("aa.1\taa.2\tcost\n", "empty"),
    ("aa.1\tcost\nA\tC\t1\n", "header"),
])
def test_malformed_tables_are_rejected(tmp_path, body, message):
    path = tmp_path / "bad.tsv"
    path.write_text(body)
    with pytest.raises(ValueError, match=message):
        load_cost_table(path)
