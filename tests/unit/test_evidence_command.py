"""The public evidence command shares annotation's model and preserves legacy recipes."""
from vdjmatch.cli import neighbours
from vdjmatch.cli.__main__ import main


def test_evidence_defaults_and_legacy_compatibility(monkeypatch):
    seen = []
    monkeypatch.setattr(neighbours, 'main', lambda args: seen.append(args) or 0)
    options = ['input.tsv', '--vdjdb', 'reference.tsv', '--locus', 'TRB',
               '--output-prefix', 'output', '--threads', '4']
    assert main(['evidence', *options]) == 0
    assert main(['neighbours', *options]) == 0
    evidence, legacy = seen
    assert (evidence.stem_depth, evidence.distance_matrix, evidence.junction_ends,
            evidence.neighbour_weighting, evidence.background, evidence.kernel_scorer) == (
            4, 'blosum62', 'trim3', 'kernel', 'vdjdb-other', 'v-bounded')
    assert legacy.stem_depth == 0
    assert legacy.neighbour_weighting == 'ball'
