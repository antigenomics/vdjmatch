"""Exercise released model probabilities through the VDJmatch batch bridge."""
import pytest

from vdjmatch.precursor.model import load_model, pgen
from vdjtools.model import native


@pytest.mark.parametrize('organism,locus,seq', [
    ('human', 'TRA', 'CAVRDSNYQLIW'), ('human', 'TRB', 'CASSIRSSYEQYF'),
    ('mouse', 'TRA', 'CAVRDSNYQLIW'), ('mouse', 'TRB', 'CASSIRSSYEQYF')])
def test_released_gene_aliases_scalar_and_batch(organism, locus, seq):
    model = load_model(locus, source='isalgo', organism=organism)
    v = model.genomic['genes_v']['gene'][0].split('*')[0]
    j = model.genomic['genes_j']['gene'][0].split('*')[0]
    expected = native.pgen_aa(model, seq, v=v, j=j)
    assert pgen(model, [seq, seq], v=[v, v+'*99'], j=[j, j+'*99'], threads=2) == [expected]*2
    nt = 'TGTGCCGTGCGTGATAGCAATTATCAGCTGATTTGG'
    nt_expected = native.pgen_nt(model, nt, v=v, j=j)
    assert native.pgen_nt_batch(model, [nt, nt], v=[v, v+'*99'], j=[j, j+'*99'], threads=2) == [nt_expected]*2
    assert pgen(model, [seq], threads=1)[0] > 0
    with pytest.raises(KeyError):
        pgen(model, [seq], v=['UNKNOWN'], threads=1)


@pytest.mark.parametrize('alias', ['learned', 'arda'])
def test_current_alias_probabilities(alias):
    seq = 'CASSIRSSYEQYF'
    assert pgen(load_model('TRB', source=alias), [seq], threads=1) == pgen(
        load_model('TRB', source='isalgo'), [seq], threads=1)


def test_zero_mass_warning_describes_selected_model(tmp_path, monkeypatch, capsys):
    import polars as pl
    import vdjmatch.precursor as precursor
    from vdjmatch.cli.__main__ import main

    path = tmp_path / 'junctions.tsv'
    pl.DataFrame({'cdr3': ['CASSIRSSYEQYF']}).write_csv(path, separator='\t')
    monkeypatch.setattr(precursor, 'summarise', lambda *a, **kw:
        pl.DataFrame({'n_zero_pgen': [1], 'n_junctions': [1]}))
    assert main(['precursor', str(path), '--source', 'isalgo', '-o', str(tmp_path/'out.tsv')]) == 0
    warning = capsys.readouterr().err
    assert '--source isalgo' in warning
    assert "This is OLGA's" not in warning
    assert 'does not inherit it' not in warning
