"""Saved-contact rendering: sequence indices must map to author numbers and insertion codes."""
from pathlib import Path
import importlib.util
import polars as pl
import pytest

pytest.importorskip('tcren')
spec = importlib.util.spec_from_file_location('nlv_panels', Path(__file__).with_name('structures_nlv.py'))
nlv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nlv)


def test_saved_contact_author_mapping(tmp_path):
    atoms = []
    residues = [('B', 40, '', 'ALA'), ('B', 42, 'A', 'VAL'), ('B', 47, '', 'GLY')]
    three = {v: k for k, v in nlv.THREE_TO_ONE.items()}
    residues += [('P', 10+i, '', three[aa]) for i, aa in enumerate(nlv.PEPTIDE)]
    for serial, (chain, number, insertion, aa) in enumerate(residues, 1):
        atoms.append(f'ATOM  {serial:5d}  CA  {aa:3s} {chain}{number:4d}{insertion:1s}   '
                     f'{float(serial):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}           C  ')
    pdb = tmp_path/'saved.pdb'
    pdb.write_text('\n'.join(atoms)+'\nEND\n')
    contacts = pl.DataFrame({'residue.index.from':[1], 'residue.aa.from':['V'],
        'residue.index.to':[4], 'residue.aa.to':['M'], 'dist':[4.0], 'tcren':[-.2],
        'atom.from':['CA'], 'atom.to':['CA']})
    mapped, context = nlv.map_contact_residues(pdb, 'B', 'P', contacts)
    assert mapped['author_residue.from'].to_list() == ['42A']
    assert mapped['author_residue.to'].to_list() == ['14']
    assert mapped['peptide_position'].to_list() == [5]
    assert mapped.select(contacts.columns).equals(contacts)
    assert context == '40+42A+47'
    nlv.WORK = tmp_path
    script = nlv.write_pymol_script(pdb, 'B', 'P', mapped, tmp_path/'panel.png', beta_context=context).read_text()
    assert 'resi 42A and name CA' in script and 'resi 14 and name CA' in script
    assert 'resi 1 and' not in script
    with pytest.raises(ValueError, match='residue identity'):
        nlv.map_contact_residues(pdb, 'B', 'P', contacts.with_columns(pl.lit('G').alias('residue.aa.from')))
