"""Run-local normalized references preserve independently loaded population scores."""
from pathlib import Path

import polars as pl
from polars.testing import assert_frame_equal
import pytest

from vdjmatch import db
from vdjmatch.cli import neighbours
from vdjmatch.runner.annotation import annotate_repertoires


@pytest.mark.parametrize('population', ['curated-benchmark', 'all-valid'])
@pytest.mark.parametrize('mixed', [False, True])
def test_mixed_populations_load_once_and_preserve_selection(tmp_path, monkeypatch, population, mixed):
    root = Path(__file__).resolve().parents[2]
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    rows = pl.DataFrame(dict(sequence_id=['a', 'b'], pair_id=['p', 'p'],
        locus=['TRA', 'TRB'], v_call=['TRAV12-1', 'TRBV19'],
        j_call=['TRAJ33', 'TRBJ2-7'], junction_aa=['CAVRDSNYQLIW', 'CASSIRSSYEQYF']))
    rows.write_csv(inputs / 'paired.tsv', separator='\t')
    if mixed:
        rows.drop('pair_id').write_csv(inputs / 'single.tsv', separator='\t')
    calls = []
    load = db.load

    def observed(*args, **kwargs):
        calls.append(kwargs)
        return load(*args, **kwargs)

    monkeypatch.setattr(db, 'load', observed)
    options = dict(reference=root / 'examples/inputs/reference.tsv', threads=1,
                   reference_population=population)
    actual = annotate_repertoires(inputs, output_prefix=tmp_path / 'shared', **options)
    assert len(calls) == 1
    assert actual['manifest']['population_service_calls'] == (3 if mixed else 1)
    service = neighbours.main

    def independent(args, **kwargs):
        kwargs.pop('reference_loader', None)
        return service(args, **kwargs)

    monkeypatch.setattr(neighbours, 'main', independent)
    calls.clear()
    expected = annotate_repertoires(inputs, output_prefix=tmp_path / 'independent', **options)
    assert len(calls) > 1
    assert_frame_equal(actual['selection'], expected['selection'], check_exact=True)
