"""Count aliases from one source must not invalidate other cohort sources."""
from pathlib import Path

import polars as pl
from polars.testing import assert_frame_equal
import pytest

from vdjmatch.runner.annotation import annotate_repertoires
from vdjmatch.cli import neighbours


@pytest.mark.parametrize('alias', ['duplicate_count', 'readcount', 'templates',
                                  'number of reads', 'count (templates/reads)', '#count'])
def test_counted_and_uncounted_sources_match_separate_runs(tmp_path, alias, monkeypatch):
    root = Path(__file__).resolve().parents[2]
    rows = pl.read_csv(root / 'examples/inputs/cohort/sample_a.tsv', separator='\t')
    counted = tmp_path / 'counted.tsv'
    uncounted = tmp_path / 'uncounted.tsv'
    rows.with_columns(pl.Series(alias, [3, 7])).write_csv(counted, separator='\t')
    rows.write_csv(uncounted, separator='\t')
    options = dict(reference=root / 'examples/inputs/reference.tsv', threads=1)
    service = neighbours.main
    service_counts = []

    def observed(args, **kwargs):
        frame = pl.read_csv(args.sample, separator='\t')
        assert alias not in frame.columns
        service_counts.append(frame['count'].to_list())
        return service(args, **kwargs)

    monkeypatch.setattr(neighbours, 'main', observed)
    together = annotate_repertoires([counted, uncounted], output_prefix=tmp_path / 'union', **options)
    separate = [annotate_repertoires(path, output_prefix=tmp_path / ('single-' + path.stem), **options)
                for path in (counted, uncounted)]
    fields = ['sample_id', 'source_query_id', 'kernel_score', 'selection_status', 'selected_epitope']
    expected = pl.concat([result['selection'].select(fields) for result in separate])
    assert_frame_equal(together['selection'].select(fields), expected, check_exact=True)
    assert service_counts == [[3, 7, 1, 1], [3, 7], [1, 1]]
    assert together['queries'][alias].cast(pl.Int64).to_list() == [3, 7, None, None]
    rows.with_columns(pl.lit(None, dtype=pl.Int64).alias(alias)).write_csv(counted, separator='\t')
    with pytest.raises(ValueError, match='without missing values'):
        annotate_repertoires([counted, uncounted], output_prefix=tmp_path / 'invalid', **options)
