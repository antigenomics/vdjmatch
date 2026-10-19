"""Synthetic cohort collection and identity restoration; no scientific computation."""
import json
from pathlib import Path

import polars as pl
import pytest

from vdjmatch.runner.annotation import annotate_repertoires


def table(path, *, locus="TRB", paired=False):
    frame = (pl.DataFrame({"clone_id": ["007", "007"], "clonotype_id": ["type-a", "type-b"],
        "query_id": ["001", "002"], "sequence_id": ["seq-a", "seq-b"], "locus": ["TRA", "TRB"],
        "v": ["TRAV1", "TRBV1"], "j": ["TRAJ1", "TRBJ1"], "junction_aa": ["CAVF", "CASSF"]}) if paired else
        pl.DataFrame({"query_id": ["001"], "sequence_id": ["seq-a"], "clonotype_id": ["type-a"],
            "locus": [locus], "v": ["TRAV1" if locus == "TRA" else "TRBV1"], "junction_aa": ["CAVF"]}))
    frame.write_csv(path, separator="\t")
    return path


@pytest.fixture
def core(monkeypatch):
    from vdjmatch.cli import neighbours as geometry
    calls = []
    def run(a, *, v_loop_distances, marginal_tables=None):
        assert v_loop_distances['species']=='human'
        raw = pl.read_csv(a.sample, separator="\t", infer_schema_length=0)
        calls.append((a, raw))
        if a.locus == "paired":
            q = raw.select("clone_id").unique(maintain_order=True).with_row_index("query_id")
        else:
            q = raw.select("query_id", "sequence_id")
        q.with_columns(pl.lit("matched").alias("status")).write_csv(a.output_prefix + ".queries.tsv", separator="\t")
        q.with_columns(pl.lit("supported").alias("selection_status")).write_csv(a.output_prefix + ".selection.tsv", separator="\t")
        Path(a.output_prefix + ".manifest.json").write_text(json.dumps({"locus": a.locus, "ingestion": {"input_rows": raw.height}}))
        return 0
    monkeypatch.setattr(geometry, "main", run)
    return calls


def test_directory_batches_once_per_population_and_restores_ids(tmp_path, core):
    inputs = tmp_path / "inputs"; inputs.mkdir()
    table(inputs / "b.tsv"); table(inputs / "a.tsv", locus="TRA")
    out = annotate_repertoires(inputs, reference="unused", output_prefix=tmp_path / "out", threads=2)
    assert [a.locus for a, _ in core] == ["TRA", "TRB"]
    assert len({i for _, raw in core for i in raw["query_id"]}) == 2
    assert out["queries"]["sample_id"].to_list() == ["a", "b"]
    assert out["queries"]["source_query_id"].to_list() == ["001", "001"]
    assert out["queries"]["source_sequence_id"].to_list() == ["seq-a", "seq-a"]
    assert out["queries"]["clonotype_id"].to_list() == ["type-a", "type-a"]
    assert all(a.selection_only and a.stem_depth == 4 and a.threads == 2 for a, _ in core)
    assert all(a.reference_population == 'curated-benchmark' for a, _ in core)
    assert Path(str(tmp_path / "out") + ".queries.tsv").exists()
    assert Path(str(tmp_path / "out") + ".selection.tsv").exists()
    assert out["manifest"]["native_threads"] == 2
    assert all(not Path(a.sample).exists() for a, _ in core)


def test_all_valid_baseline_remains_explicit(tmp_path, core):
    source = table(tmp_path / 'input.tsv')
    annotate_repertoires(source, reference='unused', output_prefix=tmp_path / 'baseline',
                         reference_population='all-valid')
    assert core[0][0].reference_population == 'all-valid'
    assert core[0][0].kernel_scorer == 'v-bounded'
    assert not core[0][0].exclude_exact


def test_pairs_are_scoped_to_samples_and_keep_both_ids(tmp_path, core):
    a = table(tmp_path / "one.tsv", paired=True); b = table(tmp_path / "two.tsv", paired=True)
    out = annotate_repertoires([a, b], reference="unused", output_prefix=tmp_path / "out")
    assert len(core) == 1 and core[0][0].locus == "paired"
    assert core[0][1]["clone_id"].n_unique() == 2
    q = out["queries"]
    assert q["clone_id"].to_list() == ["007", "007"]
    assert q["source_query_ida"].to_list() == ["001", "001"]
    assert q["source_query_idb"].to_list() == ["002", "002"]
    assert q["clonotype_ida"].to_list() == ["type-a", "type-a"]
    assert q["clonotype_idb"].to_list() == ["type-b", "type-b"]


def test_metadata_paths_relative_to_metadata_and_options(tmp_path, core):
    table(tmp_path / "input.tsv")
    pl.DataFrame({"path": ["input.tsv"], "sample_id": ["007"], "condition": ["x"]}).write_csv(tmp_path / "metadata.tsv", separator="\t")
    out = annotate_repertoires(metadata=tmp_path / "metadata.tsv", reference="unused",
                              output_prefix=tmp_path / "out", stem_weight=.5, selection_only=False)
    assert out["queries"]["sample_id"].to_list() == ["007"]
    assert out["queries"]["condition"].to_list() == ["x"]
    assert core[0][0].stem_weight == .5 and not core[0][0].selection_only


def test_single_clone_ids_do_not_invent_pairing(tmp_path, core):
    path = table(tmp_path / "single.tsv")
    pl.read_csv(path, separator="\t", infer_schema_length=0).with_columns(
        pl.lit("007").alias("clone_id")).write_csv(path, separator="\t")
    out = annotate_repertoires(path, reference="unused", output_prefix=tmp_path / "out")
    assert len(core) == 1 and core[0][0].locus == "TRB"
    assert out["queries"]["clone_id"].to_list() == ["007"]


@pytest.mark.parametrize("threads", [True, 0, -1, 1.5])
def test_invalid_budget_precedes_scoring(tmp_path, core, threads):
    with pytest.raises(ValueError, match="threads"):
        annotate_repertoires(table(tmp_path / "one.tsv"), reference="unused", threads=threads)
    assert not core


def test_duplicate_names_paths_missing_files_and_unknown_options_fail_first(tmp_path, core):
    a = table(tmp_path / "one.tsv")
    other = tmp_path / "other"; other.mkdir(); b = table(other / "one.tsv")
    for inputs in [[a, a], [a, b], [a, tmp_path / "missing.tsv"]]:
        with pytest.raises((ValueError, FileNotFoundError)):
            annotate_repertoires(inputs, reference="unused", output_prefix=tmp_path / "out")
    with pytest.raises(ValueError, match="unknown"):
        annotate_repertoires(a, reference="unused", undefined_option=True)
    assert not core


def test_empty_input_has_typed_standard_outputs(tmp_path, core):
    path = table(tmp_path / "empty.tsv")
    pl.read_csv(path, separator="\t", infer_schema_length=0).head(0).write_csv(path, separator="\t")
    out = annotate_repertoires(path, reference="unused", output_prefix=tmp_path / "out")
    assert not core
    assert out["queries"].height == out["selection"].height == 0
    assert set(out["output_paths"]) == {"queries", "selection"}


def test_actual_wide_core_preserves_pairs_and_marginal_identities(tmp_path):
    meta = dict(species="HomoSapiens", epitope="E", mhc_a="HLA-A*02:01", mhc_b="B2M", mhc_class="MHCI")
    ref = tmp_path / "reference.tsv"
    pl.DataFrame([dict(gene=locus, cdr3=seq, v=v, complex_id=1, **meta)
        for locus, seq, v in [("TRA", "CAVRDSNYQLIW", "TRAV12-1"),
                              ("TRB", "CASSLGQAYEQYF", "TRBV19")]]).write_csv(ref, separator="\t")
    path = tmp_path / "wide.tsv"
    pl.DataFrame({"clone_id": ["007", "008"], "clonotype_id": ["type-a", "type-b"],
        "query_id": ["q007", "q008"], "junction_aa_alpha": ["CAVRDSNYQLIW", None],
        "junction_aa_beta": ["CASSLGQAYEQYF"] * 2, "v_alpha": ["TRAV12-1", None],
        "v_beta": ["TRBV19"] * 2}).write_csv(path, separator="\t")
    out = annotate_repertoires(path, reference=ref, output_prefix=tmp_path / "out")
    assert out["queries"].height == out["selection"].height == 2
    assert out["queries"]["clone_id"].to_list() == ["007", "008"]
    assert out["queries"]["source_query_ida"].to_list() == ["q007", "q008"]
    assert out["queries"]["source_query_idb"].to_list() == ["q007", "q008"]
    assert {"TRA.selection", "TRB.selection"} <= out.keys()
    for lane in ["TRA", "TRB"]:
        assert set(out[lane + ".selection"]["source_query_id"]) == {"q007", "q008"}


@pytest.mark.parametrize("complete", [True, False])
def test_cell_linkage_precedes_shared_clone_identity(tmp_path, core, complete):
    loci = ["TRA", "TRB", "TRA", "TRB"] if complete else ["TRA", "TRB"]
    cells = ["cell-a", "cell-a", "cell-b", "cell-b"] if complete else ["cell-a", "cell-b"]
    path = tmp_path / "cells.tsv"
    pl.DataFrame({"cell_id": cells, "clone_id": ["clone-x"] * len(loci),
        "clonotype_id": ["type-y"] * len(loci), "sequence_id": [str(i) for i in range(len(loci))],
        "locus": loci, "junction_aa": ["CAVF" if x == "TRA" else "CASSF" for x in loci]}).write_csv(path, separator="\t")
    out = annotate_repertoires(path, reference="unused", output_prefix=tmp_path / "out")
    assert [a.locus for a, _ in core] == (["paired"] if complete else ["TRA", "TRB"])
    assert out["queries"].height == 2
    assert out["queries"]["clone_id"].to_list() == ["clone-x"] * 2
    assert out["queries"]["pair_id"].to_list() == ["cell-a", "cell-b"]
    if complete:
        assert core[0][1]["clone_id"].n_unique() == 2


def test_output_paths_cannot_overwrite_inputs(tmp_path, core):
    path = table(tmp_path / 'out.queries.tsv')
    before = path.read_bytes()
    with pytest.raises(ValueError, match='output.*input'):
        annotate_repertoires(path, reference='unused', output_prefix=tmp_path / 'out')
    assert path.read_bytes() == before and not core


@pytest.mark.parametrize('column', ['source_query_id', 'source_sequence_id'])
def test_reserved_source_identity_is_rejected_before_service(tmp_path, core, column):
    path = table(tmp_path / 'one.tsv')
    pl.read_csv(path, separator='\t', infer_schema_length=0).with_columns(
        pl.lit('pretend-source').alias(column)).write_csv(path, separator='\t')
    with pytest.raises(ValueError, match='reserved'):
        annotate_repertoires(path, reference='unused', output_prefix=tmp_path / 'out')
    assert not core


def test_source_generated_metadata_columns_are_preserved(tmp_path, core):
    path = table(tmp_path / 'one.tsv')
    pl.read_csv(path, separator='\t', infer_schema_length=0).with_columns(
        pl.lit('raw-status').alias('status')).write_csv(path, separator='\t')
    result = annotate_repertoires(path, reference='unused', output_prefix=tmp_path / 'out')
    assert result['queries']['status'].to_list() == ['matched']
    assert result['queries']['input_status'].to_list() == ['raw-status']


def test_initialized_table_pool_exceeding_library_budget_fails(tmp_path):
    import os
    import subprocess
    import sys
    code = """import polars as pl
assert pl.thread_pool_size()==4
from vdjmatch import annotate_repertoires
try:
    annotate_repertoires('missing.tsv', reference='unused', output_prefix='unwritten', threads=1)
except ValueError as exc:
    assert 'POLARS_MAX_THREADS' in str(exc), str(exc)
else:
    raise AssertionError('oversized pool accepted')
"""
    subprocess.run([sys.executable, '-c', code], cwd=tmp_path,
        env=dict(os.environ, POLARS_MAX_THREADS='4'), timeout=20, check=True)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize('role', ['reference', 'metadata', 'symlink'])
def test_all_protected_paths_precede_service(tmp_path, core, role):
    path = table(tmp_path / 'raw.tsv')
    collision = tmp_path / 'out.selection.tsv'
    kwargs = dict(reference='unused')
    if role == 'reference':
        table(collision)
        kwargs['reference'] = collision
    elif role == 'metadata':
        pl.DataFrame({'path': [str(path)]}).write_csv(collision, separator='\t')
        kwargs['metadata'] = collision
        path = None
    else:
        collision.symlink_to(path)
    with pytest.raises(ValueError, match='output.*input'):
        annotate_repertoires(path, output_prefix=tmp_path / 'out', **kwargs)
    assert not core


def test_metadata_reserved_identity_rejected(tmp_path, core):
    path = table(tmp_path / 'one.tsv')
    metadata = pl.DataFrame({'path': [str(path)], 'source_pair_id': ['pretend']})
    with pytest.raises(ValueError, match='reserved identity'):
        annotate_repertoires(metadata=metadata, reference='unused', output_prefix=tmp_path / 'out')
    assert not core


@pytest.mark.parametrize('value', ['0', '-1', 'bad', '1.5'])
def test_invalid_explicit_table_budget_rejected(tmp_path, core, monkeypatch, value):
    monkeypatch.setenv('POLARS_MAX_THREADS', value)
    with pytest.raises(ValueError, match='POLARS_MAX_THREADS'):
        annotate_repertoires('missing.tsv', reference='unused', output_prefix=tmp_path / 'out')
    assert not core
