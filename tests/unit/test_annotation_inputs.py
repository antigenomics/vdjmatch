"""Synthetic cohort collection and identity restoration; no scientific computation."""
import json
import os
from pathlib import Path
import subprocess
import sys

import polars as pl
import pytest

from vdjmatch.runner.annotation import annotate_repertoires


def test_public_api_starts_without_thread_environment(tmp_path):
    root = Path(__file__).resolve().parents[2]
    code = """
import sys
from vdjmatch import annotate_repertoires
assert 'polars' not in sys.modules
result = annotate_repertoires(sys.argv[1], reference=sys.argv[2], output_prefix=sys.argv[3])
assert result['queries'].height == 2
assert result['manifest']['reference_population'] == 'curated-benchmark'
assert result['manifest']['score_models'] == ['vx']
assert 'call_status' not in result['selection'].columns
from vdjmatch.runner import annotate_sample
assert callable(annotate_sample)
"""
    env = {k: v for k, v in os.environ.items() if k not in {
        'POLARS_MAX_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS',
        'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'}}
    result = subprocess.run([sys.executable, '-c', code,
        str(root / 'examples/inputs/cohort/sample_a.tsv'),
        str(root / 'examples/inputs/reference.tsv'), str(tmp_path / 'fresh')],
        env=env, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr


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
    def run(a, *, v_loop_distances, marginal_tables=None, reference_loader=None):
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
    out = annotate_repertoires(inputs, reference="unused", score_model="kernel", output_prefix=tmp_path / "out", threads=2)
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
    annotate_repertoires(source, reference='unused', score_model="kernel", output_prefix=tmp_path / 'baseline',
                         reference_population='all-valid')
    assert core[0][0].reference_population == 'all-valid'
    assert core[0][0].kernel_scorer == 'v-bounded'
    assert not core[0][0].exclude_exact


def test_pairs_are_scoped_to_samples_and_keep_both_ids(tmp_path, core):
    a = table(tmp_path / "one.tsv", paired=True); b = table(tmp_path / "two.tsv", paired=True)
    out = annotate_repertoires([a, b], reference="unused", score_model="kernel", output_prefix=tmp_path / "out")
    assert len(core) == 1 and core[0][0].locus == "paired"
    assert core[0][1]["clone_id"].n_unique() == 2
    q = out["queries"]
    assert q["clone_id"].to_list() == ["007", "007"]
    assert q["source_query_ida"].to_list() == ["001", "001"]
    assert q["source_query_idb"].to_list() == ["002", "002"]
    assert q["clonotype_ida"].to_list() == ["type-a", "type-a"]
    assert q["clonotype_idb"].to_list() == ["type-b", "type-b"]


def test_ordinary_paired_mouse_annotation_defaults_to_vx(tmp_path):
    from types import SimpleNamespace
    from zipfile import ZipFile

    from vdjmatch.cli.paired_vx_base import main as replay_paired
    from vdjmatch.match.geometry import load_v_loops

    model = load_v_loops("mouse")
    alpha, beta = [next(v for v, entry in model.items() if entry["locus"] == lane)
                   for lane in ("TRA", "TRB")]
    reference = tmp_path / "reference.tsv"
    pl.DataFrame({
        "gene": ["TRA", "TRB"], "cdr3": ["CASAAAAAFF", "CASSRRRRFF"],
        "v": [alpha, beta], "j": ["TRAJ1*01", "TRBJ1-1*01"], "epitope": ["PEPTIDE"] * 2,
        "mhc_a": ["H2-Kb"] * 2, "mhc_b": ["B2M"] * 2, "mhc_class": ["MHCI"] * 2,
        "species": ["MusMusculus"] * 2, "complex_id": [1, 1],
    }).write_csv(reference, separator="\t")
    pairs = tmp_path / "pairs.tsv"
    pl.DataFrame({"pair_id": ["007", "008"], "junction_aa_alpha": [" casaaaaaff ", "CASAAAAAFF"],
                  "junction_aa_beta": ["cassrrrrff", "CASSRRRRFF"], "v_alpha": [alpha] * 2, "v_beta": [beta] * 2,
                  "species": ["mouse", "human"]}
                 ).write_csv(pairs, separator="\t")
    out = annotate_repertoires(pairs, reference=reference, species="MusMusculus",
                               output_prefix=tmp_path / "paired")
    assert out["manifest"]["score_models"] == ["vx"]
    assert out["manifest"]["population_service_calls"] == 3
    assert out["manifest"]["vx_marginal_chain_rows"] == {"TRA": 2, "TRB": 2}
    assert out["selection"].select("pair_id", "score", "selected_epitope").rows() == [
        ("007", 0., "PEPTIDE"), ("008", None, None)]
    assert out["queries"]["status"].to_list() == ["matched", "unselected_species"]
    unpaired = tmp_path / "unpaired.tsv"
    pl.read_csv(reference, separator="\t").head(1).write_csv(unpaired, separator="\t")
    out = annotate_repertoires(pairs, reference=unpaired, species="MusMusculus",
                               output_prefix=tmp_path / "empty-reference")
    assert out["queries"]["status"].to_list() == ["no_reference", "unselected_species"]
    assert out["selection"]["score"].null_count() == 2
    assert out["TRA.selection"].select("source_pair_id", "score", "selected_epitope").rows() == [
        ("007", 0., "PEPTIDE"), ("008", None, None)]
    assert "call_status" not in out["TRA.selection"].columns

    release = tmp_path / "release.zip"
    with ZipFile(release, "w") as archive:
        archive.write(reference, "vdjdb.txt")
    null_query = tmp_path / "null.tsv"
    pl.read_csv(pairs, separator="\t", schema_overrides={"pair_id": pl.String}).head(1).drop(
        "species").with_columns(pl.lit("CASAAAAAFF").alias("junction_aa_alpha"),
                                pl.lit("CASSRRRRFF").alias("junction_aa_beta"),
                                pl.col("pair_id").alias("query_id")).write_csv(null_query, separator="\t")
    null_dir = tmp_path / "null"
    assert replay_paired(SimpleNamespace(queries=null_query, vdjdb=release, null_source="generated",
                                        output_dir=null_dir, exclude_exact=False, withhold_peptides=None,
                                        threads=1)) == 0
    import json

    from vdjmatch.cli import competitive
    from vdjmatch.cli.paired_vx_evalue import _source
    from vdjmatch.db.cache import sha256
    from vdjmatch.evalue import competitive as reducer

    ranking_dir = tmp_path / "raw-ranking"
    ranking_args = SimpleNamespace(queries=null_query, vdjdb=release, exclude_exact=False,
                                   output_dir=ranking_dir, threads=1)
    assert competitive.mouse_paired_main(ranking_args) == 0
    with pytest.raises(ValueError, match="fresh output directory"):
        competitive.mouse_paired_main(ranking_args)
    _, binding, *_ = _source(ranking_dir)
    _, null_binding, *_ = _source(null_dir)
    assert binding == null_binding
    ranks = pl.read_csv(ranking_dir / "ranking.tsv", separator="\t")
    assert ranks.height == 1 and ranks["score"].item() == 0
    rank_manifest = json.loads((ranking_dir / "ranking.manifest.json").read_text())
    assert rank_manifest["ranking_sha256"] == sha256(ranking_dir / "ranking.tsv")
    assert rank_manifest["base_manifest_sha256"] == sha256(ranking_dir / "manifest.json")
    assert rank_manifest["motif_reducer_source_sha256"] == sha256(Path(reducer.__file__))

    out = annotate_repertoires(pairs, reference=release, species="mouse", vx_null_dir=null_dir,
                               output_prefix=tmp_path / "calibrated")
    assert out["selection"].select("score", "E_null", "null_tail_count", "call_status").rows() == [
        (0., 1., 1, "no_supported_match"), (None, 0., 0, "no_supported_match")]
    assert out["TRA.selection"]["score"].to_list() == [0., None]
    assert "call_status" not in out["TRA.selection"].columns
    exact_excluded = annotate_repertoires(pairs, reference=release, species="mouse", exclude_exact=True,
                                         output_prefix=tmp_path / "exact-excluded")
    assert exact_excluded["selection"]["score"].to_list() == [None, None]
    assert exact_excluded["TRA.selection"]["score"].to_list() == [None, None]
    with pytest.raises(ValueError, match="exact-hit policy"):
        annotate_repertoires(pairs, reference=release, species="mouse", vx_null_dir=null_dir,
                             output_prefix=tmp_path / "incompatible", exclude_exact=True)
    (tmp_path / "collision.selection.tsv").hardlink_to(null_dir / "maxima.tsv")
    with pytest.raises(ValueError, match="output path conflicts"):
        annotate_repertoires(pairs, reference=release, species="mouse", vx_null_dir=null_dir,
                             output_prefix=tmp_path / "collision")

    single = tmp_path / "single.tsv"
    pl.DataFrame({"query_id": ["01", "02", "03"], "locus": ["TRA"] * 3,
                  "junction_aa": [" casaaaaaff ", "ZZ", "CASAAAAAFF"], "v_call": [alpha] * 3,
                  "species": ["mouse", "mouse", "human"]}).write_csv(single, separator="\t")
    out = annotate_repertoires(single, reference=release, species="mouse",
                               output_prefix=tmp_path / "single-vx")
    assert out["queries"]["status"].to_list() == ["matched", "invalid_query_or_gene", "unselected_species"]
    assert out["selection"].select("source_query_id", "score", "selected_epitope").rows() == [
        ("01", 0., "PEPTIDE"), ("02", None, None), ("03", None, None)]
    ranked = annotate_repertoires([pairs, single], reference=release, species="mouse",
                                  output_prefix=tmp_path / "mixed-ranking")
    assert ranked["manifest"]["score_models"] == ["vx"]
    assert ranked["selection"]["query_id"].n_unique() == ranked["selection"].height == 5
    assert "call_status" not in ranked["selection"].columns
    with pytest.raises(ValueError, match="null locus"):
        annotate_repertoires(single, reference=release, species="mouse", score_model="vx", vx_null_dir=null_dir,
                             output_prefix=tmp_path / "wrong-locus")
    single_null = tmp_path / "single-null.tsv"
    first = pl.read_csv(single, separator="\t", schema_overrides={"query_id": pl.String}).head(1)
    pl.concat([first, first.with_columns(pl.lit("04").alias("query_id"))]).write_csv(single_null, separator="\t")
    single_null_dir = tmp_path / "single-null"
    assert replay_paired(SimpleNamespace(queries=single_null, vdjdb=release, null_source="generated",
                                        output_dir=single_null_dir, exclude_exact=False, withhold_peptides=None,
                                        locus="TRA", threads=1)) == 0
    with pytest.raises(ValueError, match="one matching null directory"):
        annotate_repertoires([pairs, single], reference=release, species="mouse", score_model="vx",
                             vx_null_dir=null_dir, output_prefix=tmp_path / "missing-bank")
    mixed = annotate_repertoires([pairs, single], reference=release, species="mouse",
                                 vx_null_dir=[single_null_dir, null_dir],
                                 output_prefix=tmp_path / "mixed-vx")
    assert mixed["selection"]["query_id"].n_unique() == mixed["selection"].height == 5
    assert set(mixed["selection"]["locus"]) == {"TRA", "paired"}
    assert mixed["selection"].filter(pl.col("locus") == "TRA")["null_queries"].to_list() == [2] * 3
    assert mixed["selection"].filter(pl.col("locus") == "paired")["null_queries"].to_list() == [1] * 2
    assert mixed["selection"]["E_null"].sum() == 2.
    assert mixed["manifest"]["population_service_calls"] == 3
    assert mixed["TRA.selection"]["source_pair_id"].to_list() == ["007", "008"]
    assert "call_status" not in mixed["TRA.selection"].columns
    explicit = annotate_repertoires([pairs, single], reference=release, species="mouse", score_model="vx",
                                    vx_tie_policy="abstain", vx_null_dir=[null_dir, single_null_dir],
                                    output_prefix=tmp_path / "mixed-explicit")
    assert mixed["selection"].equals(explicit["selection"])
    with pytest.raises(ValueError, match="requires the vX annotation score"):
        annotate_repertoires(single, reference=release, species="mouse", score_model="kernel",
                             vx_null_dir=single_null_dir, output_prefix=tmp_path / "kernel-with-bank")
    with pytest.raises(ValueError, match="one matching null directory"):
        annotate_repertoires(single, reference=release, species="mouse", vx_null_dir=[],
                             output_prefix=tmp_path / "empty-bank")
    with pytest.raises(ValueError, match="ordinary annotation options"):
        annotate_repertoires(single, reference=release, species="mouse", vx_null_dir=single_null_dir,
                             output_prefix=tmp_path / "advanced-with-bank", selection_only=False)

    broad = pl.read_csv(reference, separator="\t")
    narrow = broad.with_columns(pl.lit(2).alias("complex_id"), pl.lit("H2-Kb:01").alias("mhc_a"))
    literal_reference = tmp_path / "restriction-hierarchy.zip"
    with ZipFile(literal_reference, "w") as archive:
        archive.writestr("vdjdb.txt", pl.concat([broad, narrow], how="vertical_relaxed").write_csv(separator="\t"))
    literal_null = tmp_path / "literal-null"
    assert replay_paired(SimpleNamespace(queries=null_query, vdjdb=literal_reference, null_source="generated",
                                        output_dir=literal_null, exclude_exact=False, withhold_peptides=None,
                                        maxima_only=True, tie_policy="literal", threads=1)) == 0
    out = annotate_repertoires(pairs, reference=literal_reference, species="mouse",
                               vx_null_dir=literal_null,
                               output_prefix=tmp_path / "literal")
    assert out["selection"].select("tied", "tie_resolved", "selected_mhc_a", "E_null").rows() == [
        (True, True, "H2-Kb:01", 1.), (False, False, None, 0.)]
    explicit_literal = annotate_repertoires(pairs, reference=literal_reference, species="mouse",
                                            score_model="vx", vx_null_dir=literal_null, vx_tie_policy="literal",
                                            output_prefix=tmp_path / "literal-explicit")
    assert out["selection"].equals(explicit_literal["selection"])
    with pytest.raises(ValueError, match="null tie policy"):
        annotate_repertoires(pairs, reference=literal_reference, species="mouse",
                             vx_null_dir=literal_null, vx_tie_policy="abstain",
                             output_prefix=tmp_path / "literal-with-abstention")

    with pytest.raises(ValueError, match="ordinary annotation options"):
        annotate_repertoires(pairs, reference=reference, species="MusMusculus",
                             output_prefix=tmp_path / "modified", score_model="vx", distance_matrix="clipped")


def test_vx_cli_custom_reference_preserves_exact_hit_policy(tmp_path):
    from vdjmatch.cli.__main__ import main

    inputs = Path(__file__).resolve().parents[2] / "examples" / "inputs"
    retained = tmp_path / "retained"
    excluded = tmp_path / "excluded"
    command = [
        "annotate", str(inputs / "cohort" / "sample_a.tsv"), "--locus", "TRB",
        "--vdjdb", str(inputs / "reference.tsv"), "--threads", "1",
    ]
    assert main([*command, "-o", str(retained)]) == 0
    assert main([*command, "-o", str(excluded), "--exclude-exact"]) == 0
    old = pl.read_csv(retained.with_suffix(".selection.tsv"), separator="\t")
    new = pl.read_csv(excluded.with_suffix(".selection.tsv"), separator="\t")
    assert old["score_model"].to_list() == new["score_model"].to_list() == ["vx", "vx"]
    assert old["selected_epitope"].to_list() == new["selected_epitope"].to_list() == ["DEMO_A", "DEMO_A"]
    assert old["score"].to_list() == pytest.approx([0.01875, 0.01875])
    assert new["score"].to_list() == [-6.0, -6.0]


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
    out = annotate_repertoires(path, reference="unused", score_model="kernel", output_prefix=tmp_path / "out")
    assert len(core) == 1 and core[0][0].locus == "TRB"
    assert out["queries"]["clone_id"].to_list() == ["007"]


@pytest.mark.parametrize("threads", [True, -1, 1.5])
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
    out = annotate_repertoires(path, reference=ref, score_model="kernel", output_prefix=tmp_path / "out")
    assert out["queries"].height == out["selection"].height == 2
    assert out["queries"]["clone_id"].to_list() == ["007", "008"]
    assert out["queries"]["source_query_ida"].to_list() == ["q007", "q008"]
    assert out["queries"]["source_query_idb"].to_list() == ["q007", "q008"]
    assert {"TRA.selection", "TRB.selection"} <= out.keys()
    for lane in ["TRA", "TRB"]:
        assert set(out[lane + ".selection"]["source_query_id"]) == {"q007", "q008"}
    vx = annotate_repertoires(path, reference=ref, output_prefix=tmp_path / "vx")
    assert vx["queries"].select("source_pair_id", "status").rows() == [
        ("007", "matched"), ("008", "invalid_alpha_query")]
    assert vx["TRB.selection"].select("source_pair_id", "selected_epitope").rows() == [
        ("007", "E"), ("008", "E")]


@pytest.mark.parametrize("complete", [True, False])
def test_cell_linkage_precedes_shared_clone_identity(tmp_path, core, complete):
    loci = ["TRA", "TRB", "TRA", "TRB"] if complete else ["TRA", "TRB"]
    cells = ["cell-a", "cell-a", "cell-b", "cell-b"] if complete else ["cell-a", "cell-b"]
    path = tmp_path / "cells.tsv"
    pl.DataFrame({"cell_id": cells, "clone_id": ["clone-x"] * len(loci),
        "clonotype_id": ["type-y"] * len(loci), "sequence_id": [str(i) for i in range(len(loci))],
        "locus": loci, "junction_aa": ["CAVF" if x == "TRA" else "CASSF" for x in loci]}).write_csv(path, separator="\t")
    out = annotate_repertoires(path, reference="unused", score_model="kernel", output_prefix=tmp_path / "out")
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
    result = annotate_repertoires(path, reference='unused', score_model="kernel", output_prefix=tmp_path / 'out')
    assert result['queries']['status'].to_list() == ['matched']
    assert result['queries']['input_status'].to_list() == ['raw-status']


def test_initialized_table_pool_within_default_allocation(tmp_path):
    root = Path(__file__).resolve().parents[2]
    code = """
import sys
import polars as pl
assert pl.thread_pool_size() == 4
from unittest.mock import patch
from vdjmatch import annotate_repertoires
with patch('vdjtools.cores.available_cores', return_value=4), patch('vdjtools.cores.physical_cores', return_value=4):
    result = annotate_repertoires(sys.argv[1], reference=sys.argv[2], output_prefix=sys.argv[3])
assert result['queries'].height == 2
assert result['manifest']['native_threads'] == 4
assert result['manifest']['table_threads'] == 4
"""
    result = subprocess.run([sys.executable, '-c', code,
        str(root / 'examples/inputs/cohort/sample_a.tsv'),
        str(root / 'examples/inputs/reference.tsv'), str(tmp_path / 'initialized')],
        env=dict(os.environ, POLARS_MAX_THREADS='4'), capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr


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
