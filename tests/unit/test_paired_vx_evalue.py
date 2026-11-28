import json

import polars as pl
import pytest

from vdjmatch.cli.paired_vx_evalue import _source, main
from vdjmatch.db.cache import sha256


def score_directory(directory, scores, *, null_source=None, excluded=False):
    directory.mkdir()
    frame = pl.DataFrame({
        "query_id": [f"{i:02d}" for i in range(len(scores))],
        "target": ["pMHC" if score is not None else None for score in scores],
        "maximum": scores, "score": scores, "tied": [False] * len(scores),
    })
    if null_source:
        frame = frame.with_columns(outcome=pl.lit("null"), null_source=pl.lit(null_source))
    else:
        # A finite diagnostic maximum on an abstaining query contributes no call.
        frame = frame.with_columns(tied=pl.col("score").is_null(),
                                   maximum=pl.col("score").fill_null(9.0))
    frame.write_csv(directory / "maxima.tsv", separator="\t")
    (directory / "targets.tsv").write_text("target\tmhc_class\npMHC\tMHCI\n")
    manifest = dict(
        inputs_sha256=dict(reference="reference", queries=directory.name),
        outputs_sha256={name: sha256(directory / name) for name in ("maxima.tsv", "targets.tsv")},
        species="mouse", exact_hits_excluded=excluded, reference_keys=1,
        v_loop_model_sha256="model", source_sha256="scorer", region_weights={"cdr1": 1},
        competitive_rule={"margin": 6}, software={"vdjmatch": "test"}, queries=len(scores),
        null_source=null_source,
    )
    (directory / "manifest.json").write_text(json.dumps(manifest))


def test_bound_evalue_round_trip_and_rejects_exact_policy_mismatch(tmp_path):
    from types import SimpleNamespace

    nulls, queries, excluded = (tmp_path / name for name in ("nulls", "queries", "excluded"))
    score_directory(nulls, [10, 9, 8, 8, 7, None, None, None, None, None], null_source="generated")
    score_directory(queries, [9, 8, None])
    score_directory(excluded, [9, 8, None], excluded=True)
    args = SimpleNamespace(queries=queries, null_dir=nulls, output_prefix=tmp_path / "out",
                           target_fpr=.5, confidence=.75)
    assert main(args) == 0
    result = pl.read_csv(tmp_path / "out.calls.tsv", separator="\t", schema_overrides={"query_id": pl.String})
    assert result["query_id"].to_list() == ["00", "01", "02"]
    assert result["E_null"].to_list() == [.2, .4, 0.]
    assert result["null_tail_count"].to_list() == [2, 4, 0]
    assert result["call_score"].to_list() == [9., None, None]
    manifest = json.loads((tmp_path / "out.manifest.json").read_text())
    assert manifest["calibration"]["empirical_E_cutoff"] == .3
    assert manifest["accepted"] == 1
    assert not manifest["calibration_bank_replayed"]
    args.queries = excluded
    with pytest.raises(ValueError, match="exact-hit policy"):
        main(args)


def test_bound_evalue_rejects_tampered_export_and_fold_input(tmp_path):
    source = tmp_path / "sample"
    score_directory(source, [9, 8, None])
    manifest = source / "manifest.json"
    original = json.loads(manifest.read_text())
    manifest.write_text(json.dumps(original | {"inputs_sha256": {"reference": "reference"}}))
    with pytest.raises(ValueError, match="incomplete paired-vx-base manifest"):
        _source(source)
    manifest.write_text(json.dumps(original))
    with (source / "maxima.tsv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(ValueError, match="maxima.tsv differs"):
        _source(source)
    manifest.write_text(json.dumps(json.loads(manifest.read_text()) | {"withholding_folds": 1}))
    with pytest.raises(ValueError, match="ordinary mouse"):
        _source(source)
