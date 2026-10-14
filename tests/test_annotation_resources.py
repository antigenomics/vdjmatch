"""Explicit cohort release gate: native budgets, parity and summed process-tree RSS.

Run with test/bench extras: python -m pytest tests/test_annotation_resources.py -q.
Synthetic junctions only; this does not recompute a scientific benchmark.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import polars as pl
import pytest


def _measured(command, report, threads):
    psutil = pytest.importorskip("psutil", reason="resource gate needs vdjmatch[bench]")
    started = time.monotonic()
    peak = children = 0
    env = dict(os.environ, POLARS_MAX_THREADS=str(threads), OPENBLAS_NUM_THREADS="1",
               OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with report.with_suffix(".log").open("w") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=env)
        try:
            while process.poll() is None:
                try:
                    root = psutil.Process(process.pid)
                    descendants = root.children(recursive=True)
                    children = max(children, len(descendants))
                    rss = 0
                    for member in [root, *descendants]:
                        try:
                            rss += member.memory_info().rss
                        except psutil.NoSuchProcess:
                            pass
                    peak = max(peak, rss)
                except psutil.NoSuchProcess:
                    pass
                assert time.monotonic() - started < 120, "cohort resource check exceeded120s"
                assert peak < 2_000_000_000, "synthetic compact cohort exceeded2GB summedRSS"
                time.sleep(.02)
        finally:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=5)
    result = dict(wall_seconds=time.monotonic()-started,
                  peak_process_tree_rss_bytes=peak, max_child_processes=children,
                  requested_native_threads=threads, exit_code=process.returncode)
    report.write_text(json.dumps(result, indent=2)+"\n")
    assert process.returncode == 0, report.with_suffix(".log").read_text()
    assert children == 0, "annotation must use batched native threads, not query worker processes"


@pytest.mark.slow
@pytest.mark.parametrize("input_format", ["simple", "airr_single", "wide_paired", "long_paired", "mixed"])
def test_directory_library_metadata_cli_parallel_resource_parity(tmp_path, input_format):
    samples = tmp_path / "samples"
    samples.mkdir()
    rows_per_sample = 1024
    paired = input_format.endswith("paired") or input_format == 'mixed'
    for name in ("one", "two", "three", "four"):
        identity = {"query_id": [f"q{i:04}" for i in range(rows_per_sample)],
            "clonotype_id": [f"type{i:04}" for i in range(rows_per_sample)],
            "clone_id": [f"clone{i:04}" for i in range(rows_per_sample)]}
        beta = ["CASSLGQAYEQYF", "CASSLGQAYEQWF", "invalid", None] * (rows_per_sample//4)
        if paired and (input_format != 'mixed' or name in ('one','two')):
            alpha = ["CAVRDSNYQLIW", None, "invalid", None] * (rows_per_sample//4)
            beta = ["CASSLGQAYEQYF", "CASSLGQAYEQWF", None, None] * (rows_per_sample//4)
            if input_format in ("wide_paired", "mixed"):
                frame = pl.DataFrame({**identity, "junction_aa_alpha": alpha, "junction_aa_beta": beta,
                    "v_alpha": ["TRAV12-1"] * rows_per_sample, "j_alpha": ["TRAJ33"] * rows_per_sample,
                    "v_beta": ["TRBV19"] * rows_per_sample, "j_beta": ["TRBJ2-7"] * rows_per_sample})
            else:
                lanes = []
                for lane, locus, sequences, v, j in [("a", "TRA", alpha, "TRAV12-1", "TRAJ33"),
                                                    ("b", "TRB", beta, "TRBV19", "TRBJ2-7")]:
                    lanes.append(pl.DataFrame({**identity, "query_id": [q+lane for q in identity["query_id"]],
                        "sequence_id": [q+lane for q in identity["query_id"]],
                        "cell_id": ["cell"+str(i) for i in range(rows_per_sample)],
                        "locus": [locus]*rows_per_sample, "junction_aa": sequences,
                        "v_call": [v]*rows_per_sample, "j_call": [j]*rows_per_sample}))
                frame = pl.concat(lanes)
        else:
            frame = pl.DataFrame({**identity, "locus": ["TRB"]*rows_per_sample,
                "v": ["TRBV19"]*rows_per_sample, "j": ["TRBJ2-7"]*rows_per_sample,
                "junction_aa": beta})
            if input_format == "airr_single":
                frame = frame.rename({"v": "v_call", "j": "j_call"}).with_columns(
                    pl.col("query_id").alias("sequence_id"))
        frame.write_csv(samples / (name+".tsv"), separator="\t")
    metadata = tmp_path / "samples.tsv"
    names = sorted(p.stem for p in samples.iterdir())
    pl.DataFrame({"sample_id": names, "path": [f"samples/{name}.tsv" for name in names]}).write_csv(metadata, separator="\t")
    reference = tmp_path / "reference.tsv"
    ref_rows = []
    for complex_id, epitope in [(1, "E"), (2, "other"), (0, "E")]:
        for locus, sequence, v in [("TRA", "CAVRDSNYQLIW", "TRAV12-1"),
                                   ("TRB", "CASSLGQAYEQYF", "TRBV19")]:
            ref_rows.append(dict(cdr3=sequence, v=v, gene=locus, complex_id=complex_id,
                species="HomoSapiens", epitope=epitope, mhc_a="HLA-A*02:01", mhc_b="B2M", mhc_class="MHCI"))
    pl.DataFrame(ref_rows).write_csv(reference, separator="\t")
    driver = tmp_path / "annotation_probe.py"
    driver.write_text('''import json,sys
from functools import wraps
from seqtree import gapblock
calls={'grouped_count_batch':[], 'grouped_paired_sum_count_batch':[]}
for name in calls:
    original=getattr(gapblock,name)
    def wrap(original,name):
        @wraps(original)
        def counted(*args,**kwargs):
            calls[name].append(kwargs.get('threads'))
            return original(*args,**kwargs)
        return counted
    setattr(gapblock,name,wrap(original,name))
mode,source,reference,prefix,threads,paired=sys.argv[1:]
threads=int(threads)
if mode=='library':
    from vdjmatch import annotate_repertoires
    annotate_repertoires(source,reference=reference,output_prefix=prefix,threads=threads)
else:
    from vdjmatch.cli.__main__ import main
    assert main(['annotate','--metadata',source,'--vdjdb',reference,'--threads',str(threads),'-o',prefix])==0
assert calls['grouped_count_batch'] and set(calls['grouped_count_batch'])=={threads}, calls
if paired in ('1','2'):
    assert calls['grouped_paired_sum_count_batch'] and set(calls['grouped_paired_sum_count_batch'])=={threads}, calls
if paired=='2':
    assert len(calls['grouped_count_batch'])==2 and len(calls['grouped_paired_sum_count_batch'])==1, calls
print(json.dumps({'native_batches':{k:len(v) for k,v in calls.items()},
    'native_thread_arguments':{k:sorted(set(v)) for k,v in calls.items()}}))
''')
    baseline = None
    for mode in ("library", "cli"):
        for threads in (1, 4):
            prefix = tmp_path / f"{mode}{threads}"
            command = [sys.executable, str(driver), mode, str(samples if mode == "library" else metadata),
                       str(reference), str(prefix), str(threads), '2' if input_format=='mixed' else str(int(paired))]
            _measured(command, prefix.with_suffix(".resources.json"), threads)
            manifest=json.loads(Path(str(prefix)+'.manifest.json').read_text())
            assert manifest['population_service_calls']==(3 if input_format=='mixed' else 1)
            assert manifest['input_samples']==len(names)
            assert manifest['input_chain_rows']==rows_per_sample*(6 if input_format=='mixed' else len(names)*(2 if paired else 1))
            assert manifest['native_threads']==threads
            assert manifest['table_threads']<=threads
            tables = {}
            for kind in ("queries", "selection"):
                frame = pl.read_csv(str(prefix)+f".{kind}.tsv", separator="\t", infer_schema_length=0)
                assert frame.height == rows_per_sample*len(names)
                assert frame["query_id"].n_unique() == frame.height
                assert set(frame["sample_id"]) == set(names)
                assert frame.group_by("sample_id").len()["len"].to_list() == [rows_per_sample]*len(names)
                assert not frame['source_clone_id'].null_count()
                assert set(frame['source_clone_id']) == {f"clone{i:04}" for i in range(rows_per_sample)}
                source_columns = (['source_query_id','source_query_ida','source_query_idb'] if input_format=='mixed' else
                                  ['source_query_ida', 'source_query_idb'] if paired else ['source_query_id'])
                for column in source_columns:
                    if input_format!='mixed':assert not frame[column].null_count()
                    suffix = column[-1] if input_format == 'long_paired' else ''
                    assert set(frame[column].drop_nulls()) == {f'q{i:04}'+suffix for i in range(rows_per_sample)}
                clonotype_columns = (['clonotype_id','clonotype_ida','clonotype_idb'] if input_format=='mixed' else
                                     ['clonotype_ida', 'clonotype_idb'] if paired else ['clonotype_id'])
                for column in clonotype_columns:
                    assert set(frame[column].drop_nulls()) == {f'type{i:04}' for i in range(rows_per_sample)}
                tables[kind] = frame.sort("sample_id", "query_id")
            if baseline is None:
                baseline = tables
            else:
                for kind in tables:
                    assert tables[kind].equals(baseline[kind]), (mode, threads, kind)
