"""Rebuild assay examples from raw AIRR tables through the VDJmatch CLI.

The LLW-versus-LLL sorted-population contrast is kept separate from the mirror's
Pgen-matched repertoire-control tasks. No classifier or threshold is fitted.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys


def prepare(benchmark):
    import polars as pl

    parts = []
    for number in (1, 2):
        raw = pl.read_csv(benchmark / f"vdjmatch/airr/dataset_{number}_TRB.tsv.gz", separator="\t")
        parts.append(raw.with_columns(pl.lit(number).alias("dataset"),
            pl.col("epitope").alias("source_epitope"),
            pl.col("binder").alias("positive"),
            (pl.col("epitope") + (" sorted positive / negative" if number == 1 else
                " / matched repertoire controls")).alias("task")))
    sorted_cells = parts[1].filter((pl.col("binder") == 1) &
        pl.col("source_epitope").is_in(["LLWNGPMAV", "LLLGIGILV"]))
    contrast = sorted_cells.with_columns(
        pl.lit("LLW / LLL sorted populations").alias("task"),
        (pl.col("source_epitope") == "LLWNGPMAV").cast(pl.Int64).alias("positive"),
        pl.lit("LLWNGPMAV").alias("epitope"))
    return pl.concat([*parts, contrast]).with_row_index("query_id")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, required=True, help="raw airr_benchmark mirror")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="new output directory")
    parser.add_argument("--threads", type=int, default=4, help="CPU worker budget")
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be positive")
    # CLI owns its own worker budget; post-analysis uses one numerical helper pool.
    for key in ("POLARS_MAX_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[key] = "1"
    os.environ["MPLBACKEND"] = "Agg"
    import numpy as np
    import polars as pl
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve

    raw = prepare(args.benchmark)
    args.output.mkdir(parents=True, exist_ok=False)
    raw.write_csv(args.output / "queries.tsv", separator="\t")
    command = [sys.executable, "-m", "vdjmatch.cli", "evidence", str(args.output / "queries.tsv"),
        "--vdjdb", str(args.reference), "--locus", "TRB", "--targets-from-sample",
        "--reference-population", "curated-benchmark", "--exclude-exact",
        "--threads", str(args.threads), "--output-prefix", str(args.output / "evidence")]
    with (args.output / "cli.log").open("w") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, timeout=300)
    keys = ["query_id", "epitope", "mhc_a", "mhc_b", "mhc_class"]
    evidence = pl.read_csv(args.output / "evidence.enrichment.tsv", separator="\t")
    scored = raw.join(evidence, on=keys, how="left", validate="1:1", nulls_equal=True)
    assert scored.height == raw.height == evidence.height
    scored.write_csv(args.output / "observations.tsv", separator="\t")
    summaries, curves, counts = [], {}, []
    for (task,), part in scored.group_by("task", maintain_order=True):
        y = part["positive"].to_numpy()
        score = part["kernel_score"].to_numpy()
        finite = np.isfinite(score)
        score = np.where(finite, score, np.nextafter(score[finite].min(), -np.inf) if finite.any() else -1.)
        fpr, tpr, thresholds = roc_curve(y, score, drop_intermediate=False)
        precision, recall, _ = precision_recall_curve(y, score)
        positives, controls = int(y.sum()), int((y == 0).sum())
        summaries.append(dict(task=task, observations=len(y), positives=positives, controls=controls,
            missing_scores=int((~finite).sum()), ROC_AUC=roc_auc_score(y, score),
            average_precision=average_precision_score(y, score)))
        curves[task] = (fpr, tpr, recall, precision)
        for threshold, tp, fp in zip(thresholds, np.rint(tpr * positives), np.rint(fpr * controls)):
            counts.append(dict(task=task, threshold=float(threshold), TP=int(tp), FP=int(fp),
                FN=positives-int(tp), TN=controls-int(fp)))
    summary = pl.DataFrame(summaries)
    summary.write_csv(args.output / "accuracy.tsv", separator="\t")
    pl.DataFrame(counts).write_csv(args.output / "operating-counts.tsv", separator="\t")
    plt.rcParams.update({"font.size":8, "axes.spines.top":False, "axes.spines.right":False,
        "pdf.fonttype":42, "svg.fonttype":"none"})
    tasks = ["NLVPMVATV sorted positive / negative", "LLW / LLL sorted populations"]
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.15), layout="constrained")
    for task, label, colour, style in zip(tasks, ["NLV sorted + / −", "LLW / LLL sorted"],
                                         ["#0072B2", "#D55E00"], ["-", "--"]):
        fpr, tpr, recall, precision = curves[task]
        row = summary.filter(pl.col("task") == task).row(0, named=True)
        axes[0].plot(fpr, tpr, color=colour, ls=style, label=f"{label}; AUC {row['ROC_AUC']:.3f}")
        axes[1].step(recall, precision, where="post", color=colour, ls=style,
            label=f"{label}; AP {row['average_precision']:.3f}")
        axes[1].axhline(row['positives']/row['observations'], color=colour, lw=.7, ls=":")
    axes[0].plot([0,1], [0,1], color="0.6", ls=":", lw=.8)
    axes[0].set(xlabel="False-positive rate", ylabel="True-positive rate", title="a  Sorted-population discrimination")
    axes[1].set(xlabel="Recall", ylabel="Precision", title="b  Recovery of target-specific receptors")
    for ax in axes:
        ax.set(xlim=(0,1), ylim=(0,1.02))
        ax.legend(frameon=False, fontsize=7)
    for extension in ("png", "svg", "pdf"):
        fig.savefig(args.output / f"sorted-populations.{extension}", dpi=300, bbox_inches="tight")
    print(summary.write_csv(separator="\t"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
