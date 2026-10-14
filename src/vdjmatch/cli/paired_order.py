"""Opt-in ordering of complete cohorts of linked-pair evidence."""
from pathlib import Path


def register(subparsers):
    parser = subparsers.add_parser(
        "paired-evidence-order", help="population-relative three-block paired evidence order",
        description="Rank a complete cohort per species/full pMHC; scores are not comparable across batches or calibrated probabilities.",
    )
    parser.add_argument("input", help="flat paired evidence TSV (one row per query/cohort/full pMHC)")
    parser.add_argument("--output", required=True, help="output TSV, preserving all input columns")
    parser.add_argument("--threads", type=int, default=1, help="Polars thread budget (default: 1)")
    parser.set_defaults(func=main)
    return parser


def main(args):
    import polars as pl
    from ..aggregate.paired_order import paired_evidence_order, GROUP, SCORES, COUNTS

    frame = pl.read_csv(args.input, separator="\t", schema_overrides={
        **{name: pl.String for name in ["query_id", *GROUP]},
        **{name: pl.Float64 for name in SCORES},
        **{name: pl.UInt32 for name in COUNTS}, "paired_available": pl.Boolean,
    })
    out = paired_evidence_order(frame)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.output, separator="\t")
    return 0
