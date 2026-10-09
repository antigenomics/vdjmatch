"""vdjmatch command-line interface.

Subcommands:
  update     fetch/cache the latest VDJdb release
  match      annotate query sample(s) against VDJdb (E-values + ranked hits + epitope summary)
  precursor  T-cell precursor frequency and unseen-junction diversity for a set of TCRs

``match`` writes hit, candidate, call, descriptive summary and ingestion TSV tables,
plus a JSON run manifest. Calls preserve query rows and identify ambiguous/no-hit cases.
Paired mode reports same-complex candidate evidence and descriptive pair counts.

``precursor`` writes one TSV, one row per group, and needs the optional extra:
  pip install 'vdjmatch[precursor]'
"""

from __future__ import annotations

import argparse
import resource
import sys
import time
from pathlib import Path

from ..config import Params


def __getattr__(name):
    """Keep numerical imports lazy until CLI resource budgets are selected."""
    import importlib

    if name in globals():
        return globals()[name]
    if name in {"db", "match"}:
        value = importlib.import_module(f"vdjmatch.{name}")
    elif name == "annotate_sample":
        value = importlib.import_module("vdjmatch.runner.multisample").annotate_sample
    else:
        raise AttributeError(name)
    globals()[name] = value
    return value


def _flat_table(frame):
    """Encode nested reference metadata as JSON cells for a flat TSV."""
    import polars as pl

    return frame.with_columns(
        *[
            pl.struct(c).struct.json_encode().alias(c)
            for c, t in frame.schema.items()
            if t.is_nested()
        ]
    )


def _cmd_update(a: argparse.Namespace) -> int:
    db = __getattr__("db")
    path = db.fetch_latest(asset=a.asset, pin=a.pin, force=a.force, cache=a.cache)
    print(f"VDJdb {a.asset} cached at {path}")
    return 0


def _peak_rss_gb() -> float:
    """Peak resident set in GB (macOS reports ru_maxrss in bytes, Linux in KiB)."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_bytes = rss if sys.platform == "darwin" else rss * 1024
    return rss_bytes / 1e9


# match flags whose JSON-config value is overridable from the CLI; sentinel default -> "not passed".
_PARAM_FLAGS = (
    "scope",
    "matrix",
    "min_score",
    "match_v",
    "match_j",
    "evalue",
    "align",
    "species",
    "threads",
)


def _resolve_params(a: argparse.Namespace) -> Params:
    """Start from defaults or ``--config`` JSON, then apply any explicitly-passed CLI flag."""
    params = Params.from_json(a.config) if a.config else Params.defaults()
    for field in _PARAM_FLAGS:
        v = getattr(a, field, None)
        if v is not None:  # argparse.SUPPRESS -> absent unless passed
            setattr(params, field, v)
    params.validate()
    return params


def _cmd_match(a: argparse.Namespace) -> int:
    import dataclasses
    import importlib.metadata
    import json
    import polars as pl
    from ..db.cache import sha256

    db, match = __getattr__("db"), __getattr__("match")
    annotate_sample = __getattr__("annotate_sample")
    p = _resolve_params(a)
    names = []
    for sample in a.samples:
        name = Path(sample).name
        if name.endswith(".gz"):
            name = name[:-3]
        name = Path(name).stem
        names.append(name)
    if len(set(names)) != len(names):
        raise ValueError("sample output names collide; use distinct input basenames")
    reference = a.vdjdb
    if reference is None:
        reference = db.fetch_latest(
            asset=a.asset, pin=a.pin, cache=getattr(a, "cache", None)
        )
    reference_species = None if p.species == "any" else p.species
    control_species = getattr(a, "control_species", None) or reference_species
    if p.evalue and control_species is None:
        raise ValueError(
            "--species any with calibration requires --control-species human or mouse"
        )
    reference_filters = {
        k: getattr(a, k, None)
        for k in (
            "epitope",
            "mhc_a",
            "mhc_b",
            "reference_id",
            "exclude_reference_ids",
            "evidence_type",
        )
    }
    vdj = db.load(
        reference,
        asset=a.asset,
        species=reference_species,
        min_score=p.min_score,
        **reference_filters,
    )
    index = match.VdjdbIndex.build(vdj)
    print(
        f"VDJdb: {vdj.height:,} observations; indexed loci {index.genes}",
        file=sys.stderr,
    )
    reference_provenance = db.provenance(reference)
    matrix = match.load_vdjam() if p.matrix == "vdjam" else ""
    Path(a.output_prefix).parent.mkdir(parents=True, exist_ok=True)
    p.to_json(f"{a.output_prefix}.params.json")
    for sample, name in zip(a.samples, names):
        t0 = time.perf_counter()
        res = annotate_sample(
            index,
            sample,
            scope=p.scope,
            matrix=matrix or None,
            species=control_species,
            with_evalue=p.evalue,
            match_v=p.match_v,
            match_j=p.match_j,
            align=p.align,
            threads=p.threads,
            progress=a.verbose,
            source=getattr(a, "input_format", "auto"),
            sequence_convention=getattr(a, "sequence_convention", None),
            paired=getattr(a, "paired", False),
            link=getattr(a, "link", None),
        )
        for kind, frame in res.items():
            _flat_table(frame).write_csv(
                f"{a.output_prefix}.{name}.{kind}.txt", separator="\t"
            )
        manifest = {
            "reference": reference_provenance,
            "sample": {"sha256": sha256(Path(sample))},
            "software": {
                pkg: importlib.metadata.version(pkg)
                for pkg in ["vdjmatch", "seqtree", "vdjtools", "polars"]
            },
            "parameters": dataclasses.asdict(p),
            "reference_filters": reference_filters,
            "control_species": control_species,
            "input_format": getattr(a, "input_format", "auto"),
            "sequence_convention": getattr(a, "sequence_convention", None)
            or "producer-defined junction",
            "paired": getattr(a, "paired", False),
            "link": getattr(a, "link", None),
            "resources": {
                "native_threads": p.threads,
                "polars_pool_threads": pl.thread_pool_size(),
                "parallel_layers": "serial samples, sequential native and table stages",
            },
            "wall_seconds": time.perf_counter() - t0,
            "process_peak_rss_gb": _peak_rss_gb(),
            "ingestion": res["ingestion"].to_dicts()[0]
            if "ingestion" in res and res["ingestion"].height
            else {},
            "reference_counts": {
                "rows": vdj.height,
                "searchable_rows": int(vdj["reference_valid"].sum())
                if "reference_valid" in vdj.columns
                else vdj.height,
            },
            "outputs": {kind: frame.height for kind, frame in res.items()},
        }
        Path(f"{a.output_prefix}.{name}.manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n"
        )
        print(
            f"[{name}] {res['hits'].height} hit rows -> {a.output_prefix}.{name}.*",
            file=sys.stderr,
        )
    return 0


#: VDJdb species -> the recombination-model organism that scores it. Scoring one species' junctions
#: against another's model does not error and does not produce zeros -- mouse VDJdb against the
#: human model returns a finite `F` for every group, a median 0.157x of the right answer with an
#: 11x spread across epitopes, so it corrupts the ranking and not just the scale. Silent and
#: plausible is the worst failure mode available, hence the guard.
_MODEL_ORGANISM = {
    "HomoSapiens": "human",
    "MusMusculus": "mouse",
    "MacacaMulatta": "monkey",
}


def _cmd_precursor(a: argparse.Namespace) -> int:
    import polars as pl

    db = __getattr__("db")

    from ..precursor import summarise

    # `olga` ships human only; mouse lives under `arda`. Switch rather than fail, since there is no
    # mouse model in the default set for the user to have meant.
    if a.organism != "human" and a.source == "olga":
        a.source = "arda"
        print(
            f"note: --organism {a.organism} has no model in 'olga'; using --source arda",
            file=sys.stderr,
        )

    if a.vdjdb or not a.samples:
        want = _MODEL_ORGANISM.get(a.species)
        if want and want != a.organism:
            raise SystemExit(
                f"--species {a.species} needs --organism {want}, but --organism is "
                f"{a.organism!r}. Scoring {a.species} junctions against the {a.organism} "
                f"recombination model runs without error and returns a plausible number that is "
                f"wrong by a factor which varies between epitopes, so it is refused rather than "
                f"warned about.\n"
                f"  run: vdjmatch precursor --vdjdb --species {a.species} --organism {want}"
            )
        frame = db.load(a.table, asset="slim", species=a.species, pin=a.pin)
        if a.mhc_class:
            frame = frame.filter(pl.col("mhc_class") == a.mhc_class)
        if not frame.height:
            raise SystemExit(
                f"VDJdb returned 0 records for species={a.species!r} "
                f"mhc_class={a.mhc_class!r} -- note species is 'HomoSapiens', "
                "not 'human'"
            )
        junction_col, group_col = "cdr3", a.group_by or "epitope"
        chain_col = a.chain_col or "gene"
        capture_col = None if a.no_unseen else (a.capture_col or "reference_id")
        label = f"VDJdb ({frame.height:,} records)"
    else:
        sep = "," if a.samples[0].endswith(".csv") else "\t"
        frame = pl.concat(
            [pl.read_csv(s, separator=sep, infer_schema_length=0) for s in a.samples],
            how="vertical_relaxed",
        )
        junction_col, group_col = a.junction_col, a.group_by
        chain_col, capture_col = a.chain_col, (None if a.no_unseen else a.capture_col)
        label = f"{len(a.samples)} file(s) ({frame.height:,} rows)"
    print(f"scoring {label}", file=sys.stderr)

    t0 = time.perf_counter()
    out = summarise(
        frame,
        junction_col=junction_col,
        group_col=group_col,
        chain_col=chain_col,
        capture_col=capture_col,
        locus=a.locus,
        source=a.source,
        organism=a.organism,
        r=a.radius,
        alpha=a.alpha,
        q=a.q,
        n_cells=a.n_cells,
        compartment=a.compartment,
        n_eff=a.n_eff,
        selection=("auto" if a.selection == "auto" else float(a.selection)),
        min_junctions=a.min_junctions,
        threads=a.threads,
        progress=a.verbose,
    )
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(a.output, separator="\t")
    print(
        f"{out.height} group(s) -> {a.output} in {time.perf_counter() - t0:.1f}s "
        f"(peak RSS {_peak_rss_gb():.2f} GB)",
        file=sys.stderr,
    )
    if a.q == 1.0:
        print(
            "note: q=1 -> F is an uncalibrated model mass; only its ranking is meaningful. "
            "Pass --q to apply a selection constant.",
            file=sys.stderr,
        )

    # A junction that passes the anchor check and still scores exactly 0 leaves the mass without a
    # signal, so say so. This is a property of the model set: `olga` faithfully reproduces OLGA's
    # deletion-bin grid, on which an allele shorter than the grid carries probability on trims it
    # cannot reach, and 5.7% of human TRA junctions in VDJdb score exactly zero as a result. The
    # refit sets do not inherit it (`learned` 0.5%, `arda` 0.0%).
    if "n_zero_pgen" in out.columns and "n_junctions" in out.columns:
        z, n = int(out["n_zero_pgen"].sum()), int(out["n_junctions"].sum())
        if z and n and z / n >= 0.01:
            print(
                f"note: {z:,}/{n:,} junctions ({100 * z / n:.1f}%) score Pgen exactly 0 under "
                f"--source {a.source} and contribute nothing to any mass. This is OLGA's "
                f"deletion-grid quirk, not your input; --source arda (or learned) does not "
                f"inherit it. Per-group counts are in the n_zero_pgen column.",
                file=sys.stderr,
            )
    return 0


_EXAMPLES = """\
examples:
  vdjmatch update                                  # cache the latest VDJdb (slim)
  vdjmatch match sample.tsv                         # annotate one AIRR sample, default reference
  vdjmatch match -o run/out --match-v *.tsv         # match V gene too, write under run/
  vdjmatch match --scope 2,1,1,2 --threads 8 s.tsv  # wider search budget, 8 threads
  vdjmatch precursor --vdjdb -o pre.txt             # precursor frequency per VDJdb epitope
  vdjmatch precursor tcrs.tsv --group-by epitope    # ... or for your own grouped TCR table
"""


_PRECURSOR_EPILOG = """\
examples:
  vdjmatch precursor --vdjdb -o pre.txt                    # every VDJdb epitope, both chains
  vdjmatch precursor --vdjdb --min-junctions 10 -r 2       # well-sampled epitopes, radius 2
  vdjmatch precursor tcrs.tsv --group-by epitope --locus TRA
  vdjmatch precursor --vdjdb --q 9.41 --n-eff 1e8          # calibrated F, plus P(>=k precursors)
  vdjmatch precursor --vdjdb --n-eff 1e8 --selection auto   # seen/unseen clonotypes, measured Q

notes:
  Sequences must be JUNCTIONS (Cys104..Phe/Trp118 inclusive), not IMGT CDR3s -- VDJdb's column is
  named `cdr3` but holds junctions. Anchor-stripped input is dropped and counted in `n_dropped`.
  `--q` is the selection constant carrying a generation probability to a repertoire frequency;
  left at 1 the reported F is a raw model mass whose ranking, not scale, is meaningful.
"""


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="vdjmatch",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_EXAMPLES,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    up = sub.add_parser("update", help="fetch/cache the latest VDJdb release")
    up.add_argument(
        "--asset",
        default="default",
        choices=["slim", "full", "default", "primary", "legacy", "airr"],
        help="bundle/table role (default: primary, historical long fallback)",
    )
    up.add_argument("--cache", default=None, help="reference download directory")
    up.add_argument(
        "--pin", default=None, help="pin a specific release tag (default: latest)"
    )
    up.add_argument(
        "--force", action="store_true", help="re-download even if already cached"
    )
    up.set_defaults(func=_cmd_update)

    m = sub.add_parser(
        "match",
        help="annotate sample(s) against VDJdb",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_EXAMPLES,
    )
    m.add_argument("samples", nargs="+", help="AIRR rearrangement sample file(s) (TSV)")
    m.add_argument(
        "-o",
        "--output-prefix",
        default="vdjmatch_out",
        help="output path prefix; files are <prefix>.<sample>.{hits,calls,summary}.txt",
    )
    m.add_argument(
        "--vdjdb", default=None, help="custom VDJdb table path (default: fetch latest)"
    )
    m.add_argument(
        "--asset",
        default="default",
        choices=["slim", "full", "default", "primary", "legacy", "airr"],
        help="bundle/table role (default: primary, historical long fallback)",
    )
    m.add_argument("--cache", default=None, help="reference download directory")
    m.add_argument(
        "--input-format",
        default="auto",
        choices=["auto", "airr", "legacy", "custom"],
        help="query sequence producer/convention; repertoire formats are detected by vdjtools",
    )
    m.add_argument(
        "--sequence-convention",
        default=None,
        choices=["junction"],
        help="explicitly declare anchor-inclusive sequence for ambiguous custom tables",
    )
    m.add_argument(
        "--paired", action="store_true", help="read linked TRA/TRB cell rows"
    )
    m.add_argument("--link", default=None, help="cell linkage column for --paired")
    m.add_argument(
        "--pin", default=None, help="pin a specific VDJdb release tag (default: latest)"
    )
    m.add_argument(
        "--config",
        default=None,
        help="JSON params file (see Params); explicit CLI flags below override its values",
    )
    m.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="show native search stages and print per-sample time / queries-per-s / peak RSS",
    )
    for selector in (
        "epitope",
        "mhc_a",
        "mhc_b",
        "reference_id",
        "exclude_reference_ids",
        "evidence_type",
    ):
        m.add_argument(
            "--" + selector.replace("_", "-"),
            action="append",
            help="reference selector (repeatable)",
        )
    m.add_argument(
        "--control-species",
        choices=["human", "mouse"],
        help="background organism; required with --species any and calibration",
    )
    # --- params (overridable from --config; SUPPRESS default => only override JSON when passed) ---
    m.add_argument(
        "--species",
        default=argparse.SUPPRESS,
        help="species filter (default: HomoSapiens); any omits reference filtering",
    )
    m.add_argument(
        "--scope",
        default=argparse.SUPPRESS,
        help="search budget as substitutions,insertions,deletions,total-edits (default: 1,0,0,1)",
    )
    m.add_argument(
        "--matrix",
        default=argparse.SUPPRESS,
        choices=["vdjam", "none"],
        help="substitution matrix: vdjam (TCR-specific, bundled) or none (unit edit cost)",
    )
    m.add_argument(
        "--min-score",
        dest="min_score",
        type=int,
        default=argparse.SUPPRESS,
        help="minimum VDJdb confidence score (default: 0)",
    )
    m.add_argument(
        "--match-v",
        dest="match_v",
        action="store_const",
        const=True,
        default=argparse.SUPPRESS,
        help="require the V gene to match as well as the CDR3",
    )
    m.add_argument(
        "--match-j",
        dest="match_j",
        action="store_const",
        const=True,
        default=argparse.SUPPRESS,
        help="require the J gene to match as well as the CDR3",
    )
    m.add_argument(
        "--no-evalue",
        dest="evalue",
        action="store_const",
        const=False,
        default=argparse.SUPPRESS,
        help="skip the control-calibrated E-value",
    )
    m.add_argument(
        "--no-align",
        dest="align",
        action="store_const",
        const=False,
        default=argparse.SUPPRESS,
        help="skip the per-hit CIGAR/alignment output",
    )
    m.add_argument(
        "--threads",
        type=int,
        default=argparse.SUPPRESS,
        help="native search threads (default: 1; 0 explicitly selects automatic)",
    )
    m.set_defaults(func=_cmd_match)

    pc = sub.add_parser(
        "precursor",
        help="precursor frequency + unseen diversity for a TCR set",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_PRECURSOR_EPILOG,
    )
    pc.add_argument(
        "samples", nargs="*", help="TCR table(s) (TSV/CSV); omit with --vdjdb"
    )
    pc.add_argument(
        "--vdjdb", action="store_true", help="score VDJdb itself instead of a file"
    )
    pc.add_argument(
        "--table", default=None, help="custom VDJdb table path (default: fetch latest)"
    )
    pc.add_argument("--pin", default=None, help="pin a specific VDJdb release tag")
    pc.add_argument(
        "--species", default="HomoSapiens", help="VDJdb species (default: HomoSapiens)"
    )
    pc.add_argument(
        "--mhc-class",
        dest="mhc_class",
        default=None,
        choices=["MHCI", "MHCII"],
        help="restrict VDJdb to one MHC class",
    )
    pc.add_argument(
        "-o", "--output", default="vdjmatch_precursor.txt", help="output TSV path"
    )
    # --- input columns ---
    pc.add_argument(
        "--junction-col",
        dest="junction_col",
        default="cdr3",
        help="junction column (Cys104..Phe/Trp118 inclusive; default: cdr3)",
    )
    pc.add_argument(
        "--group-by",
        dest="group_by",
        default=None,
        help="column to group by, one row of output per group (e.g. epitope)",
    )
    pc.add_argument(
        "--chain-col",
        dest="chain_col",
        default=None,
        help="locus column; each chain is scored with its own model",
    )
    pc.add_argument(
        "--capture-col",
        dest="capture_col",
        default=None,
        help="capture-unit column for the unseen-species fields (e.g. reference_id)",
    )
    pc.add_argument(
        "--no-unseen",
        dest="no_unseen",
        action="store_true",
        help="skip the Horvitz-Thompson unseen-species estimate",
    )
    # --- model + estimator ---
    pc.add_argument(
        "--locus", default="TRB", help="model locus when --chain-col is absent"
    )
    pc.add_argument(
        "--source", default="olga", help="recombination model set (default: olga)"
    )
    pc.add_argument(
        "--organism", default="human", help="model organism (mouse: --source arda)"
    )
    pc.add_argument(
        "-r",
        "--radius",
        type=int,
        default=1,
        help="neighbourhood radius in substitutions (default: 1, part of the estimator)",
    )
    pc.add_argument(
        "--alpha",
        type=float,
        default=0.1,
        help="cognacy retention per edit (default: 0.1, Mayer & Callan 2023)",
    )
    pc.add_argument(
        "--q",
        type=float,
        default=1.0,
        help="selection constant; 1 = uncalibrated raw mass (ALICE's TRB value: 9.41)",
    )
    pc.add_argument(
        "--min-junctions",
        dest="min_junctions",
        type=int,
        default=1,
        help="skip groups with fewer distinct junctions",
    )
    # --- absolute numbers ---
    pc.add_argument(
        "--n-cells",
        dest="n_cells",
        type=float,
        default=1e11,
        help="total T cells, for the expected precursor count (default: 1e11)",
    )
    pc.add_argument(
        "--compartment",
        type=float,
        default=1.0,
        help="fraction of the pool the restriction addresses (e.g. 0.3 for CD8)",
    )
    pc.add_argument(
        "--n-eff",
        dest="n_eff",
        type=float,
        default=None,
        help="independent rearrangements, for P(>=k precursors) and the seen/unseen "
        "clonotype counts",
    )
    pc.add_argument(
        "--selection",
        default="1.0",
        help="depth selection factor for the occupancy counts: a float, or 'auto' for "
        "the measured per-chain values (TRB 4.62, TRA 1.07)",
    )
    pc.add_argument(
        "--threads", type=int, default=0, help="worker threads (0 = all cores)"
    )
    pc.add_argument(
        "-v", "--verbose", action="store_true", help="per-group progress bar"
    )
    pc.set_defaults(func=_cmd_precursor)

    a = p.parse_args(argv)
    import os

    previous_threads = os.environ.get("POLARS_MAX_THREADS")
    if a.cmd in {"match", "update"}:
        import os
        import json

        threads = getattr(a, "threads", 1)
        if getattr(a, "config", None) and not hasattr(a, "threads"):
            threads = json.loads(Path(a.config).read_text()).get("threads", 1)
        if isinstance(threads, bool) or not isinstance(threads, int) or threads < 0:
            p.error("threads must be a nonnegative integer")
        if threads > 0:
            os.environ["POLARS_MAX_THREADS"] = str(threads)
    try:
        return a.func(a)
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        p.error(str(exc))
    finally:
        if previous_threads is None:
            os.environ.pop("POLARS_MAX_THREADS", None)
        else:
            os.environ["POLARS_MAX_THREADS"] = previous_threads


if __name__ == "__main__":
    raise SystemExit(main())
