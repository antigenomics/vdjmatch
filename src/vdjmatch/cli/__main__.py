"""vdjmatch command-line interface.

Subcommands:
  update     fetch/cache the latest VDJdb release
  match      annotate query sample(s) against VDJdb (E-values + ranked hits + epitope summary)
  precursor  T-cell precursor frequency and unseen-junction diversity for a set of TCRs

``match`` writes three TSV files per sample, prefixed with ``--output-prefix`` and the sample name:
  <prefix>.<sample>.hits.txt      per-hit CDR3 alignment (CIGAR, edits, score) for every VDJdb match
  <prefix>.<sample>.calls.txt     one predicted epitope per clonotype with its control-calibrated E-value
  <prefix>.<sample>.summary.txt   epitope-level enrichment (unique clonotypes, reads)

``precursor`` writes one TSV, one row per group, and needs the optional extra:
  pip install 'vdjmatch[precursor]'
"""
from __future__ import annotations

import argparse
import resource
import sys
import time
from pathlib import Path

from .. import db, match
from ..config import Params
from ..runner.multisample import annotate_sample


def _cmd_update(a: argparse.Namespace) -> int:
    path = db.fetch_latest(asset=a.asset, pin=a.pin, force=a.force)
    print(f"VDJdb {a.asset} cached at {path}")
    return 0


def _peak_rss_gb() -> float:
    """Peak resident set in GB (macOS reports ru_maxrss in bytes, Linux in KiB)."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_bytes = rss if sys.platform == "darwin" else rss * 1024
    return rss_bytes / 1e9


# match flags whose JSON-config value is overridable from the CLI; sentinel default -> "not passed".
_PARAM_FLAGS = ("scope", "matrix", "min_score", "match_v", "match_j", "evalue", "align",
                "species", "threads")


def _resolve_params(a: argparse.Namespace) -> Params:
    """Start from defaults or ``--config`` JSON, then apply any explicitly-passed CLI flag."""
    params = Params.from_json(a.config) if a.config else Params.defaults()
    for field in _PARAM_FLAGS:
        v = getattr(a, field, None)
        if v is not None:                              # argparse.SUPPRESS -> absent unless passed
            setattr(params, field, v)
    return params


def _cmd_match(a: argparse.Namespace) -> int:
    p = _resolve_params(a)
    vdj = db.load(a.vdjdb, asset=a.asset, species=p.species, min_score=p.min_score, pin=a.pin)
    print(f"VDJdb: {vdj.height:,} records "
          f"(species={p.species}, min_score={p.min_score})", file=sys.stderr)
    index = match.VdjdbIndex.build(vdj, species=p.species)
    print(f"indexed genes: {index.genes}", file=sys.stderr)
    matrix = match.load_vdjam() if p.matrix == "vdjam" else ""
    outdir = Path(a.output_prefix).parent
    outdir.mkdir(parents=True, exist_ok=True)
    p.to_json(f"{a.output_prefix}.params.json")        # reproducibility: effective params per run
    for sample in a.samples:
        name = Path(sample).name.split(".")[0]
        print(f"[{name}] annotating ...", file=sys.stderr)
        t0 = time.perf_counter()
        res = annotate_sample(index, sample, scope=p.scope, matrix=matrix or None,
                              species=p.species, with_evalue=p.evalue,
                              match_v=p.match_v, match_j=p.match_j, align=p.align,
                              threads=p.threads, progress=a.verbose)
        for kind, frame in res.items():
            out = f"{a.output_prefix}.{name}.{kind}.txt"
            frame.write_csv(out, separator="\t")
        print(f"[{name}] {res['hits'].height} hit rows, "
              f"{res['summary'].height} epitopes -> {a.output_prefix}.{name}.*.txt", file=sys.stderr)
        if a.verbose:
            dt = time.perf_counter() - t0
            nq = res["calls"].height
            print(f"[{name}] {dt:.1f}s, {nq / dt:,.0f} queries/s, "
                  f"peak RSS {_peak_rss_gb():.2f} GB", file=sys.stderr)
    return 0


def _cmd_precursor(a: argparse.Namespace) -> int:
    import polars as pl

    from ..precursor import summarise

    if a.vdjdb or not a.samples:
        frame = db.load(a.table, asset="slim", species=a.species, pin=a.pin)
        if a.mhc_class:
            frame = frame.filter(pl.col("mhc_class") == a.mhc_class)
        if not frame.height:
            raise SystemExit(f"VDJdb returned 0 records for species={a.species!r} "
                             f"mhc_class={a.mhc_class!r} -- note species is 'HomoSapiens', "
                             "not 'human'")
        junction_col, group_col = "cdr3", a.group_by or "epitope"
        chain_col = a.chain_col or "gene"
        capture_col = None if a.no_unseen else (a.capture_col or "reference_id")
        label = f"VDJdb ({frame.height:,} records)"
    else:
        sep = "," if a.samples[0].endswith(".csv") else "\t"
        frame = pl.concat([pl.read_csv(s, separator=sep, infer_schema_length=0)
                           for s in a.samples], how="vertical_relaxed")
        junction_col, group_col = a.junction_col, a.group_by
        chain_col, capture_col = a.chain_col, (None if a.no_unseen else a.capture_col)
        label = f"{len(a.samples)} file(s) ({frame.height:,} rows)"
    print(f"scoring {label}", file=sys.stderr)

    t0 = time.perf_counter()
    out = summarise(frame, junction_col=junction_col, group_col=group_col, chain_col=chain_col,
                    capture_col=capture_col, locus=a.locus, source=a.source, organism=a.organism,
                    r=a.radius, alpha=a.alpha, q=a.q, n_cells=a.n_cells,
                    compartment=a.compartment, n_eff=a.n_eff, min_junctions=a.min_junctions,
                    threads=a.threads, progress=a.verbose)
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(a.output, separator="\t")
    print(f"{out.height} group(s) -> {a.output} in {time.perf_counter() - t0:.1f}s "
          f"(peak RSS {_peak_rss_gb():.2f} GB)", file=sys.stderr)
    if a.q == 1.0:
        print("note: q=1 -> F is an uncalibrated model mass; only its ranking is meaningful. "
              "Pass --q to apply a selection constant.", file=sys.stderr)
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

notes:
  Sequences must be JUNCTIONS (Cys104..Phe/Trp118 inclusive), not IMGT CDR3s -- VDJdb's column is
  named `cdr3` but holds junctions. Anchor-stripped input is dropped and counted in `n_dropped`.
  `--q` is the selection constant carrying a generation probability to a repertoire frequency;
  left at 1 the reported F is a raw model mass whose ranking, not scale, is meaningful.
"""


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="vdjmatch", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=_EXAMPLES)
    sub = p.add_subparsers(dest="cmd", required=True)

    up = sub.add_parser("update", help="fetch/cache the latest VDJdb release")
    up.add_argument("--asset", default="slim", choices=["slim", "full", "default"],
                    help="which VDJdb table to fetch (default: slim)")
    up.add_argument("--pin", default=None, help="pin a specific release tag (default: latest)")
    up.add_argument("--force", action="store_true", help="re-download even if already cached")
    up.set_defaults(func=_cmd_update)

    m = sub.add_parser("match", help="annotate sample(s) against VDJdb",
                       formatter_class=argparse.RawDescriptionHelpFormatter, epilog=_EXAMPLES)
    m.add_argument("samples", nargs="+", help="AIRR rearrangement sample file(s) (TSV)")
    m.add_argument("-o", "--output-prefix", default="vdjmatch_out",
                   help="output path prefix; files are <prefix>.<sample>.{hits,calls,summary}.txt")
    m.add_argument("--vdjdb", default=None, help="custom VDJdb table path (default: fetch latest)")
    m.add_argument("--asset", default="full", choices=["slim", "full", "default"],
                   help="VDJdb table to match against (default: full — carries epitope/MHC detail)")
    m.add_argument("--pin", default=None, help="pin a specific VDJdb release tag (default: latest)")
    m.add_argument("--config", default=None,
                   help="JSON params file (see Params); explicit CLI flags below override its values")
    m.add_argument("-v", "--verbose", action="store_true",
                   help="show a per-sample progress bar and print wall time / queries-per-s / peak RSS")
    # --- params (overridable from --config; SUPPRESS default => only override JSON when passed) ---
    m.add_argument("--species", default=argparse.SUPPRESS, help="species filter (default: HomoSapiens)")
    m.add_argument("--scope", default=argparse.SUPPRESS,
                   help="search budget as substitutions,insertions,deletions,total-edits (default: 1,0,0,1)")
    m.add_argument("--matrix", default=argparse.SUPPRESS, choices=["vdjam", "none"],
                   help="substitution matrix: vdjam (TCR-specific, bundled) or none (unit edit cost)")
    m.add_argument("--min-score", dest="min_score", type=int, default=argparse.SUPPRESS,
                   help="minimum VDJdb confidence score (default: 0)")
    m.add_argument("--match-v", dest="match_v", action="store_const", const=True, default=argparse.SUPPRESS,
                   help="require the V gene to match as well as the CDR3")
    m.add_argument("--match-j", dest="match_j", action="store_const", const=True, default=argparse.SUPPRESS,
                   help="require the J gene to match as well as the CDR3")
    m.add_argument("--no-evalue", dest="evalue", action="store_const", const=False, default=argparse.SUPPRESS,
                   help="skip the control-calibrated E-value")
    m.add_argument("--no-align", dest="align", action="store_const", const=False, default=argparse.SUPPRESS,
                   help="skip the per-hit CIGAR/alignment output")
    m.add_argument("--threads", type=int, default=argparse.SUPPRESS,
                   help="worker threads (0 = all available cores)")
    m.set_defaults(func=_cmd_match)

    pc = sub.add_parser("precursor", help="precursor frequency + unseen diversity for a TCR set",
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        epilog=_PRECURSOR_EPILOG)
    pc.add_argument("samples", nargs="*", help="TCR table(s) (TSV/CSV); omit with --vdjdb")
    pc.add_argument("--vdjdb", action="store_true", help="score VDJdb itself instead of a file")
    pc.add_argument("--table", default=None, help="custom VDJdb table path (default: fetch latest)")
    pc.add_argument("--pin", default=None, help="pin a specific VDJdb release tag")
    pc.add_argument("--species", default="HomoSapiens", help="VDJdb species (default: HomoSapiens)")
    pc.add_argument("--mhc-class", dest="mhc_class", default=None, choices=["MHCI", "MHCII"],
                    help="restrict VDJdb to one MHC class")
    pc.add_argument("-o", "--output", default="vdjmatch_precursor.txt", help="output TSV path")
    # --- input columns ---
    pc.add_argument("--junction-col", dest="junction_col", default="cdr3",
                    help="junction column (Cys104..Phe/Trp118 inclusive; default: cdr3)")
    pc.add_argument("--group-by", dest="group_by", default=None,
                    help="column to group by, one row of output per group (e.g. epitope)")
    pc.add_argument("--chain-col", dest="chain_col", default=None,
                    help="locus column; each chain is scored with its own model")
    pc.add_argument("--capture-col", dest="capture_col", default=None,
                    help="capture-unit column for the unseen-species fields (e.g. reference_id)")
    pc.add_argument("--no-unseen", dest="no_unseen", action="store_true",
                    help="skip the Horvitz-Thompson unseen-species estimate")
    # --- model + estimator ---
    pc.add_argument("--locus", default="TRB", help="model locus when --chain-col is absent")
    pc.add_argument("--source", default="olga", help="recombination model set (default: olga)")
    pc.add_argument("--organism", default="human", help="model organism (mouse: --source arda)")
    pc.add_argument("-r", "--radius", type=int, default=1,
                    help="neighbourhood radius in substitutions (default: 1, part of the estimator)")
    pc.add_argument("--alpha", type=float, default=0.1,
                    help="cognacy retention per edit (default: 0.1, Mayer & Callan 2023)")
    pc.add_argument("--q", type=float, default=1.0,
                    help="selection constant; 1 = uncalibrated raw mass (ALICE's TRB value: 9.41)")
    pc.add_argument("--min-junctions", dest="min_junctions", type=int, default=1,
                    help="skip groups with fewer distinct junctions")
    # --- absolute numbers ---
    pc.add_argument("--n-cells", dest="n_cells", type=float, default=1e11,
                    help="total T cells, for the expected precursor count (default: 1e11)")
    pc.add_argument("--compartment", type=float, default=1.0,
                    help="fraction of the pool the restriction addresses (e.g. 0.3 for CD8)")
    pc.add_argument("--n-eff", dest="n_eff", type=float, default=None,
                    help="independent rearrangements, for P(>=k precursors)")
    pc.add_argument("--threads", type=int, default=0, help="worker threads (0 = all cores)")
    pc.add_argument("-v", "--verbose", action="store_true", help="per-group progress bar")
    pc.set_defaults(func=_cmd_precursor)

    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
