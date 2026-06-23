"""vdjmatch command-line interface.

Subcommands:
  update   fetch/cache the latest VDJdb release
  match    annotate query sample(s) against VDJdb (E-values + ranked hits + epitope summary)

``match`` writes three TSV files per sample, prefixed with ``--output-prefix`` and the sample name:
  <prefix>.<sample>.hits.txt      per-hit CDR3 alignment (CIGAR, edits, score) for every VDJdb match
  <prefix>.<sample>.calls.txt     one predicted epitope per clonotype with its control-calibrated E-value
  <prefix>.<sample>.summary.txt   epitope-level enrichment (unique clonotypes, reads)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .. import db, match
from ..runner.multisample import annotate_sample


def _cmd_update(a: argparse.Namespace) -> int:
    path = db.fetch_latest(asset=a.asset, pin=a.pin, force=a.force)
    print(f"VDJdb {a.asset} cached at {path}")
    return 0


def _cmd_match(a: argparse.Namespace) -> int:
    vdj = db.load(a.vdjdb, asset=a.asset, species=a.species, min_score=a.min_score, pin=a.pin)
    print(f"VDJdb: {vdj.height:,} records "
          f"(species={a.species}, min_score={a.min_score})", file=sys.stderr)
    index = match.VdjdbIndex.build(vdj, species=a.species)
    print(f"indexed genes: {index.genes}", file=sys.stderr)
    matrix = match.load_vdjam() if a.matrix == "vdjam" else ""
    outdir = Path(a.output_prefix).parent
    outdir.mkdir(parents=True, exist_ok=True)
    for sample in a.samples:
        name = Path(sample).name.split(".")[0]
        print(f"[{name}] annotating ...", file=sys.stderr)
        res = annotate_sample(index, sample, scope=a.scope, matrix=matrix or None,
                              species=a.species, with_evalue=not a.no_evalue,
                              match_v=a.match_v, match_j=a.match_j, align=not a.no_align,
                              threads=a.threads)
        for kind, frame in res.items():
            out = f"{a.output_prefix}.{name}.{kind}.txt"
            frame.write_csv(out, separator="\t")
        print(f"[{name}] {res['hits'].height} hit rows, "
              f"{res['summary'].height} epitopes -> {a.output_prefix}.{name}.*.txt", file=sys.stderr)
    return 0


_EXAMPLES = """\
examples:
  vdjmatch update                                  # cache the latest VDJdb (slim)
  vdjmatch match sample.tsv                         # annotate one AIRR sample, default reference
  vdjmatch match -o run/out --match-v *.tsv         # match V gene too, write under run/
  vdjmatch match --scope 2,1,1,2 --threads 8 s.tsv  # wider search budget, 8 threads
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
    m.add_argument("--species", default="HomoSapiens", help="species filter (default: HomoSapiens)")
    m.add_argument("--scope", default="1,0,0,1",
                   help="search budget as substitutions,insertions,deletions,total-edits (default: 1,0,0,1)")
    m.add_argument("--matrix", default="vdjam", choices=["vdjam", "none"],
                   help="substitution matrix: vdjam (TCR-specific, bundled) or none (unit edit cost)")
    m.add_argument("--min-score", type=int, default=0, help="minimum VDJdb confidence score (default: 0)")
    m.add_argument("--match-v", action="store_true", help="require the V gene to match as well as the CDR3")
    m.add_argument("--match-j", action="store_true", help="require the J gene to match as well as the CDR3")
    m.add_argument("--no-evalue", action="store_true", help="skip the control-calibrated E-value")
    m.add_argument("--no-align", action="store_true", help="skip the per-hit CIGAR/alignment output")
    m.add_argument("--threads", type=int, default=0, help="worker threads (0 = all available cores)")
    m.set_defaults(func=_cmd_match)

    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
