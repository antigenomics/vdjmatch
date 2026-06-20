"""vdjmatch command-line interface.

Subcommands:
  update   fetch/cache the latest VDJdb release
  match    annotate query sample(s) against VDJdb (E-values + ranked hits + epitope summary)
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


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="vdjmatch", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    up = sub.add_parser("update", help="fetch/cache the latest VDJdb release")
    up.add_argument("--asset", default="slim", choices=["slim", "full", "default"])
    up.add_argument("--pin", default=None, help="pin a release tag")
    up.add_argument("--force", action="store_true")
    up.set_defaults(func=_cmd_update)

    m = sub.add_parser("match", help="annotate sample(s) against VDJdb")
    m.add_argument("samples", nargs="+", help="AIRR rearrangement sample file(s)")
    m.add_argument("-o", "--output-prefix", default="vdjmatch_out")
    m.add_argument("--vdjdb", default=None, help="VDJdb table path (default: fetch latest)")
    m.add_argument("--asset", default="full", choices=["slim", "full", "default"])
    m.add_argument("--pin", default=None)
    m.add_argument("--species", default="HomoSapiens")
    m.add_argument("--scope", default="1,0,0,1", help="s,i,d,t search budget")
    m.add_argument("--matrix", default="vdjam", choices=["vdjam", "none"])
    m.add_argument("--min-score", type=int, default=0, help="min VDJdb confidence")
    m.add_argument("--match-v", action="store_true")
    m.add_argument("--match-j", action="store_true")
    m.add_argument("--no-evalue", action="store_true")
    m.add_argument("--no-align", action="store_true")
    m.add_argument("--threads", type=int, default=0)
    m.set_defaults(func=_cmd_match)

    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
