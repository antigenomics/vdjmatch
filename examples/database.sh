#!/bin/sh
# Synthetic queries and null draws -> database presence and nearest-distance tails.
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
OUT=${1:-database-demo}
CLI=${VDJMATCH:-vdjmatch}
command -v "$CLI" >/dev/null
if [ -e "$OUT" ]; then echo "Choose a fresh output directory: $OUT" >&2; exit 1; fi
mkdir -p "$OUT"
"$CLI" database-evidence "$HERE/inputs/database_queries.tsv" \
  --vdjdb "$HERE/inputs/reference.tsv" --null-sample "$HERE/inputs/database_null.tsv" \
  --species human --locus TRB --max-edits 1 --threads 4 --output-prefix "$OUT/presence"
printf 'Fresh database diagnostics written under %s\n' "$OUT"
