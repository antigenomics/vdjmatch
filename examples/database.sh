#!/bin/sh
# Synthetic queries and one raw control -> fresh Pgen and both directional tails.
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
OUT=${1:-database-demo}
CLI=${VDJMATCH:-vdjmatch}
command -v "$CLI" >/dev/null
if [ -e "$OUT" ]; then echo "Choose a fresh output directory: $OUT" >&2; exit 1; fi
mkdir -p "$OUT"
"$CLI" database-evidence "$HERE/inputs/database_queries.tsv" \
  --vdjdb "$HERE/inputs/reference.tsv" --control "$HERE/inputs/database_control.tsv" \
  --species human --locus TRB --scope 1,0,0,1 --mock-replicates 199 --seed 0 \
  --threads 4 --output-prefix "$OUT/presence"
printf 'Fresh database diagnostics written under %s\n' "$OUT"
