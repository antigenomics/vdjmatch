#!/bin/sh
# Raw synthetic cohort and paired inputs -> fresh CLI annotation tables.
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
OUT=${1:-annotation-demo}
CLI=${VDJMATCH:-vdjmatch}
REFERENCE=${2:-"$HERE/inputs/reference.tsv"}
command -v "$CLI" >/dev/null
if [ -e "$OUT" ]; then echo "Choose a fresh output directory: $OUT" >&2; exit 1; fi
mkdir -p "$OUT"
"$CLI" annotate --metadata "$HERE/inputs/samples.tsv" --vdjdb "$REFERENCE" \
  --threads 4 --output-prefix "$OUT/cohort"
"$CLI" annotate "$HERE/inputs/paired.tsv" --locus paired --vdjdb "$REFERENCE" \
  --threads 4 --output-prefix "$OUT/paired"
printf 'Fresh annotation tables written under %s\n' "$OUT"
