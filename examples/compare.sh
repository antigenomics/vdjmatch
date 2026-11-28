#!/bin/sh
# Compare exact-hit policies with the same vX scorer and reference.
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
OUT=${1:-comparison-demo}
CLI=${VDJMATCH:-vdjmatch}
command -v "$CLI" >/dev/null
if [ -e "$OUT" ]; then echo "Choose a fresh output directory: $OUT" >&2; exit 1; fi
mkdir -p "$OUT"
"$CLI" annotate "$HERE/inputs/cohort/sample_a.tsv" --vdjdb "$HERE/inputs/reference.tsv" \
  --threads 4 --output-prefix "$OUT/exact-retained"
"$CLI" annotate "$HERE/inputs/cohort/sample_a.tsv" --vdjdb "$HERE/inputs/reference.tsv" \
  --threads 4 --exclude-exact --output-prefix "$OUT/exact-excluded"
printf 'Fresh vX comparisons written under %s\n' "$OUT"
