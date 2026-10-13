#!/bin/sh
# Raw synthetic queries -> two fresh detailed CLI geometries, no fitted models.
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
OUT=${1:-comparison-demo}
CLI=${VDJMATCH:-vdjmatch}
command -v "$CLI" >/dev/null
if [ -e "$OUT" ]; then echo "Choose a fresh output directory: $OUT" >&2; exit 1; fi
mkdir -p "$OUT"
"$CLI" annotate "$HERE/inputs/cohort/sample_a.tsv" --vdjdb "$HERE/inputs/reference.tsv" \
  --threads 1 --exclude-exact --details --output-prefix "$OUT/stem4"
"$CLI" annotate "$HERE/inputs/cohort/sample_a.tsv" --vdjdb "$HERE/inputs/reference.tsv" \
  --threads 1 --exclude-exact --details --v-regions cdr12 --stem-depth 0 \
  --distance-matrix tcrdist --junction-ends trim3 --output-prefix "$OUT/vcdr12"
printf 'Fresh comparison tables written under %s\n' "$OUT"
