#!/bin/sh
# Render accepted final enrichment exports; this does not refit tests or draw mocks.
set -eu
vdjmatch repertoire-report "$1" "$2" --metric "${4:-curve_mean}" --output-prefix "$3"
