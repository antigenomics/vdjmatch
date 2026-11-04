#!/bin/sh
# Purpose: reproduce the project-local development environment for vdjmatch.
set -eu
cd "$(dirname "$0")"
DO_TESTS=0
for arg in "$@"; do
    case "$arg" in
        --test|--tests) DO_TESTS=1 ;;
        *) echo "Unknown flag: $arg" >&2; exit 2 ;;
    esac
done
if [ ! -x .venv/bin/python ]; then
    uv venv .venv
fi
extras='.[test,control]'
if [ "$DO_TESTS" -eq 1 ]; then
    extras='.[test,control,bench,diagnostics]'
fi
uv pip install --python .venv/bin/python -e "$extras"
if [ "$DO_TESTS" -eq 1 ]; then
    .venv/bin/python -c 'import vdjtools, pytest, sys; raise SystemExit(pytest.main(sys.argv[1:]))' tests/unit -q -m "not slow"
fi
