#!/bin/sh
# Purpose: reproduce the project-local development environment for vdjmatch.
set -eu
cd "$(dirname "$0")"
if [ ! -x .venv/bin/python ]; then
    uv venv .venv
fi
uv pip install --python .venv/bin/python -e '.[test,control]'
if [ "${1:-}" = "--test" ]; then
    .venv/bin/python -c 'import subprocess; subprocess.run([".venv/bin/python", "-m", "pytest", "tests/unit", "-q"], timeout=120, check=True)'
fi
