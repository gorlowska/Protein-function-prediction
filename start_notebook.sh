#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if [[ ! -x .venv/bin/python ]]; then
  printf 'Create the environment first: bash scripts/bootstrap.sh --skip-mmseqs\n' >&2
  exit 1
fi
exec .venv/bin/python -m jupyter lab "${1:-meeting_02_data_preparation.ipynb}"
