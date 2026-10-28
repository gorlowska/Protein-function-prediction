#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
PYTHON_BIN="${PYTHON_BIN:-python3}"
"$PYTHON_BIN" -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python scripts/setup_project.py "$@"
.venv/bin/python -m ipykernel install --sys-prefix --name bioinformatics \
  --display-name "Python (BioInformatics)"
.venv/bin/python scripts/verify_dataset.py
printf '\nReady. Open the notebook with: bash start_notebook.sh\n'
