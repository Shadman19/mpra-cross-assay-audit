#!/usr/bin/env bash
# Run from an OnDemand terminal in the extracted package. Use an ODU-supported
# Python 3.11/3.12 interpreter; never install over your previous experiment env.
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONNOUSERSITE=1
unset PYTHONPATH
AUDIT_BASE_PYTHON="${CARDIAC_PYTHON:-python3}"
AUDIT_BASE_PYTHON="$(command -v "$AUDIT_BASE_PYTHON")"
"$AUDIT_BASE_PYTHON" -c 'import sys; assert (3,11)<=sys.version_info[:2]<=(3,12), "Use Python 3.11 or 3.12"'
"$AUDIT_BASE_PYTHON" -m venv .venv
export CARDIAC_PYTHON="$PWD/.venv/bin/python"
"$CARDIAC_PYTHON" -m pip install --upgrade pip
# CUDA 11.8 wheels work with a suitably recent NVIDIA driver without compiling.
"$CARDIAC_PYTHON" -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu118
"$CARDIAC_PYTHON" -m pip install -r requirements.txt
"$CARDIAC_PYTHON" code/download_reference.py
"$CARDIAC_PYTHON" code/prepare_data.py
"$CARDIAC_PYTHON" -m unittest discover -s tests -v
"$CARDIAC_PYTHON" code/preflight.py --metadata-only
printf 'Setup complete. In your submitting shell run: export CARDIAC_PYTHON="%s/.venv/bin/python"\n' "$PWD"
