#!/usr/bin/env bash
# One-command CPU reproduction of all exploratory diagnostics and a check of every
# manuscript number against machine-generated tables. Primary and follow-up audit
# results are read from the committed result files; GPU scoring is described in
# reliability_audit_v2/README.md and is not needed here.
set -euo pipefail
cd "$(dirname "$0")"
python3 -m venv .venv_cpu
.venv_cpu/bin/pip install -q -r exploratory/requirements.txt
.venv_cpu/bin/python -W ignore exploratory/diagnostics.py
.venv_cpu/bin/python exploratory/check_paper_numbers.py
