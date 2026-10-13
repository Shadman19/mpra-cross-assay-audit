#!/usr/bin/env bash
#SBATCH --job-name=subgroup-audit
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/subgroup-%j.out
#SBATCH --error=logs/subgroup-%j.err
# No gpu partition, GPU allocation, foundation-model inference or API access.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the extracted Subgroup_Audit directory}"
: "${CARDIAC_PYTHON:?Export your existing v2 Python interpreter}"
: "${RELIABILITY_BASE_DIR:?Export the path to the completed Reliability_Audit_v2}"
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH
mkdir -p results
exec 9>results/subgroup_audit.lock
flock -n 9 || { echo 'Subgroup audit already running'; exit 1; }
"$CARDIAC_PYTHON" -m unittest discover -s tests -v
"$CARDIAC_PYTHON" code/run_audit.py --base "$RELIABILITY_BASE_DIR"
"$CARDIAC_PYTHON" code/assay_descriptives.py
"$CARDIAC_PYTHON" code/update_outline.py
if command -v pdflatex >/dev/null 2>&1; then
  bash conference/compile.sh || echo "CPU audit complete; inspect the optional TeX compile log before using the PDF."
else
  echo "CPU audit complete and manuscript sources updated; compile the conference folder on Overleaf."
fi
