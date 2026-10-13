#!/usr/bin/env bash
#SBATCH --job-name=reliability-audit
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/audit-%j.out
#SBATCH --error=logs/audit-%j.err
# Uses ODU's default CPU partition, not an invented partition name.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export CARDIAC_PYTHON="${CARDIAC_PYTHON:-$PWD/.venv/bin/python}"
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH
mkdir -p runtime/locks
exec 9>runtime/locks/primary-audit.lock
flock -n 9 || { echo "Primary audit already running"; exit 1; }
"$CARDIAC_PYTHON" code/audit.py --shards 8 --bootstrap 10000
"$CARDIAC_PYTHON" code/build_preprint.py
