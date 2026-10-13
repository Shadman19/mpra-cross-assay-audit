#!/usr/bin/env bash
# Safe to submit again: complete shards exit without inference; partial shards resume.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
export CARDIAC_PYTHON="${CARDIAC_PYTHON:-$PWD/.venv/bin/python}"
AUDIT_JOB_ID="$(sbatch --parsable --export=ALL hpc/score_gpu.slurm)"
AUDIT_JOB_ID="${AUDIT_JOB_ID%%;*}"
printf 'GPU array %s submitted. Inspect squeue and logs. Exit 75 means resumable, not complete.\n' "$AUDIT_JOB_ID"
sbatch --dependency="afterok:$AUDIT_JOB_ID" --export=ALL hpc/audit_cpu.slurm
