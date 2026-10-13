#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
: "${CARDIAC_PYTHON:?Set your existing v2 interpreter}"
: "${RELIABILITY_BASE_DIR:?Set your completed v2 package path}"
sbatch --export=ALL hpc/subgroup_cpu.slurm
