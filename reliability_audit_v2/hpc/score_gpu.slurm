#!/usr/bin/env bash
#SBATCH --job-name=reliability-v2
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=05:59:00
#SBATCH --array=0-15%2
#SBATCH --signal=B:USR1@300
#SBATCH --output=logs/score-%A_%a.out
#SBATCH --error=logs/score-%A_%a.err
# Each task uses one GPU and one model. No minimum GPU-memory request is made:
# preflight reports actual hardware and stops visibly on insufficient resources.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the package root}"
export CARDIAC_PYTHON="${CARDIAC_PYTHON:-$PWD/.venv/bin/python}"
export PYTHONNOUSERSITE=1 CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH
: "${SLURM_ARRAY_TASK_ID:?Use the provided array submission}"
AUDIT_SHARD=$(( SLURM_ARRAY_TASK_ID % 8 ))
AUDIT_MODEL=enformer
if (( SLURM_ARRAY_TASK_ID >= 8 )); then AUDIT_MODEL=borzoi; fi
mkdir -p runtime/locks
# No two jobs may append to the same checkpoint or verify the same model at once.
exec 9>"runtime/locks/${AUDIT_MODEL}-shard${AUDIT_SHARD}.lock"
flock -n 9 || { echo 'This shard is already running'; exit 1; }
(
  flock 8
  "$CARDIAC_PYTHON" code/preflight.py --gpu --model "$AUDIT_MODEL"
) 8>"runtime/locks/${AUDIT_MODEL}-preflight.lock"
"$CARDIAC_PYTHON" code/score_open_models.py --model "$AUDIT_MODEL" --shard "$AUDIT_SHARD" --shards 8 --max-seconds 19200 &
AUDIT_CHILD=$!
trap 'kill -USR1 "$AUDIT_CHILD" 2>/dev/null || true' USR1 TERM
set +e
wait "$AUDIT_CHILD"
AUDIT_EXIT=$?
# A trapped signal can interrupt wait before Python finishes committing its row.
if kill -0 "$AUDIT_CHILD" 2>/dev/null; then wait "$AUDIT_CHILD"; AUDIT_EXIT=$?; fi
exit "$AUDIT_EXIT"
