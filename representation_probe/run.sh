#!/usr/bin/env bash
# Representation probe on Wahab. Run from ~/cardiac_r21:  bash representation_probe/run.sh
# Safe to rerun: finished context groups are skipped, so a rerun resumes timed-out shards.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
BASE="$HOME/cardiac_r21/audit2/Reliability_Audit_v2"
REPO="$HOME/cardiac_r21/repo/mpra-cross-assay-audit"
ROWS="$REPO/subgroup_audit/results/annotated_audit_variants.csv.gz"
OUT="$HOME/cardiac_r21/probe_features"; RES="$HOME/cardiac_r21/probe_results"
CPU_PY="$HOME/cardiac_r21/env/bin/python"; GPU_PY="${CARDIAC_PYTHON:-$BASE/.venv/bin/python}"; [ -x "$GPU_PY" ] || GPU_PY="$CPU_PY"
mkdir -p "$OUT" "$RES" "$HERE/logs"
for f in "$BASE/runtime/hg38.2bit" "$ROWS" "$GPU_PY" "$CPU_PY"; do [ -e "$f" ] || { echo "Missing $f"; exit 1; }; done
GPU=$(sbatch --parsable --partition=gpu --gres=gpu:1 --cpus-per-task=4 --mem=64G --time=05:59:00 \
  --array=0-7%4 --signal=B:USR1@300 --job-name=probe-gpu --output="$HERE/logs/gpu-%A_%a.out" \
  --wrap="cd $BASE && export PYTHONNOUSERSITE=1 CUBLAS_WORKSPACE_CONFIG=:4096:8 OMP_NUM_THREADS=4; unset PYTHONPATH; \
          $GPU_PY $HERE/score_features.py --base $BASE --rows $ROWS --out $OUT --shard \$SLURM_ARRAY_TASK_ID --shards 8")
CPU=$(sbatch --parsable --dependency=afterok:$GPU --cpus-per-task=8 --mem=48G --time=02:00:00 --job-name=probe-cpu \
  --output="$HERE/logs/cpu-%j.out" \
  --wrap="export OMP_NUM_THREADS=8; $CPU_PY -W ignore $HERE/probe.py --features $OUT --rows $ROWS --audit-code $REPO/subgroup_audit/code --out $RES")
echo "Submitted GPU array $GPU and analysis $CPU."
echo "When done:  tail -25 $HERE/logs/cpu-$CPU.out"
