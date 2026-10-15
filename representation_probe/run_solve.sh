#!/bin/bash
#SBATCH -J solve -c 16 --mem=96G -t 8:00:00 -o /home/spath004/cardiac_r21/representation_probe/logs/solve-%j.out
cd ~/cardiac_r21
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
~/cardiac_r21/env/bin/python representation_probe/solve.py --features ~/cardiac_r21/probe_features --rows ~/cardiac_r21/repo/mpra-cross-assay-audit/subgroup_audit/results/annotated_audit_variants.csv.gz --out ~/cardiac_r21/probe_results/solve
