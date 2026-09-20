#!/bin/bash
#SBATCH --job-name=sk2b_factorial
#SBATCH --array=0-23%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/factorial/factorial_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/factorial/factorial_%A_%a.err

set -euo pipefail
HPC=/lustre1/g/aos_shihuang/sk2bgrow-hpc
SRC=$HPC/src
BENCH=$HPC/bench/factorial
mkdir -p "$BENCH/results" "$HPC/logs/factorial"
export PYTHONPATH="$SRC/python${PYTHONPATH:+:$PYTHONPATH}"
PY=$HPC/micromamba/envs/pilea/bin/python
"$PY" "$SRC/benches/factorial/simulate_factorial.py" \
  --task-id "$SLURM_ARRAY_TASK_ID" \
  --reps 10 \
  --out-dir "$BENCH/results"
