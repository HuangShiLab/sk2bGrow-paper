#!/bin/bash
set -euo pipefail
HPC=/lustre1/g/aos_shihuang/sk2bgrow-hpc
SRC=$HPC/src
BENCH=$HPC/bench/factorial
mkdir -p "$BENCH/results" "$HPC/logs/factorial"
chmod 0755 "$SRC/benches/factorial/simulate_factorial.py" \
  "$SRC/benches/factorial/aggregate_factorial.py"
JOB=$(sbatch --parsable "$SRC/benches/factorial/run_factorial.sbatch")
echo "$JOB" > "$BENCH/latest_array_job_id"
sbatch --dependency=afterok:"$JOB" --job-name=sk2b_factorial_agg \
  --cpus-per-task=1 --mem=2G --time=00:20:00 \
  --output="$HPC/logs/factorial/aggregate_%j.out" \
  --error="$HPC/logs/factorial/aggregate_%j.err" \
  --wrap="export PYTHONPATH=$SRC/python:\$PYTHONPATH; $HPC/micromamba/envs/pilea/bin/python $SRC/benches/factorial/aggregate_factorial.py --results-dir $BENCH/results --long-output $BENCH/factorial_long.tsv --contrast-output $BENCH/factorial_contrasts.tsv"
echo "submitted factorial array $JOB and aggregate dependency job"
