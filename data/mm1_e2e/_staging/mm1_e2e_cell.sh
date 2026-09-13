#!/usr/bin/env bash
#SBATCH --job-name=mm1_e2e
#SBATCH --partition=amd
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/mm1_e2e/cell.%A_%a.out
# Task 2 (review R1): full-sample mm=1 end-to-end benchmark, C5 MG samples.
# Methodology mirrors baseline 01_cell.sh exactly (two timed stages), with
# --max-mismatch 1 added to the count stage:
#   stage 1: sk2bgrow profile --no-stats --threads 8 --quiet  (== baseline count)
#   stage 2: python -m sk2bgrow.cli profile ... --use-rust-windows --count-model ztp
#            (== baseline stats)
# Baseline count (mm=2, src build) covered stage 1 only; e2e = stage1 + stage2,
# same decomposition as c5_cost_table.py "sk2bgrow_total_min".
# The db is opened read-only; counts are kept (not rm'd) for anchor-retention.
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
C5=$BASE/bench/C5
E2E=$BASE/bench/mm1_e2e
BIN=$BASE/src_m4/target/release/sk2bgrow
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
export PYTHONPATH=$BASE/src_m4/python:$PYTHONPATH

t=$SLURM_ARRAY_TASK_ID
run=$(awk -F"\t" -v task="$t" '$1 == task {print $2}' "$C5/samples.tsv")
[ -n "$run" ] || { echo "no sample for task $t"; exit 1; }
fq1=$BASE/fastq/PRJNA974210/${run}_1.fastq.gz
fq2=$BASE/fastq/PRJNA974210/${run}_2.fastq.gz
[ -s "$fq1" ] && [ -s "$fq2" ] || { echo "missing fastq for $run"; exit 1; }

out=$E2E/res/A_${run}
mkdir -p "$out"

# ---- stage 1: count at mm=1 (same口径 as baseline count stage) ----
if [ -f "$out/time.txt" ] && grep -qa "Exit status: 0" "$out/time.txt"; then
  echo "skip count $run (time.txt ok)"
else
  env SK2B_COUNT_TIMING=1 /usr/bin/time -v "$BIN" profile \
    "$fq1" "$fq2" -d "$C5/db" -o "$out" \
    --no-stats --threads 8 --quiet --max-mismatch 1 \
    2> "$out/time.txt"
fi
echo "--- count $run"; grep -a "phase-timing\|Elapsed\|Maximum resident\|Exit status" "$out/time.txt" || true

# ---- stage 2: stats (identical invocation to baseline 01_cell.sh) ----
if [ -f "$out/output.tsv" ]; then
  echo "skip stats $run (output.tsv exists)"
else
  /usr/bin/time -v "$PY" -m sk2bgrow.cli profile "$out"/*.counts.tsv \
    --db "$C5/db" --output "$out" --windows "$out/windows.tsv" \
    --use-rust-windows --count-model ztp \
    > "$out/stats.log" 2> "$out/stats.time.txt"
fi
echo "--- stats $run"; grep -a "Elapsed\|Maximum resident\|Exit status" "$out/stats.time.txt" || true

# ---- stage 3: per-genome anchor-retention analysis (mm1 vs mm2 baseline) ----
"$PY" "$E2E/scripts/analyze_sample.py" "$run" || echo "ANALYSIS_FAILED $run"

echo "MM1_E2E_DONE task=$t run=$run"
