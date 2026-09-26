#!/usr/bin/env bash
set -euo pipefail
# SLURM array worker for public real 2bRAD route-B profiling.
# Manifest lines: sample_id<TAB>fastq_path<TAB>output_dir
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
DB=$BASE/bench/ecc_routeB/db_gtdb_bcgI_sp
BIN=$BASE/src_m4/target/release/sk2bgrow
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
export PYTHONPATH=$BASE/src/python:$PYTHONPATH
MANIFEST=${MANIFEST:-$BASE/bench/public_2brad/profile_manifest.tsv}
SCRATCH=${SCRATCH:-/tmp/$USER/sk2bgrow_public_routeB}
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "Set SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")
IFS=$'\t' read -r SAMPLE FASTQ OUTDIR <<< "$LINE"
[[ -s "$FASTQ" ]] || { echo "Missing FASTQ: $FASTQ" >&2; exit 2; }
COUNT_DIR=$SCRATCH/${SLURM_ARRAY_TASK_ID}_${SAMPLE}
mkdir -p "$COUNT_DIR" "$OUTDIR"
echo "[public-routeB] sample=$SAMPLE fastq=$FASTQ threads=${SLURM_CPUS_PER_TASK:-4}"
"$BIN" profile "$FASTQ" \
  --db "$DB" --output "$COUNT_DIR" --mode 2brad --max-mismatch 0 \
  --threads "${SLURM_CPUS_PER_TASK:-4}" --no-stats --quiet
COUNT_FILE=$(find "$COUNT_DIR" -maxdepth 1 -type f -name '*.counts.tsv' | head -1)
[[ -n "$COUNT_FILE" ]] || { echo "No count table produced for $SAMPLE" >&2; exit 2; }
"$PY" -m sk2bgrow.cli profile "$COUNT_FILE" \
  --db "$DB" --output "$OUTDIR/noGC_glm" --no-gc-correct --method glm \
  >"$OUTDIR/stats.log" 2>&1
# Count tables are the largest intermediate and are reproducible from raw FASTQ.
find "$COUNT_DIR" -type f \( -name '*.counts.tsv' -o -name 'windows.tsv' -o -name '*.stats.json' \) -delete
rmdir "$COUNT_DIR" 2>/dev/null || true
echo "[public-routeB] complete sample=$SAMPLE output=$OUTDIR/noGC_glm/output.tsv"
