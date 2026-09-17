#!/usr/bin/env bash
# multiseed five-arm pipeline for one instance seed (85 cells):
#   A: sk2bgrow count (PE) -> stats --windows
#   B: same stats + --method sorted
#   E: armE_counts_pe.py (scale 104, pooled PE) -> stats (no --windows)
#   C: pilea default / relaxed (-x 0 -z 0 -c 0)
# All stages skip-if-exists; new files only under multiseed/seed<SEED>/.
set -uo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
RC=$BASE/bench/repro_check
MS=$RC/multiseed
S=$1  # instance seed
ROOT=$MS/seed$S
BIN=$BASE/src/target/release/sk2bgrow
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
PILEA=$BASE/bench/C1/vendor/env-pilea138/bin/pilea
PILEA_PY=$BASE/bench/C1/vendor/env-pilea138/bin/python
DB=$BASE/bench/C1/db
PILEADB=$BASE/bench/C1/pileadb
CACHE=$BASE/bench/F1/cache_s104.pkl
DRIVER=$RC/armE_counts_pe.py
export BASE ROOT BIN PY PILEA PILEA_PY DB PILEADB CACHE DRIVER
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
MEDIA="M6 M3 M2 RUN_OUT M27 M25 M24 M22 M23 M19 M18 M17 M13 M12 M1 M10 M4"
DEPTHS="0.5 1 2 5 10"

count_a() { # depth medium
  local d=$1 m=$2
  local srr; srr=$(awk -v m="$m" '$2==m{print $1; exit}' "$BASE/bench/C1/run_medium.tsv")
  local out="$ROOT/counts/$srr.${d}x"
  mkdir -p "$out"
  [ -s "$out/$srr.${d}x.counts.tsv" ] && { echo "$d $m skip"; return 0; }
  "$BIN" profile "$ROOT/fq/$srr.${d}x_1.fq.gz" "$ROOT/fq/$srr.${d}x_2.fq.gz" \
    -d "$DB" -o "$out" --no-stats --threads 4 --quiet 2> "$out/time.txt" \
    && echo "$d $m rc=0" || echo "$d $m FAIL"
}
stats_ab() { # arm(A|B) depth medium
  local arm=$1 d=$2 m=$3
  local srr; srr=$(awk -v m="$m" '$2==m{print $1; exit}' "$BASE/bench/C1/run_medium.tsv")
  local cell="$ROOT/counts/$srr.${d}x"
  local out="$ROOT/out_${arm}/$srr.${d}x"
  mkdir -p "$out"
  [ -s "$out/output.tsv" ] && { echo "$arm $d $m skip"; return 0; }
  if [ "$arm" = A ]; then
    "$PY" -m sk2bgrow.cli profile "$cell/$srr.${d}x.counts.tsv" --db "$DB" \
      --output "$out" --windows "$cell/windows.tsv" > "$out/stats.log" 2>&1
  else
    "$PY" -m sk2bgrow.cli profile "$cell/$srr.${d}x.counts.tsv" --db "$DB" \
      --output "$out" --windows "$cell/windows.tsv" --method sorted > "$out/stats.log" 2>&1
  fi
  echo "$arm $d $m rc=$?"
}
sketch_e() { # depth medium
  local d=$1 m=$2
  local srr; srr=$(awk -v m="$m" '$2==m{print $1; exit}' "$BASE/bench/C1/run_medium.tsv")
  local name="$srr.${d}x_pe.fq"
  [ -s "$ROOT/sketch/$name.counts.tsv" ] && { echo "$d $m skip"; return 0; }
  "$PILEA_PY" "$DRIVER" -s 104 --cache "$CACHE" -o "$ROOT/sketch" --name "$name" \
    "$ROOT/fq/$srr.${d}x_1.fq.gz" "$ROOT/fq/$srr.${d}x_2.fq.gz" > /dev/null 2>&1
  echo "$d $m rc=$?"
}
stats_e() { # depth medium
  local d=$1 m=$2
  local srr; srr=$(awk -v m="$m" '$2==m{print $1; exit}' "$BASE/bench/C1/run_medium.tsv")
  local counts="$ROOT/sketch/$srr.${d}x_pe.fq.counts.tsv"
  local out="$ROOT/out_E/$srr.${d}x"
  mkdir -p "$out"
  [ -s "$out/output.tsv" ] && { echo "$d $m skip"; return 0; }
  [ -s "$counts" ] || { echo "$d $m NO_COUNTS"; return 1; }
  "$PY" -m sk2bgrow.cli profile "$counts" --db "$DB" --output "$out" \
    > "$out/stats.log" 2>&1
  echo "$d $m rc=$?"
}
pilea_c() { # mode depth medium
  local mode=$1 d=$2 m=$3
  local srr; srr=$(awk -v m="$m" '$2==m{print $1; exit}' "$BASE/bench/C1/run_medium.tsv")
  local out="$ROOT/pilea/$mode/$srr.${d}x"
  mkdir -p "$out"
  [ -s "$out/output.tsv" ] && { echo "$mode $d $m skip"; return 0; }
  local extra=()
  [ "$mode" = relaxed ] && extra=(-x 0 -z 0 -c 0)
  "$PILEA" profile "$ROOT/fq/$srr.${d}x_1.fq.gz" "$ROOT/fq/$srr.${d}x_2.fq.gz" \
    -d "$PILEADB" -o "$out" -t 4 "${extra[@]}" > "$out/pilea.log" 2>&1
  echo "$mode $d $m rc=$?"
}
export -f count_a stats_ab sketch_e stats_e pilea_c

jobs() { for d in $DEPTHS; do for m in $MEDIA; do echo "$d $m"; done; done; }

{
  echo "start pipeline seed=$S $(date)"
  echo "-- stage count_A"; jobs | xargs -P 8 -n 2 bash -c 'count_a "$@"' _
  echo "-- stage stats_A"; jobs | xargs -P 8 -n 2 bash -c 'stats_ab A "$@"' _
  echo "-- stage stats_B"; jobs | xargs -P 8 -n 2 bash -c 'stats_ab B "$@"' _
  echo "-- stage sketch_E"; mkdir -p "$ROOT/sketch"; jobs | xargs -P 8 -n 2 bash -c 'sketch_e "$@"' _
  echo "-- stage stats_E"; jobs | xargs -P 8 -n 2 bash -c 'stats_e "$@"' _
  echo "-- stage pilea_C"; jobs | xargs -P 6 -n 2 bash -c 'pilea_c default "$@"' _
  jobs | xargs -P 6 -n 2 bash -c 'pilea_c relaxed "$@"' _
  echo "ALL_DONE pipeline seed=$S $(date)"
} > "$MS/pipeline_$S.log" 2>&1
touch "$MS/PIPELINE_$S.Done"
