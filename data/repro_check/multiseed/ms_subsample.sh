#!/usr/bin/env bash
# multiseed subsampling: replicates of the C1 Zheng grid at new seeds.
# Same calibration as 01_subsample.sh: frac = depth*GLEN/base_count (manifest
# col3, PE total); seqkit sample -s <INSTANCE_SEED> on both mates (same seed
# -> mates stay in sync, as C1 verified). One instance seed per run of this
# script. Info tsv per cell records exact read counts + mate-sync check.
set -uo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
SEQKIT=$BASE/tools/seqkit
MANI=$BASE/manifests/PRJNA615952_counts.tsv
GLEN=4641652
FQDIR=$BASE/fastq/PRJNA615952
SEED=$1   # instance seed, e.g. 20260912
OUT=$BASE/bench/repro_check/multiseed/seed$SEED/fq
mkdir -p "$OUT"
DEPTHS=(0.5 1 2 5 10)
CODES=(50 100 200 500 1000)
export BASE SEQKIT MANI GLEN FQDIR SEED OUT

sub_one() { # run depth
  local run=$1 depth=$2
  local name="${run}.${depth}x"
  local bases frac
  bases=$(awk -v r="$run" '$1==r{print $3}' "$MANI")
  frac=$(awk -v d=$depth -v g=$GLEN -v b=$bases 'BEGIN{printf "%.10f", d*g/b}')
  for mate in 1 2; do
    local o="$OUT/${name}_${mate}.fq.gz"
    [ -s "$o" ] && continue
    local tmp="${o%.fq.gz}.tmp.fq.gz"
    "$SEQKIT" sample -j 4 -p "$frac" -s "$SEED" "$FQDIR/${run}_${mate}.fastq.gz" -o "$tmp" \
      && mv "$tmp" "$o" || { echo "FAIL sample $name $mate"; return 1; }
  done
  local n1 n2 h1 h2 sync
  n1=$("$SEQKIT" stats -T "$OUT/${name}_1.fq.gz" | awk 'NR==2{print $4}')
  n2=$("$SEQKIT" stats -T "$OUT/${name}_2.fq.gz" | awk 'NR==2{print $4}')
  h1=$("$SEQKIT" seq -n -i "$OUT/${name}_1.fq.gz" | md5sum | cut -d" " -f1)
  h2=$("$SEQKIT" seq -n -i "$OUT/${name}_2.fq.gz" | md5sum | cut -d" " -f1)
  sync=FAIL; [ "$n1" = "$n2" ] && [ "$h1" = "$h2" ] && sync=OK
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$run" "$SEED" "$depth" "-" "$SEED" "$frac" "$n1" "$n2" "$sync" \
    > "$OUT/${name}.info.tsv"
  echo "$name n1=$n1 n2=$n2 sync=$sync"
}
export -f sub_one

MEDIA="M6 M3 M2 RUN_OUT M27 M25 M24 M22 M23 M19 M18 M17 M13 M12 M1 M10 M4"
{
  echo "start seed=$SEED $(date)"
  for m in $MEDIA; do
    srr=$(awk -v m="$m" '$2==m{print $1; exit}' "$BASE/bench/C1/run_medium.tsv")
    for d in "${DEPTHS[@]}"; do echo "$srr $d"; done
  done | xargs -P 6 -n 2 bash -c 'sub_one "$@"' _
  echo "ALL_DONE seed=$SEED $(date)"
} > "$BASE/bench/repro_check/multiseed/subsample_$SEED.log" 2>&1
touch "$BASE/bench/repro_check/multiseed/SUBSAMPLE_$SEED.Done"
