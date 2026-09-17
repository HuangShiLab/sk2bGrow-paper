#!/usr/bin/env bash
# F3 step 3: run the sk2bgrow counting layer at --max-mismatch 0/1/2 on the
# F3 simulated reads (error rates 0 / 0.001 / 0.01), enzyme panel k=16, against
# the existing F2 index databases. Idempotent: existing stats.json are skipped.
#SBATCH --job-name=F3_count
#SBATCH --partition=intel
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/F3/count.%j.out
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
F2=$BASE/bench/F2
F3=$BASE/bench/F3
BIN=$BASE/src/target/release/sk2bgrow
mkdir -p "$F3/cells" "$BASE/logs/F3"

for fq in "$F3"/fq/*.fq.gz; do
  name=$(basename "$fq" .fq.gz)          # {gid}.{errtag}
  gid=${name%%.*}
  errtag=${name##*.}
  for mm in 0 1 2; do
    out="$F3/cells/$gid/$errtag/mm$mm"
    mkdir -p "$out"
    [ -f "$out/$name.stats.json" ] && { echo "skip $name mm$mm"; continue; }
    "$BIN" profile "$fq" -d "$F2/db/${gid}_k16" -o "$out" \
      --no-stats --threads 4 --max-mismatch "$mm" --quiet \
      && echo "ok $name mm$mm" || echo "FAIL $name mm$mm"
  done
done
echo F3_COUNT_DONE
