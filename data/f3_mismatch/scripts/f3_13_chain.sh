#!/usr/bin/env bash
# F3 chain driver: simulate reads, then Rust counts, then replica attribution.
#SBATCH --job-name=F3_chain
#SBATCH --partition=intel
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/F3/chain.%j.out
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
F2=$BASE/bench/F2
F3=$BASE/bench/F3
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
BIN=$BASE/src/target/release/sk2bgrow
mkdir -p "$F3/fq" "$F3/cells" "$F3/parts" "$BASE/logs/F3"

echo "== simreads =="
"$PY" "$F3/scripts/f3_11_simreads.py"

echo "== rust counts =="
for fq in "$F3"/fq/*.fq.gz; do
  name=$(basename "$fq" .fq.gz)
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

echo "== replica =="
tail -n +2 "$F2/genomes.tsv" | cut -f1 | \
  xargs -P 8 -I{} "$PY" "$F3/scripts/f3_13_replica.py" {}

echo "== aggregate =="
"$PY" "$F3/scripts/f3_14_aggregate.py"
echo F3_CHAIN_DONE
