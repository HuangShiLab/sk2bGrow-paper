#!/usr/bin/env bash
# F3 smoke test on one genome (g08 E. coli) before launching the full run.
#SBATCH --job-name=F3_smoke
#SBATCH --partition=intel
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/F3/smoke.%j.out
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
F2=$BASE/bench/F2
F3=$BASE/bench/F3
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
BIN=$BASE/src/target/release/sk2bgrow
GID=g08_Escherichia_coli
mkdir -p "$F3/fq" "$F3/cells" "$F3/parts" "$BASE/logs/F3"

echo "== census =="
"$PY" "$F3/scripts/f3_10_nearneighbor.py" "$GID"

echo "== simreads =="
"$PY" "$F3/scripts/f3_11_simreads.py" --only "$GID"

echo "== rust counts =="
for errtag in e0 e0001 e01; do
  fq="$F3/fq/$GID.$errtag.fq.gz"
  for mm in 0 1 2; do
    out="$F3/cells/$GID/$errtag/mm$mm"
    mkdir -p "$out"
    "$BIN" profile "$fq" -d "$F2/db/${GID}_k16" -o "$out" \
      --no-stats --threads 4 --max-mismatch "$mm" --quiet \
      && echo "ok $errtag mm$mm"
  done
done

echo "== replica =="
"$PY" "$F3/scripts/f3_13_replica.py" "$GID"

echo "== aggregate =="
"$PY" "$F3/scripts/f3_14_aggregate.py"
echo F3_SMOKE_DONE
