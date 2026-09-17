#!/usr/bin/env bash
# F3 census driver: near-collision census for all 18 genomes (xargs -P 8).
#SBATCH --job-name=F3_nn
#SBATCH --partition=intel
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/F3/nn.%j.out
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
F3=$BASE/bench/F3
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
mkdir -p "$F3/parts" "$BASE/logs/F3"

tail -n +2 "$BASE/bench/F2/genomes.tsv" | cut -f1 | \
  xargs -P 8 -I{} "$PY" "$F3/scripts/f3_10_nearneighbor.py" {}
echo F3_NN_DONE
