#!/usr/bin/env bash
#SBATCH --job-name=dedup_mech
#SBATCH --partition=amd
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/dedup_mech/dedup_mech.%j.out
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
D=$BASE/bench/dedup_mech
$PY $D/dedup_mech.py $D/dedup_mechanism
echo "DEDUP_MECH_DONE"
