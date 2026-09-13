#!/usr/bin/env bash
#SBATCH --job-name=P2_capdedup
#SBATCH --partition=amd
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/logs/P2/capdedup.%j.out
set -euo pipefail
/lustre1/g/aos_shihuang/sk2bgrow-hpc/micromamba/envs/sk2bgrow/bin/python -u \
  /lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/P2_sigma/scripts/cap_dedup_sens.py
