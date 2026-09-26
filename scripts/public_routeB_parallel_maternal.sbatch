#!/usr/bin/env bash
#SBATCH -p amd
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 8:00:00
#SBATCH -J pub-routeB-par
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
ROOT=$BASE/bench/public_2brad
export SK2BGROW_BIN=$BASE/src_m4/target/release/sk2bgrow
export SK2BGROW_PY=$BASE/micromamba/envs/sk2bgrow/bin/python
export PYTHONPATH=$BASE/src/python:${PYTHONPATH:-}
"$SK2BGROW_PY" "$BASE/scripts/public_routeB_profile_parallel.py" \
  --manifest "$ROOT/profile_manifest_maternal.tsv" \
  --db "$BASE/bench/ecc_routeB/db_gtdb_bcgI_sp" \
  --scratch "/tmp/$USER/sk2bgrow_public_routeB_parallel" \
  --threads-per-sample 2 --workers 8 \
  --status "$ROOT/logs/maternal_status.tsv"
