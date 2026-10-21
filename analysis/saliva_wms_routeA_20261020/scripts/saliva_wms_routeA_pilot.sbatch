#!/usr/bin/env bash
#SBATCH -p amd
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --job-name=salivaWMS-rA-pilot
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/real_routeB_20261012/logs/saliva_wms_routeA_pilot_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/real_routeB_20261012/logs/saliva_wms_routeA_pilot_%j.err
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
ROOT=$BASE/bench/real_routeB_20261012
export SK2BGROW_BIN=$BASE/src_20261012/target/release/sk2bgrow
export SK2BGROW_PY=$BASE/micromamba/envs/sk2bgrow/bin/python
export PYTHONPATH=$BASE/src/python:${PYTHONPATH:-}
"$SK2BGROW_PY" "$BASE/scripts/real_routeB_20261012/public_routeA_parallel_wms.py" \
  --manifest "$ROOT/analysis/saliva_wms/manifests/saliva_wms_pilot.tsv" \
  --db "$BASE/bench/ecc_routeB/db_gtdb_bcgI_sp" \
  --scratch "$ROOT/analysis/saliva_wms/scratch_pilot" \
  --threads-per-sample 16 --workers 1 \
  --status "$ROOT/logs/saliva_wms_routeA_pilot_status.tsv"
