#!/usr/bin/env bash
#SBATCH -p amd
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH --job-name=salivaWMS-sum
#SBATCH --output=/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/real_routeB_20261012/logs/saliva_wms_summarize_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/real_routeB_20261012/logs/saliva_wms_summarize_%j.err
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/sk2bgrow-hpc
ROOT=$BASE/bench/real_routeB_20261012
PY=$BASE/micromamba/envs/sk2bgrow/bin/python
"$PY" "$BASE/scripts/real_routeB_20261012/saliva_wms_summarize.py" \
  --root "$ROOT/analysis/saliva_wms" \
  --routeB-long "$ROOT/analysis/all_output_long.tsv.gz" \
  --routeB-metadata "$ROOT/manifests/all_routeB_metadata.tsv" \
  --outdir "$ROOT/analysis/saliva_wms/paired_analysis"
