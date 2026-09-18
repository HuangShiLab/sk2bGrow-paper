# C5 explicit sorted-policy refusion provenance

## Purpose

This sensitivity arm asks how much of the legacy C5 QC recall depended on
deliberately using sorted-rank regression, while holding the rest of the
statistics layer at the current review-final policy. It is **not** the default
and **not** an exact reconstruction of legacy C5.

## Runs

- HPC branch: `review-final`
- Code commit: `929f4c2` (`/lustre1/g/aos_shihuang/sk2bgrow-hpc/src`)
- Initial sorted array: SLURM `4077173`; tasks 0–7 completed.
- Task 8 retry after a compute-node communication failure: SLURM `4077222`;
  completed.
- Dependent aggregation job: SLURM `4077223`; completed.
- Aggregation log: `hpc_agg_4077223.log`.

## Method

For each of the nine C5 samples, retained `windows.rates.tsv` rows with finite
positive `log2_se` were fitted with `fit_windows(method='sorted')`, fused and
QC-scored with current code. Counting was not rerun. Because current fusion,
r² exclusion, signed estimates and QC differ from the original C5 statistics,
these outputs quantify the sorted fallback under current policy rather than
reproducing the legacy arm exactly.

## Result

Explicit sorted refusion passed QC for 108 of 4,698 genome×sample observations
(2–30 per sample), compared with 26 under current `auto` and 484 in legacy C5.
The legacy negative association between QC and contig count was not reproduced
(raw/partial Spearman ρ for contig count: +0.043/+0.078). Relative to legacy
C5, the median absolute change in log₂PTR was 1.70 and 514 of 4,698 QC
decisions changed.

## Files

- `review/c5_sorted_results.tsv`: merged genome×sample outputs.
- `review/c5_sorted_sample_summary.tsv`: per-sample QC and fit summaries.
- `review/c5_sorted_qc_x_quality.tsv`: row-matched MAG-quality data.
- `review/c5_sorted_coverage_control.tsv` and `.txt`: associations.
- `review/c5_sorted_vs_previous.tsv` and `.txt`: direct legacy comparison.
