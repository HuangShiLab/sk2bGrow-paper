# C5 current-policy refusion provenance

## Purpose

The original C5 result tables predated the signed fixed-origin statistics and
the later conservative treatment of fragmented references. These runs reuse the
retained per-window rate tables and rerun only the statistics/fusion stage, so
all changes below isolate the stats-stage policy from read counting and cost.

## Runs

- HPC branch: `review-final`
- Code commit: `929f4c2` (`/lustre1/g/aos_shihuang/sk2bgrow-hpc/src`)
- Current-policy refusion array: SLURM `4076617`, 9/9 COMPLETED, exit 0.
- Current-policy aggregation job: SLURM `4076831`, COMPLETED, exit 0.
- Explicit sorted-policy sensitivity array: SLURM `4077173`; tasks 0–7 completed.
  Task 8 was retried as SLURM `4077222` after a compute-node communication failure
  and completed successfully. Sorted-policy aggregation job: SLURM `4077223`,
  dependent on `afterok:4077222`, completed successfully.

## Method

For each of the nine C5 metagenome samples, `windows.rates.tsv` was filtered to
rows with finite positive `log2_se`, fitted with `fit_windows`, fused with the
current fusion rules, and assembled with the current report/QC code. The
current-policy arm uses `method='auto'`; therefore fragmented MAGs no longer
fall back to sorted-rank fitting. The sensitivity arm explicitly uses
`method='sorted'` to quantify how much of the legacy C5 recall depended on that
fallback.

Coverage, dispersion, fraction, and containment fields in the refusion outputs
are inherited from the original C5 count-stage records because the refusion-only
job does not recompute count-stage statistics. PTR estimates and QC decisions
are recomputed.

## Current-policy result

The current `auto` policy returns rows for all 522 MAGs, but passes QC for only
26 genome×sample observations across the nine samples (range 1–6 per sample),
compared with 484 in the legacy C5 arm. This is not a compute failure: it is the
intended consequence of refusing unsafe sorted-rank fallback on fragmented
references. Consequently, the legacy C5 QC/fragmentation correlations cannot be
silently presented as a current-default result.

## Files

- `aggregate_qc_np.py`: current-policy aggregation script.
- `hpc_review/c5_current_results.tsv`: current-policy genome×sample outputs.
- `hpc_review/c5_current_sample_summary.tsv`: per-sample counts.
- `hpc_review/c5_current_vs_previous.tsv`: row-matched comparison with legacy C5.
- `hpc_review/c5_current_coverage_control.tsv`: raw and coverage-partial
  Spearman associations with MAG quality.
- `hpc_agg_4076831.log`: HPC aggregation report.
- Sorted-policy outputs and provenance are in `../c5_refuse_sorted/`.
