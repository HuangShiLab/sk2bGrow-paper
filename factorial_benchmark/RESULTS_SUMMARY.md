# Factorial benchmark status

## Job status

- Final array: **4088324**, tasks 0-23: **24/24 COMPLETED**.
- Aggregate dependency: **4088325**: **COMPLETED**.
- Wall time: 13 s-1 min 58 s per array task; aggregate <1 min.
- Peak RSS: 76-81 MiB per task.
- Requested: 1 CPU, 2 GiB, 2 h per task.
- Persistent output: 196 KiB parts + 104 KiB long table + 92 KiB contrast
  table. No FASTQ or per-anchor intermediates were retained.

## Design covered

- strains: 4/8/16/32
- per-strain depth: 0.5/1/2/4/8x
- shared-anchor fraction: 0/0.05/0.15
- abundance: even/10:1/100:1
- landmark source: coordinate-grid 2bRAD-like/random FracMinHash-like
- estimator: coordinate V-fit/sorted-rank
- 10 replicates; 54,000 strain-level estimates per estimator

## Headline observations

These are **mechanistic count-level simulation** results, not an end-to-end
replacement for the read-based Zheng/Pilea benchmark.

1. **Estimator.** Under the no-shared-anchor control, coordinate V-fit was much
   more accurate than sorted-rank at shallow depth (mean RMSE 0.115 versus
   0.641 log2 units at 0.5x). Both improved with depth, but sorted-rank remained
   worse over the tested range.
2. **Shared-anchor confusion.** Coordinate V-fit was nearly unbiased when
   anchors were private, but became negatively biased as shared anchors
   increased. At 8x and 15% shared anchors, mean bias reached -1.32 log2 units.
   This supports treating shared/ambiguous anchors as a distinct accuracy
   axis, rather than collapsing it into ordinary depth.
3. **Depth and confusion are not interchangeable.** Increasing depth reduced
   sampling error in the no-shared control but amplified the shared-anchor
   bias term in this simulator. Thus a deeper wrong or ambiguous observation
   profile does not solve multi-strain assignment.
4. **Abundance imbalance.** For coordinate V-fit, the mean RMSE increased from
   0.504 (even) to 0.557 (10:1) and 0.593 log2 units (100:1). Sorted-rank was
   more strongly affected: 0.304 to 0.380 and 0.624 log2 units, respectively.
5. **Landmark placement.** With no shared anchors and matched density, the
   coordinate-grid and random-sketch arms were equivalent (paired RMSE delta
   0.002 log2 units for coordinate V-fit; -0.014 for sorted-rank). Differences
   appeared only once shared anchors were introduced; this isolates shared-anchor
   handling rather than landmark density as the driver.

## Caution

The 2bRAD arm is an idealized coordinate grid, and the FracMinHash arm is an
idealized random sketch at equal density. Real enzyme motifs, GC bias, sequence
error, fragmentation, annotation error, and Pilea's shipped gates are not
modeled. Do not use these results alone to claim universal 2bRAD superiority.
The defensible conclusion is that shallow-depth performance and multi-strain
shared-anchor robustness are separable failure modes and should be reported
separately.

## HPC cleanup on 2026-09-20

- Removed completed Python `__pycache__` directories.
- Removed 1.17 GB of micromamba package tarballs/index cache.
- Compressed 57,637 old sketch-benchmark log files (123 MiB) to a verified
  1.8 MiB archive in `sk2bgrow-hpc/logs_archives/`, then removed originals.
- Current group usage was 46.13 TiB of 50 TiB after cleanup; the small increase
  relative to the earlier audit is from concurrent jobs/other group activity.
