# C5 fast-mode benchmark (k8 + mismatch 1)

## Aim

Test whether the C5 MAG runtime can be reduced without changing the current
conservative reporting policy, by using the top-8 enzyme panel and mismatch 1.
This is a shotgun (`--mode wms`) runtime benchmark. For real 2bRAD libraries,
mismatch 1 is not recommended because F3 showed that tolerant matching can open
a wrong-coordinate rescue channel.

## Design

- Dataset: C5 RBC metagenome, PRJNA974210, 522 bacterial MAG references.
- Full-depth sample: `SRR28338156`.
- Reads: paired-end shotgun; the profile-level read total was 67,423,986.
- Current baseline: 16 enzymes, mismatch 2, no screen, current conservative
  fragmented-reference policy (`C5_refuse_current`).
- Fast arm: top-8 ranked enzymes
  (`CjePI,CjeI,AlfI,Hin4I,HaeIV,BcgI,Bsp24I,BsaXI`), mismatch 1, no screen.
- Threads: 8.
- Job: SLURM 4095386 (`C5_k8_noscr`).
- The full-depth count TSV was deleted after successful statistics to conserve
  Lustre space; the per-genome output and summaries are retained.

## Runtime and footprint

See `full_sample_runtime.tsv`. The fast arm reduced wall time from 21.17 h to
1.26 h (**16.75x**) on one full-depth C5 sample. Peak RSS fell from 15.37 GB to
10.38 GB. Pilea defaults remain faster (5.72 min for this sample; the fast arm
is 13.25x slower than Pilea).

The current-policy runtime is a two-stage construction: the original
k16/mismatch-2 count stage plus the retained-window refit under current policy.
The fast arm is end-to-end (counting plus statistics).

## Count conservation and stage timing on a 1M-pair subset

See `subset_1M_k8_mm1_vs_k16_mm1_counts.tsv` and
`subset_1M_stage_timing_corrected.tsv`. Relative to 16-enzyme mismatch-1
counting, the 8-enzyme arm retained 76.58% of total assigned count (median
per-genome ratio 76.44%) and lost no detectable genomes. The log10 count totals
were highly correlated (r=0.9972). The database had 18,820,392 rather than
24,109,804 anchors (78.06%).

The earlier mixed-stage comparison has been superseded by a stage-level table.
For the 1M-pair subset, k16/mm1 used 111.25 s counting and 683.01 s statistics
(794.26 s end to end), whereas k8/mm1 used 45.14 s counting and 721.57 s
statistics (766.71 s end to end). Thus k8 mainly reduced counting, not the
Python statistics stage.

## Current-policy output agreement

See `full_sample_k8_mm1_vs_k16_current_contrast.tsv` and
`full_sample_output_summary.tsv`. Under the current fragmented-reference policy,
both arms returned eight finite PTRs. The current 16-enzyme/mismatch-2 arm
passed QC for four genomes and the fast arm for three; all three fast QC calls
were also QC calls in the current arm. Median absolute log2 PTR difference on
those three calls was 0.0194. Across all eight finite calls, Pearson r was 0.291
and median absolute difference was 0.0924. Thus the runtime reduction did not
create extra QC-pass calls, but the number of comparable calls is small.

## Screen result

A containment screen was tested separately on a 1M-pair subset
(`screen_equivalence_summary.tsv` and `screen_equivalence_per_genome.tsv`). It
retained 457/522 genomes with no extra genomes and 99.35% of count mass, but it
was not per-anchor neutral: 14,182 nonzero anchors were lost, 4,096 were gained,
and counts changed for 412/522 genomes. The full-depth screen arm took slightly
longer than no-screen (4,721 s versus 4,551 s) and gave the same three QC calls,
so the screen is not useful for this 522-MAG C5 benchmark. It remains relevant
for much larger GTDB-scale databases, where its measured lookup speedup was
21.1x and it reduced absent-reference false positives from 4,048/4,700 to zero.

## Interpretation

For shotgun users, k8 + mismatch 1 is a practical fast configuration, not yet a
validated shipped default: it is 16.75x faster than the current sk2bGrow C5
configuration and agrees on the few QC-pass calls. Confirmation is needed across
all nine C5 samples. For 2bRAD users, k8 rather than k16 is the relevant cost
reduction, but mismatch should remain 0 unless explicitly validated; the Zheng
panel-size benchmark already shows that 8 enzymes are sufficient (Table 7).
