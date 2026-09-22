# Statistics review scope

- Input reviewed: latest `manuscript/manuscript.md`, Table 11 and Table 12,
  `factorial_benchmark/factorial_long.tsv`, C5 fast-mode summaries and
  `data/mm1_e2e/mm1_full_sample_cost.tsv`.
- Boundary: no raw count-level factorial replicates were retained, and no raw C5
  FASTQs were reprocessed for this audit.
- Study design readout: controlled Zheng attribution grid; post-hoc count-level
  factorial; C5 mismatch-1 nine-sample counting benchmark; one-sample C5
  end-to-end fast-mode case study.
- Independent unit and replication readout: Zheng media are the bootstrap unit;
  factorial strain-level estimates are pooled with cell-level bootstrap of
  factorial cells; C5 samples are the unit for runtime and mismatch benchmarks.

# Major statistical issues

## P0. Factorial text and Table 11 used different summaries

Evidence: Results 4.1 originally reported pooled landmark arms, whereas Table 11
reported the regular-coordinate arm.

Why it matters: readers could infer that the large random-arm effects were
represented by the much smaller regular-arm table values.

Fix applied: Results now states that values are n-weighted pooled summaries.
Table 11 caption states that it shows the regular-coordinate arm only. New
files `tables/table11_pooled_sources_summary_ci.tsv` and
`tables/table11_by_source_summary_ci.tsv` provide pooled, arm-specific and
bootstrap summaries.

## P1. Factorial estimates lacked uncertainty

Evidence: ten replicates per cell were summarized as point means.

Fix applied: pooled factorial summaries now use cell-level bootstrap 95%
intervals from 10,000 resamples. These intervals describe factorial-cell
variability, not read-level pipeline validation.

## P1. C5 fast-mode inference is sparse

Evidence: one full-depth sample yielded eight finite calls, four current QC
passes, three fast QC passes and three shared QC calls.

Fix applied: Table 12 now reports all-finite Pearson r, all-finite median
absolute difference and shared-QC median absolute difference. The text labels
the result a case study and requires all-nine-sample confirmation before a
default change.

## P2. Screen equivalence was initially aggregate only

Evidence: earlier reporting used genome counts and total count mass.

Fix applied: new anchor- and genome-level tables report 99.44% nonzero-anchor
retention, 14,182 lost anchors, 4,096 gained anchors, 412 changed genomes and
65 dropped genomes. The manuscript now says the screen is not per-anchor neutral
on this reference set.

# Ready-to-paste revision

The revised Results text now reports n-weighted pooled factorial summaries with
bootstrap intervals, explicitly separates regular-coordinate from random-sketch
effects, and labels C5 fast-mode results as a one-sample deployment case study.

# AUTHOR_INPUT_NEEDED

- None for the current descriptive statements.
- All-nine-sample C5 resources are needed before making a shipped-default claim.

# Reviewer-risk note

A statistical reviewer may still ask for a prespecified factorial model with
strain, depth, shared-anchor fraction, abundance ratio and landmark-source
terms, and for paired per-genome C5 uncertainty once all nine samples are
available.
