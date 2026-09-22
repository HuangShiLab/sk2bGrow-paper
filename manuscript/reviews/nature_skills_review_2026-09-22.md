# Nature-skills review package, 2026-09-22

## Execution record

- `nature-reviewer` was run as three isolated reviewer agents. Reviewer 1
  emphasized methodological validity and statistical inference. Reviewer 2
  emphasized computational reproducibility, runtime and scalability. Reviewer 3
  emphasized relevance and readability for Microbiome readers.
- `nature-statistics` guidance was used to audit experimental units, effect
  summaries, uncertainty and the factorial reporting mismatch.
- `nature-polishing` guidance was used to tighten the abstract and claim
  hierarchy.
- `nature-data` guidance was used to add the C5 fast-mode benchmark artifacts to
  the data availability and repository record.
- No reviewer saw another completed report before submitting its own report.

# Reviewer 1

## Readiness summary

The manuscript is promising and unusually transparent about negative controls,
code-policy changes, survivorship and previous claim corrections. However, the
primary real-data titration uses one bacterial strain and one complete
reference. The factorial and C5 fast-mode analyses lack complete inferential
reporting. Several claims extend beyond the sampled systems. The C5 fast-mode
conclusion is based on one full-depth sample and three shared QC-pass calls.

## Who would be interested in the results, and why

Developers of metagenomic growth-rate methods, microbial ecologists seeking
replication signals from low-depth data, and bioinformaticians working on
reference-guided counting would be interested. The manuscript is also relevant
to readers concerned with policy-dependent recall.

## Major strengths

1. Pilea is evaluated with shipped gates and gates disabled.
2. The landmark-source by estimator design directly tests attribution.
3. Measured growth rate and predicted log2PTR are separated.
4. Negative results and retired claims are reported.
5. The factorial separates depth from shared-anchor ambiguity.
6. The C5 fast-mode benchmark reports runtime, output and count-retention
   evidence.

## Major Concerns

### R1-M1

Severity Major. Blocking No. Axis technical soundness and generalization.

Claim pointer. Abstract and Results 2.1.
Evidence pointer. Methods 3.1.
Concern. The primary real-data titration uses one E. coli strain and one
complete reference, but parts of the abstract imply broader estimator behavior.
Why it matters. Chromosome architecture, GC, plasmids and fragmented MAGs are
not represented in that titration.
Resolution test. Bound the claim to the Zheng system or add independent
multi-genome and community benchmarks.

### R1-M2

Severity Major. Blocking Yes before the factorial claims are finalized. Axis
statistical reporting and internal consistency.

Claim pointer. Results 4.1.
Evidence pointer. Methods 6.6, Table 11 and factorial summary tables.
Concern. The original text reported pooled values that did not match the
regular-coordinate values in Table 11. The table omitted random-arm, strain,
abundance and uncertainty summaries.
Resolution test. Reconcile pooled and arm-specific values, state the weighting,
report bootstrap intervals, and publish complete factorial summaries. This has
been addressed in the current revision by n-weighted summaries, bootstrap
intervals and supplementary tables.

### R1-M3

Severity Major. Blocking No if labelled a case study; Yes for a default-change
claim. Axis evidence adequacy.

Claim pointer. C5 fast-mode paragraph and Table 12.
Evidence pointer. `data/c5_fast_bench` and Table 12.
Concern. The full-depth C5 comparison uses one sample and three shared QC calls,
with no interval.
Resolution test. Run all nine C5 samples under current policy and fast policy, or
restrict the claim to a feasibility case study.

### R1-M4

Severity Major. Blocking No. Axis validity of QC policy.

Claim pointer. C5 current and legacy policy results.
Concern. C5 has no external PTR truth. Low QC recall cannot by itself
distinguish conservative rejection from over-aggressive rejection.
Resolution test. Validate QC against simulated known truth, monocultures,
replicates or independent biological evidence.

### R1-M5

Severity Major. Blocking No. Axis uncertainty.

Claim pointer. F6 and C5 fast-mode results.
Concern. Some newer conclusions originally lacked uncertainty summaries.
Resolution test. Add cell-level bootstrap intervals and complete all-sample C5
resources. The factorial bootstrap intervals have now been added.

### R1-M6

Severity Major. Blocking No. Axis claim boundary.

Claim pointer. Abstract and Results 4.
Concern. The read-level attribution and count-level simulator support related
but non-identical claims.
Resolution test. Label them complementary and avoid implying one data-generating
process.

## Minor Comments

R1-m1. Distinguish the nine-sample counting speedup from the one-sample
end-to-end fast-mode result.

R1-m2. Explain the mixed-stage 1M timing table. This was caused by comparing a
count-only k16 time with an end-to-end k8 time and has been corrected with a
stage-level table.

R1-m3. State that Table 12 is descriptive and has no intervals for three shared
QC calls.

R1-m4. Add a table of the three shared QC calls and the fourth current-only QC
call.

R1-m5. Explain the relationship between the factorial shared-anchor assignment
model and production EM/mapping.

R1-m6. Add sensitivity of the simulator to genome size, landmark layout and PTR
distribution.

R1-m7. Define genome loss operationally.

R1-m8. Define common returned sets for comparative accuracy.

R1-m9. Report survivor denominators beside 0.5x sketch estimates.

R1-m10. Clarify the M4 build label.

# Reviewer 2

## Readiness summary

The manuscript is transparent and provides useful C5 fast-mode provenance,
including sample accession, read count, enzyme panel, mismatch, threads, wall
time, RSS, outputs, QC calls and SLURM job 4095386. The factorial also has strong
job-level provenance. However, the full-depth fast-mode result is one sample.
The initial 1M subset runtime table mixed a count-only k16 time with an
end-to-end k8 time. The containment screen evidence was aggregate rather than
per-anchor. Reproducibility metadata should be strengthened.

## Who would be interested in the results, and why

Microbiome researchers interested in low-depth growth profiling, MAG panels and
scalable reference indexing would be interested.

## Major strengths

1. Legacy and current C5 policies are separated.
2. Runtime, RSS, threads, outputs and QC calls are reported for the fast mode.
3. The factorial has explicit SLURM provenance.
4. One-sample and all-nine-sample caveats are present.
5. The containment screen is tested for absent-genome false positives.
6. Storage and deletion decisions are disclosed.

## Major Concerns

### R2-M1

Severity Major. Blocking No for a case study, Yes for a default change. Axis
scalability.

Claim pointer. C5 fast-mode paragraph and Table 12.
Evidence pointer. `full_sample_runtime.tsv` and Table 12.
Concern. The full-depth fast-mode result uses one sample.
Resolution test. Run all nine C5 samples or an independent low/medium/high
presence set, with per-sample runtime, RSS, outputs and agreement.

### R2-M2

Severity Major. Blocking Yes until resolved. Axis reproducibility.

Claim pointer. Initial 1M subset timing table.
Evidence pointer. `subset_1M_k8_mm1_vs_k16_mm1_counts.tsv`.
Concern. The original table compared count-only k16 time with end-to-end k8
time.
Resolution test. Replace it with stage-level timing. This has now been done in
`subset_1M_stage_timing_corrected.tsv`: k16/mm1 used 111.25 s counting plus
683.01 s statistics, whereas k8/mm1 used 45.14 s counting plus 721.57 s
statistics.

### R2-M3

Severity Major. Blocking No for C5, Yes for database-scale screen claims. Axis
benchmark equivalence.

Claim pointer. Containment-screen discussion.
Evidence pointer. `screen_equivalence_summary.tsv` and per-genome table.
Concern. Aggregate count mass does not establish per-anchor neutrality.
Resolution test. Report anchor-level and genome-level retention. The new audit
reports 99.44% nonzero-anchor retention, 14,182 lost anchors, 4,096 gained
anchors, 412 changed genomes and 65 dropped genomes.

### R2-M4

Severity Major. Blocking No because the nine-sample table exists, but Yes until
explicitly cited. Axis provenance.

Claim pointer. Abstract mismatch-1 counting speedup.
Evidence pointer. `data/mm1_e2e/mm1_full_sample_cost.tsv`.
Concern. The nine-sample resource evidence should be cited directly.
Resolution test. Cite the machine-readable table in Results and availability.
This has now been added.

### R2-M5

Severity Major. Blocking No after clarification. Axis benchmark design.

Claim pointer. Table 12.
Concern. The current-policy baseline combines an original count run with a
retained-window refit.
Resolution test. State the stage construction explicitly. This has now been
added to the Results text and table caption.

### R2-M6

Severity Major. Blocking No. Axis statistical reporting.

Claim pointer. C5 output comparison.
Concern. Three shared QC calls are too few, and all-finite agreement is weaker
than shared-QC agreement.
Resolution test. Report all-finite correlation and difference alongside the
shared-QC difference. This has now been added.

### R2-M7

Severity Major. Blocking No. Axis reproducibility.

Claim pointer. Computational environment and availability.
Concern. Job IDs and commits are useful but not a full manifest.
Resolution test. Add a reproducibility manifest with commands, jobs, inputs,
outputs and retained artifacts.

## Minor Comments

R2-m1. Separate counting-only from end-to-end speedups.

R2-m2. Replace "roughly one order of magnitude" with measured ranges.

R2-m3. Explain NA fields in output summaries.

R2-m4. Add the full-depth screen runtime to machine-readable output.

R2-m5. Clarify that M4 is a code/build label, not hardware.

R2-m6. Retitle Table 12 as a paired current-versus-fast contrast.

R2-m7. State that the C5 fast-mode run omitted the screen.

R2-m8. Report genome-level count-retention quantiles as well as aggregate
retention.

R2-m9. Preserve simulator seed and deterministic metadata.

R2-m10. State retained and deleted C5 artifacts.

R2-m11. Keep one-sample and nine-sample ratios distinct in all tables.

R2-m12. Report screen scale, threshold, memory and build cost.

# Reviewer 3

## Readiness summary

The paper is candid and data rich, but currently reads partly as an internal
benchmark dossier. The core distinction from Pilea is well supported and
important. The manuscript needs clearer claim hierarchy, consolidated practical
guidance, and more accessible terminology for Microbiome readers.

## Who would be interested in the results, and why

Shotgun metagenome researchers, multi-strain community analysts, 2bRAD users and
developers of sketch-based PTR methods would be interested in coordinate-aware
shallow-depth fitting, shared-anchor ambiguity and practical deployment limits.

## Major strengths

1. The 2 by 2 attribution design separates estimator geometry from landmark
   source.
2. The factorial separates depth from shared-anchor ambiguity.
3. Negative results, compression, runtime and QC loss are transparent.
4. Fecal concordance and shallow subsampling connect the work to real microbiome
   data.
5. The exact-deduplication result is practical and mechanistically supported.

## Major Concerns

### R3-M1

Severity Major. Blocking No. Axis readability and inference.

Claim pointer. Abstract and Zheng results.
Evidence pointer. Gate-ablation results.
Concern. Reported fraction, gates-off sensitivity and conditional accuracy can be
misread.
Resolution test. State these three quantities separately in the abstract and
primary results.

### R3-M2

Severity Major. Blocking No. Axis readability.

Claim pointer. Results 4 and Discussion.
Evidence pointer. Attribution experiment.
Concern. The Pilea comparison is distributed across internal labels and dense
benchmark prose.
Resolution test. Add an early plain-language comparison and define all internal
labels.

### R3-M3

Severity Major. Blocking No. Axis practical relevance.

Claim pointer. Results and Conclusions.
Evidence pointer. Panel-size, mismatch, deduplication, GC and C5 analyses.
Concern. Practical guidance is scattered.
Resolution test. Consolidate data-type-specific recommendations, evidence level
and limitations.

### R3-M4

Severity Major. Blocking No. Axis claim boundary.

Claim pointer. Fecal and C5 sections.
Evidence pointer. Three fecal samples and 26/4,698 current C5 QC observations.
Concern. Real-community validation and MAG-scale utility are limited.
Resolution test. State sample-size and current-policy yield limits in the
abstract and Discussion.

### R3-M5

Severity Major. Blocking No. Axis simulation boundary.

Claim pointer. F6 factorial.
Evidence pointer. Methods 6.6 and factorial design.
Concern. Exact factorial magnitudes depend on a simplified assignment model.
Resolution test. Label the magnitudes as mechanistic summaries and state what
cannot be inferred about production pipelines.

### R3-M6

Severity Major. Blocking No. Axis readability.

Claim pointer. Abstract.
Evidence pointer. Current abstract.
Concern. The abstract covers too many findings and dilutes the primary claim.
Resolution test. Focus on one primary claim, secondary deployment findings and
key limitations.

## Minor Comments

R3-m1. Define usable PTR quantitatively in the abstract.
R3-m2. Consider whether the title should foreground coordinate-aware fitting.
R3-m3. Remove internal revision provenance from the journal title page.
R3-m4. Move most repository paths to Data and Code Availability.
R3-m5. Add a terminology box for F1-F6, C5, arms and policies.
R3-m6. Pair "usable 1-2x" immediately with compression caveats.
R3-m7. Consider a compact gate-behaviour figure.
R3-m8. Explain the non-monotonic fecal subsampling correlation.
R3-m9. Warn earlier about sparse current MAG-scale QC yield.
R3-m10. Define genome loss precisely.
R3-m11. Describe the low-GC boundary jointly with genome size.
R3-m12. Contextualize the low per-anchor R2 in the deduplication mechanism.
R3-m13. Consider a prospective-user decision workflow figure.
R3-m14. Break up paragraphs longer than approximately 300 words.

# Cross-review synthesis

## Consensus strengths

All three reviewers found the attribution design, transparency and practical
C5 runtime experiment compelling. All identified coordinate-aware fitting, not
landmark determinism, as the paper's central mechanistic finding.

## Consensus blocking concerns

1. Resolve the original mismatch between pooled factorial text and the
   regular-arm Table 11. Addressed with n-weighted summaries, source-specific
   summaries and bootstrap intervals.
2. Resolve the mixed-stage 1M timing comparison. Addressed with
   `subset_1M_stage_timing_corrected.tsv`.
3. Keep the full-depth k8/mm1 C5 result explicitly a one-sample case study. The
   manuscript now does so.

## Other consensus major concerns

1. Complete k8/mm1 validation across all nine C5 samples before changing the
   shipped default.
2. Keep real 2bRAD route-B at mismatch 0 until separately validated.
3. Add screen equivalence beyond aggregate mass. Addressed with anchor- and
   genome-level audit tables.
4. Add reproducibility manifests for C5 benchmark commands and artifacts.
5. Make current C5 QC yield a central limitation rather than a late caveat.
6. Tighten the abstract and distinguish accuracy, recall and gate behavior.

## Minor revision checklist

1. Remove internal revision provenance from the title page.
2. Cite Table 12 in the C5 fast-mode paragraph.
3. Report all-finite as well as shared-QC C5 differences.
4. Explain the original 1M timing discrepancy.
5. Add screen-equivalence distributions.
6. Report factorial uncertainty summaries.
7. Define genome loss and comparative denominators in final production editing.
8. Add a terminology table and consolidated deployment guidance during the next
   full revision.

## Post-review manuscript actions

- The abstract was tightened from approximately 331 to 227 words.
- The internal revision-provenance line was removed from the journal title page.
- Table 12 was retitled as a paired current-versus-fast contrast.
- Table 12 now reports median QC-pass log2PTR, all-finite agreement and
  shared-QC agreement.
- The C5 screen discussion now reports anchor-level and genome-level equivalence.
- The 1M-pair subset timing table was replaced with a stage-level comparison.
- Factorial text now uses n-weighted summaries, source-specific summaries and
  10,000-cell bootstrap intervals.

## Broad-interest / significance readout

The paper is relevant to Microbiome readers interested in shallow metagenomes,
growth dynamics, MAG-scale computation and reduced-representation assays. The
most important contribution is the controlled demonstration that coordinate-aware
fitting, rather than 2bRAD identity, enables shallow-depth PTR inference.

## Risk / unsupported claims

- Do not claim universal 2bRAD superiority.
- Do not present factorial magnitudes as production predictions.
- Do not treat the one-sample C5 fast-mode run as a validated default.
- Do not interpret current C5 QC loss as proof that all rejected estimates were
  wrong.
- Do not treat containment screening as per-anchor neutral on the C5 reference
  set.
