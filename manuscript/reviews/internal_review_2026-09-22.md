# Internal review — 2026-09-22

## Scope

Fresh referee-style review of `manuscript/manuscript.md` after addition of:

1. F6, the count-level factorial decomposition of estimator, depth and
   shared-anchor ambiguity; and
2. the post-hoc C5 fast-mode benchmark using the top-8 enzyme panel and
   mismatch 1.

The review also checked Markdown-to-DOCX packaging, figure/table numbering,
data availability, numerical claims added in the latest revisions, and whether
the manuscript avoids overclaiming relative to the new analyses.

## Overall assessment

The paper is much stronger than the earlier version because it now separates
three claims that were previously entangled:

1. shallow-depth PTR extraction is driven by coordinate-aware V-fitting;
2. 2bRAD supplies a wet-lab-realizable, QC-capable landmark system rather than
   universal statistical superiority; and
3. shared-anchor ambiguity is a distinct multi-strain failure mode.

The C5 fast-mode experiment is a useful addition for shotgun users: k8 plus
mismatch 1 reduced one full-depth sample from 21.17 h to 1.26 h while preserving
the current conservative policy. The manuscript correctly labels this as a
one-sample deployment benchmark rather than a new primary grid.

Recommendation: **acceptable after minor/substantive-polish revision**, with the
blocking submission items listed below.

## Major points

### 1. Validate k8/mismatch-1 beyond one C5 sample

The new C5 result is compelling but uses one full-depth sample. It reduced
runtime from 21.17 h to 1.26 h and preserved three of the four current-policy
QC calls, with median absolute log2 PTR difference 0.0194 on shared QC calls.
This is enough to justify a candidate fast mode, not yet a shipped default. The
manuscript now calls it a “candidate default”; that wording should be retained.

Recommended action: run the same paired k16/mm2 versus k8/mm1 comparison across
all nine C5 samples before release, or explicitly list this as pending
validation.

### 2. Do not imply that 8-enzyme validation transfers automatically to real 2bRAD

The 8-enzyme result is strongest for Zheng shotgun simulations and C5 shotgun
runtime. The real-2bRAD data are BcgI libraries and do not test a wet-lab
8-enzyme 2bRAD protocol. In addition, mismatch 1 was benchmarked on shotgun
data; route-B reads should remain mismatch 0 until separately validated. The
manuscript states this, but the cover letter and abstract should not imply that
real multi-enzyme 2bRAD PTR accuracy has already been demonstrated.

### 3. Separate statistical attribution from deployment engineering

The Zheng attribution and F6 factorial answer a mechanistic question: the
coordinate-aware estimator, not landmark determinism, drives the shallow-depth
gain. The C5 fast-mode benchmark answers an engineering question: a smaller
panel and lower mismatch can make MAG-scale shotgun runs practical. These should
remain distinct in claims. The current Discussion mostly does this, but later
editing should preserve the distinction.

### 4. Conservative C5 QC remains a scientific limitation

Under current policy only 26 of 4,698 C5 observations pass QC. The fast-mode
result does not solve this; it reduces runtime without inflating QC yield. The
manuscript should continue to avoid describing this as “poor sensitivity” alone.
The QC loss is an intended consequence of refusing unsafe coordinate fits on
fragmented references. Prospective work still needs scaffold-first or
high-quality-reference strategies to increase valid yield.

### 5. Literature coverage is too narrow

There are only eight references. For *Microbiome*, additional PTR methods and
applications should be cited, including iRep, CoPTR, DEMIC/PtrRC or equivalent
approaches, and disease-cohort PTR applications. The absence of these methods is
the most significant scholarly gap.

## Minor points and mechanical audit

1. Table 12 had been duplicated in Markdown; fixed in this revision. The final
   manuscript contains exactly one Table 12 block.
2. Figure numbering is now 1–9 in narrative order and all image files exist.
3. Table numbering is now 1–12; all tables are declared in the Tables section.
4. The abstract is ~331 words; check the current *Microbiome* limit during final
   formatting.
5. The manuscript is ~24,200 words including tables. This may be acceptable for
   a methods paper but should be tightened before final submission.
6. The title page, affiliations, ORCID, corresponding author, funding, and
   CRediT statements remain placeholders and are submission blockers.
7. The References section is BibTeX-free manual numbered text; this is
   acceptable for review but should be checked against journal style.
8. The factorial analysis is count-level and uses a deliberately conservative
   shared-anchor assignment model. The manuscript should continue to avoid
   describing it as read-level validation.
9. The C5 screen is not useful at the 522-MAG presence rate in C5 but remains
   relevant for GTDB-scale databases. The manuscript states this conditional
   interpretation correctly.
10. All nine figures and twelve tables were present in the regenerated DOCX.

## Numerical spot checks

- F6: 360 non-estimator cells × 10 replicates × mean 15 strains = 54,000
  genome-level estimates per estimator. This matches the Methods text.
- C5 fast mode: 76,218.63 / 4,551.00 = 16.75× runtime reduction.
- C5 fast mode versus Pilea defaults: 4,551 / 343.46 = 13.25×.
- C5 k8 database: 18,820,392 / 24,109,804 = 78.06% of current anchors.
- Matched 1M-pair subset: k8 retains 76.58% of total count, median per-genome
  ratio 76.44%, no genomes lost, log10 count Pearson r = 0.9972.
- Current-policy comparison: current k16/mm2 returns 8 finite PTRs and 4 QC
  passes; k8/mm1 returns 8 finite PTRs and 3 QC passes; all 3 fast QC calls are
  shared; median absolute log2 PTR difference 0.0194.

## Packaging audit after fixes

- Figures: 9/9 present and correctly ordered.
- Tables: 12/12 present; duplicate Table 12 removed.
- Images in DOCX: 9.
- Tables recognized in DOCX: 24 blocks, because several tables have multiple
  lettered sub-blocks.
- No stale Fig. 5–9 filenames remain.
- No `git diff --check` whitespace errors.
