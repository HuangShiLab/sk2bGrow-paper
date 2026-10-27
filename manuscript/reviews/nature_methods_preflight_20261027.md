# Nature Methods preflight review and revision log (2026-10-27)

## Scope

This review covers `manuscript/manuscript.md`, the compact Nature Methods
submission file, and `manuscript/manuscript_full_benchmark_20261027.md`, the
complete benchmark compendium retained for provenance, supplementary tables and
detailed methods.

## Recommendation

**Conditional accept for submission preparation.** The scientific narrative,
data record, benchmark evidence and revised claims are now in a compact Nature
Methods Brief Communication format. The manuscript is not yet literally
submittable because author names, affiliations, corresponding author, funding,
CRediT contributions and the saliva ethics/consent statement require human
entry. These fields must not be invented.

## Key changes made

1. Changed the target journal from *Microbiome* to *Nature Methods*.
2. Reorganised `manuscript.md` into a compact submission file. The main text is
   approximately 1,635 words including Summary, Introduction, Results,
   Discussion and Conclusions, but excluding Methods, Declarations, References
   and Supplementary information. The previous full benchmark compendium is
   preserved in `manuscript_full_benchmark_20261027.md`.
3. Reduced the Summary to 154 words.
4. Changed Results from a long benchmark log into five claim-linked sections:
   isolate accuracy, attribution, fragmentation, real metagenomes and runtime.
5. Added the matched saliva WMS/2bRAD-M cohort as Methods section 3.8 and
   Results discussion, with explicit WMS and 2bRAD parameters.
6. Corrected the saliva sample mapping: `S11_1` means 11 AM for participant 1.
7. Separated Data availability and Code availability.
8. Added PRJNA1131785 and the matched saliva analysis path to the
   availability statements.
9. Added Table 14 for the WMS/2bRAD concordance analysis.
10. Revised the abstract, introduction, discussion and conclusions to avoid
    claiming absolute PTR validation from the saliva cohort.

## Numerical consistency checks

- Summary word count: 154.
- Main text excluding Methods, Declarations, References and Supplementary
  information: approximately 1,635 words.
- E. coli benchmark: r = 0.92 and 0.96 at 1x and 2x.
- Saliva WMS GLM: 13,123 finite estimates, 9,112 exploratory robust estimates,
  361 strict-QC estimates.
- Matched 2bRAD GLM: 14,040 finite estimates, 8,785 exploratory robust
  estimates, 26 strict-QC estimates.
- Shared finite WMS-GLM/2bRAD-GLM pairs: 11,343; Pearson r = 0.460;
  Spearman rho = 0.436.
- Both exploratory robust pairs: 6,434; Pearson r = 0.418;
  Spearman rho = 0.403; median absolute difference 0.361 log2PTR.
- Both strict QC pairs: 2.

## Remaining human-input blockers

These are the only reasons the file cannot yet be uploaded verbatim:

1. Complete author list, affiliations, corresponding author and ORCID.
2. Funding statement.
3. CRediT author contributions.
4. Saliva cohort ethics approval body, protocol identifier and consent
   statement.

## Editorial-risk assessment

The strongest claim is now coordinate-aware fitting rather than landmark
choice. The saliva experiment is described as exploratory assay concordance,
not growth-rate truth. The WMS-auto policy is shown to be weak, and the strict
QC overlap is explicitly stated. These wording changes prevent the new data
from overreaching its validation role.
