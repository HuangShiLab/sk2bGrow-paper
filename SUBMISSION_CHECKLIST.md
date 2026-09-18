# Submission checklist

## Scientific and statistical claims

- [x] Primary Zheng benchmark uses paired-end reads, signed fixed-origin fitting, and single-pass GC correction.
- [x] Residual two-pass GC correction is excluded from the primary benchmark and retained only as a diagnostic.
- [x] Matched-density landmark-source comparisons are described as bootstrap non-significance rather than equivalence.
- [x] Title and framing identify coordinate-aware V-fitting, not landmark determinism, as the source of the shallow-depth gain.
- [ ] Update the C5 section to distinguish legacy C5 outputs from the current-policy refusion.
- [ ] Add the explicit sorted-estimator C5 sensitivity arm once HPC array 4077173 finishes and is aggregated.

## Numerical and consistency checks

- [ ] Verify every abstract number against committed tables.
- [ ] Verify Results, figure captions, and table captions use the same primary-grid instance counts.
- [ ] Check that C5 cost, recall, and QC analyses state the exact code policy used.
- [ ] Regenerate all figures and tables after the final C5 update.
- [ ] Render manuscript.md to DOCX/PDF and inspect page layout.

## Manuscript packaging

- [x] Replace the draft header with a revision provenance note.
- [x] Add a formal References section.
- [ ] Verify all author names, affiliations, ORCID IDs, and corresponding author details.
- [ ] Complete journal-specific title-page requirements.
- [ ] Finalize data and code availability text after repository tagging/release.
- [ ] Complete a figure/table caption audit.

## Code and data reproducibility

- [x] Copy current-policy C5 refusion outputs and aggregate statistics from HPC.
- [ ] Copy explicit sorted-policy C5 sensitivity outputs from HPC.
- [ ] Record HPC SLURM job IDs and code commits for every final analysis.
- [ ] Reconcile the local, HPC, and GitHub sk2bGrow code branches.
- [ ] Tag and archive the exact code release.
- [ ] Push the final paper repository.

## Target and cover letter

- [x] Add a modular journal cover-letter draft.
- [ ] Select target journal and tailor formatting and word/figure limits.
- [ ] Confirm whether real-community disease/application data will be added before submission or reserved for follow-up.

