# Submission checklist

## Target journal

- [x] Target journal selected: *Microbiome* (BMC/Springer Nature).
- [x] Manuscript reordered to BMC Research-article flow: Background, Methods, Results, Discussion, Conclusions.
- [x] BMC-style Declarations added.
- [x] Cover letter targeted to *Microbiome*.
- [x] Reviewer-suggestion worksheet added; at least three candidates included.
- [ ] Confirm author names, affiliations, ORCID IDs, and corresponding-author details.
- [ ] Confirm funding statement.
- [ ] Confirm reviewer institutional emails and conflict checks.

## Scientific and statistical claims

- [x] Primary Zheng benchmark uses paired-end reads, signed fixed-origin fitting, and single-pass GC correction.
- [x] Residual two-pass GC correction is excluded from the primary benchmark and retained only as a diagnostic.
- [x] Matched-density landmark-source comparisons are described as bootstrap non-significance rather than equivalence.
- [x] Title and framing identify coordinate-aware V-fitting, not landmark determinism, as the source of the shallow-depth gain.
- [x] Count-level factorial is labelled post hoc and is not represented as read-level validation.
- [x] C5 legacy outputs are distinguished from current-policy refusion.
- [x] Explicit sorted-estimator C5 sensitivity arm is reported.

## Numerical and consistency checks

- [x] Verify every abstract number against committed tables.
- [x] Verify Results, figure captions, and table captions use the same primary-grid instance counts.
- [x] Check that C5 cost, recall, and QC analyses state the exact code policy used.
- [x] Regenerate all figures and tables after the final C5 update.
- [x] Render a *Microbiome* review draft to DOCX.
- [ ] Complete final visual QA after embedding all tables and verifying every page.
- [ ] Complete a figure/table caption audit against journal length limits.

## Manuscript packaging

- [x] Add a title-page scaffold.
- [x] Add keywords and a list of abbreviations.
- [x] Move data/code availability into BMC Declarations.
- [ ] Finalize author metadata and CRediT contributions.
- [ ] Finalize acknowledgements and funding.
- [ ] Apply final BMC production formatting at revision if requested.

## Code and data reproducibility

- [x] Copy current-policy C5 refusion outputs and aggregate statistics from HPC.
- [x] Copy explicit sorted-policy C5 sensitivity outputs from HPC.
- [x] Record HPC SLURM job IDs and code commits for the primary Zheng grid and C5 refusion arms.
- [ ] Reconcile the local, HPC, and GitHub sk2bGrow code branches.
- [ ] Tag and archive the exact code release.
- [ ] Push the final paper repository after author metadata and final QA.

## Remaining decisions before submission

- [ ] Confirm whether real-community disease/application data will be added before submission or reserved for follow-up.
- [ ] Decide whether figures and tables must be embedded at first mention for the initial submission or supplied in the current end-of-document layout, which BMC accepts as a flexible review draft.
